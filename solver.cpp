#include "solver.h"
#include "utils.h"
#include <cmath>
#include <time.h>
#include <limits>
#include <algorithm>
#include <bitset>

namespace cplex_tap {

    struct compare{
        int key;
        compare(int const &i): key(i) { }
        bool operator()(int const &i)
        {
            return (i == key);
        }};

    ILOMIPINFOCALLBACK1(MyCallback, IloInt, num) {
        IloEnv env = getEnv();
        double this_z = getIncumbentObjValue();
        double best_z = getBestObjValue();
        //TODO only print if we have better solution ....
        if (1)
            env.out() << "[INFO MIP Callback]" << " CLK " << clock() << " Z " << this_z << std::endl;
    }

    Solution Solver::solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production, bool seed, string warmStart) const {
        std::cout << "CLK_RATE " << CLOCKS_PER_SEC << std::endl;
        std::cout << "Starting Model generation ..." << std::endl;

        // Init CPLEX environment and model objects
        IloEnv env;
        IloModel model(env);

        const uint64_t n = tap.size();

        // Variables
        IloArray<IloNumVarArray> x(env, n + 2u);
        IloNumVarArray s(env, n);
        IloNumVarArray u(env, n);

        // Init variables
        init_vars(env, n, x, s, u);

        // Init Constraints
        build_constraints(dist_bound, time_bound, env, model, n, x, s);
        if (!progressive)
            build_subtour_const_all(env, model, n, x, u);

        //
        // --- Objective (1) ---
        //
        IloExpr expr(env);
        for (auto i = 0u; i < n; ++i) {
            expr += tap.interest(i) * s[i];
        }
        IloObjective obj(env, expr, IloObjective::Maximize);
        model.add(obj);
        std::cout << "Added Objective to model\n";

        // Free some memory 
        expr.end();

        //Init solver
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::TimeLimit, 1800);
        //cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap	, 0.00001);
        cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 16000);
        if (!production)
            cplex.setParam(IloCplex::Param::Threads, 1);
        else
            cplex.setParam(IloCplex::Param::Threads, 8);
        // Export model to file
        if (debug) {
            cplex.exportModel("tap_debug_model.lp");
            IloCplex::Callback mycallback = cplex.use(MyCallback(env, 10));
        }

        if (seed) {
            warm_start(warmStart, env, n, x, s, cplex);
        }

        bool solved = false;
        time_t start, end;
        start = clock();
        std::cout << "CLK_START " << start << std::endl;
        try {
            solved = cplex.solve();
        }
        catch (const IloException& e) {
            std::cerr << "\n\n--- CPLEX Exception ---\n";
            std::cerr << e << "\n";
            env.end();
            throw;
        }
        end = clock();
        double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;

        if (solved) {
            // If CPLEX successfully solved the model, print the results
            std::cout << "\n--- Solver success ---\n";
            std::cout << "    Status: " << cplex.getStatus() << "\n";
            std::cout << "    Objective: " << cplex.getObjValue() << "\n";
            if (debug)
                dump(cplex, x, env, s, u);
            print_solution(cplex, x);

            if (cplex.getStatus() == IloAlgorithm::Feasible && time_to_sol < 3595){
                return Solution(-time_to_sol, cplex.getObjValue(), get_solution(cplex, x));
            }

            if (progressive) {
                cout << "Looking for subtours in solution" << endl;
                vector<vector<int>> subtours = getSubtours(n, x, cplex);
                int p_iters = 0;
                while (!subtours.empty()) {
                    vector<int> nodes = flatten(subtours);
                    IloExpr exp(env);
                    IloModel toUpdate = cplex.getModel();
                    stringstream vname;
                    cout << "Adding constraints" << endl;
                    for (auto k = 0u; k < nodes.size(); ++k) {
                        for (auto l = 0u; l < nodes.size(); ++l) {
                            int i = nodes.at(k), j = nodes.at(l);
                            exp = ((n - 1) * (1 - x[i + 1][j + 1])) - u[i] + u[j]; // >= 1
                            vname << "mtz_" << i << "_" << j;
                            toUpdate.add(IloRange(env, 1, exp, IloInfinity, vname.str().c_str()));
                            vname.str("");
                            exp.clear();
                        }
                    }
                    start = clock();
                    cplex.solve();
                    end = clock();
                    time_to_sol += (double) (end - start) / (double) CLOCKS_PER_SEC;
                    cout << "Looking for subtours in solution" << endl;
                    subtours = getSubtours(n, x, cplex);
                    p_iters++;
                }
                print_solution(cplex, x);
                cout << "Progressive Mode : iterations = " << p_iters << endl;
            }
        }
        else {
            std::cerr << "\n--- Solver Error ---\n";
            std::cerr << "    Status: " << cplex.getStatus() << "\n";
            std::cerr << "    Error details: " << cplex.getCplexStatus() << "\n";
        }
        Solution result = Solution(time_to_sol, cplex.getObjValue(), get_solution(cplex, x));
        env.end();
        return result;
    }

    //query indexes start at 0
    vector<vector<int>> Solver::getSubtours(const uint64_t n, const IloArray<IloNumVarArray> &x,
                                            const IloCplex &cplex) const {
        cout << "Subtour routine started" << endl;
        vector<vector<int>> subtours;
        vector<int> solution = get_solution(cplex, x);
        int** mutableX = getX_center(cplex, x);
        //prune solution
        for (int i = 0; i < solution.size()-1; ++i) {
            mutableX[solution.at(i)-1][solution.at(i+1)-1] = 0;
        }
        cout << "  prunning" << endl;
        //check sum if > 0 must be subtours
        cout << "|Offending edges| = " << getSumX(n, mutableX) << endl;
        // find the subtour(s)

        while (getSumX(n, mutableX) > 0){
            //find an edge
            int ei = 0, ej = 0;
            for (int i = 0; i < n; ++i) {
                if (ej + ei > 0)
                    break;
                for (int j = 0; j < n; ++j) {
                    if (mutableX[i][j] == 1){
                        ei = i;
                        ej = j;
                        break;
                    }

                }
            }
            //For the edge find cycle
            cout << "Found an edge between " << ei << " and " << ej << endl << "Subtour is ";
            vector<int> subtour;
            int p_size = 0;
            int current_pos = ei;
            bool init = false;
            while (true) {
                if (init){
                    p_size = subtour.size();
                    if (std::none_of(subtour.begin(), subtour.end(), compare(current_pos))){
                        subtour.push_back(current_pos);
                        cout << current_pos << " ";
                    }
                } else{
                    subtour.push_back(current_pos);
                    init = true;
                    cout << current_pos << " ";
                }
                for (int i = 0; i < n; ++i) {
                    if (mutableX[current_pos][i] > 0) {
                        current_pos = i;
                        break;
                    }
                }
                if (p_size == subtour.size())
                    break;
            } //while (p_size != subtour.size());
            cout << endl;
            //Prune from matrix
            cout << "Pruning ...";
            for (int i = 0; i < subtour.size()-1; ++i) {
                //should catch the loop ?
                mutableX[subtour.at(i)][subtour.at(i+1)] = 0;
            }
            mutableX[subtour.at(subtour.size()-1)][subtour.at(0)] = 0;
            cout << " done." << endl;
            //Add to collection
            subtours.push_back(subtour);
        }
        //cleanup after ourselves
        for (int i = 0; i < n; ++i) {
            delete mutableX[i];
        }
        delete mutableX;
        cout << "Subtour routine is done" << endl;
        return subtours;
    }

    int Solver::getSumX(const uint64_t n, int *const *mutableX) const {
        int sum = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                sum += mutableX[i][j];
            }
        }
        return sum;
    }

    void Solver::build_subtour_const_all(const IloEnv &env, const IloModel &model, const uint64_t n,
                                         const IloArray<IloNumVarArray> &x, const IloNumVarArray &u) const {
        IloExpr expr(env);
        stringstream vname;

        // Create subtour elimination constraints (8)
        IloArray<IloRangeArray> mtz(env, n);
        for (auto i = 1u; i <= n; ++i) {
            mtz[i-1] = IloRangeArray(env, n);
            for (auto j = 1u; j <= n; ++j) {
                expr = ((n - 1) * (1 - x[i][j])) - u[i - 1] + u[j - 1]; // >= 1
                vname << "mtz_" << i << "_" << j;
                mtz[i - 1][j - 1] = IloRange(env, 1, expr, IloInfinity, vname.str().c_str());
                vname.str("");
                expr.clear();
            }
            // Add constraints (8)_i to the model
            model.add(mtz[i-1]);
        }
        cout << "Added (8) to model\nConstraint building complete.\n";
    }

    void
    Solver::build_constraints(int dist_bound, int time_bound, const IloEnv &env, const IloModel &model,
                              const uint64_t n, const IloArray<IloNumVarArray> &x, const IloNumVarArray &s) const {

        IloExpr expr(env);
        std::stringstream vname;
        IloRangeArray inbound_arcs(env, n);  // Constraints (5)
        IloRangeArray outbound_arcs(env, n); // Constraints (6)

        // Create constraint (2)
        for (auto i = 1; i <= n; ++i) {
            for(auto j = 1; j <= n; ++j) {
                expr += x[i][j] * (int) tap.dist(i - 1, j - 1);
            }
        }
        IloRange interestingness(env, -IloInfinity, expr, dist_bound, "distance_epsilon");
        expr.clear();
        // Add constraints (2) to the model
        model.add(interestingness);
        cout << "Added (2) to model\n";

        // Create constraint (4)
        for (auto i = 0; i < n; ++i) {
            expr += s[i] * (int) tap.time(i);
        }
        IloRange exec_time(env, 0, expr, time_bound, "time_epsilon");
        expr.clear();
        // Add constraints (4) to the model
        model.add(exec_time);
        cout << "Added (4) to model\n";

        // Create constraints (5)
        // quelques soit j de 1 a n
        for (auto j = 1u; j <= n; ++j) {
            for (auto i = 0u; i <= n; ++i) {
                if (j != i)
                    expr += x[i][j];
            }
            expr -= s[j - 1];

            vname << "inbound_" << j;
            inbound_arcs[j-1] = IloRange(env, 0, expr, 0, vname.str().c_str()); // = 0
            vname.str("");
            expr.clear();
        }
        // Add constraints (5) to the model
        model.add(inbound_arcs);
        cout << "Added (5) to model\n";

        // Create constraints (6)
        // quelques soit i de 1 a n
        for (auto i = 1u; i <= n; ++i) {
            for (auto j = 1u; j <= n+1; ++j) {
                if (j != i)
                    expr += x[i][j];
            }
            expr += -s[i - 1];
            vname << "outbound_" << i;
            outbound_arcs[i-1] = IloRange(env, 0, expr, 0, vname.str().c_str());
            vname.str("");
            expr.clear();
        }
        // Add constraints (6) to the model
        model.add(outbound_arcs);
        cout << "Added (6) to model\n";

        // Create constraint (7s)
        for (auto i = 1u; i <= n; ++i) {
            expr += x[0][i];
        }
        IloRange path_start(env, 1, 1, "path_start");
        path_start.setExpr(expr);
        expr.clear();
        // Add constraints (7s) to the model
        model.add(path_start);
        cout << "Added (7s) to model\n";

        // Create constraint (7e)
        for (auto i = 1u; i <= n; ++i) {
            expr += x[i][n+1u];
        }
        IloRange path_end(env, 1, 1, "path_end"); // Constraint (7)
        path_end.setExpr(expr);
        expr.clear();
        // Add constraints (7e) to the model
        model.add(path_end);
        cout << "Added (7e) to model\n";

        //Forbidden links (start to end, end to start ...)
        for (auto i = 0; i <= n + 1; ++i) {
            expr += x[i][0];
            expr += x[n+1][i];
        }
        expr += x[0][n + 1];
        IloRange path_form(env, 0, expr, 0, "path_structure");
        expr.clear();
        model.add(path_form);

        expr.end();
    }

    void Solver::init_vars(const IloEnv &env, const uint64_t n, IloArray<IloNumVarArray> &x, IloNumVarArray &s,
                           IloNumVarArray &u) const {
        // Init variables x for arcs
        std::stringstream vname;
        for (auto i = 0; i <= n+1; ++i) {
            x[i] = IloNumVarArray(env, n + 2u);
            for (auto j = 0u; j <= n + 1; ++j) {
                vname << "x_" << i << "_" << j;
                if (i == j )
                    x[i][j] = IloNumVar(env, 0, 0, IloNumVar::Bool, vname.str().c_str());
                else
                    x[i][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, vname.str().c_str());
                vname.str("");
            }
        }

        // Init variables for selection
        for (auto i = 0; i < n; ++i) {
            vname << "s_" << i;
            s[i] = IloNumVar(env, 0, 1, IloNumVar::Bool, vname.str().c_str());
            vname.str("");
        }

        // Init variables for MTZ subtour elimination and enforce part of (8)
        for (auto i = 1u; i <= n; ++i) {
            vname << "u_" << i;
            u[i-1] = IloNumVar(env, 2, n, IloNumVar::Int, vname.str().c_str());
            vname.str("");
        }

        cout << "Variables Init. complete\n";
    }



}
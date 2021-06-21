#include "SolverVPLS.h"
#include "random"

namespace cplex_tap {

    ILOMIPINFOCALLBACK1(MyCallbackVPLS, IloInt, num) {
        IloEnv env = getEnv();
        double this_z = getIncumbentObjValue();
        double best_z = getBestObjValue();
        //TODO only print if we have better solution ....
        if (1)
        env.out() << "  [INFO MIP Callback]" << " CLK " << clock() << " Z " << this_z << std::endl;
    }

    double
    SolverVPLS::solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production, bool seed, string warmStart) const{
        std::cout << "CLK_RATE " << CLOCKS_PER_SEC << std::endl;
        std::cout << "Starting Model generation ...\n";

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
        cplex.setParam(IloCplex::Param::TimeLimit, max_init_time);
        IloCplex::Callback mycallback = cplex.use(MyCallbackVPLS(env, 10));
        //cplex.setParam(IloCplex::IntSolLim, 4);// remove
        if (!production)
            cplex.setParam(IloCplex::Param::Threads, 1);
        else
            cplex.setParam(IloCplex::Param::Threads, 8);

        if (seed) {
            warm_start(warmStart, env, n, x, s, cplex);
        }

        bool solved = false;
        time_t start, end;
        start = clock();
        try {
            solved = cplex.solve();
        }
        catch (const IloException &e) {
            std::cerr << "\n\n--- CPLEX Exception ---\n";
            std::cerr << e << "\n";
            env.end();
            throw;
        }
        end = clock();
        double time_to_sol = (double) (end - start) / (double) CLOCKS_PER_SEC;

        if (solved) {
            // If CPLEX successfully solved the model, print the results
            std::cout << "\n--- Solver success ---\n";
            std::cout << "    Status: " << cplex.getStatus() << "\n";
            std::cout << "    Objective: " << cplex.getObjValue() << "\n";
            if (debug)
                dump(cplex, x, env, s, u);
            print_solution(cplex, x);

            // use rd() instead of seed for non determinism
            //std::random_device rd;
            std::mt19937 mt(42);
            std::vector<int> current_fixed;
            cplex.setParam(IloCplex::Param::TimeLimit, max_epoch_time);
            //cplex.setParam(IloCplex::Param::Emphasis::MIP, 4);
            //cplex.setParam(IloCplex::IntSolLim, 9223372036800000000);

            vector<double> zvalues;

            time_t clk = clock();
            std::cout << "  CLK_START " << start << std::endl;

            cout << "Starting MPLS heurisitc max iterations " << max_iter << endl;
            for (auto iter = 0; iter < max_iter; ++iter){
                cout << "  Starting iteration " << iter << endl;
                time_t clk = clock();
                std::cout << "  CLK_START_ITER " << start << std::endl;
                vector<int> solution = get_solution(cplex, x);
                //if (solution.)
                std::uniform_int_distribution<int> dist(1, solution.size() - h);
                int wstart = dist(mt);
                int wend = wstart + h;
                cout << "  window=[" << wstart << "," << wend << "]" << endl;

                IloNumArray vals_s(env);
                try {
                    cout << "  Status:" << cplex.getStatus() << endl;
                    cplex.getValues(vals_s, s);
                }
                catch (const IloException &e) {
                    std::cerr << "\n\n--- CPLEX Exception (VPLS) ---\n";
                    std::cerr << e << "\n";
                    env.end();
                    throw;
                }

                if(iter > 0){
                    for (const auto& f: current_fixed) {
                        s[f].setBounds(0,1);
                    }
                    current_fixed.clear();
                    cout << "  clear ok" << endl;
                }

                for (auto j = 0u; j < solution.size(); ++j) {
                    if (!(j >= wstart && j <= wend)) {
                        int value = vals_s[solution.at(j)-1] > 0.5;
                        current_fixed.push_back(solution.at(j)-1);
                        s[solution.at(j)-1].setBounds(value, value);

                    }
                }
                cout << "  set ok" << endl;
                //Report
                start = clock();
                solved = cplex.solve();
                end = clock();
                time_to_sol += (double) (end - start) / (double) CLOCKS_PER_SEC;
                print_solution(cplex, x);
                zvalues.push_back(cplex.getObjValue());
                cout << "  Z=" << cplex.getObjValue() << endl;

                //Check convergence criterion, no significant change in Z in the last 5 iterations
                if (iter >= 4){
                    bool check = true;
                    for (int i = iter; i > iter - 4 ; --i) {
                        check = check && (zvalues.at(i) - zvalues.at(i-1) < 0.001);

                    }
                    if (check) {
                        cout << "  VPLS Converged at iteration " << iter << endl;
                        break;
                    }
                }
            }



        } else {
            std::cerr << "\n--- Solver Error ---\n";
            std::cerr << "    Status: " << cplex.getStatus() << "\n";
            std::cerr << "    Error details: " << cplex.getCplexStatus() << "\n";
        }

        env.end();
        return time_to_sol;
    }
}
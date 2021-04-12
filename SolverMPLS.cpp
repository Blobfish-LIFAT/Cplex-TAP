#include "SolverMPLS.h"
#include "random"

namespace cplex_tap {

    double
    SolverMPLS::solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production) const{
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
        cplex.setParam(IloCplex::Param::TimeLimit, 60);
        cplex.setParam(IloCplex::IntSolLim, 5);
        if (!production)
            cplex.setParam(IloCplex::Param::Threads, 1);
        else
            cplex.setParam(IloCplex::Param::Threads, 8);


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

            //TODO move parameters
            int h = 10;
            int max_iter = 20;
            // use rd() instead of seed for non determinism
            //std::random_device rd;
            std::mt19937 mt(42);
            std::uniform_int_distribution<int> dist(1, n-1 - h);
            std::vector<IloRange> current_fixed;

            cplex.setParam(IloCplex::Param::TimeLimit, 600);
            cplex.setParam(IloCplex::IntSolLim, 9223372036800000000);

            cout << "Starting MPLS heurisitc max iterations " << max_iter << endl;
            for (auto iter = 0; iter < max_iter; ++iter){
                cout << "  Starting iteration " << iter << endl;
                int wstart = dist(mt);
                int wend = wstart + h;
                cout << "  window=[" << wstart << "," << wend << "]" << endl;
                if(iter > 0){
                    for (const auto& f: current_fixed) {
                        cplex.getModel().remove(f);
                    }
                    current_fixed.clear();
                    cout << "  clear ok" <<endl;
                }

                IloNumArray vals_s(env);
                cplex.getValues(vals_s, s);
                for (auto j = 0u; j < n; ++j) {
                    if (!(j >= wstart && j <= wend)) {
                        int value = vals_s[j] > 0.5;
                        IloRange fixed(env, value, s[j], value);
                        model.add(fixed);
                        current_fixed.push_back(fixed);
                    }
                }

                cplex.solve();
                print_solution(cplex, x);
                cout << "  Z=" << cplex.getObjValue() << endl;
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
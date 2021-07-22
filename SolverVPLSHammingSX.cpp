#include "SolverVPLSHammingSX.h"

namespace cplex_tap {

    ILOMIPINFOCALLBACK1(VPLS_HAMMINGSX_INFO_CBK, IloInt, num) {
        IloEnv env = getEnv();
        double this_z = getIncumbentObjValue();
        //TODO only print if we have better solution ....
        env.out() << "  [INFO MIP Callback]" << " CLK " << clock() << " Z " << this_z << std::endl;
    }

    Solution
    SolverVPLSHammingSX::solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production,
                                       bool seed, string warmStart) const {
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
        if (debug) IloCplex::Callback mycallback = cplex.use(VPLS_HAMMINGSX_INFO_CBK(env, 10));
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

        // If CPLEX successfully solved the model, print the results
        if (solved) {
            std::cout << "\n--- Solver success ---\n" << "    Status: " << cplex.getStatus() << "\n"
                      << "    Objective: " << cplex.getObjValue() << "\n";
            print_solution(cplex, x);

            cplex.setParam(IloCplex::Param::TimeLimit, max_epoch_time);
            // constraint
            IloRange hamming_cstr;
            vector<double> zvalues;

            time_t clk = clock();
            std::cout << "  CLK_START " << start << "\nStarting MPLS heurisitc max iterations " << max_iter << endl;
            for (auto iter = 0; iter < max_iter; ++iter) {
                // Begin the iteration
                time_t clk = clock();
                cout << "  Starting iteration " << iter << "\n  CLK_START_ITER " << start << std::endl;

                // Fetch current solution
                // Has to do it before touching previous constraint to avoid 1217
                int** xvalues = getX_center(cplex, x);
                std::vector<int> solution = get_solution(cplex, x);
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

                if (iter > 0) {
                    model.remove(hamming_cstr);
                    cout << "  clear ok" << endl;
                }

                IloExpr hamming_expr(env);
                // S variables
                for (auto j = 0u; j < n; ++j) {
                    int value = vals_s[j] > 0.5;
                    if (value == 1){
                        hamming_expr += 1 - s[j];
                    } else{
                        hamming_expr += 1 - (1 - s[j]);
                    }
                }
                // X variables
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        if (i != j){
                            // indexes on x off by one from values
                            if (xvalues[i][j] == 1){
                                hamming_expr += 1 - x[i+1][j+1];
                            } else{
                                hamming_expr += 1 - (1 - x[i+1][j+1]);
                            }
                        }
                    }
                }
                cout << "  building bound" << endl;
                hamming_cstr = IloRange(env, 0, hamming_expr, h, "Hamming_bound");
                cout << "  adding bound" << endl;
                model.add(hamming_cstr);
                cout << "  set ok" << endl;

                //Report
                start = clock();
                solved = cplex.solve();
                end = clock();
                time_to_sol += (double) (end - start) / (double) CLOCKS_PER_SEC;
                print_solution(cplex, x);
                cout << "  Z=" << cplex.getObjValue() << endl;
                zvalues.push_back(cplex.getObjValue());

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
        Solution result = Solution(time_to_sol, cplex.getObjValue(), get_solution(cplex, x));
        env.end();
        return result;
    }
}
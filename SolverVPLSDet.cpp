//
// Created by chanson on 6/22/2021.
//

#include "SolverVPLSDet.h"

namespace cplex_tap {

    ILOMIPINFOCALLBACK1(VPLS_DET_INFO_CBK, IloInt, num) {
    IloEnv env = getEnv();
    double this_z = getIncumbentObjValue();
    //TODO only print if we have better solution ....
    env.out() << "  [INFO MIP Callback]" << " CLK " << clock() << " Z " << this_z << std::endl;
    }

double
SolverVPLSDet::solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production, bool seed, string warmStart) const{
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
    IloCplex::Callback mycallback = cplex.use(VPLS_DET_INFO_CBK(env, 10));
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
        std::cout << "\n--- Solver success ---\n" << "    Status: " << cplex.getStatus() << "\n" << "    Objective: " << cplex.getObjValue() << "\n";
        print_solution(cplex, x);

        std::vector<int> current_fixed;
        cplex.setParam(IloCplex::Param::TimeLimit, max_epoch_time);
        vector<double> zvalues;
        // variables for positioning of the window
        bool prev_success = true;
        bool prev_final_subseq = false;
        int prev_start_pos = -1;

        time_t clk = clock();
        std::cout << "  CLK_START " << start << "\nStarting MPLS heurisitc max iterations " << max_iter << endl;
        for (auto iter = 0; iter < max_iter; ++iter){
            // if we are at the end of the sequence and we successfully solved the sub-problem
            if (prev_final_subseq && prev_success){
                cout << "  VPLS Converged at iteration " << iter << endl;
                break; // Stop
            }

            // Begin the iteration
            time_t clk = clock();
            cout << "  Starting iteration " << iter << "\n  CLK_START_ITER " << start << std::endl;

            // Fetch current solution
            vector<int> solution = get_solution(cplex, x);

            // Draw window start
            int window_start = 0;
            if (!prev_success) {
                window_start = prev_start_pos + (h / 2);
            }
            prev_start_pos = window_start;
            int wend = std::min(window_start + h, (int) solution.size() - 1);
            if (wend ==  solution.size() - 1)
                prev_final_subseq = true;
            else
                prev_final_subseq = false;
            cout << "  window=[" << window_start << "," << wend << "]" << endl;

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
                if (!(j >= window_start && j <= wend)) {
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

            if (cplex.getCplexStatus())

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
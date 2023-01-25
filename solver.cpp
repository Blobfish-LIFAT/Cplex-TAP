#include "solver.h"
#include "utils.h"
#include <cmath>
#include <time.h>
#include <limits>
#include <algorithm>

namespace cplex_tap {

    ILOMIPINFOCALLBACK1(MyCallback, IloInt, num) {
        IloEnv env = getEnv();
        double this_z = getIncumbentObjValue();
        double best_z = getBestObjValue();
        //TODO only print if we have better solution ....
        if (1)
            env.out() << "[INFO MIP Callback]" << " CLK " << clock() << " Z " << this_z << std::endl;
    }

    Solution Solver::solve(int dist_bound, int time_bound, bool seed, string warmStart) const {
        //std::cout << "CLK_RATE " << CLOCKS_PER_SEC << std::endl;
        if (debug) std::cout << "Starting Model generation ..." << std::endl;

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
        if (debug) cout << "Variables Init. complete\n";

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
        if (debug) std::cout << "Added Objective to model\n";

        // Free some memory 
        expr.end();

        //Init solver
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::TimeLimit, timeout);
        cplex.setParam(IloCplex::Param::Threads, 1);
        //cplex.setParam(IloCplex::Param::MIP::Display, 0);
        //cplex.setParam(IloCplex::Param::Simplex::Display, 0);
        //cplex.setOut(env.getNullStream());

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
        if (debug) std::cout << "CLK_START " << start << std::endl;
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
            if (!shutUp) {
                std::cout << "\n--- Solver success ---\n";
                std::cout << "    Status: " << cplex.getStatus() << "\n";
                std::cout << "    Objective: " << cplex.getObjValue() << "\n";
                print_solution(cplex, x);
            }
            if (debug) {
                dump(cplex, x, env, s, u);
            }


        }
        else {
            std::cerr << "\n--- Solver Error ---\n";
            std::cerr << "    Status: " << cplex.getStatus() << "\n";
            std::cerr << "    Error details: " << cplex.getCplexStatus() << "\n";
        }

        Solution result = Solution(cplex.getStatus() == IloAlgorithm::Optimal, time_to_sol, cplex.getObjValue(), get_solution(cplex, x), cplex.getNnodes());
        cplex.end();
        env.end();
        return result;
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
        if (debug) cout << "Added (8) to model\nConstraint building complete.\n";
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
        if (debug) cout << "Added (2) to model\n";

        // Create constraint (4)
        for (auto i = 0; i < n; ++i) {
            expr += s[i] * (int) tap.time(i);
        }
        IloRange exec_time(env, 0, expr, time_bound, "time_epsilon");
        expr.clear();
        // Add constraints (4) to the model
        model.add(exec_time);
        if (debug) cout << "Added (4) to model\n";

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
        if (debug) cout << "Added (5) to model\n";

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
        if (debug)  cout << "Added (6) to model\n";

        // Create constraint (7s)
        for (auto i = 1u; i <= n; ++i) {
            expr += x[0][i];
        }
        IloRange path_start(env, 1, 1, "path_start");
        path_start.setExpr(expr);
        expr.clear();
        // Add constraints (7s) to the model
        model.add(path_start);
        if (debug) cout << "Added (7s) to model\n";

        // Create constraint (7e)
        for (auto i = 1u; i <= n; ++i) {
            expr += x[i][n+1u];
        }
        IloRange path_end(env, 1, 1, "path_end"); // Constraint (7)
        path_end.setExpr(expr);
        expr.clear();
        // Add constraints (7e) to the model
        model.add(path_end);
        if (debug) cout << "Added (7e) to model\n";

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
                           IloNumVarArray &u) {
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
    }

    int Solver::getTimeout() const {
        return timeout;
    }

    void Solver::setTimeout(int timeout) {
        Solver::timeout = timeout;
    }


}
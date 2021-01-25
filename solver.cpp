#include "solver.h"
#include <cmath>
#include <limits>

namespace cplex_tap {
    void Solver::solve_and_print() const {
        std::cout << "Starting Solver\n";
        //TODO parameterize this
        int int_bound = 50;
        int time_bound = 25;

        // Init CPLEX environment and model objects
        IloEnv env;
        IloModel model(env);

        const auto n = tap.size();

        // Variables
        IloArray<IloNumVarArray> x(env, n + 1);
        IloNumVarArray s(env, n);
        IloNumVarArray u(env, n);

        // String stream trick from stackoverflow
        std::stringstream vname;

        // Init variables x for arcs
        for (auto i = 0; i < n+1; ++i) {
            x[i] = IloNumVarArray(env, n + 1);
            for (auto j = 0; j < n+1; ++j) {
                vname << "x_" << i << "_" << j;
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
        for (auto i = 1; i < n; ++i) {
            vname << "u_" << i;
            u[i] = IloNumVar(env, 2, n, IloNumVar::Int, vname.str().c_str());
            vname.str("");
        }

        std::cout << "Variables Init. complete\n";

        // Constraints
        IloRange interestingness(env, int_bound, IloInfinity, "interestingness_epsilon");
        IloRange exec_time(env, 0, time_bound, "time_epsilon");
        IloRangeArray inbound_arcs(env, n);  // Constraints (5)
        IloRangeArray outbound_arcs(env, n); // Constraints (6)
        IloRange path_start(env, 1, 1, "path_start"); // Constraint (7)
        IloRange path_end(env, 1, 1, "path_end"); // Constraint (7)
        IloArray<IloRangeArray> mtz(env, n); // Constraints (8)

        // Expresion object for constraints/objective
        IloExpr expr(env);

        // Create constraint (1)
        for (auto i = 0; i < n; ++i) {
           expr += s[i] * tap.interest(i);  
        }
        interestingness.setExpr(expr);
        expr.clear();
        // Add constraints (1) to the model
        model.add(interestingness);
        std::cout << "Added (1) to model\n";

        // Create constraint (4)
        for (auto i = 1; i < n; ++i) {
            expr += s[i] * (int) tap.time(i);
        }
        exec_time.setExpr(expr);
        expr.clear();
        // Add constraints (4) to the model
        model.add(exec_time);
        std::cout << "Added (4) to model\n";

        // Create constraints (5)
        for (auto i = 0u; i < n; ++i) {
            for (auto j = 0u; j < n+1; ++j) {
                expr += x[i][j] - s[i];
            }

            vname << "inbound_" << i;
            inbound_arcs[i] = IloRange(env, 0, expr, 0, vname.str().c_str()); // = 0
            vname.str("");
            expr.clear();
        }
        // Add constraints (5) to the model
        model.add(inbound_arcs);
        std::cout << "Added (5) to model\n";

        // Create constraints (6)
        for (auto i = 0u; i < n; ++i) {
            for (auto j = 0u; j < n+1; ++j) {
                expr += x[j][i] - s[i];
            }

            vname << "outbound_" << i;
            outbound_arcs[i] = IloRange(env, 0, expr, 0, vname.str().c_str());
            vname.str("");
            expr.clear();
        }
        // Add constraints (6) to the model
        model.add(outbound_arcs);
        std::cout << "Added (6) to model\n";

        // Create constraint (7s)
        for (auto i = 1; i < n; ++i) {
            expr += x[0][i];          
        }
        path_start.setExpr(expr);
        expr.clear();
        // Add constraints (7s) to the model
        model.add(path_start);

        // Create constraint (7e)
        for (auto i = 1; i < n; ++i) {
            expr += x[i][n+1];
        }
        path_end.setExpr(expr);
        expr.clear();
        // Add constraints (7e) to the model
        model.add(path_end);
        std::cout << "Added (7) to model\n";

        // Create subtour elimination constraints (8)
        mtz[0] = IloRangeArray(env);
        for (auto i = 1u; i < n; ++i) {
            mtz[i] = IloRangeArray(env, n);
            for (auto j = 1u; j < n; ++j) {
                expr = (n-1) * (1 - x[i][j]) - u[i] + u[j]; // >= 1

                vname << "mtz_" << i << "_" << j;
                mtz[i][j] = IloRange(env, 1, expr, IloInfinity, vname.str().c_str());
                vname.str("");
                expr.clear();
            }
            // Add constraints (8)_i to the model
            model.add(mtz[i]);
        }
        std::cout << "Added (8) to model\nConstraint building complete.";

        //
        // --- Objective (2) ---
        //
        for (auto i = 0u; i < n; ++i) {
            for (auto j = 0u; j < n; ++j) {
                expr += (int) tap.dist(i, j) * x[i][j];
            }
        }
        IloObjective obj(env, expr, IloObjective::Minimize);
        model.add(obj);

        // Free some memory 
        expr.end();

        //Init solver
        IloCplex cplex(model);
        // Export model to file
        cplex.exportModel("tap_debug_model.lp");

        bool solved = false;
        try {
            solved = cplex.solve();
        }
        catch (const IloException& e) {
            std::cerr << "\n\n--- CPLEX Exception ---\n";
            std::cerr << e << "\n";
            env.end();
            throw;
        }

        if (solved) {
            // If CPLEX successfully solved the model, print the results
            std::cout << "\n\n--- Solver success ---\n";
            std::cout << "\tStatus: " << cplex.getStatus() << "\n";
            std::cout << "\tObjective: " << cplex.getObjValue() << "\n";
        }
        else {
            std::cerr << "\n\n--- Solver Error ---\n";
            std::cerr << "\tStatus: " << cplex.getStatus() << "\n";
            std::cerr << "\tError details: " << cplex.getCplexStatus() << "\n";
        }

        env.end();
    }


}
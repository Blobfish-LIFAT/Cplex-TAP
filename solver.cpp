#include "solver.h"
#include <cmath>
#include <time.h>
#include <limits>

namespace cplex_tap {
    double Solver::solve_and_print(int dist_bound, int time_bound) const {
        bool const debug = false;
        std::cout << "Starting Solver\n";

        // Init CPLEX environment and model objects
        IloEnv env;
        IloModel model(env);

        const uint64_t n = tap.size();

        // Variables
        IloArray<IloNumVarArray> x(env, n + 2u);
        IloNumVarArray s(env, n);
        IloNumVarArray u(env, n);

        // String stream trick from stackoverflow
        std::stringstream vname;

        // Init variables x for arcs
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

        std::cout << "Variables Init. complete\n";

        // Constraints
        IloRangeArray inbound_arcs(env, n);  // Constraints (5)
        IloRangeArray outbound_arcs(env, n); // Constraints (6)
        
        

        // Expresion object for constraints/objective
        IloExpr expr(env);

        // Create constraint (2)
        for (auto i = 1; i <= n; ++i) {
            for(auto j = 1; j <= n; ++j)
                expr += x[i][j] * (int) tap.dist(i-1, j-1);  
        }
        IloRange interestingness(env, -IloInfinity, expr, dist_bound, "distance_epsilon");
        expr.clear();
        // Add constraints (2) to the model
        model.add(interestingness);
        std::cout << "Added (2) to model\n";

        // Create constraint (4)
        for (auto i = 0; i < n; ++i) {
            expr += s[i] * (int) tap.time(i);
        }
        IloRange exec_time(env, 0, expr, time_bound, "time_epsilon");
        expr.clear();
        // Add constraints (4) to the model
        model.add(exec_time);
        std::cout << "Added (4) to model\n";

        // Create constraints (5)
        // queqlues soit j de 1 a n
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
        std::cout << "Added (5) to model\n";

        // Create constraints (6)
        // queqlues soit i de 1 a n
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
        std::cout << "Added (6) to model\n";

        // Create constraint (7s)
        for (auto i = 1u; i <= n; ++i) {
            expr += x[0][i];          
        }
        IloRange path_start(env, 1, 1, "path_start");
        path_start.setExpr(expr);
        expr.clear();
        // Add constraints (7s) to the model
        model.add(path_start);
        std::cout << "Added (7s) to model\n";

        // Create constraint (7e)
        for (auto i = 1u; i <= n; ++i) {
            expr += x[i][n+1u];
        }
        IloRange path_end(env, 1, 1, "path_end"); // Constraint (7)
        path_end.setExpr(expr);
        expr.clear();
        // Add constraints (7e) to the model
        model.add(path_end);
        std::cout << "Added (7e) to model\n";

        //Forbidden links (start to end, end to start ...)
        for (auto i = 0; i <= n + 1; ++i) {
            expr += x[i][0];
            expr += x[n+1][i];
        }
        expr += x[0][n + 1];
        IloRange path_form(env, 0, expr, 0, "path_structure");
        expr.clear();
        model.add(path_form);
        
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
        std::cout << "Added (8) to model\nConstraint building complete.\n";

        //
        // --- Objective (1) ---
        //
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
        cplex.setParam(IloCplex::TiLim, 3600);
        // Export model to file
        cplex.exportModel("tap_debug_model.lp");

        bool solved = false;
        time_t start, end;
        start = clock();
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
        double time_to_sol = (double)(end - start) / (double)CLK_TCK;

        if (solved) {
            // If CPLEX successfully solved the model, print the results
            std::cout << "\n--- Solver success ---\n";
            std::cout << "    Status: " << cplex.getStatus() << "\n";
            std::cout << "    Objective: " << cplex.getObjValue() << "\n";

            if (debug) {
                IloNumArray vals_s(env);
                cplex.getValues(vals_s, s);
                env.out() << "s = " << vals_s << endl;

                IloNumArray vals_u(env);
                cplex.getValues(vals_u, u);
                env.out() << "u = " << vals_u << endl;

                print_X(cplex, x);
            }
            //print_solution(cplex, x);

        }
        else {
            std::cerr << "\n--- Solver Error ---\n";
            std::cerr << "    Status: " << cplex.getStatus() << "\n";
            std::cerr << "    Error details: " << cplex.getCplexStatus() << "\n";
        }

        env.end();
        return time_to_sol;
    }

    void Solver::print_X(const IloCplex& cplex, const IloArray<IloNumVarArray>& x) const {
        const uint64_t n = tap.size();

        for (IloInt i = 0; i <= n + 1; ++i) {
            for (IloInt j = 0; j <= n + 1; ++j) {
                if (j == i) {
                    cout << "X" << ", ";
                }
                else {
                    cout << cplex.getValue(x[i][j]) << ", ";
                }
            }
            cout << endl;
        }
    }

    void Solver::print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x)  const {
        const uint64_t n = tap.size();
        const auto start_node = 0u;
        bool first = true;
        auto current_node = start_node;
        

        while (current_node != n + 1) {
            if (first)
                first = false;
            else
                std::cout << current_node << " ";
            for (auto i = 0u; i <= n+1; ++i) {
                if (cplex.getValue(x[current_node][i]) > 0) {
                    current_node = i;
                    break;
                }
            }
        }
    
    }
}
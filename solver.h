#pragma once

#include "instance.h"
#include <iostream>

#ifndef IL_STD
#define IL_STD
#endif

#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

namespace cplex_tap {
    class Solver {
        // The TAP instance
        const Instance& tap;


        vector<int> get_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x)  const {
            const uint64_t n = tap.size();
            const auto start_node = 0u;
            bool first = true;
            auto current_node = start_node;
            vector<int> path;

            while (current_node != n + 1) {
                if (first)
                    first = false;
                else
                    path.push_back(current_node);
                for (auto i = 0u; i <= n+1; ++i) {
                    if (i!= current_node && cplex.getValue(x[current_node][i]) > 0) {
                        current_node = i;
                        break;
                    }
                }
            }
            return path;
        }

        void print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x)  const {
            cout << "SOLUTION: ";
            vector<int> path = get_solution(cplex, x);
            for (int i : path)
                std::cout << i << ' ';
            std::cout << std::endl;
        }

        // dump decision variables
        void dump(const IloCplex& cplex, const IloArray<IloNumVarArray>& x, const IloEnv& env, const IloNumVarArray &s, const IloNumVarArray &u) const {
            cout << "--- DUMP STARTED --" << endl;
            const uint64_t n = tap.size();
            IloNumArray vals_s(env);
            cplex.getValues(vals_s, s);
            env.out() << "s = " << vals_s << endl;
            IloNumArray vals_u(env);
            try {
                cplex.getValues(vals_u, u);
                env.out() << "u = " << vals_u << endl;
            } catch (IloException& e) {
                cerr << "Concert Exception: " << e << endl;
            }
            for (IloInt i = 0; i <= n + 1; ++i) {
                for (IloInt j = 0; j <= n + 1; ++j) {
                    if (j == i) {
                        cout << "X" << " ";
                    }
                    else {
                        try{
                            cout << abs(cplex.getValue(x[i][j])) << " ";
                        }catch (IloException& e) {
                            cout << "N" << " ";
                        }

                    }
                }
                cout << endl;
            }
            cout << "--- DUMP COMPLETE --" << endl;
        }

        int** getX_center(const IloCplex& cplex, const IloArray<IloNumVarArray>& x) const {
            const uint64_t n = tap.size();
            int** out = new int*[n];
            for (IloInt i = 1; i <= n; ++i) {
                out[i-1] = new int[n];
                for (IloInt j = 1; j <= n; ++j) {
                    if (i != j){
                        int temp = 0;
                        try {
                            temp = ((int) cplex.getValue(x[i][j]) > 0.000005);
                        } catch (IloException& e) {
                            cerr << "Concert Exception: " << e << endl;
                        }
                        out[i-1][j-1] = temp;
                    }
                    else
                        out[i-1][j-1] = 0;
                    //cout << out[i-1][j-1] << " ";
                }
                //cout << endl;
            }
            return out;
        }

    public:

        // Builds a solver the specified instance
        explicit Solver(const Instance& tap) : tap{ tap } {}

        // Run solver and dump result to stdout
        double solve_and_print(int dist_bound, int time_bound, bool debug) const;

        void init_vars(const IloEnv &env, const uint64_t n, IloArray<IloNumVarArray> &x, IloNumVarArray &s,
                       IloNumVarArray &u) const;

        void
        build_constraints(int dist_bound, int time_bound, const IloEnv &env, const IloModel &model, const uint64_t n,
                          const IloArray<IloNumVarArray> &x, const IloNumVarArray &s) const;

        void
        build_subtour_const_all(const IloEnv &env, const IloModel &model, const uint64_t n,
                                const IloArray<IloNumVarArray> &x,
                                const IloNumVarArray &u) const;

        int getSumX(const uint64_t n, int *const *mutableX) const;

        //find subtours in a relaxed solution
        vector<vector<int>> getSubtours(const uint64_t n, const IloArray<IloNumVarArray> &x, const IloCplex &cplex) const;
    };
}


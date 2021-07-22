#pragma once

#include "instance.h"
#include <iostream>
#include <algorithm>

#ifndef IL_STD
#define IL_STD
#endif

#include <cstring>
#include <ilcplex/ilocplex.h>
#include "Solution.h"
ILOSTLBEGIN

namespace cplex_tap {
    class Solver {
    protected:
        // The TAP instance
        const Instance& tap;

        //virtual IloCplex::CallbackI* IloCplex::CallbackI::duplicateCallback() const = 0;

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
                            temp = ((int) cplex.getValue(x[i][j]) > 0.005);
                        } catch (IloException& e) {
                            cerr << "Concert Exception: " << e << endl;
                            exit(-1);
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

    public:

        // Builds a solver the specified instance
        explicit Solver(const Instance& tap) : tap{ tap } {}

        // Run solver and dump result to stdout
        virtual Solution solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production, bool seed, string warmStart) const;

        void warm_start(std::string warm_file, IloEnv &env, uint64_t n,  IloArray <IloNumVarArray> &x, IloNumVarArray &s, IloCplex &cplex) const{
            cout << "loading warm start from " << warm_file << endl;
            //Read file
            vector<int> ssol;
            std::ifstream ifs(warm_file);
            std::string line;
            std::getline(ifs, line);
            cout << "  Warm solution :" << line << endl;
            ifs.close();
            // Load starting solution
            string token;
            size_t pos = 0;
            while ((pos = line.find(" ")) != string::npos) {
                token = line.substr(0, pos);
                ssol.push_back(std::stoi(token, nullptr));
                line.erase(0, pos + 1);
            }
            // Build MIP Start
            IloNumVarArray startVar(env);
            IloNumArray startVal(env);
            // Build S vector
            for (int i = 0; i < n; ++i) {
                startVar.add(s[i]);
                startVal.add(std::find(ssol.begin(), ssol.end(), i) != ssol.cend());
            }
            // Build X matrix
            // First line - start node
            int sn = ssol.at(0) + 1;
            for (int j = 1; j <= n; ++j) {
                startVar.add(x[0][j]);
                startVal.add(sn == j);
            }
            // Last column - finish node
            int en = ssol.at(ssol.size()-1) + 1;
            for (int j = 1; j <= n; ++j) {
                startVar.add(x[j][n+1]);
                startVal.add(en == j);
            }
            for (int i = 1; i <= n; ++i) {
                for (int j = 1; j <= n; ++j) {
                    if (i != j){
                        startVar.add(x[i][j]);
                        std::vector<int>::iterator jit = std::find(ssol.begin(), ssol.end(), j-1);
                        // first node or not in solution
                        if (jit == ssol.begin() || jit == ssol.end()){
                            startVal.add(0);
                        }
                        else{
                            // if we are on the right line
                            if(*std::prev(jit) == (i-1)){
                                startVal.add(1);
                            } else {
                                startVal.add(0);
                            }
                        }
                    }
                }
            }

            cplex.addMIPStart(startVar, startVal);
        }

        // Ids from 1
        vector<int> get_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x)  const {
            const uint64_t n = tap.size();
            const auto start_node = 0u;
            bool first = true;
            auto current_node = start_node;
            vector<int> path;

            uint64_t iter = 0;
            uint64_t max_iter = n+2;

            while (current_node != n + 1) {
                //std::cout << current_node << "|" << n+1 << std::endl;
                if (first)
                    first = false;
                else
                    path.push_back(current_node);
                for (auto i = 0u; i <= n+1; ++i) {
                    if (i!= current_node && cplex.getValue(x[current_node][i]) > 0.5) {
                        current_node = i;
                        break;
                    }
                }
                iter++;
                //no valid solution: eg cycle
                if (iter >= max_iter) {
                    return vector<int>();
                }
            }
            return path;
        }

        void warm_start(std::vector<int> warm_sol, IloEnv &env, uint64_t n,  IloArray <IloNumVarArray> &x, IloNumVarArray &s, IloCplex &cplex) const{
            // Build MIP Start
            IloNumVarArray startVar(env);
            IloNumArray startVal(env);
            // Build S vector
            for (int i = 0; i < n; ++i) {
                startVar.add(s[i]);
                startVal.add(std::find(warm_sol.begin(), warm_sol.end(), i) != warm_sol.cend());
            }
            cout << "build s" << endl;
            // Build X matrix
            // First line - start node
            int sn = warm_sol.at(0) + 1;
            for (int j = 1; j <= n; ++j) {
                startVar.add(x[0][j]);
                startVal.add(sn == j);
            }
            cout << "build start" << endl;
            // Last column - finish node
            int en = warm_sol.at(warm_sol.size()-1) + 1;
            for (int j = 1; j <= n; ++j) {
                startVar.add(x[j][n+1]);
                startVal.add(en == j);
            }
            cout << "build end" << endl;
            for (int i = 1; i <= n; ++i) {
                for (int j = 1; j <= n; ++j) {
                    if (i != j){
                        startVar.add(x[i][j]);
                        std::vector<int>::iterator jit = std::find(warm_sol.begin(), warm_sol.end(), j-1);
                        // first node or not in solution
                        if (jit == warm_sol.begin() || jit == warm_sol.end()){
                            startVal.add(0);
                        }
                        else{
                            // if we are on the right line
                            if(*std::prev(jit) == (i-1)){
                                startVal.add(1);
                            } else {
                                startVal.add(0);
                            }
                        }
                    }
                }
            }

            cplex.addMIPStart(startVar, startVal);
        }

        // Ids from 1
        void print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x)  const {
            cout << "SOLUTION: ";
            vector<int> path = get_solution(cplex, x);
            for (int i : path)
                std::cout << i << ' ';
            std::cout << std::endl;
        }

    };
}


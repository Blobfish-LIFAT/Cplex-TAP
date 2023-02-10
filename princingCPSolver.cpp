//
// Created by alex on 08/02/23.
//

#include "princingCPSolver.h"
#include "Query.h"
#include <random>
#include "JVMAdapter.h"
#include <numeric>
#include "solver.h"
#include <ilcp/cp.h>

namespace cplex_tap {
    Solution princingCPSolver::solve() const {
        std::cout << "[INFO] CLK_RATE " << CLOCKS_PER_SEC << std::endl;
        int global_t = 0;

        vector<Query> rmpQSet;
        rmpQSet.insert(rmpQSet.end(), extStarting.begin(), extStarting.end());

        const auto timeWeights = pricingIST.getDimWeights();


        // Build a rmp problem instances
        Instance rmpIST = buildRMPInstance(rmpQSet);

        // Init CPLEX environment and model objects
        IloEnv cplex;
        IloModel pricing(cplex);

        /*
         *  Data about queries in the RMP
         */
        vector<vector<bool>> S;
        S.reserve(rmpQSet.size());
        for (int i = 0; i < rmpQSet.size(); ++i) {
            vector<bool> tmp;
            for (int j = 0; j < pricingIST.getNbDims(); ++j) {
                bool appears = false;
                const string dimName = pricingIST.getDimName(j);
                for (int k = 0; k < rmpQSet[i].getLeftPredicate().size(); ++k) {
                    appears |= rmpQSet[i].getLeftPredicate()[k].first == dimName;
                }
                for (int k = 0; k < rmpQSet[i].getRightPredicate().size(); ++k) {
                    appears |= rmpQSet[i].getRightPredicate()[k].first == dimName;
                }
                tmp.emplace_back(appears);
            }
            S.emplace_back(tmp);
        }

        /*
         * Variables Declaration
         */
        IloNumVarArray cpLeftMeasure(cplex, pricingIST.getNbMeasures());
        IloNumVarArray cpRightMeasure(cplex, pricingIST.getNbMeasures());
        IloNumVarArray cpGroupBy(cplex, pricingIST.getNbDims());
        IloNumVarArray cpSelection(cplex, pricingIST.getNbDims());
        IloArray<IloNumVarArray> cpLeftSel(cplex, pricingIST.getNbDims());
        IloArray<IloNumVarArray> cpRightSel(cplex, pricingIST.getNbDims());
        //IloArray<IloNumVarArray> cpSelExclusive(cplex, pricingIST.getNbDims());
        // original model vars
        IloArray<IloNumVarArray> tap_x(cplex, rmpQSet.size() + 3u);
        IloNumVarArray tap_s(cplex, rmpQSet.size() + 1);
        //IloNumVarArray tap_u(cplex, rmpQSet.size() + 1);
        // HV
        auto HV3 = 10 * (2 * (pricingIST.getNbDims() + pricingIST.getNbMeasures() + 1) + 1);
        // Linearization variables
        IloNumVarArray lin_D_in(cplex, rmpQSet.size() + 1, 0, HV3, IloNumVar::Float);
        IloNumVarArray lin_D_out(cplex, rmpQSet.size() + 1, 0, HV3, IloNumVar::Float);


        /*
         * Variables Initialisation
         */
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            cpGroupBy[i] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, ("Gamma_" + std::to_string(i)).c_str());
            cpSelection[i] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, ("psi_" + std::to_string(i)).c_str());
            cpLeftSel[i] = IloNumVarArray(cplex, pricingIST.getAdSize(i));
            cpRightSel[i] = IloNumVarArray(cplex, pricingIST.getAdSize(i));
            //cpSelExclusive[i] = IloNumVarArray(cplex, pricingIST.getAdSize(i));
            for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                cpLeftSel[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,
                                            ("l_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                cpRightSel[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,
                                             ("r_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                //cpSelExclusive[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,("c_" + std::to_string(i) + "," + std::to_string(j)).c_str());
            }
        }
        for (int i = 0; i < pricingIST.getNbMeasures(); ++i) {
            cpRightMeasure[i] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, ("alpha_" + std::to_string(i)).c_str());
            cpLeftMeasure[i] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, ("beta_" + std::to_string(i)).c_str());
        }
        // Init variables x for arcs
        std::stringstream vname;
        for (auto i = 0; i < rmpQSet.size() + 3; ++i) {
            tap_x[i] = IloNumVarArray(cplex, rmpQSet.size() + 3);
            for (auto j = 0u; j < rmpQSet.size() + 3; ++j) {
                if (i == j)
                    tap_x[i][j] = IloNumVar(cplex, 0, 0, IloNumVar::Bool,
                                            ("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                else
                    tap_x[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,
                                            ("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
            }
        }
        // Init variables for (query) selection
        for (auto i = 0; i < rmpQSet.size() + 1; ++i) {
            tap_s[i] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, ("s_" + std::to_string(i)).c_str());
        }
        // Take the new query;
        tap_s[rmpQSet.size()].setLB(1);

        // Init variables for MTZ subtour elimination and enforce part of (8)
        /*for (auto i = 1u; i <= rmpQSet.size() + 1; ++i) {
            vname << "u_" << i;
            tap_u[i - 1] = IloNumVar(cplex, 2, rmpQSet.size() + 1, TYPE_VAR_TAP, vname.str().c_str());
            vname.str("");
        }*/


        //
        // --- Pricing objective ---
        //
        IloExpr expr(cplex);
        for (auto i = 0u; i < rmpQSet.size(); ++i)
            expr += rmpIST.interest(i) * tap_s[i];
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            expr += cpSelection[i] * pricingIST.getDimWeight(i);
        }
        IloObjective obj(cplex, expr, IloObjective::Maximize);
        pricing.add(obj);
        expr.clear();
        if (debug) std::cout << "[INFO]Added Objective to pricing model\n";

        /*
         *  --- Pricing Constraints ---
         */
        // One measure per side
        for (int i = 0; i < pricingIST.getNbMeasures(); ++i)
            expr += cpLeftMeasure[i];
        pricing.add(IloRange(cplex, 1, expr, 1, "left_one_measure"));
        expr.clear();
        for (int i = 0; i < pricingIST.getNbMeasures(); ++i)
            expr += cpRightMeasure[i];
        pricing.add(IloRange(cplex, 1, expr, 1, "right_one_measure"));
        expr.clear();
        // One group key
        for (int i = 0; i < pricingIST.getNbDims(); ++i)
            expr += cpGroupBy[i];
        pricing.add(IloRange(cplex, 1, expr, 1, "one_gk"));
        expr.clear();
        // One or more for selection
        for (int i = 0; i < pricingIST.getNbDims(); ++i)
            expr += cpSelection[i];
        //TODO pricingIST.getNbDims()
        pricing.add(IloRange(cplex, 1, expr, 1, "one_more_selection"));
        expr.clear();
        // No overlap selection / group by
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            pricing.add(
                    IloRange(cplex, 0, cpSelection[i] + cpGroupBy[i], 1, ("sel_gk_" + std::to_string(i)).c_str()));
        }
        //Allow for similar predicates in left and right column if measure is different only
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            IloExpr lexp(cplex);
            IloExpr rexp(cplex);
            for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                lexp += cpLeftSel[i][j];
                rexp += cpRightSel[i][j];
            }
            pricing.add(IloRange(cplex, 0, lexp - cpSelection[i], 0));
            pricing.add(IloRange(cplex, 0, rexp - cpSelection[i], 0));
        }

        /*
         *  --- Symmetry breaking ---
         */

        for (auto k = 0u; k < pricingIST.getNbDims(); ++k) {
            for (auto i = 0u; i < k; ++i) {
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    expr += cpLeftSel[i][j];
                }
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    expr -= cpRightSel[i][j];;
                }
            }
            pricing.add(IloRange(cplex, -IloInfinity, expr, 0, ("sym_brk_series1_" + std::to_string(k)).c_str()));
            expr.clear();
        }

        /*
         *  --- Original Model Constraints ---
         */
        for (auto j = 1u; j <= rmpQSet.size() + 1; ++j) {
            for (auto i = 0u; i <= rmpQSet.size() + 1; ++i) {
                if (j != i)
                    expr += tap_x[i][j];
            }
            expr -= tap_s[j - 1];
            pricing.add(IloRange(cplex, 0, expr, 0, ("inbound_" + std::to_string(j)).c_str())); // = 0
            expr.clear();
        }

        for (auto i = 1u; i <= rmpQSet.size() + 1; ++i) {
            for (auto j = 1u; j <= rmpQSet.size() + 1 + 1; ++j) {
                if (j != i)
                    expr += tap_x[i][j];
            }
            expr += -tap_s[i - 1];
            pricing.add(IloRange(cplex, 0, expr, 0, ("outbound_" + std::to_string(i)).c_str()));
            expr.clear();
        }

        // start
        for (auto i = 1u; i <= rmpQSet.size() + 1; ++i)
            expr += tap_x[0][i];
        pricing.add(IloRange(cplex, 1, expr, 1, "path_start"));
        expr.clear();

        //end
        for (auto i = 1u; i <= rmpQSet.size() + 1; ++i)
            expr += tap_x[i][rmpQSet.size() + 1 + 1u];
        pricing.add(IloRange(cplex, 1, expr, 1, "path_end")); // Constraint (7)
        expr.clear();

        //Forbidden links (start to end, end to start ...)
        for (auto i = 0; i <= rmpQSet.size() + 1 + 1; ++i) {
            expr += tap_x[i][0];
            expr += tap_x[rmpQSet.size() + 1 + 1][i];
        }
        expr += tap_x[0][rmpQSet.size() + 1 + 1];
        pricing.add(IloRange(cplex, 0, expr, 0, "path_structure"));
        expr.clear();

        //subtour elimination constraints
        /*
        IloArray<IloRangeArray> mtz(cplex, rmpQSet.size() + 1);
        for (auto i = 1u; i <= rmpQSet.size() + 1; ++i) {
            mtz[i - 1] = IloRangeArray(cplex, rmpQSet.size() + 1);
            for (auto j = 1u; j <= rmpQSet.size() + 1; ++j) {
                expr = ((rmpQSet.size() + 1 - 1) * (1 - tap_x[i][j])) - tap_u[i - 1] + tap_u[j - 1]; // >= 1
                vname << "mtz_" << i << "_" << j;
                mtz[i - 1][j - 1] = IloRange(cplex, 1, expr, IloInfinity, vname.str().c_str());
                vname.str("");
                expr.clear();
            }
            pricing.add(mtz[i - 1]);
        }*/

        // Time
        for (auto i = 0; i < rmpQSet.size(); ++i)
            expr += tap_s[i] * (int) rmpIST.time(i);
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            expr += cpSelection[i] * timeWeights[i];
        }
        //expr + lin_T;
        pricing.add(IloRange(cplex, 0, expr, time_bound, "time_epsilon"));
        expr.clear();
        /*for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            expr += cpSelection[i] * timeWeights[i];
        }
        expr -= (1 - tap_s[rmpQSet.size()]) * HV2;
        expr -= lin_T;
        pricing.add(IloRange(cplex, -IloInfinity, expr, 0, "time_linearization"));
        expr.clear();
        pricing.add(IloRange(cplex, 0, (tap_s[rmpQSet.size()]*HV2)-lin_T));*/

        //Distance
        for (auto i = 1; i <= rmpQSet.size(); ++i) {
            for (auto j = 1; j <= rmpQSet.size(); ++j) {
                expr += tap_x[i][j] * (int) rmpIST.dist(i - 1, j - 1);
            }
        }
        for (auto i = 0; i < rmpQSet.size() + 1; ++i)
            expr += lin_D_in[i] + lin_D_out[i];
        pricing.add(IloRange(cplex, -IloInfinity, expr, dist_bound, "distance_epsilon"));
        expr.clear();
        for (int i = 0; i < rmpQSet.size(); ++i) {
            for (int k = 0; k < pricingIST.getNbDims(); ++k) {
                expr += cpSelection[k] * (1 - S[i][k]) + S[i][k] * (1 - cpSelection[k]);
            }
            expr -= (1 - tap_x[i][rmpQSet.size()+1]) * HV3;
            expr -= lin_D_in[i];
            pricing.add(IloRange(cplex, -IloInfinity, expr, 0));
            expr.clear();
        }
        for (int i = 0; i < rmpQSet.size(); ++i) {
            for (int k = 0; k < pricingIST.getNbDims(); ++k) {
                expr += cpSelection[k] * (1 - S[i][k]) + S[i][k] * (1 - cpSelection[k]);
            }
            expr -= (1 - tap_x[rmpQSet.size()+1][i]) * HV3;
            expr -= lin_D_out[i];
            pricing.add(IloRange(cplex, -IloInfinity, expr, 0));
            expr.clear();
        }
        for (int i = 0; i < rmpQSet.size(); ++i) {
            pricing.add(IloRange(cplex, 0, (tap_x[i][rmpQSet.size()+1]*HV3) - lin_D_in[i]));
            pricing.add(IloRange(cplex, 0, (tap_x[rmpQSet.size()+1][i]*HV3) - lin_D_out[i]));
        }

        /*
         *  --- Constraints forbidding having same query as existing one ---
         */
        for (auto q : rmpQSet) {
            int var_cnt = 0;
            // GB Key
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                int keyID = pricingIST.getDimId(q.getGbAttribute());
                if (i == keyID)
                    expr += cpGroupBy[i];
                else
                    expr += (1 - cpGroupBy[i]);
                var_cnt++;
            }
            // Measure - Left
            for (int i = 0; i < pricingIST.getNbMeasures(); ++i) {
                int lMeasureID = pricingIST.getMeasureId(q.getMeasureLeft());
                if (i == lMeasureID)
                    expr += cpLeftMeasure[i];
                else
                    expr += (1 - cpLeftMeasure[i]);
                var_cnt++;
            }
            // Measure - Right
            for (int i = 0; i < pricingIST.getNbMeasures(); ++i) {
                int rMeasureID = pricingIST.getMeasureId(q.getMeasureRight());
                if (i == rMeasureID)
                    expr += cpRightMeasure[i];
                else
                    expr += (1 - cpRightMeasure[i]);
                var_cnt++;
            }
            // Selection predicate - Left
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                bool dimPresent = false;
                int valueIdx = -1;
                for (auto p : q.getLeftPredicate()){
                    if (p.first == pricingIST.getDimName(i)) {
                        dimPresent = true;
                        valueIdx = p.second;
                    }
                }
                if (dimPresent){
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        if (valueIdx == j)
                            expr += cpLeftSel[i][j];
                        else
                            expr += (1 - cpLeftSel[i][j]);
                        var_cnt++;
                    }
                } else {
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        expr += (1 - cpLeftSel[i][j]);
                        var_cnt++;
                    }
                }
            }
            // Selection predicate - Right
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                bool dimPresent = false;
                int valueIdx = -1;
                for (auto p : q.getRightPredicate()){
                    if (p.first == pricingIST.getDimName(i)) {
                        dimPresent = true;
                        valueIdx = p.second;
                    }
                }
                if (dimPresent){
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        if (valueIdx == j)
                            expr += cpRightSel[i][j];
                        else
                            expr += (1 - cpRightSel[i][j]);
                        var_cnt++;
                    }
                } else {
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        expr += (1 - cpRightSel[i][j]);
                        var_cnt++;
                    }
                }
            }

            pricing.add(IloRange(cplex, 0, expr, var_cnt - 1));
            expr.clear();
        }

        /*
         *  --- Constraints forbidding having symmetries of an existing query ---
         */
        for (auto q : rmpQSet) {
            int var_cnt = 0;
            // GB Key
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                int keyID = pricingIST.getDimId(q.getGbAttribute());
                if (i == keyID)
                    expr += cpGroupBy[i];
                else
                    expr += (1 - cpGroupBy[i]);
                var_cnt++;
            }
            // Measure - Left
            for (int i = 0; i < pricingIST.getNbMeasures(); ++i) {
                int lMeasureID = pricingIST.getMeasureId(q.getMeasureLeft());
                if (i == lMeasureID)
                    expr += cpLeftMeasure[i];
                else
                    expr += (1 - cpLeftMeasure[i]);
                var_cnt++;
            }
            // Measure - Right
            for (int i = 0; i < pricingIST.getNbMeasures(); ++i) {
                int rMeasureID = pricingIST.getMeasureId(q.getMeasureRight());
                if (i == rMeasureID)
                    expr += cpRightMeasure[i];
                else
                    expr += (1 - cpRightMeasure[i]);
                var_cnt++;
            }
            // Selection predicate - Left
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                bool dimPresent = false;
                int valueIdx = -1;
                for (auto p : q.getLeftPredicate()){
                    if (p.first == pricingIST.getDimName(i)) {
                        dimPresent = true;
                        valueIdx = p.second;
                    }
                }
                if (dimPresent){
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        if (valueIdx == j)
                            expr += cpRightSel[i][j];
                        else
                            expr += (1 - cpRightSel[i][j]);
                        var_cnt++;
                    }
                } else {
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        expr += (1 - cpRightSel[i][j]);
                        var_cnt++;
                    }
                }
            }
            // Selection predicate - Right
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                bool dimPresent = false;
                int valueIdx = -1;
                for (auto p : q.getRightPredicate()){
                    if (p.first == pricingIST.getDimName(i)) {
                        dimPresent = true;
                        valueIdx = p.second;
                    }
                }
                if (dimPresent){
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        if (valueIdx == j)
                            expr += cpLeftSel[i][j];
                        else
                            expr += (1 - cpLeftSel[i][j]);
                        var_cnt++;
                    }
                } else {
                    for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                        expr += (1 - cpLeftSel[i][j]);
                        var_cnt++;
                    }
                }
            }

            pricing.add(IloRange(cplex, 0, expr, var_cnt - 1));
            expr.clear();
        }

        //Init solver
        IloCP cp(pricing);

        cp.solve();

        /*
             *  --- Recover query from solution ---
             */
        //IloNumArray solSelection(cplex);
        IloNumArray solGBKey(cplex);
        IloNumArray solLeftMeasure(cplex);
        IloNumArray solRightMeasure(cplex);
        IloArray<IloNumArray> solLeftSelection(cplex, pricingIST.getNbDims());
        IloArray<IloNumArray> solRightSelection(cplex, pricingIST.getNbDims());

        //cplex_solver.getValues(solSelection, cpSelection);
        cp.getValues(cpGroupBy, solGBKey);
        cp.getValues(cpLeftMeasure, solLeftMeasure);
        cp.getValues(cpRightMeasure, solRightMeasure);
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            IloNumArray line(cplex);
            cp.getValues(cpLeftSel[i], line);
            solLeftSelection[i] = line;
        }
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            IloNumArray line(cplex);
            cp.getValues(cpRightSel[i], line);
            solRightSelection[i] = line;
        }

        int gbKeyIdx = 0;
        while (gbKeyIdx < pricingIST.getNbDims()) {
            if (abs(solGBKey[gbKeyIdx] - 1) < 10e-6)
                break;
            else
                gbKeyIdx++;
        }
        int lmIdx = 0;
        while (lmIdx < pricingIST.getNbMeasures()) {
            if (abs(solLeftMeasure[lmIdx] - 1) < 10e-6)
                break;
            else
                lmIdx++;
        }
        int rmIdx = 0;
        while (rmIdx < pricingIST.getNbMeasures()) {
            if (abs(solRightMeasure[rmIdx] - 1) < 10e-6)
                break;
            else
                rmIdx++;
        }
        std::vector<std::pair<std::string, int>> lPredicate;
        std::vector<std::pair<std::string, int>> rPredicate;
        for (int i = 0; i < pricingIST.getNbDims(); ++i) {
            for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                if (abs(solLeftSelection[i][j] - 1) < 10e-6)
                    lPredicate.emplace_back(make_pair(pricingIST.getDimName(i), j));
                if (abs(solRightSelection[i][j] - 1) < 10e-6)
                    rPredicate.emplace_back(make_pair(pricingIST.getDimName(i), j));
            }
        }

        Query picked = Query(pricingIST.getTableName(), "sum", pricingIST.getDimName(gbKeyIdx),
                             pricingIST.getMeasureName(lmIdx), pricingIST.getMeasureName(rmIdx),
                             lPredicate, rPredicate);

        cout << "[Pricing Query] - " << picked << endl;


        cp.end();
        cplex.end();



        std::vector<int> empty;
        return Solution(false, 0, 0, empty, 0);
    }

    Instance princingCPSolver::buildRMPInstance(vector<Query>& queries) const {
        vector<double> interest = JVMAdapter::getInterest(queries, pricingIST);
        vector<int> time = JVMAdapter::getTime(queries, pricingIST);
        vector<vector<int>> distMatrix;
        distMatrix.reserve(queries.size());
        for (int i = 0; i < queries.size(); ++i) {
            vector<int> line;
            line.reserve(queries.size());
            for (int j = 0; j < queries.size(); ++j) {
                line.emplace_back(queries[i].dist(queries[j]));
            }
            distMatrix.emplace_back(line);
        }
        return {static_cast<int>(queries.size()), interest, time, distMatrix};
    }
} // cplex_tap
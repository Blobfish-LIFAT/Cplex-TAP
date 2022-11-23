//
// Created by alex on 21/07/22.
//

#include "pricingSolver.h"
#include "Query.h"
#include <random>
#include "JVMAdapter.h"
#include "solver.h"
#include <numeric>


namespace cplex_tap {

    Solution pricingSolver::solve() const {
        std::cout << "[INFO] CLK_RATE " << CLOCKS_PER_SEC << std::endl;

        vector<Query> rmpQSet;

        if (extStarting.size() == 0) {
            int starting_count = 50;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::cout << "[STEP] Building Initial RMP query set" << std::endl;
            std::uniform_int_distribution<> rdAttr(0, pricingIST.getNbDims() - 1);
            for (; rmpQSet.size() < starting_count;) {
                int lAttrID = rdAttr(gen);
                //int rAttrID = rdAttr(gen);
                int measureID = 0;
                int gbAttr = rdAttr(gen);
                while (gbAttr == lAttrID) {// || gbAttr == rAttrID
                    gbAttr = rdAttr(gen);
                }
                std::uniform_int_distribution<> rdValLeft(0, pricingIST.getAdSize(lAttrID) - 1);
                //std::uniform_int_distribution<> rdValRight(0, pricingIST.getAdSize(rAttrID)-1);

                std::vector<std::pair<string, int> > lPredicate = {{pricingIST.getDimName(lAttrID), rdValLeft(gen)}};
                std::vector<std::pair<string, int> > rPredicate = {{pricingIST.getDimName(lAttrID), rdValLeft(gen)}};
                //std::vector<std::pair<string, int> > rPredicate = { {pricingIST.getDimName(rAttrID), rdValRight(gen)}};
                Query rdQ = Query(pricingIST.getTableName(), "sum", pricingIST.getDimName(gbAttr),
                                  pricingIST.getMeasureName(measureID), pricingIST.getMeasureName(measureID),
                                  lPredicate, rPredicate);

                bool duplicate = false;
                for (Query &q : rmpQSet)
                    duplicate |= q == rdQ;
                if (!duplicate )rmpQSet.emplace_back(rdQ);
            }
            std::cout << "[STEP][END] Building Initial RMP query set " << rmpQSet.size() << std::endl;
        } else{
            rmpQSet.insert(rmpQSet.end(), extStarting.begin(), extStarting.end());
        }

        double prevRmpObj = 0;
        vector<double> objValues;
        bool isNewQuerySelected = true;
        Solution rmpSol(false, 0, 0, std::vector<int>());

        int it = 0;
        time_t start;
        start = clock();
        while (it++ < 300) {

            std::cout << "[STEP] Building RMP model" << std::endl;
            Instance rmpIST = buildRMPInstance(rmpQSet);
            Solver tapSolver = Solver(rmpIST);
            rmpSol = tapSolver.solve(dist_bound, time_bound, false, "");
            std::cout << "[STEP][END] Building RMP model - z*=" << std::to_string(rmpSol.z) << " (prev=" << prevRmpObj << ")" << std::endl;
            prevRmpObj = rmpSol.z;
            objValues.emplace_back(rmpSol.z);

            cout << "[Solution DUMP]";
            for (int i = 0; i < rmpSol.sequence.size(); ++i) {
                cout << rmpQSet[rmpSol.sequence[i]-1];
                if (i < rmpSol.sequence.size() - 1)
                    cout << ";";
            }
            cout << endl;

            if (assessConvergence(objValues))
                break;

            // Init CPLEX environment and model objects
            IloEnv cplex;
            IloModel pricing(cplex);

            /*
             *  Data about queries in the RMP
             */
            vector<vector<bool>> S;
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
            IloArray<IloNumVarArray> cpSelExclusive(cplex, pricingIST.getNbDims());
            // original model vars
            IloArray<IloNumVarArray> tap_x(cplex, rmpQSet.size() + 3u);
            IloNumVarArray tap_s(cplex, rmpQSet.size() + 1);
            IloNumVarArray tap_u(cplex, rmpQSet.size() + 1);
            // HV1 sum of all weights
            auto dimWeights = pricingIST.getDimWeights();
            auto HV1 = std::accumulate(dimWeights.begin(), dimWeights.end(), decltype(dimWeights)::value_type(0));
            // HV1 sum of all time weights
            auto timeWeights = pricingIST.getDimWeights();
            auto HV2 = std::accumulate(timeWeights.begin(), timeWeights.end(), decltype(timeWeights)::value_type(0));
            auto HV3 = 10 * (2 * (pricingIST.getNbDims() + pricingIST.getNbMeasures() + 1) + 1);
            // Linearization variables
            IloNumVar lin_I(cplex, 0, HV1, IloNumVar::Float, "I");
            IloNumVar lin_T(cplex, 0, HV2, IloNumVar::Float, "T");
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
                cpSelExclusive[i] = IloNumVarArray(cplex, pricingIST.getAdSize(i));
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    cpLeftSel[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,
                                                ("l_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                    cpRightSel[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,
                                                 ("r_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                    cpSelExclusive[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,
                                                     ("c_" + std::to_string(i) + "," + std::to_string(j)).c_str());
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
                        tap_x[i][j] = IloNumVar(cplex, 0, 0, IloNumVar::Float,
                                                ("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                    else
                        tap_x[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Float,
                                                ("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                }
            }
            // Init variables for (query) selection
            for (auto i = 0; i < rmpQSet.size(); ++i) {
                tap_s[i] = IloNumVar(cplex, 0, 1, IloNumVar::Float, ("s_" + std::to_string(i)).c_str());
            }
            // Last one from the pricing is binary taken or not
            tap_s[rmpQSet.size()] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, "s_new");

            // Init variables for MTZ subtour elimination and enforce part of (8)
            for (auto i = 1u; i <= rmpQSet.size() + 1; ++i) {
                vname << "u_" << i;
                tap_u[i - 1] = IloNumVar(cplex, 2, rmpQSet.size() + 1, IloNumVar::Float, vname.str().c_str());
                vname.str("");
            }


            //
            // --- Pricing objective ---
            //
            IloExpr expr(cplex);
            for (auto i = 0u; i < rmpQSet.size(); ++i)
                expr += rmpIST.interest(i) * tap_s[i];
            expr += lin_I;
            IloObjective obj(cplex, expr, IloObjective::Maximize);
            pricing.add(obj);
            expr.clear();
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                expr += cpSelection[i] * pricingIST.getDimWeight(i);
            }
            expr += (1 - tap_s[rmpQSet.size()]) * HV1;
            expr -= lin_I;
            pricing.add(IloRange(cplex, 0, expr, IloInfinity, "lin_I"));
            expr.clear();
            pricing.add(IloRange(cplex, 0, (tap_s[rmpQSet.size()]*HV1) - lin_I));
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
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    pricing.add(IloRange(cplex, 0, cpLeftSel[i][j] + cpRightSel[i][j], 1)); // (12)
                    // (13)
                    pricing.add(IloRange(cplex, 0, ((cpLeftSel[i][j] + cpRightSel[i][j]) / 2.0) - cpSelExclusive[i][j],
                                         IloInfinity));
                    expr += cpSelection[i];
                    for (int k = 0; k < pricingIST.getAdSize(i); ++k) {
                        if (k != j) {
                            expr -= cpLeftSel[i][k];
                            expr -= cpRightSel[i][k];
                        }
                    }
                    expr -= cpSelExclusive[i][j];
                    pricing.add(IloRange(cplex, -IloInfinity, expr, 0));
                    expr.clear();
                    // (14)
                    for (int k = 0; k < pricingIST.getNbMeasures(); ++k) {
                        pricing.add(IloRange(cplex, 0, 2 - cpSelExclusive[i][j] - cpLeftMeasure[k] - cpRightMeasure[k],
                                             IloInfinity));
                    }
                }
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
            }

            // Time
            for (auto i = 0; i < rmpQSet.size(); ++i)
                expr += tap_s[i] * (int) rmpIST.time(i);
            expr + lin_T;
            pricing.add(IloRange(cplex, 0, expr, time_bound, "time_epsilon"));
            expr.clear();
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                expr += cpSelection[i] * timeWeights[i];
            }
            expr -= (1 - tap_s[rmpQSet.size()]) * HV2;
            expr -= lin_T;
            pricing.add(IloRange(cplex, -IloInfinity, expr, 0, "time_linearization"));
            expr.clear();
            pricing.add(IloRange(cplex, 0, (tap_s[rmpQSet.size()]*HV2)-lin_T));

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

            //Init solver
            IloCplex cplex_solver(pricing);
            cplex_solver.setParam(IloCplex::Param::TimeLimit, 3600);
            cplex_solver.setParam(IloCplex::Param::Threads, 1);
            cplex_solver.setOut(cplex.getNullStream());

            bool solved = false;
            try {
                solved = cplex_solver.solve();
            }
            catch (const IloException &e) {
                std::cerr << "\n\n--- CPLEX Exception (Pricing) ---\n";
                std::cerr << e << "\n";
                cplex.end();
                throw;
            }
            if (solved) {
                //std::cout << "\n--- Solver success ---\n";
                //std::cout << "    Status: " << cplex_solver.getStatus() << "\n";
                //std::cout << "    Objective: " << cplex_solver.getObjValue() << "\n";
            } else {
                std::cerr << "\n--- Solver Error (Pricing) ---\n";
                std::cerr << "    Status: " << cplex_solver.getStatus() << "\n";
                std::cerr << "    Error details: " << cplex_solver.getCplexStatus() << "\n";
            }

            /*
             *  --- Recover query from solution ---
             */
            //IloNumArray solSelection(cplex);
            IloNumArray solGBKey(cplex);
            IloNumArray solLeftMeasure(cplex);
            IloNumArray solRightMeasure(cplex);
            IloArray<IloNumArray> solLeftSelection(cplex, pricingIST.getNbDims());
            IloArray<IloNumArray> solRightSelection(cplex, pricingIST.getNbDims());
            isNewQuerySelected = cplex_solver.getValue(tap_s[rmpQSet.size()]);

            //cplex_solver.getValues(solSelection, cpSelection);
            cplex_solver.getValues(solGBKey, cpGroupBy);
            cplex_solver.getValues(solLeftMeasure, cpLeftMeasure);
            cplex_solver.getValues(solRightMeasure, cpRightMeasure);
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                IloNumArray line(cplex);
                cplex_solver.getValues(line, cpLeftSel[i]);
                solLeftSelection[i] = line;
            }
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                IloNumArray line(cplex);
                cplex_solver.getValues(line, cpRightSel[i]);
                solRightSelection[i] = line;
            }

            int gbKeyIdx = 0;
            for (; gbKeyIdx < pricingIST.getNbDims(); ++gbKeyIdx) {
                if (solGBKey[gbKeyIdx] == 1)
                    break;
            }
            int lmIdx = 0;
            for (; lmIdx < pricingIST.getNbMeasures(); ++lmIdx) {
                if (solLeftMeasure[lmIdx] == 1)
                    break;
            }
            int rmIdx = 0;
            for (; rmIdx < pricingIST.getNbMeasures(); ++rmIdx) {
                if (solRightMeasure[rmIdx] == 1)
                    break;
            }
            std::vector<std::pair<std::string, int>> lPredicate;
            std::vector<std::pair<std::string, int>> rPredicate;
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    if (solLeftSelection[i][j] == 1)
                        lPredicate.emplace_back(make_pair(pricingIST.getDimName(i), j));
                    if (solRightSelection[i][j] == 1)
                        rPredicate.emplace_back(make_pair(pricingIST.getDimName(i), j));
                }
            }

            Query picked = Query(pricingIST.getTableName(), "sum", pricingIST.getDimName(gbKeyIdx),
                                 pricingIST.getMeasureName(lmIdx), pricingIST.getMeasureName(rmIdx),
                                 lPredicate, rPredicate);

            cout << "[Pricing Query] " << isNewQuerySelected << " - " << picked << endl;
            time_t end = clock();
            cout << "[Pricing z value] " << cplex_solver.getObjValue() << endl;
            double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
            cout << "[TIME][ITER][s] " << time_to_sol << endl;
            rmpQSet.emplace_back(picked);
            cplex_solver.end();
            cplex.end();
        }

        /*for (int i = 0; i < rmpQSet.size(); ++i) {
            cout << i << " - " << rmpQSet[i] << endl;
        }
        for (int i = 0; i < rmpSol.sequence.size(); ++i) {
            cout << rmpSol.sequence[i] << " ";
        }
        cout << endl;*/
        cout << "[OBJ]";
        for (auto it = objValues.begin(); it != objValues.end() ; ++it) {
            cout << *it;
            if (it != objValues.end()-1)
                cout << std::string(";");

        }
        cout <<endl;
        cout << "[INFO] iterations " << objValues.size() << endl;

        for (auto i : rmpSol.sequence) {
            cout << rmpQSet[i] << endl;
        }
        return rmpSol;
    }

    bool pricingSolver::assessConvergence(vector<double> objValues){
        int depth = 150;
        double epsilon = 10e-8;

        if (objValues.size() < depth)
            return false;

        bool change = false;
        for (auto it = objValues.rbegin();*it != objValues[objValues.size()-depth];it++){
            change |= *it - *(it + 1) > epsilon;
        }

        return !change;

    }

    Instance pricingSolver::buildRMPInstance(vector<Query> queries) const {
        vector<double> interest = JVMAdapter::getInterest(queries, pricingIST);
        vector<int> time = JVMAdapter::getTime(queries, pricingIST);
        vector<vector<int>> distMatrix;
        for (int i = 0; i < queries.size(); ++i) {
            vector<int> line;
            for (int j = 0; j < queries.size(); ++j) {
                line.emplace_back(queries[i].dist(queries[j]));
            }
            distMatrix.emplace_back(line);
        }
        return {static_cast<int>(queries.size()), interest, time, distMatrix};
    }
} // cplex_tap
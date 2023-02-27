//
// Created by alex on 21/07/22.
//

#include "pricingSolver.h"
#include "Query.h"
#include <random>
#include "JVMAdapter.h"
#include "solver.h"
#include <numeric>
#include "SolverVPLSHammingSX.h"

//#define TYPE_VAR_TAP ILOINT
#define TYPE_VAR_TAP ILOFLOAT
#define NO_PRINT false

namespace cplex_tap {

    ILOMIPINFOCALLBACK0(print_z) {
        IloEnv env = getEnv();
        env.out() << "[INFO MIP Callback]" << " CLK " << getCplexTime() << " Z " << getBestObjValue() << " ? " << hasIncumbent() << std::endl;
    }

    Solution pricingSolver::solve() const {
        std::cout << "[INFO] CLK_RATE " << CLOCKS_PER_SEC << std::endl;
        int global_t = 0;

        vector<Query> rmpQSet;
        rmpQSet.insert(rmpQSet.end(), extStarting.begin(), extStarting.end());

        double prevRmpObj = 0;
        vector<double> objValues;
        bool isNewQuerySelected = true;
        //Solution rmpSol(false, 0, 0, std::vector<int>());
        //Solution prevRMPSol(false, 0, 0, std::vector<int>());

        int it = 0;
        time_t start;
        start = clock();
        while (it++ < 100) {

            if (!NO_PRINT) std::cout << "[STEP] Building RMP model" << std::endl;
            Instance rmpIST = buildRMPInstance(rmpQSet);
            /*
            Solver tapSolver = Solver(rmpIST);
            tapSolver.setTimeout(master_it_timeout);
            prevRMPSol = rmpSol;
            rmpSol = tapSolver.solve(dist_bound, time_bound, false, "");
            if (!NO_PRINT) std::cout << "[STEP][END] Building RMP model - z*=" << std::to_string(rmpSol.z) << std::endl;
            prevRmpObj = rmpSol.z;
            objValues.emplace_back(rmpSol.z);

            if (debug) {
                cout << "[Solution DUMP]";
                for (int i = 0; i < rmpSol.sequence.size(); ++i) {
                    cout << rmpQSet[rmpSol.sequence[i] - 1];
                    if (i < rmpSol.sequence.size() - 1)
                        cout << ";";
                }
                cout << endl;
            }

            // Check convergence criterion
            if (assessConvergence(objValues)){
                cout << "[BREAK] Reason: convergence" << endl;
                break;
            }
            if (!rmpSol.optimal){
                //cout << "[BREAK] Reason: rmp timeout" << endl;
                //break;
            }*/

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
            //IloNumVarArray tap_u(cplex, rmpQSet.size() + 1);
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
                        tap_x[i][j] = IloNumVar(cplex, 0, 0, TYPE_VAR_TAP,
                                                ("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                    else
                        tap_x[i][j] = IloNumVar(cplex, 0, 1, TYPE_VAR_TAP,
                                                ("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                }
            }
            // Init variables for (query) selection
            for (auto i = 0; i < rmpQSet.size(); ++i) {
                tap_s[i] = IloNumVar(cplex, 0, 1, TYPE_VAR_TAP, ("s_" + std::to_string(i)).c_str());
            }
            // Last one from the pricing is binary taken or not
            //TODO remove me to allow selection or not
            tap_s[rmpQSet.size()] = IloNumVar(cplex, 1, 1, IloNumVar::Bool, "s_new");

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
            //expr += lin_I;
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                expr += cpSelection[i] * pricingIST.getDimWeight(i);
            }
            IloObjective obj(cplex, expr, IloObjective::Maximize);
            pricing.add(obj);
            expr.clear();

            //expr += (1 - tap_s[rmpQSet.size()]) * HV1;
            //expr -= lin_I;
            //pricing.add(IloRange(cplex, 0, expr, IloInfinity, "lin_I"));
            //expr.clear();
            //pricing.add(IloRange(cplex, 0, (tap_s[rmpQSet.size()]*HV1) - lin_I));
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
            }/*

            for (auto i = 0u; i < pricingIST.getNbDims(); ++i) {
                for (int k = 0; k < pricingIST.getAdSize(i); ++k) {
                    for (int j = 0; j < k; ++j) {
                        expr += cpLeftSel[i][j];
                    }
                    for (int j = 0; j < k; ++j) {
                        expr -= cpRightSel[i][j];;
                    }
                    expr -= 1;
                    for (int l = 0; l < pricingIST.getAdSize(i); ++l) {
                        expr += cpRightSel[i][l];
                    }
                    pricing.add(IloRange(cplex, -IloInfinity, expr, 0,
                                         ("sym_brk_series2_" + std::to_string(k) + "_" + std::to_string(i)).c_str()));
                    expr.clear();
                }
            }*/

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
            IloCplex cplex_solver(pricing);
            cplex_solver.setParam(IloCplex::Param::TimeLimit, pricing_it_timeout);
            cplex_solver.setParam(IloCplex::Param::Threads, 1);
            cplex_solver.setParam(IloCplex::Param::Preprocessing::Symmetry, cplex_sym);
            cplex_solver.setParam(IloCplex::Param::Preprocessing::Aggregator, 0);
            cplex_solver.setParam(IloCplex::Param::Preprocessing::Presolve	, 0);
            cplex_solver.setParam(IloCplex::Param::MIP::Display, 0); //5 is max
            cplex_solver.setParam(IloCplex::Param::Simplex::Display, 0);
            cplex_solver.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort, 0);
            cplex_solver.setParam(IloCplex::Param::Emphasis::MIP	, 2);
            cplex_solver.setParam(IloCplex::Param::MIP::SubMIP::NodeLimit, 50);
            //cplex_solver.setOut(cplex.getNullStream());

            int allowed_restarts = 29;

            //IloCplex::Callback mycallback = cplex_solver.use(print_z(cplex));
            std::cout << "[INFO] CLK_START " << clock() << std::endl;

            bool solved = false;
            try {
                solved = cplex_solver.solve();
            }
            catch (const IloException &e) {
                std::cout << "\n\n--- CPLEX Exception (Pricing) ---\n";
                std::cout << e << "\n";
                cplex.end();
                throw;
            }
            if (!solved) {
                if (cplex_solver.getCplexStatus() == IloCplex::AbortTimeLim) {
                    while (cplex_solver.getCplexStatus() == IloCplex::AbortTimeLim && allowed_restarts > 0 &&
                           !(cplex_solver.getStatus() == IloAlgorithm::Optimal ||
                             cplex_solver.getStatus() == IloAlgorithm::Feasible)) {
                        std::cout << "[WARN] restarting solver after time limit" << endl;
                        allowed_restarts--;
                        cplex_solver.solve();
                    }
                } else {
                    std::cout << "\n--- Solver Error (Pricing) ---\n";
                    std::cout << "    Status: " << cplex_solver.getStatus() << "\n";
                    std::cout << "    Error details: " << cplex_solver.getCplexStatus() << "\n";
                }
            } else {
            }

            if (!(cplex_solver.getStatus() == IloAlgorithm::Optimal || cplex_solver.getStatus() == IloAlgorithm::Feasible)){
                cout << "[Info] pricing timeout without feasible solution" << endl;
                break;
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

            if (!NO_PRINT) cout << "[Pricing Query] " << isNewQuerySelected << " - " << picked << endl;
            time_t end = clock();
            if (!NO_PRINT) cout << "[Pricing z value] " << cplex_solver.getObjValue() << endl;
            double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
            if (!NO_PRINT) cout << "[TIME][ITER][s] " << time_to_sol << endl;
            rmpQSet.emplace_back(picked);

            cplex_solver.end();
            cplex.end();

            if (time_to_sol > global_timeout){
                cout << "[BREAK] Reason: global timeout" << endl;
                break;
            }
        }

        cout << "[OBJ]";
        for (auto it = objValues.begin(); it != objValues.end() ; ++it) {
            cout << *it;
            if (it != objValues.end()-1)
                cout << std::string(";");

        }
        cout <<endl;
        cout << "[INFO] iterations " << objValues.size() << endl;

        /*
        // If we time out on last mip we keep the previous solution
        if (rmpSol.optimal){
            for (auto i : rmpSol.sequence) {
                cout << rmpQSet[i] << endl;
            }
            return rmpSol;
        } else{
            for (auto i : prevRMPSol.sequence) {
                cout << rmpQSet[i] << endl;
            }
            return prevRMPSol;
        }*/

        Instance rmpIST = buildRMPInstance(rmpQSet);
        auto tapSolver = Solver(rmpIST);
        tapSolver.setTimeout(master_it_timeout);
        auto final_sol = tapSolver.solve(dist_bound, time_bound, false, "");
        cout << "[MASTER] " << final_sol.z << "|" << final_sol.optimal << endl;

        for ( auto qid : final_sol.sequence){
            cout << rmpQSet[qid] << endl;
        }

        auto matheuristic = SolverVPLSHammingSX(rmpIST, 15, 15, 30, 20);
        auto mathsol = matheuristic.solve(dist_bound, time_bound, false, "");
        cout << "[MASTER][VPLS] " << mathsol.z << " | " << final_sol.optimal <<  " | " << mathsol.time << endl;

        return final_sol;

}

bool pricingSolver::assessConvergence(vector<double> objValues){
    int depth = 50;
    double epsilon = 10e-8;

    if (objValues.size() < depth)
        return false;

    bool change = false;
    for (auto it = objValues.rbegin();*it != objValues[objValues.size()-depth];it++){
        change |= *it - *(it + 1) > epsilon;
    }

    return !change;

}

Instance pricingSolver::buildRMPInstance(vector<Query>& queries) const {
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

int pricingSolver::getPricingItTimeout() const {
    return pricing_it_timeout;
}

void pricingSolver::setPricingItTimeout(int pricingItTimeout) {
    pricing_it_timeout = pricingItTimeout;
}

int pricingSolver::getMasterItTimeout() const {
    return master_it_timeout;
}

void pricingSolver::setMasterItTimeout(int masterItTimeout) {
    master_it_timeout = masterItTimeout;
}

void pricingSolver::setCplexSym(int cplexSym) {
    cplex_sym = cplexSym;
}
} // cplex_tap
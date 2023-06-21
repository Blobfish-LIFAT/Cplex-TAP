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
#include "KnapsackSolver.h"

#define NO_PRINT true

namespace cplex_tap {

    pricingSolver::ConfMap pricingSolver::confMap_ = {
            //Symetry breaking -> Constraints (21) and (22) : 1 - Constraints (23) : 2 - CPLEX : 3
            //CPLEX Presolve -> ON : 1 - OFF : 0
            //CPLEX SubMIP Node limit -> 500 - 250 - 50
            //CPLEX MIP Emphasis -> balanced : 0 - feasibility : 1 - optimality : 2
            {"best", { 1, 0, 250, 1}},
            {"default", { 3, 1, 500, 0}},
            {"cstr_2122", { 1, 1, 500, 0}},
            {"cstr_23", { 2, 1, 500, 0}},
            {"presolve_off", { 3, 0, 500, 0}},
            {"node_lim_50", { 3, 1, 50, 0}},
            {"node_lim_250", { 3, 1, 250, 0}},
            {"emp_feas", { 3, 1, 500, 1}},
            {"emp_opt", { 3, 1, 500, 2}}

    };

    ILOMIPINFOCALLBACK0(print_z) {
        IloEnv env = getEnv();
        env.out() << "[INFO MIP Callback]" << " CLK " << getCplexTime() << " Z " << getBestObjValue() << " ? " << hasIncumbent() << std::endl;
    }

    Solution pricingSolver::solve() const {
        std::cout << "[INFO] CLK_RATE " << CLOCKS_PER_SEC << std::endl;
        int global_t = 0;

        const int N_Dims = pricingIST.getNbDims();
        vector<int> timing_flat;
        for (int i = 0; i < N_Dims; ++i) {
            for (int j = i + 1; j < N_Dims; ++j) {
                timing_flat.push_back(pricingIST.getTiming()[i][j]);
            }
        }

        vector<Query> rmpQSet;
        rmpQSet.insert(rmpQSet.end(), extStarting.begin(), extStarting.end());

        for (int qid = 0; qid < rmpQSet.size(); ++qid) {
            std::cout << rmpQSet[qid] << endl;
        }

        double prevRmpObj = 0;
        vector<double> objValues;
        bool isNewQuerySelected = true;

        vector<vector<double>> I;

        UserProfile* sing2 = UserProfile::GetInstance();
        std::unordered_map<std::string, std::unordered_map<std::string,int>> imap = sing2->value();
        ActiveDomains* adSingleton = ActiveDomains::GetInstance();
        std::unordered_map<std::string,std::vector<std::string>> ad = adSingleton->value();

        for (int dimId = 0; dimId < N_Dims; ++dimId) {
            double tableRows = pricingIST.getNbRows();
            string dim = pricingIST.getDimName(dimId);
            vector<double> intCoefs;
            intCoefs.reserve(pricingIST.getAdSize(dimId));
            for (int i = 0; i < pricingIST.getAdSize(dimId); ++i) {
                intCoefs.push_back(imap[dim][ad[dim][i]] / tableRows);
            }
            I.push_back(intCoefs);
        }

        int it = 0;
        time_t start;
        start = clock();
        while (it++ < 100) {

            if (!NO_PRINT) std::cout << "[STEP] Building RMP model" << std::endl;
            Instance rmpIST = buildRMPInstance(rmpQSet);

            // Init CPLEX environment and model objects
            IloEnv cplex;
            IloModel pricing(cplex);

            /*
             *  Data about queries in the RMP
             */
            vector<vector<bool>> S; // selection attribute
            vector<vector<bool>> G; // Group by attribute
            for (int i = 0; i < rmpQSet.size(); ++i) {
                vector<bool> tmp;
                vector<bool> tmp2;
                for (int j = 0; j < N_Dims; ++j) {
                    bool appears = false;
                    const string dimName = pricingIST.getDimName(j);
                    for (int k = 0; k < rmpQSet[i].getLeftPredicate().size(); ++k) {
                        appears |= rmpQSet[i].getLeftPredicate()[k].first == dimName;
                    }
                    for (int k = 0; k < rmpQSet[i].getRightPredicate().size(); ++k) {
                        appears |= rmpQSet[i].getRightPredicate()[k].first == dimName;
                    }
                    tmp.emplace_back(appears);
                    tmp2.emplace_back(rmpQSet[i].getGbAttribute() == dimName);
                }
                S.emplace_back(tmp);
                G.emplace_back(tmp2);
            }

            /*
             * Variables Declaration
             */
            IloNumVarArray cpLeftMeasure(cplex, pricingIST.getNbMeasures());
            IloNumVarArray cpRightMeasure(cplex, pricingIST.getNbMeasures());
            IloBoolVarArray cpGroupBy(cplex, N_Dims);
            IloBoolVarArray cpSelection(cplex, N_Dims);
            IloArray<IloNumVarArray> cpLeftSel(cplex, N_Dims);
            IloArray<IloNumVarArray> cpRightSel(cplex, N_Dims);
            IloBoolVarArray cpCombiAB(cplex, (N_Dims * (N_Dims - 1) ) / 2);
            IloBoolVarArray cpCombiBA(cplex, (N_Dims * (N_Dims - 1) ) / 2);
            //IloArray<IloNumVarArray> cpSelExclusive(cplex, pricingIST.getNbDims());
            // original model vars
            IloArray<IloNumVarArray> tap_x(cplex, rmpQSet.size() + 3u);
            IloNumVarArray tap_s(cplex, rmpQSet.size() + 1);
            //IloNumVarArray tap_u(cplex, rmpQSet.size() + 1);
            // HV1 sum of all weights
            //auto dimWeights = pricingIST.getDimWeights();
            //auto HV1 = std::accumulate(dimWeights.begin(), dimWeights.end(), decltype(dimWeights)::value_type(0));
            // HV1 sum of all time weights
            //auto timeWeights = pricingIST.getDimWeights();
            //auto HV2 = std::accumulate(timeWeights.begin(), timeWeights.end(), decltype(timeWeights)::value_type(0));
            auto HV3 = 10 * (2 * (N_Dims + pricingIST.getNbMeasures() + 1) + 1);
            // Linearization variables
            //IloNumVar lin_I(cplex, 0, HV1, IloNumVar::Float, "I");
            //IloNumVar lin_T(cplex, 0, HV2, IloNumVar::Float, "T");
            IloNumVarArray lin_D_in(cplex, rmpQSet.size() + 1, 0, HV3, IloNumVar::Float);
            IloNumVarArray lin_D_out(cplex, rmpQSet.size() + 1, 0, HV3, IloNumVar::Float);


            /*
             * Variables Initialisation
             */
            for (int i = 0; i < N_Dims; ++i) {
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
                    //cpSelExclusive[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, ("c_" + std::to_string(i) + "," + std::to_string(j)).c_str());
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
                    if (i == j) {
                        if (useFloat)
                            tap_x[i][j] = IloNumVar(cplex, 0, 0, IloNumVar::Float,("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                        else
                            tap_x[i][j] = IloNumVar(cplex, 0, 0, IloNumVar::Bool,("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                    } else {
                        if (useFloat)
                            tap_x[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Float,("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                        else
                            tap_x[i][j] = IloNumVar(cplex, 0, 1, IloNumVar::Bool,("x_" + std::to_string(i) + "," + std::to_string(j)).c_str());
                    }
                }
            }

            // Init variables for (query) selection
            for (auto i = 0; i < rmpQSet.size(); ++i) {
                if (useFloat)
                    tap_s[i] = IloNumVar(cplex, 0, 1, IloNumVar::Float, ("s_" + std::to_string(i)).c_str());
                else
                    tap_s[i] = IloNumVar(cplex, 0, 1, IloNumVar::Bool, ("s_" + std::to_string(i)).c_str());
            }
            // Last one from the pricing is binary taken or not
            //TODO remove me to allow selection or not
            tap_s[rmpQSet.size()] = IloNumVar(cplex, 1, 1, IloNumVar::Bool, "s_new");


            //
            // --- Pricing objective ---
            //
            IloExpr expr(cplex);
            for (auto i = 0u; i < rmpQSet.size(); ++i)
                expr += rmpIST.interest(i) * tap_s[i];
            for (int i = 0; i < N_Dims; ++i) {
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    expr += I[i][j] * cpLeftSel[i][j];
                    expr += I[i][j] * cpRightSel[i][j];
                }
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
            for (int i = 0; i < N_Dims; ++i)
                expr += cpGroupBy[i];
            pricing.add(IloRange(cplex, 1, expr, 1, "one_gk"));
            expr.clear();
            // One or more for selection
            for (int i = 0; i < N_Dims; ++i)
                expr += cpSelection[i];
            //TODO pricingIST.getNbDims()
            pricing.add(IloRange(cplex, 1, expr, 1, "one_more_selection"));
            expr.clear();
            // No overlap selection / group by
            for (int i = 0; i < N_Dims; ++i) {
                pricing.add(
                        IloRange(cplex, 0, cpSelection[i] + cpGroupBy[i], 1, ("sel_gk_" + std::to_string(i)).c_str()));
            }
            //Allow for similar predicates in left and right column if measure is different only
            for (int i = 0; i < N_Dims; ++i) {
                IloExpr lexp(cplex);
                IloExpr rexp(cplex);
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    lexp += cpLeftSel[i][j];
                    rexp += cpRightSel[i][j];
                }
                pricing.add(IloRange(cplex, 0, lexp - cpSelection[i], 0));
                pricing.add(IloRange(cplex, 0, rexp - cpSelection[i], 0));
            }
            for (int i = 0; i < N_Dims; ++i) {
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    pricing.add(IloRange(cplex, 0, cpLeftSel[i][j] + cpRightSel[i][j], 1)); // (12)
                }
            }
            /*
             *  --- Symmetry breaking ---
             */
            if (confMap_[selectedConf][0] == 1) {
                for (auto k = 0u; k < N_Dims; ++k) {
                    for (auto i = 0u; i < k; ++i) {
                        for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                            expr += cpLeftSel[i][j];
                        }
                        for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                            expr -= cpRightSel[i][j];;
                        }
                    }
                    pricing.add(
                            IloRange(cplex, -IloInfinity, expr, 0, ("sym_brk_series1_" + std::to_string(k)).c_str()));
                    expr.clear();
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


            // Time
            int pos = 0;
            for (int i = 0; i < N_Dims; ++i) {
                for (int j = i + 1; j < N_Dims; ++j) {
                    pricing.add(IloRange(cplex, cpCombiAB[pos] - cpSelection[i] , 0));
                    pricing.add(IloRange(cplex, cpCombiAB[pos] - cpGroupBy[j] , 0));
                    pricing.add(IloRange(cplex, - 1 - cpCombiAB[pos] + cpGroupBy[j] + cpSelection[i], 0));

                    pricing.add(IloRange(cplex, cpCombiBA[pos] - cpSelection[j] , 0));
                    pricing.add(IloRange(cplex, cpCombiBA[pos] - cpGroupBy[i] , 0));
                    pricing.add(IloRange(cplex, - 1 - cpCombiBA[pos] + cpGroupBy[i] + cpSelection[j], 0));
                    pos ++;
                }
            }

            IloExpr time_expr(cplex);
            for (auto i = 0; i < rmpQSet.size(); ++i)
                time_expr += tap_s[i] * (int) rmpIST.time(i);
            for (int i = 0; i < (N_Dims * (N_Dims - 1))/2; ++i) {
                time_expr += cpCombiAB[i] * timing_flat[i];
                time_expr += cpCombiBA[i] * timing_flat[i];
            }
            pricing.add(IloRange(cplex, 0, time_expr , time_bound, "time_epsilon"));


            //Distance
            IloExpr distance_expr(cplex);
            for (auto i = 1; i <= rmpQSet.size(); ++i) {
                for (auto j = 1; j <= rmpQSet.size(); ++j) {
                    distance_expr += tap_x[i][j] * (int) rmpIST.dist(i - 1, j - 1);
                }
            }
            for (auto i = 0; i < rmpQSet.size() + 1; ++i)
                distance_expr += lin_D_in[i] + lin_D_out[i];
            pricing.add(IloRange(cplex, -IloInfinity, distance_expr, dist_bound, "distance_epsilon"));

            for (int i = 0; i < rmpQSet.size(); ++i) {
                for (int k = 0; k < N_Dims; ++k) {
                    expr += cpSelection[k] * (1 - S[i][k]) + S[i][k] * (1 - cpSelection[k]);
                }
                expr -= (1 - tap_x[i][rmpQSet.size()+1]) * HV3;
                expr -= lin_D_in[i];
                pricing.add(IloRange(cplex, -IloInfinity, expr, 0));
                expr.clear();
            }
            for (int i = 0; i < rmpQSet.size(); ++i) {
                for (int k = 0; k < N_Dims; ++k) {
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
                for (int i = 0; i < N_Dims; ++i) {
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
                for (int i = 0; i < N_Dims; ++i) {
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
                for (int i = 0; i < N_Dims; ++i) {
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
            if(confMap_[selectedConf][0] == 2) {
                for (auto q: rmpQSet) {
                    int var_cnt = 0;
                    // GB Key
                    for (int i = 0; i < N_Dims; ++i) {
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
                    for (int i = 0; i < N_Dims; ++i) {
                        bool dimPresent = false;
                        int valueIdx = -1;
                        for (auto p: q.getLeftPredicate()) {
                            if (p.first == pricingIST.getDimName(i)) {
                                dimPresent = true;
                                valueIdx = p.second;
                            }
                        }
                        if (dimPresent) {
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
                    for (int i = 0; i < N_Dims; ++i) {
                        bool dimPresent = false;
                        int valueIdx = -1;
                        for (auto p: q.getRightPredicate()) {
                            if (p.first == pricingIST.getDimName(i)) {
                                dimPresent = true;
                                valueIdx = p.second;
                            }
                        }
                        if (dimPresent) {
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
            }

            //Init solver
            IloCplex cplex_solver(pricing);
            cplex_solver.setParam(IloCplex::Param::TimeLimit, pricing_it_timeout);
            cplex_solver.setParam(IloCplex::Param::Threads, 1);
            if (confMap_[selectedConf][0] != 3) {
                cplex_solver.setParam(IloCplex::Param::Preprocessing::Symmetry, cplex_sym);
            }
            cplex_solver.setParam(IloCplex::Param::Preprocessing::Aggregator, confMap_[selectedConf][1]);
            cplex_solver.setParam(IloCplex::Param::Preprocessing::Presolve	, confMap_[selectedConf][1]);
            cplex_solver.setParam(IloCplex::Param::MIP::Display, 0); //5 is max
            cplex_solver.setParam(IloCplex::Param::Simplex::Display, 0);
            cplex_solver.setParam(IloCplex::Param::MIP::Strategy::HeuristicEffort, 0);
            cplex_solver.setParam(IloCplex::Param::Emphasis::MIP	, confMap_[selectedConf][3]);
            cplex_solver.setParam(IloCplex::Param::MIP::SubMIP::NodeLimit, confMap_[selectedConf][2]);
            //cplex_solver.setParam(IloCplex::Param::MIP::Tolerances::MIPGap	, 0.0);
            //cplex_solver.setParam(IloCplex::Param::Simplex::Tolerances::Optimality	, 10e-9);
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
                std::cout << "[Info] pricing timeout without feasible solution" << endl;
                break;
            }

            /*
             *  --- Recover query from solution ---
             */
            //IloNumArray solSelection(cplex);
            IloNumArray solGBKey(cplex);
            IloNumArray solLeftMeasure(cplex);
            IloNumArray solRightMeasure(cplex);
            IloArray<IloNumArray> solLeftSelection(cplex, N_Dims);
            IloArray<IloNumArray> solRightSelection(cplex, N_Dims);
            isNewQuerySelected = cplex_solver.getValue(tap_s[rmpQSet.size()]);

            //cplex_solver.getValues(solSelection, cpSelection);
            cplex_solver.getValues(solGBKey, cpGroupBy);
            cplex_solver.getValues(solLeftMeasure, cpLeftMeasure);
            cplex_solver.getValues(solRightMeasure, cpRightMeasure);
            for (int i = 0; i < N_Dims; ++i) {
                IloNumArray line(cplex);
                cplex_solver.getValues(line, cpLeftSel[i]);
                solLeftSelection[i] = line;
            }
            for (int i = 0; i < N_Dims; ++i) {
                IloNumArray line(cplex);
                cplex_solver.getValues(line, cpRightSel[i]);
                solRightSelection[i] = line;
            }

            int gbKeyIdx = 0;
            while (gbKeyIdx < N_Dims) {
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
            for (int i = 0; i < N_Dims; ++i) {
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


            if (!NO_PRINT) std::cout << "[Pricing Query] " << isNewQuerySelected << " - " << picked << endl;

            if (!NO_PRINT) std::cout << "[Pricing z value] " << cplex_solver.getObjValue() << endl;
            if (!NO_PRINT) std::cout << "[Pricing t value] " << cplex_solver.getValue(time_expr) << endl;
            if (!NO_PRINT) std::cout << "[Pricing d value] " << cplex_solver.getValue(distance_expr) << endl;

            time_t end = clock();
            double time_to_sol = (double)(end - start) / (double)CLOCKS_PER_SEC;
            if (!NO_PRINT) std::cout << "[TIME][ITER][s] " << time_to_sol << endl;
            rmpQSet.emplace_back(picked);

            cplex_solver.end();
            cplex.end();

            if (time_to_sol > global_timeout){
                std::cout << "[BREAK] Reason: global timeout" << endl;
                break;
            }

            if (!NO_PRINT){
                std::cout << "[KS Sol] ";
                KnapsackSolver ksolve(pricingIST);
                Solution ksol = ksolve.solve(rmpQSet, time_bound, dist_bound);
                for (int q : ksol.sequence) {
                    std::cout << rmpQSet[q] << " ";
                }
                std::cout << endl;
                std::cout << "[KS z value] " << ksol.z << endl;
            }
        }

        if (!NO_PRINT) {
            std::cout << "[OBJ]";
            for (auto it = objValues.begin(); it != objValues.end(); ++it) {
                std::cout << *it;
                if (it != objValues.end() - 1)
                    std::cout << std::string(";");

            }
            std::cout << endl;
        }
        std::cout << "[INFO] iterations " << objValues.size() << endl;
        std::cout << "[INFO][TIME] iterations " << (double)(::clock() - start) / (double)CLOCKS_PER_SEC << endl;


        Instance rmpIST = buildRMPInstance(rmpQSet);
        auto tapSolver = Solver(rmpIST);
        tapSolver.setTimeout(master_it_timeout);
        auto final_sol = tapSolver.solve(dist_bound, time_bound, false, "");
        std::cout << "[MASTER] " << final_sol.z << "|" << final_sol.optimal << endl;

        for ( auto qid : final_sol.sequence){
            //std::cout << rmpQSet[qid] << endl;
        }

        auto matheuristic = SolverVPLSHammingSX(rmpIST, 15, 15, 30, 20);
        auto mathsol = matheuristic.solve(dist_bound, time_bound, false, "");
        std::cout << "[MASTER][VPLS] " << mathsol.z << " | " << final_sol.optimal <<  " | " << mathsol.time << endl;

        return final_sol;

}


Instance pricingSolver::buildRMPInstance(vector<Query>& queries) const {
    vector<double> interest = JVMAdapter::getInterest(queries, pricingIST);
    vector<double> time = JVMAdapter::getTime(queries, pricingIST);
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

    void pricingSolver::setSelectedConf(const string &selectedConf) {
        if (confMap_.find(selectedConf) == confMap_.end()) {
            cout << "[WARNING] Configuration " << selectedConf << " not found. Using " << pricingSolver::selectedConf << " instead." << endl;
        } else {
            pricingSolver::selectedConf = selectedConf;
        }

    }

    void pricingSolver::setUseFloat(bool useFloat) {
        pricingSolver::useFloat = useFloat;
    }

    void pricingSolver::setGlobalTimeout(int globalTimeout) {
        global_timeout = globalTimeout;
    }
} // cplex_tap
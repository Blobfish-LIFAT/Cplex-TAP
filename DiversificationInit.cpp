//
// Created by alex on 10/11/22.
//

#include "DiversificationInit.h"
#include <utility>
#include <ilcplex/ilocplex.h>


namespace cplex_tap {

    vector<Query> DiversificationInit::build(int setSize) {
        while (baseSet.size() <= setSize){
            // Init CPLEX environment and model objects
            IloEnv cplex;
            IloModel initProblem(cplex);

            // Selection from queries
            vector<vector<bool>> S;
            for (int i = 0; i < baseSet.size(); ++i) {
                vector<bool> tmp;
                for (int j = 0; j < pricingIST.getNbDims(); ++j) {
                    bool appears = false;
                    const string dimName = pricingIST.getDimName(j);
                    for (int k = 0; k < baseSet[i].getLeftPredicate().size(); ++k) {
                        appears |= baseSet[i].getLeftPredicate()[k].first == dimName;
                    }
                    for (int k = 0; k < baseSet[i].getRightPredicate().size(); ++k) {
                        appears |= baseSet[i].getRightPredicate()[k].first == dimName;
                    }
                    tmp.emplace_back(appears);
                }
                S.emplace_back(tmp);
            }
            vector<int> time = JVMAdapter::getTime(baseSet, pricingIST);
            vector<double> interest = JVMAdapter::getInterest(baseSet, pricingIST);

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

            //
            // --- Objective ---
            //

            // Majorants
            int N_t = *max_element(std::begin(time), std::end(time));
            N_t = max(N_t, std::accumulate(pricingIST.getDimTimes().begin(), pricingIST.getDimTimes().end(), 0));
            N_t = 2 * (N_t + 1);
            int N_i = *max_element(std::begin(interest), std::end(interest));
            N_i = max(N_i, std::accumulate(pricingIST.getDimWeights().begin(), pricingIST.getDimWeights().end(), 0));
            N_i = 2 * (N_i + 1);

            IloExpr expr(cplex);
            for (int i = 0; i < baseSet.size(); ++i) {
                IloExpr sumT(cplex);
                IloExpr sumI(cplex);
                for (int j = 0; j < pricingIST.getNbDims(); ++j) {
                    expr += cpSelection[j]*(1 - S[i][j]) + (1 - cpSelection[j]) * S[i][j];
                    sumT += pricingIST.getDimTime(j) * cpSelection[j];
                    sumI += pricingIST.getDimWeight(j) * cpSelection[j];
                    // Square thing
                    //expr += (pricingIST.getDimWeight(j) * cpSelection[j] - S[i][j] * pricingIST.getDimWeight(j))*(pricingIST.getDimWeight(j) * cpSelection[j] - S[i][j] * pricingIST.getDimWeight(j));
                    //expr += (pricingIST.getDimTime(j) * cpSelection[j] - S[i][j] * pricingIST.getDimTime(j))*(pricingIST.getDimTime(j) * cpSelection[j] - S[i][j] * pricingIST.getDimTime(j));
                }
                // linearization of abs value terms
                IloNumVar Z_t(cplex, ("Zt_" + to_string(i)).c_str());
                IloNumVar Y_t(cplex, 0, 1, IloNumVar::Bool, ("Y_" + to_string(i)).c_str());
                initProblem.add(IloRange(cplex, 0, sumT - time[i] + N_t * Y_t - Z_t, IloInfinity));
                initProblem.add(IloRange(cplex, 0, - sumT + time[i] + N_t * (1 - Y_t) - Z_t, IloInfinity));
                initProblem.add(IloRange(cplex, sumT - time[i] - Z_t, 0));
                initProblem.add(IloRange(cplex, - sumT + time[i] - Z_t, 0));

                IloNumVar Z_i(cplex, ("Zt_" + to_string(i)).c_str());
                IloNumVar Y_i(cplex, 0, 1, IloNumVar::Bool, ("Y_" + to_string(i)).c_str());
                initProblem.add(IloRange(cplex, 0, sumI - interest[i] + N_i * Y_i - Z_i, IloInfinity));
                initProblem.add(IloRange(cplex, 0, - sumI + interest[i] + N_i * (1 - Y_i) - Z_i, IloInfinity));
                initProblem.add(IloRange(cplex, sumI - interest[i] - Z_i, 0));
                initProblem.add(IloRange(cplex, - sumI + interest[i] - Z_i, 0));

                expr += Z_t;
                expr += Z_i;
            }
            IloObjective obj(cplex, expr, IloObjective::Maximize);
            initProblem.add(obj);
            expr.clear();


            /*
             *  --- Pricing Constraints ---
             */
            // One measure per side
            for (int i = 0; i < pricingIST.getNbMeasures(); ++i)
                expr += cpLeftMeasure[i];
            initProblem.add(IloRange(cplex, 1, expr, 1, "left_one_measure"));
            expr.clear();
            for (int i = 0; i < pricingIST.getNbMeasures(); ++i)
                expr += cpRightMeasure[i];
            initProblem.add(IloRange(cplex, 1, expr, 1, "right_one_measure"));
            expr.clear();
            // One group key
            for (int i = 0; i < pricingIST.getNbDims(); ++i)
                expr += cpGroupBy[i];
            initProblem.add(IloRange(cplex, 1, expr, 1, "one_gk"));
            expr.clear();
            // One for selection
            for (int i = 0; i < pricingIST.getNbDims(); ++i)
                expr += cpSelection[i];
            initProblem.add(IloRange(cplex, 1, expr, 1, "one_more_selection"));
            expr.clear();
            // No overlap selection / group by
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                initProblem.add(
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
                initProblem.add(IloRange(cplex, 0, lexp - cpSelection[i], 0));
                initProblem.add(IloRange(cplex, 0, rexp - cpSelection[i], 0));
            }
            for (int i = 0; i < pricingIST.getNbDims(); ++i) {
                for (int j = 0; j < pricingIST.getAdSize(i); ++j) {
                    initProblem.add(IloRange(cplex, 0, cpLeftSel[i][j] + cpRightSel[i][j], 1)); // (12)
                    // (13)
                    initProblem.add(IloRange(cplex, 0, ((cpLeftSel[i][j] + cpRightSel[i][j]) / 2.0) - cpSelExclusive[i][j],
                                             IloInfinity));
                    expr += cpSelection[i];
                    for (int k = 0; k < pricingIST.getAdSize(i); ++k) {
                        if (k != j) {
                            expr -= cpLeftSel[i][k];
                            expr -= cpRightSel[i][k];
                        }
                    }
                    expr -= cpSelExclusive[i][j];
                    initProblem.add(IloRange(cplex, -IloInfinity, expr, 0));
                    expr.clear();
                    // (14)
                    for (int k = 0; k < pricingIST.getNbMeasures(); ++k) {
                        initProblem.add(IloRange(cplex, 0, 2 - cpSelExclusive[i][j] - cpLeftMeasure[k] - cpRightMeasure[k],
                                                 IloInfinity));
                    }
                }
            }

            /*
            *  --- Constraints forbidding having same query as existing one ---
            */
            for (auto q : baseSet) {
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

                initProblem.add(IloRange(cplex, 0, expr, var_cnt - 1));
                expr.clear();
            }

            //Init solver
            IloCplex cplex_solver(initProblem);
            cplex_solver.setParam(IloCplex::Param::TimeLimit, 300);
            cplex_solver.setParam(IloCplex::Param::Threads, 1);
            cplex_solver.setParam(IloCplex::Param::MIP::Display, 0);
            cplex_solver.setParam(IloCplex::Param::Simplex::Display, 0);
            //cplex_solver.setParam(IloCplex::Param::Preprocessing::QToLin, 0);
            if (!debug) cplex_solver.setOut(cplex.getNullStream());

            bool solved = false;
            try {
                solved = cplex_solver.solve();
            }
            catch (const IloException &e) {
                std::cerr << "\n\n--- CPLEX Exception (Init) ---\n";
                std::cerr << e << "\n";
                cplex.end();
                throw;
            }
            if (solved) {
                if (debug) {
                    std::cout << "\n--- Solver success (Init) ---\n";
                    std::cout << "    Status: " << cplex_solver.getStatus() << "\n";
                    std::cout << "    Objective: " << cplex_solver.getObjValue() << "\n";
                }
            } else {
                std::cerr << "\n--- Solver Error (Init) ---\n";
                std::cerr << "    Status: " << cplex_solver.getStatus() << "\n";
                std::cerr << "    Error details: " << cplex_solver.getCplexStatus() << "\n";
            }

            /*
            *  --- Recover query from solution ---
            */
            IloNumArray solGBKey(cplex, ILOINT);
            IloNumArray solLeftMeasure(cplex);
            IloNumArray solRightMeasure(cplex);
            IloArray<IloNumArray> solLeftSelection(cplex, pricingIST.getNbDims());
            IloArray<IloNumArray> solRightSelection(cplex, pricingIST.getNbDims());

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
            baseSet.emplace_back(picked);
        }
        return baseSet;
    }
} // cplex_tap
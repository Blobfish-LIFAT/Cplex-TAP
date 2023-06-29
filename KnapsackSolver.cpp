//
// Created by alex on 26/01/23.
//

#include <algorithm>
#include "KnapsackSolver.h"
#include "JVMAdapter.h"

namespace cplex_tap {
    Solution KnapsackSolver::solve(std::vector<Query> queries, int timeBudget, int maxDistance) {
        time_t start, end;
        start = clock();

        std::cout << "[INFO] KS Heuristic : Init |Q|=" << queries.size() << std::endl;
        int size = queries.size();

        std::vector<int> solution;
        std::vector<std::pair<int,double>> order;
        vector<double> interests = JVMAdapter::getInterest(queries, (ist));
        vector<double> times = JVMAdapter::getTime(queries, (ist));
        for (int i = 0; i < size; i++) {
            order.emplace_back(std::make_pair(i, interests[i]/times[i]) );
        }
        //Sort by ratio descending
        //cout << "[INFO] KS Heuristic : begin sort" << endl;
        sort(order.begin(), order.end(), sortbysec_rev);
        //cout << "[INFO] KS Heuristic : sort done" << endl;

        /*for (int i = 0; i < 100; ++i) {
            cout << order[i].first << "|" << order[i].second << " ";
        }
        cout << endl;*/

        double total_dist = 0;
        double total_time = 0;
        double z = 0;


        //std::cout << "[INFO] KS Heuristic : Construction Solution" << std::endl;
        for (int i = 0; i < size; i++)
        {
            int current = order[i].first;

            if (timeBudget - (total_time + times[current]) >= 0){
                double backup = total_dist;
                total_dist += insert_opt(&solution, current, &queries);
                if (total_dist > maxDistance){
                    //rollback and check next querry
                    auto removed = remove_if(solution.begin(), solution.end(), [&](int el) -> bool {return el == current;});
                    std::vector<int> rollback_sol;
                    for (auto q = solution.begin(); q != removed; ++q)
                        rollback_sol.emplace_back(*q);
                    solution = rollback_sol;
                    total_dist = backup;
                } else {
                    total_time += times[current];
                    z += interests[current];
                }
            }
        }
        //std::cout << "[INFO] KS Heuristic : eptime=" << total_time << "/" << timeBudget << std::endl;
        //std::cout << "[INFO] KS Heuristic : epdist=" << total_dist << "/" << maxDistance << std::endl;
        //std::cout << "[INFO] KS Heuristic : solution is done z=" << z << std::endl;

        end = clock();
        double time = (double)(end - start) / (double)CLOCKS_PER_SEC;
        return {false, time, z, solution};
    }

    double KnapsackSolver::insert_opt(std::vector<int> *solution, int candidate, std::vector<Query> *queries) {
        if (solution->empty()){
            solution->emplace_back(candidate);
            return 0;
        }
        double best_insert_cost = 10e50;// large enough
        int best_insert_pos = -1;
        Query candidateQuery = queries->at(candidate);
        for (int i = 0; i < solution->size() + 1; i++) {
            double new_cost = 0;
            // insert at first position
            if (i == 0){
                new_cost += candidateQuery.dist(queries->at(solution->at(0)));
            } else if (i < solution->size()){
                int current_querry = solution->at(i);
                new_cost += candidateQuery.dist(queries->at(solution->at(i-1)));
                new_cost += candidateQuery.dist(queries->at(current_querry));
                new_cost -= queries->at(solution->at(i-1)).dist(queries->at(current_querry));
            } else {
                new_cost += queries->at(solution->at(solution->size()-1)).dist(candidateQuery);
            }
            if (new_cost < best_insert_cost){
                best_insert_cost = new_cost;
                best_insert_pos = i;
            }
        }
        solution->insert(solution->begin() + best_insert_pos, candidate);
        return best_insert_cost;
    }


} // cplex_tap
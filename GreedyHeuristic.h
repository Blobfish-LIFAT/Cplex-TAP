//
// Created by chanson on 5/19/2021.
//

#ifndef CPLEX_TEST_GREEDYHEURISTIC_H
#define CPLEX_TEST_GREEDYHEURISTIC_H

#include "instance.h"
#include "Element.h"

namespace cplex_tap {

    class GreedyHeuristic {
    public:
        std::vector<int> getBasicSolution(Instance &ist, double eptime, double epdist) {
            std::vector<int> solution;

            vector<Element> order = new ArrayList<>();
            for (int i = 0; i < ist.size; i++) {
                order.push_back(Element(i, ist.interest(i)/ist.time(i)));
            }
            std::sort(order.begin(), order.end());
            std::reverse(order.begin(), order.end());

            double total_dist = 0;
            double total_time = 0;
            double z = 0;

            for (int i = 0; i < ist.size; i++)
            {
                int current = order.at(i).index;

                if (eptime - (total_time + ist.costs[current]) > 0){
                    double backup = total_dist;
                    total_dist += insert_opt(solution, current, ist.distances, total_dist);
                    if (total_dist > epdist){
                        //rollback and check next querry
                        solution.remove(Integer.valueOf(current));
                        total_dist = backup;
                        continue;
                    }
                    total_time += ist.costs[current];


                    z += ist.interest[current];
                }
            }*/
            return solution;
        }
    protected:
    double insert_opt(std::vector<int> solution, int candidate, std::vector<std::vector<std::uint32_t>> distances, double base_dist) {
            if (solution.empty()){
                solution.push_back(candidate);
                return 0;
            }
            double best_insert_cost = 10e50;// large enough
            int best_insert_pos = -1;
            for (int i = 0; i < solution.size() + 1; i++) {
                double new_cost = 0;
                // insert at first position
                if (i == 0){
                    new_cost += distances[candidate][solution.at(0)];
                } else if (i < solution.size()){
                    int current_querry = solution.at(i);
                    new_cost += distances[candidate][solution.at(i-1)];
                    new_cost += distances[candidate][current_querry];
                    new_cost -= distances[solution.at(i-1)][current_querry];
                } else {
                    new_cost += distances[solution.at(solution.size()-1)][candidate];
                }
                if (new_cost < best_insert_cost){
                    best_insert_cost = new_cost;
                    best_insert_pos = i;
                }
            }
        if (best_insert_pos != solution.size() - 1) {
            auto it = solution.begin();
            solution.insert(it + best_insert_pos + 1, candidate);
        }
        else
            solution.push_back(candidate);
        return best_insert_cost;
        }
    };

}
#endif //CPLEX_TEST_GREEDYHEURISTIC_H

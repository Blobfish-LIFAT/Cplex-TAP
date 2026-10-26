//
// Created by achan on 26/10/2026.
//

#ifndef CPLEX_TEST_KSWARMSTART_H
#define CPLEX_TEST_KSWARMSTART_H

#include <cstdint>
#include <vector>
#include <algorithm>
#include <ctime>
#include <limits>

#include "instance.h"
#include "Solution.h"
#include "Element.h"

namespace cplex_tap {

class KSWarmStart
{
protected:
    const Instance& tap;

public:
    explicit KSWarmStart(const Instance& tap) : tap(tap) {}

    // Evaluates optimal insertion position for candidate into the current solution
    static double insert_opt(std::vector<int> &solution, int candidate, const Instance &tap) {
        if (solution.empty()) {
            solution.push_back(candidate);
            return 0.0;
        }
        double best_insert_cost = 1e50;
        int best_insert_pos = -1;
        const int n = static_cast<int>(solution.size());

        for (int i = 0; i <= n; ++i) {
            double new_cost = 0.0;
            if (i == 0) {
                new_cost += tap.dist(candidate, solution[0]);
            } else if (i < n) {
                int current_query = solution[i];
                new_cost += tap.dist(candidate, solution[i - 1]);
                new_cost += tap.dist(candidate, current_query);
                new_cost -= tap.dist(solution[i - 1], current_query);
            } else {
                new_cost += tap.dist(solution[n - 1], candidate);
            }

            if (new_cost < best_insert_cost) {
                best_insert_cost = new_cost;
                best_insert_pos = i;
            }
        }

        solution.insert(solution.begin() + best_insert_pos, candidate);
        return best_insert_cost;
    }

    // Computes the total distance along the sequence of queries
    static double sequenceDistance(const std::vector<int> &tour, const Instance &tap) {
        if (tour.size() < 2) return 0.0;
        double d = 0.0;
        for (size_t i = 0; i < tour.size() - 1; ++i) {
            d += tap.dist(tour[i], tour[i + 1]);
        }
        return d;
    }

    // Returns the query in sol with the lowest interest
    static int argMinInt(const std::vector<int> &sol, const Instance &tap) {
        if (sol.empty()) return -1;
        int victim = sol.front();
        double min_val = tap.interest(victim);
        for (int q : sol) {
            double val = tap.interest(q);
            if (val < min_val) {
                min_val = val;
                victim = q;
            }
        }
        return victim;
    }

    // Solves TAP using Knapsack-style greedy heuristic with 1-opt post-processing
    static Solution solve(const Instance &tap, int dist_bound, int time_bound) {
        clock_t start = clock();

        std::vector<int> solution;
        std::vector<Element> order;
        order.reserve(tap.size());

        for (uint32_t i = 0; i < tap.size(); ++i) {
            double cost = tap.time(i);
            double val = (cost > 0) ? (tap.interest(i) / cost) : std::numeric_limits<double>::infinity();
            order.emplace_back(static_cast<int>(i), val);
        }

        std::stable_sort(order.begin(), order.end(), [](const Element &a, const Element &b) {
            return a.value > b.value;
        });

        double total_dist = 0.0;
        double total_time = 0.0;
        double z = 0.0;

        for (size_t i = 0; i < order.size(); ++i) {
            int current = order[i].index;

            if (total_time + tap.time(current) <= static_cast<double>(time_bound)) {
                double backup = total_dist;
                total_dist += insert_opt(solution, current, tap);
                if (total_dist > dist_bound) {
                    // Rollback and check next query
                    solution.erase(std::remove(solution.begin(), solution.end(), current), solution.end());
                    total_dist = backup;
                    continue;
                }
                total_time += tap.time(current);
                z += tap.interest(current);
            }
        }

        // Post-processing: try removing the least interesting query (victim)
        // and replacing it with higher-interest candidates
        while (!solution.empty()) {
            std::vector<int> prevSolution = solution;
            int victim = argMinInt(solution, tap);
            solution.erase(std::remove(solution.begin(), solution.end(), victim), solution.end());
            total_dist = sequenceDistance(solution, tap);
            total_time -= tap.time(victim);
            z -= tap.interest(victim);

            std::vector<Element> potential;
            for (uint32_t i = 0; i < tap.size(); ++i) {
                int query_idx = static_cast<int>(i);
                if (query_idx != victim &&
                    std::find(solution.begin(), solution.end(), query_idx) == solution.end() &&
                    tap.interest(victim) < tap.interest(query_idx)) {
                    potential.emplace_back(query_idx, tap.interest(query_idx));
                }
            }

            std::stable_sort(potential.begin(), potential.end(), [](const Element &a, const Element &b) {
                return a.value > b.value;
            });

            bool success = false;
            for (const auto &e : potential) {
                if (total_time + tap.time(e.index) <= static_cast<double>(time_bound)) {
                    double backup = total_dist;
                    total_dist += insert_opt(solution, e.index, tap);
                    if (total_dist > dist_bound) {
                        // Rollback and check next query
                        solution.erase(std::remove(solution.begin(), solution.end(), e.index), solution.end());
                        total_dist = backup;
                        continue;
                    }
                    total_time += tap.time(e.index);
                    z += tap.interest(e.index);
                    success = true;
                    break;
                }
            }

            if (!success) {
                solution = prevSolution;
                z = 0.0;
                for (int idx : solution) {
                    z += tap.interest(idx);
                }
                break;
            }
        }

        clock_t end = clock();
        double time_to_sol = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        return Solution(time_to_sol, z, solution);
    }

    Solution solve(int dist_bound, int time_bound) const {
        return solve(tap, dist_bound, time_bound);
    }
};

} // namespace cplex_tap

using KSWarmStart = cplex_tap::KSWarmStart;

#endif //CPLEX_TEST_KSWARMSTART_H

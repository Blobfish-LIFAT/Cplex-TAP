//
// Created by alex on 30/01/23.
//

#ifndef CPLEX_TEST_KMPPINIT_H
#define CPLEX_TEST_KMPPINIT_H

#include "CGTAPInstance.h"
#include "Query.h"
#include "RandomInit.h"
#include "KnapsackSolver.h"

namespace cplex_tap {

    class KMppInit {
    protected:
        CGTAPInstance pricingIST;
        long random_set_size = 10000;
        int ep_dist, ep_time;
    public:
        KMppInit(const CGTAPInstance &pricingIst, int epDist, int epTime) : pricingIST(pricingIst), ep_dist(epDist),
                                                                          ep_time(epTime) {}

        std::vector<Query> build(int size){
            RandomInit rdini = RandomInit(pricingIST);
            vector<Query> rands = rdini.build(random_set_size);

            vector<Query> out;

            while (out.size() < size/3) {
                KnapsackSolver ks = KnapsackSolver(&pricingIST);
                Solution sol = ks.solve(rands, ep_time, ep_dist);
                vector<Query> sol_q;
                for (int i = 0; i < sol.sequence.size(); ++i) {
                    sol_q.emplace_back(rands[sol.sequence[i]]);
                }
                for (int i = 0; i < sol.sequence.size() && out.size() < size/3; ++i) {
                    out.emplace_back(rands[sol.sequence[i]]);
                }

                auto itr = std::remove_if(rands.begin(), rands.end(),
                                          [&](Query q) { return std::find(sol_q.begin(), sol_q.end(), q) != sol_q.end(); });
                rands.erase(itr, rands.end());
            }


            // KMeans ++ style diversification
            std::random_device rd;
            std::mt19937 gen(rd());

            while (out.size() < size){
                vector<double> weights;

                for (int i = 0; i < rands.size(); ++i) {
                    // For each remaining query find the distance to the nearest query already selected
                    double closest_dist = rands[i].dist(out[0]);
                    for (int j = 1; j < out.size(); ++j) {
                        double dist = rands[i].dist(out[j]);
                        if (dist < closest_dist)
                            closest_dist = dist;
                    }
                    weights.emplace_back(closest_dist);
                }

                std::discrete_distribution<> rdQuery(weights.begin(), weights.end());
                int idx = rdQuery(gen);

                out.emplace_back(rands[idx]);
                rands.erase( rands.begin() + idx);
            }

            return out;
        }

    };

} // cplex_tap

#endif //CPLEX_TEST_KMPPINIT_H

//
// Created by alex on 30/01/23.
//

#ifndef CPLEX_TEST_KSINIT_H
#define CPLEX_TEST_KSINIT_H

#include "CGTAPInstance.h"
#include "Query.h"
#include "RandomInit.h"
#include "KnapsackSolver.h"

namespace cplex_tap {

    class KSInit {
    protected:
        CGTAPInstance pricingIST;
        long random_set_size = 50000;
        int ep_dist, ep_time;
    public:
        KSInit(const CGTAPInstance &pricingIst, int epDist, int epTime) : pricingIST(pricingIst), ep_dist(epDist),
                                                                          ep_time(epTime) {}

    std::vector<Query> build(int size){
        RandomInit rdini = RandomInit(pricingIST);
        vector<Query> rands = rdini.build(random_set_size);

        vector<Query> out;

        while (out.size() < size) {
            KnapsackSolver ks = KnapsackSolver(&pricingIST);
            Solution sol = ks.solve(rands, ep_time, ep_dist);
            vector<Query> sol_q;
            for (int i = 0; i < sol.sequence.size(); ++i) {
                sol_q.emplace_back(rands[sol.sequence[i]]);
            }
            for (int i = 0; i < sol.sequence.size() && out.size() < size; ++i) {
                out.emplace_back(rands[sol.sequence[i]]);
            }

            auto itr = std::remove_if(rands.begin(), rands.end(),
                                      [&](Query q) { return std::find(sol_q.begin(), sol_q.end(), q) != sol_q.end(); });
            rands.erase(itr, rands.end());
        }

        return out;
        }
    };

} // cplex_tap

#endif //CPLEX_TEST_KSINIT_H

//
// Created by alex on 21/07/22.
//

#include "instance.h"
#include "CGTAPInstance.h"
#include <iostream>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include "Solution.h"
#include "Query.h"

#ifndef CPLEX_TEST_PRICINGSOLVER_H
#define CPLEX_TEST_PRICINGSOLVER_H

namespace cplex_tap {

    class pricingSolver {
    protected:
        // debug prints ?
        bool debug;
        // The TAP restricted master problem instance
        //const Instance& rmpIST;
        // The TAP pricing problem instance
        const CGTAPInstance &pricingIST;

        //ep constraints
        int dist_bound;
        int time_bound;

    public:
        // Builds a solver with the specified instance
        explicit pricingSolver(const CGTAPInstance &ist, int dist_bound, int time_bound) : debug{false},
                                                                                           pricingIST{ist},
                                                                                           dist_bound{dist_bound},
                                                                                           time_bound{time_bound} {}

        explicit pricingSolver(const CGTAPInstance &ist, int dist_bound, int time_bound, bool debug) : debug(debug),
                                                                                                       pricingIST{ist},
                                                                                                       dist_bound{dist_bound},
                                                                                                       time_bound{time_bound} {}

        // Run solver
        virtual Solution solve() const;

        Instance buildRMPInstance(vector<Query> queries) const;
    };

} // cplex_tap

#endif //CPLEX_TEST_PRICINGSOLVER_H

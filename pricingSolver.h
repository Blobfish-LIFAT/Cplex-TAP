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

        // starting queries
        vector<Query> extStarting;
        int randomStartSize = 0;

    public:
        // Builds a solver with the specified instance
        explicit pricingSolver(const CGTAPInstance &ist, int dist_bound, int time_bound) : debug{false},
                                                                                           pricingIST{ist},
                                                                                           dist_bound{dist_bound},
                                                                                           time_bound{time_bound},
                                                                                           extStarting(vector<Query>()){}

        explicit pricingSolver(const CGTAPInstance &ist, int dist_bound, int time_bound, int rdStartSize) : debug{false},
                                                                                           pricingIST{ist},
                                                                                           dist_bound{dist_bound},
                                                                                           time_bound{time_bound},
                                                                                           extStarting(vector<Query>()),
                                                                                           randomStartSize(rdStartSize){}

        explicit pricingSolver(const CGTAPInstance &ist, int dist_bound, int time_bound, vector<Query> initSet, bool debug) :
                                                                                                       debug(debug),
                                                                                                       pricingIST{ist},
                                                                                                       dist_bound{dist_bound},
                                                                                                       time_bound{time_bound},
                                                                                                       extStarting(initSet){}

        explicit pricingSolver(const CGTAPInstance &ist, int dist_bound, int time_bound, vector<Query> initSet) :
                debug(false),
                pricingIST{ist},
                dist_bound{dist_bound},
                time_bound{time_bound},
                extStarting(initSet){}

        // Run solver
        virtual Solution solve() const;

        Instance buildRMPInstance(vector<Query> queries) const;

        static bool assessConvergence(vector<double> objValues);
    };

} // cplex_tap

#endif //CPLEX_TEST_PRICINGSOLVER_H

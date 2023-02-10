//
// Created by alex on 08/02/23.
//

#ifndef CPLEX_TEST_PRINCINGCPSOLVER_H
#define CPLEX_TEST_PRINCINGCPSOLVER_H

#include "Query.h"
#include "CGTAPInstance.h"
#include "Solution.h"
#include "instance.h"

namespace cplex_tap {

    class princingCPSolver {
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

        //time
        int pricing_it_timeout = 10;
        int master_it_timeout = 600;
        int global_timeout = 3600;

    public:
        // Builds a solver with the specified instance
        explicit princingCPSolver(const CGTAPInstance &ist, int dist_bound, int time_bound, vector<Query> initSet, bool debug) :
                debug(debug),
                pricingIST{ist},
                dist_bound{dist_bound},
                time_bound{time_bound},
                extStarting(initSet){}

        explicit princingCPSolver(const CGTAPInstance &ist, int dist_bound, int time_bound, vector<Query> initSet) :
                debug(false),
                pricingIST{ist},
                dist_bound{dist_bound},
                time_bound{time_bound},
                extStarting(initSet){}

        // Run solver
        virtual Solution solve() const;

        cplex_tap::Instance buildRMPInstance(vector<Query>& queries) const;
    };

} // cplex_tap

#endif //CPLEX_TEST_PRINCINGCPSOLVER_H

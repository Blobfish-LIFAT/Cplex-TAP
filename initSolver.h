//
// Created by alex on 10/11/22.
//

#ifndef CPLEX_TEST_INITSOLVER_H
#define CPLEX_TEST_INITSOLVER_H

#include "CGTAPInstance.h"
#include "Query.h"
#include "JVMAdapter.h"


namespace cplex_tap {

    class initSolver {
    protected:
        const CGTAPInstance &pricingIST;
        const int setSize;
        const bool debug;
    public:
        initSolver(const CGTAPInstance &pricingIst, int setSize);
        initSolver(const CGTAPInstance &pricingIst, int setSize, bool debug);
        vector<Query> build_starting_set();
    };

} // cplex_tap

#endif //CPLEX_TEST_INITSOLVER_H

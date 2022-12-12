//
// Created by alex on 12/12/22.
//

#ifndef CPLEX_TEST_INTENSIFICATIONINIT_H
#define CPLEX_TEST_INTENSIFICATIONINIT_H

#include "CGTAPInstance.h"
#include "Query.h"

namespace cplex_tap {

    class IntensificationInit {
    protected:
        const CGTAPInstance &pricingIST;
        const bool debug;
    public:
        IntensificationInit(const CGTAPInstance &pricingIst);
        IntensificationInit(const CGTAPInstance &pricingIst, bool debug);
        vector<Query> build(int setSize);
    };

} // cplex_tap

#endif //CPLEX_TEST_INTENSIFICATIONINIT_H

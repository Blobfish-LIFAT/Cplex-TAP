//
// Created by alex on 12/12/22.
//

#ifndef CPLEX_TEST_INTENSIFICATIONINIT_H
#define CPLEX_TEST_INTENSIFICATIONINIT_H

#include "CGTAPInstance.h"
#include "Query.h"
#include "JVMAdapter.h"

namespace cplex_tap {

    class IntensificationInit {
    protected:
        const CGTAPInstance &pricingIST;
        const bool debug;
        std::vector<Query> baseSet;
    public:
        IntensificationInit(const CGTAPInstance &pricingIst, std::vector<Query> baseSet) : pricingIST(pricingIst), debug(false) {
            this->baseSet = vector<Query>(std::move(baseSet));
        };
        IntensificationInit(const CGTAPInstance &pricingIst, std::vector<Query> baseSet, bool debug) : pricingIST(pricingIst), debug(debug) {
            this->baseSet = vector<Query>(std::move(baseSet));
        };
        vector<Query> build(int setSize);
    };

} // cplex_tap

#endif //CPLEX_TEST_INTENSIFICATIONINIT_H

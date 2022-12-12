//
// Created by alex on 10/11/22.
//

#ifndef CPLEX_TEST_DIVERSIFICATIONINIT_H
#define CPLEX_TEST_DIVERSIFICATIONINIT_H

#include "CGTAPInstance.h"
#include "Query.h"
#include "JVMAdapter.h"


namespace cplex_tap {

    class DiversificationInit {
    protected:
        const CGTAPInstance &pricingIST;
        const bool debug;
        std::vector<Query> baseSet;
    public:
        DiversificationInit(const CGTAPInstance &pricingIst, std::vector<Query> baseSet) : pricingIST(pricingIst), debug(false) {
            this->baseSet = vector<Query>(std::move(baseSet));
        };
        DiversificationInit(const CGTAPInstance &pricingIst, std::vector<Query> baseSet, bool debug) : pricingIST(pricingIst), debug(debug) {
            this->baseSet = vector<Query>(std::move(baseSet));
        };
        vector<Query> build(int setSize);
    };

} // cplex_tap

#endif //CPLEX_TEST_DIVERSIFICATIONINIT_H

//
// Created by alex on 19/01/23.
//

#ifndef CPLEX_TEST_RANDTHENSORTINIT_H
#define CPLEX_TEST_RANDTHENSORTINIT_H

#include "CGTAPInstance.h"
#include "Query.h"
#include "JVMAdapter.h"
#include "RandomInit.h"

namespace cplex_tap {

    class randThenSortInit {
    protected:
        const CGTAPInstance &pricingIST;
        const bool debug;
        const int rand_per_selected = 500; // how many queries to generate for 1 query returned
        static bool sortbysec_rev(const pair<int,double> &a, const pair<int,double> &b){
            return (a.second > b.second);
        }
    public:
        randThenSortInit(const CGTAPInstance &pricingIst) : pricingIST(pricingIst), debug(false) {};
        randThenSortInit(const CGTAPInstance &pricingIst, bool debug) : pricingIST(pricingIst), debug(debug) {};
        vector<Query> build(int setSize);

    };

} // cplex_tap

#endif //CPLEX_TEST_RANDTHENSORTINIT_H

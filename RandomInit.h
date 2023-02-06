//
// Created by alex on 12/12/22.
//

#ifndef CPLEX_TEST_RANDOMINIT_H
#define CPLEX_TEST_RANDOMINIT_H

#include <vector>
#include <random>
#include "Query.h"
#include "CGTAPInstance.h"

namespace cplex_tap {

    class RandomInit {
    protected:
        CGTAPInstance pricingIST;
        long seed = -1;
    public:
        explicit RandomInit(const CGTAPInstance &pricingIst) : pricingIST(pricingIst) {}

        RandomInit(const CGTAPInstance &pricingIst, long seed) : pricingIST(pricingIst), seed(seed) {}

    public: std::vector<Query> build(int size){
            vector<Query> rmpQSet;
            std::random_device rd;
                std::mt19937 gen(seed != -1 ? seed : rd());
                std::uniform_int_distribution<> rdAttr(0, pricingIST.getNbDims() - 1);
                for (; rmpQSet.size() < size;) {
                    int lAttrID = rdAttr(gen);
                    //int rAttrID = rdAttr(gen);
                    int measureID = 0;
                    int gbAttr = rdAttr(gen);
                    while (gbAttr == lAttrID) {// || gbAttr == rAttrID
                        gbAttr = rdAttr(gen);
                    }
                    std::uniform_int_distribution<> rdValLeft(0, pricingIST.getAdSize(lAttrID) - 1);
                    //std::uniform_int_distribution<> rdValRight(0, pricingIST.getAdSize(rAttrID)-1);
                    int laval = rdValLeft(gen);
                    int rval = rdValLeft(gen);
                    while (rval == laval){
                        rval = rdValLeft(gen);
                    }
                    std::vector<std::pair<string, int> > lPredicate = {{pricingIST.getDimName(lAttrID), laval}};
                    std::vector<std::pair<string, int> > rPredicate = {{pricingIST.getDimName(lAttrID), rval}};
                    //std::vector<std::pair<string, int> > rPredicate = { {pricingIST.getDimName(rAttrID), rdValRight(gen)}};
                    Query rdQ = Query(pricingIST.getTableName(), "sum", pricingIST.getDimName(gbAttr),
                                      pricingIST.getMeasureName(measureID), pricingIST.getMeasureName(measureID),
                                      lPredicate, rPredicate);

                    bool duplicate = false;
                    for (Query &q : rmpQSet)
                        duplicate |= q == rdQ;
                    if (!duplicate){
                        rmpQSet.emplace_back(rdQ);
                    }
                }
            return rmpQSet;
    }

    };

} // cplex_tap

#endif //CPLEX_TEST_RANDOMINIT_H

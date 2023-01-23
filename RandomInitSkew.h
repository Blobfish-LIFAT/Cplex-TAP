//
// Created by alex on 23/01/23.
//

#ifndef CPLEX_TEST_RANDOMINITSKEW_H
#define CPLEX_TEST_RANDOMINITSKEW_H

#include <vector>
#include <random>
#include "Query.h"
#include "CGTAPInstance.h"

namespace cplex_tap {

    class RandomInitSkew {
    protected:
        CGTAPInstance pricingIST;
        long seed = -1;
    public:
        explicit RandomInitSkew(const CGTAPInstance &pricingIst) : pricingIST(pricingIst) {}

        RandomInitSkew(const CGTAPInstance &pricingIst, long seed) : pricingIST(pricingIst), seed(seed), was_seeded(true) {}

    protected:
        bool was_seeded = false;

    public: std::vector<Query> build(int size){
            vector<Query> rmpQSet;
            std::random_device rd;
            std::mt19937 gen(was_seeded ? seed : rd());
            std::uniform_int_distribution<> rdAttr(0, pricingIST.getNbDims() - 1);
            for (; rmpQSet.size() < size;) {
                int lAttrID = rdAttr(gen);
                //int rAttrID = rdAttr(gen);
                int measureID = 0;
                int gbAttr = rdAttr(gen);
                while (gbAttr == lAttrID) {// || gbAttr == rAttrID
                    gbAttr = rdAttr(gen);
                }

                std::discrete_distribution<> rdAttSkewed(pricingIST.getDimWeights().begin(), pricingIST.getDimWeights().end());

                int laval = rdAttSkewed(gen);
                int rval = rdAttSkewed(gen);
                while (rval == laval){
                    rval = rdAttSkewed(gen);
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

#endif //CPLEX_TEST_RANDOMINITSKEW_H

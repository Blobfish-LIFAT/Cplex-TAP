//
// Created by alex on 26/01/23.
//

#ifndef CPLEX_TEST_KNAPSACKSOLVER_H
#define CPLEX_TEST_KNAPSACKSOLVER_H

#include <vector>
#include "Query.h"
#include "CGTAPInstance.h"

namespace cplex_tap {

    class KnapsackSolver {
    protected:
        CGTAPInstance *ist;
    public:
        explicit KnapsackSolver(CGTAPInstance *ist);

    protected:

        static bool sortbysec_rev(const std::pair<int,double> &a, const std::pair<int,double> &b){
            return (a.second > b.second);
        }
        double insert_opt(std::vector<int> *solution, int candidate, std::vector<Query> *queries);

    public:
        std::vector<int> solve(std::vector<Query> queries, int timeBudget, int maxDistance);


    };

} // cplex_tap

#endif //CPLEX_TEST_KNAPSACKSOLVER_H

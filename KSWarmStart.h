//
// Created by achan on 26/10/2026.
//

#ifndef CPLEX_TEST_KSWARMSTART_H
#define CPLEX_TEST_KSWARMSTART_H
#include "instance.h"
#include "Solution.h"


class KSWarmStart
{
protected:
    const cplex_tap::Instance& tap;
public:
    KSWarmStart();

    cplex_tap::Solution solve(int dist_bound, int time_bound) const
    {

    }
};


#endif //CPLEX_TEST_KSWARMSTART_H

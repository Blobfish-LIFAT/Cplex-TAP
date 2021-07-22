//
// Created by chanson on 6/22/2021.
//

#ifndef CPLEX_TEST_SOLVERVPLSDET_H
#define CPLEX_TEST_SOLVERVPLSDET_H

#include "solver.h"

namespace cplex_tap {
    class SolverVPLSDet : public Solver {
    public:
        // Builds a solver the specified instance
        explicit SolverVPLSDet(const Instance &tap, int maxIter, int h, int maxInitTime, int maxEpochTime)
                : Solver{tap}, max_iter(maxIter), h(h), max_epoch_time(maxEpochTime), max_init_time(maxInitTime) {}

        Solution solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production, bool seed,
                               string warmStart) const override;

    protected:
        int h;
        int max_iter;
        int max_init_time;
        int max_epoch_time;
    };

}

#endif //CPLEX_TEST_SOLVERVPLSDET_H

#ifndef CPLEX_TEST_SOLVERMPLS_H
#define CPLEX_TEST_SOLVERMPLS_H

#include "solver.h"

namespace cplex_tap {
    class SolverMPLS : public Solver {
    public:
        // Builds a solver the specified instance
        explicit SolverMPLS(const Instance& tap)
        : Solver{ tap }
        {}

        double solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production) const;
    };
}
#endif //CPLEX_TEST_SOLVERMPLS_H

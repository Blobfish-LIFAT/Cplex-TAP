#ifndef CPLEX_TEST_SOLVERMPLS_H
#define CPLEX_TEST_SOLVERMPLS_H

#include "solver.h"

namespace cplex_tap {
    class SolverMPLS : public Solver {
    public:
        // Builds a solver the specified instance
        explicit SolverMPLS(const Instance &tap, int maxIter, int h)
        : Solver{ tap }, max_iter(maxIter), h(h) {}

        double solve_and_print(int dist_bound, int time_bound, bool progressive, bool debug, bool production) const;

    protected:
        int h;
        int max_iter;
    };
}
#endif //CPLEX_TEST_SOLVERMPLS_H

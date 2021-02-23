#pragma once

#include "instance.h"
#include <iostream>

#ifndef IL_STD
#define IL_STD
#endif

#include <cstring>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

namespace cplex_tap {
    class Solver {
        // The TAP instance
        const Instance& tap;

        // Prints the solution from CPLEX to stdout
        void print_solution(const IloCplex& cplex, const IloArray<IloNumVarArray>& x) const;

        // dump decision variables
        void print_X(const IloCplex& cplex, const IloArray<IloNumVarArray>& x) const;

    public:

        // Builds a solver the specified instance
        explicit Solver(const Instance& tap) : tap{ tap } {}

        // Run solver and dump result to stdout
        double solve_and_print(int int_bound, int time_bound) const;

    };
}


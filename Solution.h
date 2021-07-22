//
// Created by chanson on 7/19/2021.
//

#ifndef CPLEX_TEST_SOLUTION_H
#define CPLEX_TEST_SOLUTION_H

namespace cplex_tap {
    class Solution {
    public:
        Solution(double time, double z, const std::vector<int> &sequence) : time(time), z(z), sequence(sequence) {}

        double time;
        double z;
        std::vector<int> sequence;

    };
}

#endif //CPLEX_TEST_SOLUTION_H

//
// Created by chanson on 7/19/2021.
//

#ifndef CPLEX_TEST_SOLUTION_H
#define CPLEX_TEST_SOLUTION_H

namespace cplex_tap {
    class Solution {
    public:
        Solution(double time, double z, const std::vector<int> &sequence, int nodes) : time(time), z(z), sequence(sequence), nodes(nodes) {}
        Solution(double time, double z, const std::vector<int> &sequence) : time(time), z(z), sequence(sequence), nodes(0) {}

        double time;
        double z;
        int nodes;
        std::vector<int> sequence;

    };
}

#endif //CPLEX_TEST_SOLUTION_H

//
// Created by chanson on 7/19/2021.
//

#ifndef CPLEX_TEST_SOLUTION_H
#define CPLEX_TEST_SOLUTION_H

namespace cplex_tap {
    class Solution {
    public:
        Solution(bool optimal, double time, double z, const std::vector<int> &sequence, int nodes) : time(time), z(z), nodes(nodes), sequence(sequence), optimal(optimal) {}
        Solution(bool optimal, double time, double z, const std::vector<int> &sequence) : time(time), z(z), nodes(0), sequence(sequence), optimal(optimal) {}

        double time;
        double z;
        int nodes;
        bool optimal;
        std::vector<int> sequence;

    };
}

#endif //CPLEX_TEST_SOLUTION_H

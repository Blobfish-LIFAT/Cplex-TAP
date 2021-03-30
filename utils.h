//
// Created by chanson on 3/30/2021.
//

#ifndef CPLEX_TEST_UTILS_H
#define CPLEX_TEST_UTILS_H
#include <vector>

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> & vec) {
    std::vector<T> result;
    for (const auto & v : vec)
        result.insert(result.end(), v.begin(), v.end());
    return result;
}
#endif //CPLEX_TEST_UTILS_H

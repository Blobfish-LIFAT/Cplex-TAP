#pragma once

#include <vector>
#include <iostream>
#include <string>

namespace cplex_tap {
    class Instance {
        // Number of queries
        std::uint32_t nbQueries;

        // Distance between queries
        std::vector<std::vector<std::uint32_t>> distances;

        // Time to run a query
        std::vector<std::uint32_t> times;

        // Interstingness of a query
        std::vector<double> interests;



    public:
        const std::vector<std::vector<std::uint32_t>> &getDistances() const;

        // Load instance from a text file
        Instance(std::string file_path);

        //Build instance from memory
        Instance(int nbQ, std::vector<double> interest, std::vector<int> time, std::vector<std::vector<int>> distance);

        // Size of the instance
        std::uint32_t size() const { return nbQueries; }

        // Getters
        std::uint32_t dist(std::uint32_t i, std::uint32_t j) const { return distances[i][j]; }
        std::uint32_t time(std::uint32_t i) const { return times[i]; }
        double interest(std::uint32_t i) const { return interests[i]; }


    };
}

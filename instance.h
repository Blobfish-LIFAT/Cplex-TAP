#pragma once

#include <vector>
#include <iostream>
#include <string>

namespace cplex_tap {
    class Instance {
        // Number of queries
        int nbQueries;

        // Distance between queries
        std::vector<std::vector<int>> distances;

        // Time to run a query
        std::vector<double> times;

        // Interstingness of a query
        std::vector<double> interests;



    public:
        [[nodiscard]] const std::vector<std::vector<int>> &getDistances() const;

        // Load instance from a text file
        Instance(std::string file_path);

        //Build instance from memory
        Instance(int nbQ, std::vector<double> interest, std::vector<double> time, std::vector<std::vector<int>> distance);

        virtual ~Instance();

        // Size of the instance
        [[nodiscard]] int size() const { return nbQueries; }

        // Getters
        [[nodiscard]] int dist(int i, int j) const { return distances[i][j]; }
        [[nodiscard]] double time(int i) const { return times[i]; }
        [[nodiscard]] double interest(int i) const { return interests[i]; }

        //write
        void write(std::string path);

    };
}

#include "instance.h"
#include <cmath>
#include <random>
#include <iterator>
#include <algorithm>
#include <fstream>


namespace cplex_tap {

    Instance::Instance(std::string file_path) {
        using namespace std;
        cout << "File path: " << file_path << endl;

        ifstream tapFile(file_path);
        string line;
        
        // Fetch query number
        getline(tapFile, line);         
        nbQueries = stoi(line, nullptr);
        cout << nbQueries << " Queries to load" << endl;
        
        // Load interestingness
        vector<double> interests_(nbQueries);
        getline(tapFile, line);
        size_t pos = 0;
        int i = 0;
        string token;
        while ((pos = line.find(" ")) != string::npos) {
            token = line.substr(0, pos);
            interests_[i++] = stod(token, nullptr);
            line.erase(0, pos + 1);
        }
        interests_[i++] = stod(line, nullptr);
        interests = interests_;
        cout << "  Interest loaded" << endl;

        // Load time to run
        vector<double> times_(nbQueries);
        getline(tapFile, line);
        pos = 0;
        i = 0;
        token = "";
        while ((pos = line.find(" ")) != string::npos) {
            token = line.substr(0, pos);
            times_[i++] = stod(token, nullptr);
            line.erase(0, pos + 1);
        }
        times_[i++] = stoi(line, nullptr);
        times = times_;
        cout << "  Time loaded" << endl;

        vector<vector<int>> distances_; //(nbQueries, vector<uint32_t>(nbQueries));
        // Load distnces
        for (auto j = 0; j < nbQueries; j++)
        {
            vector<int>  line_vec(nbQueries);
            getline(tapFile, line);
            pos = 0;
            i = 0;
            token = "";
            while ((pos = line.find(" ")) != string::npos) {
                token = line.substr(0, pos);
                line_vec[i++] = stoi(token, nullptr);
                line.erase(0, pos + 1);
            }
            line_vec[i++] = stod(line, nullptr);
            distances_.push_back(line_vec);
           
        }
        distances = distances_;
        cout << "  Distances loaded" << endl;


        tapFile.close();
    }

    const std::vector<std::vector<int>> &Instance::getDistances() const {
        return distances;
    }

    Instance::Instance(int nbQ, std::vector<double> interest, std::vector<double> time,
                       std::vector<std::vector<int>> distance) {
        nbQueries = nbQ;
        interests = interest;
        times = time;
        for (int i = 0; i < nbQ; ++i) {
            std::vector<uint32_t> line;
            for (int j = 0; j < nbQ; ++j) {
                line.emplace_back(distance[i][j]);
            }
            distances.emplace_back(line);
        }

    }

    Instance::~Instance() {
        interests.clear();
        times.clear();
        for (auto line: distances) {
            line.clear();
        }
        distances.clear();
    }

}
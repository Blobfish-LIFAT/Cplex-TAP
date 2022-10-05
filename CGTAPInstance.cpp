//
// Created by alex on 13/09/22.
//

#include "CGTAPInstance.h"
#include <fstream>
#include <iostream>

namespace cplex_tap {
    CGTAPInstance::CGTAPInstance(std::string folder_path) {
        bool DEBUG = true;
        using namespace std;

        string file_path = folder_path + "/__main.dat";
        if (DEBUG) cout << "Main instance file path: " << file_path << endl;

        ifstream mainFile(file_path);
        string line;

        // Fetch table name
        getline(mainFile, line);
        TName = string(line);

        // number of rows
        getline(mainFile, line);
        nbRows = stoi(line, nullptr);

        nbDims = 1;
        // Load dimensions names
        getline(mainFile, line);
        size_t pos = 0;
        int i = 0;
        string token;
        while ((pos = line.find(',')) != string::npos) {
            dimNames.push_back(line.substr(0, pos));
            line.erase(0, pos + 1);
            nbDims++;
        }
        dimNames.push_back(string(line));
        if (DEBUG) cout << "  Dimensions loaded" << endl;

        // Load active domain sizes
        vector<int> ad_sizes_(nbDims);
        getline(mainFile, line);
        pos = 0;i = 0;token = "";
        while ((pos = line.find(',')) != string::npos) {
            token = line.substr(0, pos);
            ad_sizes_[i++] = stoi(token, nullptr);
            line.erase(0, pos + 1);
        }
        ad_sizes_[i++] = stoi(line, nullptr);
        adSizes = ad_sizes_;
        if (DEBUG) cout << "  AD Sizes loaded" << endl;

        nbMeasures = 1;
        // load measures names
        getline(mainFile, line);
        pos = 0;i = 0;token = "";
        while ((pos = line.find(',')) != string::npos) {
            measNames.push_back(line.substr(0, pos));
            line.erase(0, pos + 1);
            nbMeasures++;
        }
        measNames.push_back(string(line));
        if (DEBUG) cout << "  Dimensions loaded" << endl;

        mainFile.close();

        //Load active domains
        for (string att : dimNames){
            string adFilePath = folder_path + "/" + att + "_ad.dat";
            ifstream adFile(adFilePath);
            vector<string> values;

            for (auto j = 0; j < getAdSize(att); j++) {
                string tmp;
                getline(adFile, tmp);
                values.push_back(string(tmp));
            }
            activeDomains.push_back(values);
            adFile.close();
        }

        //Load weights (for linear estimators)
        ifstream attWFile(folder_path + "/dim_weights.dat");
        pos = 0;i = 0;token = "";
        getline(attWFile, line);
        while ((pos = line.find(',')) != string::npos) {
            token = line.substr(0, pos);
            dimWeights.push_back(stod(token, nullptr));
            line.erase(0, pos + 1);
        }
        dimWeights.push_back(stod(line, nullptr));
        attWFile.close();

        ifstream timeWFile(folder_path + "/dim_time.dat");
        pos = 0;i = 0;token = "";
        getline(timeWFile, line);
        while ((pos = line.find(',')) != string::npos) {
            token = line.substr(0, pos);
            dimTimes.push_back(stod(token, nullptr));
            line.erase(0, pos + 1);
        }
        dimTimes.push_back(stod(line, nullptr));
        timeWFile.close();
    }

    const vector<double> &CGTAPInstance::getDimWeights() const {
        return dimWeights;
    }

    const vector<double> &CGTAPInstance::getDimTimes() const {
        return dimTimes;
    }
} // cplex_tap
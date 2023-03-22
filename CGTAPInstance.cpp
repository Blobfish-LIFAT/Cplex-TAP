//
// Created by alex on 13/09/22.
//

#include "CGTAPInstance.h"
#include <fstream>
#include <iostream>
#include "ActiveDomains.h"
#include "UserProfile.h"
#include "FakeTimeStats.h"
#include <fstream>
#include <nlohmann/json.hpp>


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
        std::unordered_map<std::string,std::vector<std::string>> ad_map;
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
            ad_map.insert({att, values});
        }
        ActiveDomains* adSingleton = ActiveDomains::GetInstance(ad_map);

        //Load user profile
        std::unordered_map<std::string, std::unordered_map<std::string,int>> up_map;
        string json_path = folder_path + "/_freq.json";
        using json = nlohmann::json;
        std::ifstream profile(json_path);
        json data = json::parse(profile);

        for (auto it = data.begin(); it != data.end(); ++it) {
            //std::cout << "key: " << it.key() << endl;
            up_map.insert({it.key(), it.value().get<std::unordered_map<std::string,int>>()});
        }
        UserProfile* adUp = UserProfile::GetInstance(up_map);

        std::unordered_map<std::string, std::unordered_map<std::string,int>> time_map;
        json_path = folder_path + "/_freq_timing.json";
        using json = nlohmann::json;
        std::ifstream profile2(json_path);
        data = json::parse(profile2);

        for (auto it = data.begin(); it != data.end(); ++it) {
            //std::cout << "key: " << it.key() << endl;
            time_map.insert({it.key(), it.value().get<std::unordered_map<std::string,int>>()});
        }
        FakeTimeStats* tmap = FakeTimeStats::GetInstance(time_map);

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
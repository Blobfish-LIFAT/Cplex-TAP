//
// Created by alex on 13/09/22.
//

#ifndef CPLEX_TEST_CGTAPINSTANCE_H
#define CPLEX_TEST_CGTAPINSTANCE_H

#include <string>
#include <vector>
#include <algorithm>

using namespace std;

namespace cplex_tap {

    class CGTAPInstance {
        // Number of tuples
        int nbRows;
        // Number of dimensions
        int nbDims;
        // Number of measures
        int nbMeasures;

        //table name
        string TName;
        // dim names
        std::vector<string> dimNames;
        // measures names
        std::vector<string> measNames;

        //Dimensions Active domains
        std::vector<int> adSizes;
        std::vector<std::vector<string>> activeDomains;

        //weights
        std::vector<double> dimWeights;
        std::vector<double> dimTimes;

    public:

        // Load instance from a folder
        CGTAPInstance(std::string folder_path);

        int getNbRows() const {
            return nbRows;
        }

        int getNbDims() const {
            return nbDims;
        }

        const vector<double> &getDimWeights() const;

        int getNbMeasures() const {
            return nbMeasures;
        }

        const string &getTableName() const {
            return TName;
        }

        int getAdSize(int n) const {
            return adSizes[n];
        }

        int getAdSize(string atName) {
            int idx = 0;
            for(auto e : dimNames){
                if (e == atName)
                    break;
                idx++;
            }
            return adSizes[idx];
        }

        double getDimWeight(int n) const {
            return dimWeights[n];
        }

        double getDimWeight(string atName) const{
            int idx = 0;
            for(auto e : dimNames){
                if (e == atName)
                    break;
                idx++;
            }
            return dimWeights[idx];
        }

        double getDimTime(int n) const {
            return dimTimes[n];
        }

        double getDimTime(string atName) const{
            int idx = 0;
            for(auto e : dimNames){
                if (e == atName)
                    break;
                idx++;
            }
            return dimTimes[idx];
        }

        string getDimName(int n) const{
            return dimNames[n];
        }

        string getMeasureName(int n) const{
            return measNames[n];
        }

        const vector<double> &getDimTimes() const;

    };

} // cplex_tap

#endif //CPLEX_TEST_CGTAPINSTANCE_H

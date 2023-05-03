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

        //Time weights
        std::vector<std::vector<int>>  timing;
    public:
        const vector<std::vector<int>> &getTiming() const;

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
            int idx = -1;
            for(int i = 0; i < dimNames.size(); ++i){
                if (dimNames[i] == atName) {
                    idx = i;
                    break;
                }
            }
            return adSizes[idx];
        }

        int getDimId(string atName) const{
            int idx = -1;
            for(int i = 0; i < dimNames.size(); ++i){
                if (dimNames[i] == atName) {
                    idx = i;
                    break;
                }
            }
            return idx;
        }

        int getMeasureId(string atName) const{
            int idx = -1;
            for(int i = 0; i < measNames.size(); ++i){
                if (measNames[i] == atName) {
                    idx = i;
                    break;
                }
            }
            return idx;
        }

        string getDimName(int n) const{
            return dimNames[n];
        }

        string getMeasureName(int n) const{
            return measNames[n];
        }


    };

} // cplex_tap

#endif //CPLEX_TEST_CGTAPINSTANCE_H

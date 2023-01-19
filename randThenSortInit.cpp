//
// Created by alex on 19/01/23.
//

#include "randThenSortInit.h"

namespace cplex_tap {
    vector<Query> randThenSortInit::build(int setSize) {
        RandomInit rdini = RandomInit(pricingIST);
        vector<Query> rands = rdini.build(setSize*rand_per_selected);
        vector<double> interests = JVMAdapter::getInterest(rands, pricingIST);
        vector<int> times = JVMAdapter::getTime(rands, pricingIST);
        vector<pair<int,double>> pos_ratio;
        for (int i = 0; i < rands.size(); ++i) {
            pos_ratio.emplace_back(make_pair(i, interests[i]/static_cast<double>(times[i])));
        }
        sort(pos_ratio.begin(), pos_ratio.end(), sortbysec_rev);
        vector<Query> pool;
        if (debug) cout << endl << "RATIOS: ";
        for (int i = 0; i < setSize; ++i) {
            pool.emplace_back(rands[pos_ratio[i].first]);
            if (debug) cout << pos_ratio[i].second << " ";
        }
        if (debug) cout << endl;
        return pool;
    }
} // cplex_tap
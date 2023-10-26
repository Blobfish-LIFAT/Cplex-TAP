//
// Created by alex on 27/09/22.
//

#ifndef CPLEX_TEST_JVMADAPTER_H
#define CPLEX_TEST_JVMADAPTER_H

#include "Query.h"
#include "CGTAPInstance.h"
#include "set"
#include "ActiveDomains.h"
#include "UserProfile.h"
#include <iostream>


class JVMAdapter {

public:
    /*static vector<double> getInterest(const std::vector<cplex_tap::Query>& qs, const cplex_tap::CGTAPInstance& ist) {
        using namespace cplex_tap;
        vector<double> interest;
        interest.reserve(qs.size());

        for (const auto& q: qs) {
            std::set<std::string> usedAttributes;
            //usedAttributes.insert(q.getGbAttribute());
            for (const auto& pair: q.getLeftPredicate())
                usedAttributes.insert(pair.first);
            for (const auto& pair: q.getRightPredicate())
                usedAttributes.insert(pair.first);

            double  tmp = 0;
            for (const std::string& item: usedAttributes) {
                tmp += ist.getDimWeight(item);
            }
            interest.emplace_back(tmp);
        }

        return interest;
    }*/

    static vector<double> getInterest(const std::vector<cplex_tap::Query>& qs, const cplex_tap::CGTAPInstance& ist) {
        UserProfile* upPT = UserProfile::GetInstance();
        std::unordered_map<std::string, std::unordered_map<std::string,int>> freqs = upPT->value();

        ActiveDomains* adSingleton = ActiveDomains::GetInstance();
        std::unordered_map<std::string,std::vector<std::string>> ad = adSingleton->value();

        using namespace cplex_tap;
        vector<double> interest;
        interest.reserve(qs.size());
        double table_card = ist.getNbRows();

        for (auto  q : qs){
            string dim = q.getLeftPredicate()[0].first;
            int lsel_id = q.getLeftPredicate()[0].second;
            int rsel_id = q.getRightPredicate()[0].second;
            string lsel = ad[dim][lsel_id];
            string rsel = ad[dim][rsel_id];
            double l_freq = freqs[dim][lsel];
            double r_freq = freqs[dim][rsel];

            interest.push_back((l_freq + r_freq)/table_card);
        }

        return interest;
    }

    static vector<double> getTime(const std::vector<cplex_tap::Query>& qs, const cplex_tap::CGTAPInstance& ist){
        using namespace cplex_tap;
        vector<double> interest;
        interest.reserve(qs.size());

        for (const auto&  q : qs){
            interest.push_back(ist.getTiming()[ist.getDimId(q.getGbAttribute())][ist.getDimId(q.getLeftPredicate()[0].first)]);
        }

        return interest;
    }



};


#endif //CPLEX_TEST_JVMADAPTER_H

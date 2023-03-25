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
#include "FakeTimeStats.h"
#include <cpr/cpr.h>
#include <iostream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

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
        FakeTimeStats* upPT = FakeTimeStats::GetInstance();
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

            interest.push_back( 10000 * ( (l_freq + r_freq)/table_card ) );
        }

        return interest;
    }

    /*
    static vector<int> getTime(const std::vector<cplex_tap::Query>& qs, const cplex_tap::CGTAPInstance& ist){

        using namespace cplex_tap;
        vector<int> times;
        times.reserve(qs.size());

        for (const auto& q: qs) {
            std::set<std::string> usedAttributes;
            //usedAttributes.insert(q.getGbAttribute());
            for (const auto& pair: q.getLeftPredicate())
                usedAttributes.insert(pair.first);
            for (const auto& pair: q.getRightPredicate())
                usedAttributes.insert(pair.first);

            double tmp = 0;
            for (const std::string& item: usedAttributes) {
                tmp += ist.getDimTime(item);
            }
            times.emplace_back(round(tmp));
        }

        return times;
        /*
        auto ad = ActiveDomains::GetInstance()->value();
        json body;

        for (auto q : qs) {
            json lp;
            for (int i = 0; i < q.getLeftPredicate().size(); ++i) {
                lp[q.getLeftPredicate()[i].first] = ad[string(
                        q.getLeftPredicate()[i].first)][q.getLeftPredicate()[i].second];
            }
            json rp;
            for (int i = 0; i < q.getRightPredicate().size(); ++i) {
                rp[q.getRightPredicate()[i].first] = ad[string(
                        q.getRightPredicate()[i].first)][q.getRightPredicate()[i].second];
            }
            json r = {{"agg",            q.getAgg()},
                      {"gb",             q.getGbAttribute()},
                      {"leftMeasure",    q.getMeasureLeft()},
                      {"rightMeasure",   q.getMeasureRight()},
                      {"leftPredicate",  lp},
                      {"rightPredicate", rp}
            };
            body.push_back(r);
        }

        cpr::Response r = cpr::Post(cpr::Url{"http://localhost:8000/time"}, cpr::Body(body.dump()),
                                    cpr::Header{{"Content-Type", "application/json"}}, cpr::VerifySsl{false});
        cout << "Query done, code " <<  r.status_code << endl;

        json response = json::parse(r.text);
        vector<int> times;
        for (auto e : response){
            //cout << e << endl;
            times.push_back(e);
        }
        return response;
    }

    /*
     *         std::cout << "Response Error: " << r.error.message << std::endl;
        std::cout << "Response Error Code: " << (int)r.error.code << std::endl;
        std::cout << "Response Status Code: " << r.status_code << std::endl;
        std::cout << "Response Text: " << std::endl;
     */

};


#endif //CPLEX_TEST_JVMADAPTER_H

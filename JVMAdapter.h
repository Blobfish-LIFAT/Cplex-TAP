//
// Created by alex on 27/09/22.
//

#ifndef CPLEX_TEST_JVMADAPTER_H
#define CPLEX_TEST_JVMADAPTER_H

#include "Query.h"
#include "CGTAPInstance.h"
#include "set"
#include "ActiveDomains.h"
#include <cpr/cpr.h>
#include <iostream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

class JVMAdapter {

public:
    static vector<double> getInterest(std::vector<cplex_tap::Query> qs, cplex_tap::CGTAPInstance ist) {
        using namespace cplex_tap;
        vector<double> interest;

        for (auto q: qs) {
            std::set<std::string> usedAttributes;
            usedAttributes.insert(q.getGbAttribute());
            for (auto pair: q.getLeftPredicate())
                usedAttributes.insert(pair.first);
            for (auto pair: q.getRightPredicate())
                usedAttributes.insert(pair.first);

            double  tmp = 0;
            for (std::string item: usedAttributes) {
                tmp += ist.getDimWeight(item);
            }
            interest.emplace_back(tmp);
        }

        return interest;
    }

    static vector<int> getTime(std::vector<cplex_tap::Query> qs, cplex_tap::CGTAPInstance ist){
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

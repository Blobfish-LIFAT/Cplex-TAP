#ifndef CPLEX_TEST_CHECKINSIGHT_H
#define CPLEX_TEST_CHECKINSIGHT_H

#include "CGTAPInstance.h"
#include "Query.h"

#include <cpr/cpr.h>
#include <nlohmann/json.hpp>
#include <iostream>

using json = nlohmann::json;

namespace cplex_tap {

    class CheckInsight {

    public :
        static bool checkForInsight(const cplex_tap::Query &q, const cplex_tap::CGTAPInstance &ist) {

            using namespace cplex_tap;


            ActiveDomains *adSingleton = ActiveDomains::GetInstance();
            std::unordered_map<std::string, std::vector<std::string>> ad = adSingleton->value();

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
            json body = {{"agg",            q.getAgg()},
                         {"gb",             q.getGbAttribute()},
                         {"leftMeasure",    q.getMeasureLeft()},
                         {"rightMeasure",   q.getMeasureRight()},
                         {"leftPredicate",  lp},
                         {"rightPredicate", rp}
            };

            cpr::Response r = cpr::Post(cpr::Url{"http://localhost:4242/insight"}, cpr::Body(body.dump()),
                                        cpr::Header{{"Content-Type", "application/json"}}, cpr::VerifySsl{false});
            std::cout << "Query done, code " << r.status_code << endl;

            json response = json::parse(r.text);

            std::cout << "Response Error: " << r.error.message << std::endl;
            std::cout << "Response Error Code: " << (int) r.error.code << std::endl;
            std::cout << "Response Status Code: " << r.status_code << std::endl;
            std::cout << "Response Text: " << std::endl;

            return response["response"];
        }


        static vector<bool>
        checkForInsights(const std::vector<cplex_tap::Query> &qs, const cplex_tap::CGTAPInstance &ist) {
            using namespace cplex_tap;
            ActiveDomains *adSingleton = ActiveDomains::GetInstance();
            std::unordered_map<std::string, std::vector<std::string>> ad = adSingleton->value();

            json body;

            for (auto q: qs) {
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

            cpr::Response r = cpr::Post(cpr::Url{"http://localhost:4242/insights"}, cpr::Body(body.dump()),
                                        cpr::Header{{"Content-Type", "application/json"}}, cpr::VerifySsl{false});
            cout << "Query done, code " << r.status_code << endl;

            json response = json::parse(r.text);
            vector<bool> times;
            for (auto e: response) {
                times.push_back(e);
            }
            return response;
        }


    };

} // cplex_tap

#endif //CPLEX_TEST_CHECKINSIGHT_H

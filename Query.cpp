//
// Created by alex on 18/09/22.
//

#include "Query.h"

namespace cplex_tap {
    const std::string &Query::getTable() const {
        return table;
    }

    const std::string &Query::getGbAttribute() const {
        return GBAttribute;
    }

    const std::string &Query::getMeasureLeft() const {
        return measureLeft;
    }

    const std::string &Query::getMeasureRight() const {
        return measureRight;
    }

    const std::vector<std::pair<std::string, int>> &Query::getLeftPredicate() const {
        return leftPredicate;
    }

    const std::vector<std::pair<std::string, int>> &Query::getRightPredicate() const {
        return rightPredicate;
    }

    double Query::getInterest() const {
        return interest;
    }

    void Query::setInterest(double interest) {
        Query::interest = interest;
    }

    int Query::getTime() const {
        return time;
    }

    void Query::setTime(int time) {
        Query::time = time;
    }

    Query::Query(const std::string &table, const std::string &agg, const std::string &gbAttribute, const std::string &measureLeft,
                 const std::string &measureRight, const std::vector<std::pair<std::string, int>> &leftPredicate,
                 const std::vector<std::pair<std::string, int>> &rightPredicate) : table(table.substr(0, 256)),
                                                                                           agg(agg),
                                                                                           GBAttribute(gbAttribute),
                                                                                           measureLeft(measureLeft),
                                                                                           measureRight(measureRight),
                                                                                           leftPredicate(leftPredicate),
                                                                                           rightPredicate(
                                                                                                   rightPredicate) {}

    int Query::dist(Query &other) const {
        int diffs = 0;
        // Agg function changed ?
        if(agg != other.agg)  diffs += 1;
        // Measure changed ?
        if(measureLeft != other.measureLeft) diffs += 1;
        if(measureRight != other.measureRight) diffs += 1;

        // predicates
        if(leftPredicate.size() != other.leftPredicate.size()) {
            diffs += 2;
        } else{
            bool changed = false;
            for (int i = 0; i < leftPredicate.size(); ++i) {
                    changed |= leftPredicate[i] != other.leftPredicate[i];
            }
            diffs += 2*changed;
        }
        if(rightPredicate.size() != other.rightPredicate.size()) {
            diffs += 2;
        } else{
            bool changed = false;
            for (int i = 0; i < rightPredicate.size(); ++i) {
                changed |= rightPredicate[i] != other.rightPredicate[i];
            }
            diffs += 2*changed;
        }

        // Group by dimension
        if (GBAttribute != other.GBAttribute){
            diffs += 5;
        }

        return diffs;
    }

    const std::string &Query::getAgg() const {
        return agg;
    }

    bool Query::operator==(const Query &rhs) const {
        return interest == rhs.interest &&
               time == rhs.time &&
               table == rhs.table &&
               GBAttribute == rhs.GBAttribute &&
               measureLeft == rhs.measureLeft &&
               measureRight == rhs.measureRight &&
               agg == rhs.agg &&
               leftPredicate == rhs.leftPredicate &&
               rightPredicate == rhs.rightPredicate;
    }

    bool Query::operator!=(const Query &rhs) const {
        return !(rhs == *this);
    }
} // cplex_tap
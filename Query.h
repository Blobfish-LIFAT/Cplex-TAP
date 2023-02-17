//
// Created by alex on 18/09/22.
//

#include <string>
#include <vector>
#include <utility>
#include <boost/functional/hash.hpp>

#ifndef CPLEX_TEST_QUERY_H
#define CPLEX_TEST_QUERY_H

namespace cplex_tap {

    class Query {
        double interest;
        int time;

        std::string table;
        std::string GBAttribute;
        std::string measureLeft;
        std::string measureRight;
        std::string agg;

        std::vector<std::pair<std::string, int>> leftPredicate;
        std::vector<std::pair<std::string, int>> rightPredicate;


    public:
        Query(const std::string &table, const std::string &agg, const std::string &gbAttribute, const std::string &measureLeft,
              const std::string &measureRight, const std::vector<std::pair<std::string, int>> &leftPredicate,
              const std::vector<std::pair<std::string, int>> &rightPredicate);

        double getInterest() const;

        void setInterest(double interest);

        int getTime() const;

        const std::string &getAgg() const;

        void setTime(int time);
        const std::string &getTable() const;

        const std::string &getGbAttribute() const;

        const std::string &getMeasureLeft() const;

        const std::string &getMeasureRight() const;

        const std::vector<std::pair<std::string, int>> &getLeftPredicate() const;

        const std::vector<std::pair<std::string, int>> &getRightPredicate() const;

        int dist(Query&  other) const;

        bool operator==(const Query &rhs) const;

        bool operator!=(const Query &rhs) const;

        friend std::ostream& operator<<(std::ostream &strm, const Query &q) {
            strm << std::string("Query(") << q.table << std::string(",");
            strm << q.agg << std::string(",");
            strm << q.measureLeft << std::string(",");
            strm << q.measureRight << std::string(",");
            strm << q.GBAttribute << std::string(",");
            for (auto i = q.leftPredicate.begin(); i != q.leftPredicate.end(); ++i) {
                strm << i->first << std::string("=") << std::to_string(i->second);
                if (i != q.leftPredicate.end()-1)
                    strm << std::string(" AND ");
            }
            strm << std::string(",");
            for (auto i = q.rightPredicate.begin(); i != q.rightPredicate.end(); ++i) {
                strm << i->first << std::string("=") << std::to_string(i->second);
                if (i != q.rightPredicate.end()-1)
                    strm << std::string(" AND ");
            }
            //strm << std::string(",");
            return strm << std::string(")");
        }
    };

} // cplex_tap

namespace std
{
    template <>
    struct hash<cplex_tap::Query>
    {
        size_t operator()(const cplex_tap::Query& k) const
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, k.getTime());
            boost::hash_combine(seed, k.getInterest());
            boost::hash_combine(seed, k.getInterest());
            boost::hash_combine(seed, k.getTable());
            boost::hash_combine(seed, k.getAgg());
            boost::hash_combine(seed, k.getGbAttribute());
            boost::hash_combine(seed, k.getMeasureLeft());
            boost::hash_combine(seed, k.getMeasureRight());
            boost::hash_combine(seed, k.getLeftPredicate());
            boost::hash_combine(seed, k.getRightPredicate());

        return seed;
        }
    };
}

#endif //CPLEX_TEST_QUERY_H

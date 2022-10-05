//
// Created by alex on 27/09/22.
//

#ifndef CPLEX_TEST_JVMADAPTER_H
#define CPLEX_TEST_JVMADAPTER_H

#include "Query.h"
#include "CGTAPInstance.h"
#include "set"

class JVMAdapter {

public:
    static double getInterest(cplex_tap::Query q, cplex_tap::CGTAPInstance ist){
        using namespace cplex_tap;
        double interest = 0;

        std::set<std::string> usedAttributes;
        usedAttributes.insert(q.getGbAttribute());
        for (auto pair : q.getLeftPredicate())
            usedAttributes.insert(pair.first);
        for (auto pair : q.getRightPredicate())
            usedAttributes.insert(pair.first);

        for (std::string item : usedAttributes) {
            interest += ist.getDimWeight(item);
        }

        return interest;
    }

    static int getTime(cplex_tap::Query q, cplex_tap::CGTAPInstance ist){
        return 1;
    }

};


#endif //CPLEX_TEST_JVMADAPTER_H

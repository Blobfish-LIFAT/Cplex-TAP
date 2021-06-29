//
// Created by chanson on 6/29/2021.
//

#ifndef CPLEX_TEST_ELEMENT_H
#define CPLEX_TEST_ELEMENT_H


class Element {
public:
    Element(int index, double value) : index(index), value(value) {}

public:
    int index;
    double value;

    bool operator<(const Element& other) const
    {
        return value < other.value;
    }

    bool operator>(const Element& other) const
    {
        return value > other.value;
    }
};


#endif //CPLEX_TEST_ELEMENT_H

//
// Created by alex on 22/03/23.
//

#ifndef CPLEX_TEST_FAKETIMESTATS_H
#define CPLEX_TEST_FAKETIMESTATS_H


#include <mutex>
#include <vector>
#include <unordered_map>

class FakeTimeStats{

private:
    static FakeTimeStats * pinstance_;
    static std::mutex mutex_;

protected:
    FakeTimeStats(const std::unordered_map<std::string, std::unordered_map<std::string,int>> value): value_(value)    {}
    ~FakeTimeStats() {}    std::unordered_map<std::string, std::unordered_map<std::string,int>> value_;

public:

    FakeTimeStats(FakeTimeStats &other) = delete;
    void operator=(const FakeTimeStats &) = delete;

    static FakeTimeStats *GetInstance(const std::unordered_map<std::string, std::unordered_map<std::string,int>>& value);
    static FakeTimeStats *GetInstance();

    std::unordered_map<std::string, std::unordered_map<std::string,int>> value() const{
        return value_;
    }
};


#endif //CPLEX_TEST_FAKETIMESTATS_H

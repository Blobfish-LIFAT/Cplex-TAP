//
// Created by alex on 16/03/23.
//

#ifndef CPLEX_TEST_USERPROFILE_H
#define CPLEX_TEST_USERPROFILE_H


#include <mutex>
#include <vector>
#include <unordered_map>

class UserProfile{

private:
    static UserProfile * pinstance_;
    static std::mutex mutex_;

protected:
    UserProfile(const std::unordered_map<std::string, std::unordered_map<std::string,int>> value): value_(value)    {}
    ~UserProfile() {}    std::unordered_map<std::string, std::unordered_map<std::string,int>> value_;

public:

    UserProfile(UserProfile &other) = delete;
    void operator=(const UserProfile &) = delete;

    static UserProfile *GetInstance(const std::unordered_map<std::string, std::unordered_map<std::string,int>>& value);
    static UserProfile *GetInstance();

    std::unordered_map<std::string, std::unordered_map<std::string,int>> value() const{
        return value_;
    }
};


#endif //CPLEX_TEST_USERPROFILE_H




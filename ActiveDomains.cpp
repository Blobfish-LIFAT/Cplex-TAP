//
// Created by alex on 07/10/22.
//
#include "ActiveDomains.h"



ActiveDomains* ActiveDomains::pinstance_{nullptr};
std::mutex ActiveDomains::mutex_;

ActiveDomains *ActiveDomains::GetInstance(const std::unordered_map<std::string,std::vector<std::string>>& value)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (pinstance_ == nullptr)
    {
        pinstance_ = new ActiveDomains(value);
    }
    return pinstance_;

}ActiveDomains *ActiveDomains::GetInstance()
{
    std::lock_guard<std::mutex> lock(mutex_);
    return pinstance_;
}

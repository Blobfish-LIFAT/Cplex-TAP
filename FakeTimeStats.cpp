//
// Created by alex on 22/03/23.
//

#include "FakeTimeStats.h"

FakeTimeStats* FakeTimeStats::pinstance_{nullptr};
std::mutex FakeTimeStats::mutex_;

FakeTimeStats *FakeTimeStats::GetInstance(const std::unordered_map<std::string, std::unordered_map<std::string,int>>& value)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (pinstance_ == nullptr)
    {
        pinstance_ = new FakeTimeStats(value);
    }
    return pinstance_;

}FakeTimeStats *FakeTimeStats::GetInstance()
{
    std::lock_guard<std::mutex> lock(mutex_);
    return pinstance_;
}
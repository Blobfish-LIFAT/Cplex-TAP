#include "UserProfile.h"

UserProfile* UserProfile::pinstance_{nullptr};
std::mutex UserProfile::mutex_;

UserProfile *UserProfile::GetInstance(const std::unordered_map<std::string, std::unordered_map<std::string,int>>& value)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (pinstance_ == nullptr)
    {
        pinstance_ = new UserProfile(value);
    }
    return pinstance_;

}UserProfile *UserProfile::GetInstance()
{
    std::lock_guard<std::mutex> lock(mutex_);
    return pinstance_;
}

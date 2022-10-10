#include <mutex>
#include <vector>
#include <unordered_map>

class ActiveDomains
{

private:
    static ActiveDomains * pinstance_;
    static std::mutex mutex_;

protected:
    ActiveDomains(const std::unordered_map<std::string,std::vector<std::string>> value): value_(value)
    {}
    ~ActiveDomains() {}
    std::unordered_map<std::string,std::vector<std::string>> value_;

public:

    ActiveDomains(ActiveDomains &other) = delete;
    void operator=(const ActiveDomains &) = delete;

    static ActiveDomains *GetInstance(const std::unordered_map<std::string,std::vector<std::string>>& value);
    static ActiveDomains *GetInstance();

    std::unordered_map<std::string,std::vector<std::string>> value() const{
        return value_;
    }
};

#ifndef FRAMEWORK_UTIL_UTIL_TIME_H__
#define FRAMEWORK_UTIL_UTIL_TIME_H__

#include "../Base.h"

class UtilTime
{
public:
    static uint32_t getNow();
    static uint64_t getNowMS();
    static uint64_t getNowUS();
};


class RunTimeCalculator
{
public:
    explicit RunTimeCalculator(const char* szFuncName);
    ~RunTimeCalculator();

private:
    std::string     _strFunName;
    clock_t         _start;
};

#endif  // FRAMEWORK_UTIL_UTIL_TIME_H__

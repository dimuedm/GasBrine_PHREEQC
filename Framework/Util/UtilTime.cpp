#include "UtilTime.h"
#ifdef WIN32
    #include <winsock2.h>
    #include <time.h>
#else
    #include <time.h>
#endif
#include "UtilLog.h"

uint32_t UtilTime::getNow()
{
#ifdef WIN32
    struct timeval tv;
    struct tm tm;
    SYSTEMTIME wtm;

    GetLocalTime(&wtm);
    tm.tm_year = wtm.wYear - 1900;
    tm.tm_mon = wtm.wMonth - 1;
    tm.tm_mday = wtm.wDay;
    tm.tm_hour = wtm.wHour;
    tm.tm_min = wtm.wMinute;
    tm.tm_sec = wtm.wSecond;
    tm.tm_isdst = -1;
    tv.tv_sec = (uint32_t)mktime(&tm);
    //tv.tv_usec = wtm.wMilliseconds * 1000;
    return ((uint32_t)tv.tv_sec);
#else
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return ((uint32_t)tv.tv_sec);
#endif
}


uint64_t UtilTime::getNowMS()
{
#ifdef WIN32
    struct timeval tv;
    struct tm tm;
    SYSTEMTIME wtm;

    GetLocalTime(&wtm);
    tm.tm_year = wtm.wYear - 1900;
    tm.tm_mon = wtm.wMonth - 1;
    tm.tm_mday = wtm.wDay;
    tm.tm_hour = wtm.wHour;
    tm.tm_min = wtm.wMinute;
    tm.tm_sec = wtm.wSecond;
    tm.tm_isdst = -1;
    tv.tv_sec = (uint32_t)mktime(&tm);
    tv.tv_usec = wtm.wMilliseconds * 1000;
    return ((uint64_t)tv.tv_sec * 1000 + (uint64_t)tv.tv_usec / 1000);
#else
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return ((uint64_t)tv.tv_sec * 1000 + (uint64_t)tv.tv_usec / 1000);
#endif
}
    
    

uint64_t UtilTime::getNowUS()
{
#ifdef WIN32
    struct timeval tv;
    struct tm tm;
    SYSTEMTIME wtm;

    GetLocalTime(&wtm);
    tm.tm_year = wtm.wYear - 1900;
    tm.tm_mon = wtm.wMonth - 1;
    tm.tm_mday = wtm.wDay;
    tm.tm_hour = wtm.wHour;
    tm.tm_min = wtm.wMinute;
    tm.tm_sec = wtm.wSecond;
    tm.tm_isdst = -1;
    tv.tv_sec = (uint32_t)mktime(&tm);
    tv.tv_usec = wtm.wMilliseconds * 1000;
    return ((uint64_t)tv.tv_sec * 1000000 + (uint64_t)tv.tv_usec);
#else
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return ((uint64_t)tv.tv_sec * 1000000 + (uint64_t)tv.tv_usec);
#endif
}


RunTimeCalculator::RunTimeCalculator(const char* szFuncName)
{
    _strFunName = szFuncName;
    _start = clock();
}


RunTimeCalculator::~RunTimeCalculator()
{
    double runT = (clock() - _start);
    LOG(TRACE) << "when " << _strFunName.c_str() << ", using time " << runT << " ms.";
}

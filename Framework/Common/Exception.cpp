#include "Exception.h"
#include <cstring>
//#include <stdio.h>
#include <cerrno>
//#include <errno.h>
//#include <cstdlib>
//using namespace std;

Exception::Exception(const string &buffer)
    : m_buffer_(buffer), m_code_(0)
{
    //    getBacktrace();
}

Exception::Exception(const string &buffer, int err)
{
    string strerr = strerror(err);
    m_buffer_ = buffer + " :" + strerr;
    m_code_ = err;
    //    getBacktrace();
}

Exception::~Exception() throw()
{
}

const char* Exception::what() const throw()
{
    return m_buffer_.c_str();
}

void Exception::getBacktrace()
{
    //void * array[64];
    //int nSize = backtrace(array, 64);
    //char ** symbols = backtrace_symbols(array, nSize);

    //for (int i = 0; i < nSize; i++)
    //{
    //    m_buffer_ += symbols[i];
    //    m_buffer_ += "\n";
    //}
    //free(symbols);
}

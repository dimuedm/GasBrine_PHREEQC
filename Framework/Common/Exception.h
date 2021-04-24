#ifndef __KRI_EXCEPTION_H__
#define __KRI_EXCEPTION_H__

#include <string>
#include <stdexcept>

using namespace std;

class Exception : public exception
{
public:
    explicit Exception(const string &buffer);

    Exception(const string &buffer, int err);

    virtual ~Exception() throw();

    virtual const char* what() const throw();

    int getErrCode() { return m_code_; }

private:
    void getBacktrace();

private:
    string          m_buffer_;
    int             m_code_;
};

#endif


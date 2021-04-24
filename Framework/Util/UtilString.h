#ifndef FRAMEWORK_UTIL_UTIL_STRING_H_
#define FRAMEWORK_UTIL_UTIL_STRING_H_

#include "../Base.h"

class UtilString
{
public:
    static string trim(const string &str);
    static string trimLeft(const string &str);
    static string trimRight(const string &str);
    static bool isTrimEmpty(const string &str);

    static void toUpperInplace(string &str);
    static void toLowerInplace(string &str);
    static string toUpper(const string &str);
    static string toLower(const string &str);

    static bool isContains(const string &s, const string &match);
    static bool isBeginsWith(const string &s, const string &match);
    static bool isEndsWith(const string &s, const string &match);
    static bool isdigit2(const string &s);
    static bool isxdigit(const string &s);
    
    static int32_t mapToInt(const string &s);

    static string replace(const string &str, const string &sFrom, const string &sTo);
    static string format(const char *pszFormat, ...);

    static string hex2Str(const string &sHex);
    static string str2Hex(const string &sStr);

    static int remove(string& str, char ch);
    static int remove(string& str, const char* pszTrim);
    static int remove(string& str, const string &strTrim);

    static int removeNot(string& str, const char* pszTrim);
   
    static int32_t compareNoCase(const string& str1, const string& str2);
    static int32_t nCompareNoCase(const string& str1, const string& str2, size_t count);
    static int32_t nCompareCase(const string& str1, const string& str2, size_t count);
    
    struct streamstr
    {
        ostringstream os;

        template <typename T>
        streamstr &operator<< (const T &t) { os << t; return *this; }

        typedef ostream &(*manip)(ostream &);
        streamstr &operator<< (manip pf) { os << pf; return *this; }

        operator string () const { return os.str(); }
    };

    template <typename T>
    static string tostr(const T &t);

    template <typename T>
    static T strto(const string &s);

    template <typename T>
    static T strto(const string &s, const T &def)
    {
        if (s.empty()) return def;
        return strto<T>(s);
    }

    template<typename T>
    static T strto(const char* pszStr, int iLen = -1, int* piRet = NULL);

    static void splitString(const string &str, const string &sep, bool bStrip, vector<string> &vResult);
    static vector<string> splitString(const string &str, const string &sep, bool bStrip = true);

    template <typename T> static void splitString(const string &str, const string &sep, bool bStrip, vector<T> &vResult);
    template <typename T> static vector<T> splitString(const string &str, const string &sep, bool bStrip = true);

    static string joinString(const vector<string> &v, const string &sep);
    template <typename T> static string joinString(const vector<T> &v, const string &sep);

    static void splitString2(const string &str, const string &sep1, const string &sep2, vector<pair<string, string> > &vResult);
    static void splitString2(const string &str, const string &sep1, const string &sep2, map<string, string> &mResult);
    static map<string, string> splitString2(const string &str, const string &sep1, const string &sep2);

    template <typename T1, typename T2>
    static void splitString2(const string &str, const string &sep1, const string &sep2, vector<pair<T1, T2> > &vResult);
    template <typename T1, typename T2>
    static void splitString2(const string &str, const string &sep1, const string &sep2, map<T1, T2> &mResult);
    template <typename T1, typename T2>
    static map<T1, T2> splitString2(const string &str, const string &sep1, const string &sep2);

    static string joinString2(const vector<pair<string, string> > &v, const string &sep1, const string &sep2);
    template <typename T1, typename T2>
    static string joinString2(const vector<pair<T1, T2> > &v, const string &sep1, const string &sep2);
    template <typename T1, typename T2>
    static string joinString2(const map<T1, T2> &v, const string &sep1, const string &sep2);

    static string repeat(const string &str, uint32_t n);
    static string repeat(const string &str, const string &sep, uint32_t n);

    static const string &getEmptyString();

    static uint64_t parseHumanReadableSize(const string &s);
    
    enum REPACE_MODE
    {
        REPLACE_MODE_ALL_MATCH,
        REPLACE_MODE_KEEP_VARIABLE_ON_MISS,
        REPLACE_MODE_EMPTY_VARIABLE_ON_MISS,
        REPLACE_MODE_EMPTY_RESULT_ON_MISS,
    };
    static string replaceVariable(const string &str, const map<string, string> &mVariable, REPACE_MODE mode);

    static string joinURLParam(const map<string, string> &mParam, const string &sep1 = "&", const string &sep2 = "=");
    static void splitURLParam(const string &sParam, map<string, string> &mParam, const string &sep1 = "&", const string &sep2 = "=");
    static bool splitIni(const string& strLine, string& strName, string& strValue, char chDel = '=');
    static void splitPair(const string& strLine, const string& sep, string& strFirst, string& strSecond);
};


namespace impl
{
    template <typename T>
    struct tostr_helper
    {
        static string tostr(const T &t)
        {
            ostringstream ss;
            ss << t;
            return ss.str();
        }
    };

    template <>
    struct tostr_helper <float>
    {
        static string tostr(float t)
        {
            ostringstream ss;
            float r = (uint32_t)1 << numeric_limits<float>::digits;
            if (t <= r && t >= -r && floor(t) == t)
            {
                ss << (int32_t)t;
            }
            else
            {
                ss << setprecision(8) << t;
            }
            return ss.str();
        }
    };

    template <>
    struct tostr_helper <double>
    {
        static string tostr(double t)
        {
            ostringstream ss;
            double r = (uint64_t)1 << numeric_limits<double>::digits;
            if (t <= r && t >= -r && floor(t) == t)
            {
                ss << (int64_t)t;
            }
            else
            {
                ss << setprecision(15) << t;
            }
            return ss.str();
        }
    };

    template <typename T>
    struct lexical_castor
    {
        static T cast(const string &s)
        {
            istringstream ss(s);
            T t;
            ss >> t;
            return t;
        }
    };

    template<> struct lexical_castor <uint64_t>
    {
        static uint64_t cast(const string &s)
        {
            return strtoull(s.c_str(), NULL, 0);
        }
    };

    template<> struct lexical_castor <bool>
    {
        static bool cast(const string &s) { return lexical_castor<uint64_t>::cast(s) ? true : false; }
    };

    template<> struct lexical_castor <int8_t>
    {
        static int8_t cast(const string &s) { return static_cast<int8_t>(lexical_castor<uint64_t>::cast(s)); }
    };

    template<> struct lexical_castor <uint8_t>
    {
        static uint8_t cast(const string &s) { return static_cast<uint8_t>(lexical_castor<uint64_t>::cast(s)); }
    };

    template<> struct lexical_castor <int16_t>
    {
        static int16_t cast(const string &s) { return static_cast<int16_t>(lexical_castor<uint64_t>::cast(s)); }
    };

    template<> struct lexical_castor <uint16_t>
    {
        static uint16_t cast(const string &s) { return static_cast<uint16_t>(lexical_castor<uint64_t>::cast(s)); }
    };

    template<> struct lexical_castor <int32_t>
    {
        static int32_t cast(const string &s) { return static_cast<int32_t>(lexical_castor<uint64_t>::cast(s)); }
    };

    template<> struct lexical_castor <uint32_t>
    {
        static uint32_t cast(const string &s) { return static_cast<uint32_t>(lexical_castor<uint64_t>::cast(s)); }
    };

    template<> struct lexical_castor <int64_t>
    {
        static int64_t cast(const string &s) { return static_cast<int64_t>(lexical_castor<uint64_t>::cast(s)); }
    };

    template<> struct lexical_castor <float>
    {
        static float cast(const string &s)
        {
            return strtof(s.c_str(), NULL);
        }
    };

    template<> struct lexical_castor <double>
    {
        static double cast(const string &s)
        {
            return strtod(s.c_str(), NULL);
        }
    };

    template<> struct lexical_castor <string>
    {
        static const string &cast(const string &s)
        {
            return s;
        }
    };

}

template <typename T>
string UtilString::tostr(const T &t)
{
    return impl::tostr_helper<T>::tostr(t);
}

template <typename T>
T UtilString::strto(const string &s)
{
    return impl::lexical_castor<T>::cast(s);
}

template<typename T>
T UtilString::strto(const char* pszStr, int iLen, int* piRet)
{
    const char* pszEnd = (iLen > 0 ? pszStr + iLen -1 : NULL);
    char* p = (char*) pszStr;
    if (p == NULL || *p == 0)
    {
        if (piRet)
        {
            *piRet = -1;
        }
        return 0;
    }

    const bool bIsUnsigned = (((T)(-1)) > 0);
    T basedata = 1;
    T makedata = (T)-1;
    const T AbsMax = bIsUnsigned ? makedata : (basedata<<(sizeof(T)*8 - 1)^makedata);
    const T Overflow = AbsMax / 10;

    while (*p == ' '&& (pszEnd ? p <= pszEnd : true))
    {
        ++p;
    }

    int sign = 0;
    if (*p == '-' || *p == '+')
    {
        if (bIsUnsigned)
        {
            if (piRet)
            {
                *piRet = -2;
            }
            return 0;
        }

        sign = (*p == '-' ? -1 : 1);
        ++p;
    }

    while (*p == '0' && (pszEnd ? p <= pszEnd : true))
    {
        ++p;
    }

    T ret = 0;
    for (; *p != 0 && (pszEnd ? p <= pszEnd : true); ++p)
    {
        if (*p >= '0' && *p <= '9')
        {
            if (ret > Overflow)
            {
                if (piRet)
                {
                    *piRet = -3;
                }
                return 0;
            }
            else if (ret == Overflow)
            {
                if (bIsUnsigned)
                {
                    if (*p > '5')
                    {
                        if (piRet)
                        {
                            *piRet = -3;
                        }
                        return 0;
                    }
                }
                else
                {
                    if ((sign < 0 && *p > '8') || (sign >= 0 && *p > '7'))
                    {
                        if (*p > '5')
                        {
                            if (piRet)
                            {
                                *piRet = -3;
                            }
                            return 0;
                        }
                    }
                }
            }

            ret = ret * 10 + (*p - '0');
        }
        else if (*p == ' ')
        {
            for (; *p != 0 && (pszEnd ? p <= pszEnd : true); ++p)
            {
                if (*p != ' ')
                {
                    if (piRet)
                    {
                        *piRet = -2;
                    }
                    return 0;
                }
            }
            break;
        }
        else
        {
            if (piRet)
            {
                *piRet = -2;
            }
            return 0;
        }
    }
    if(piRet)
    {
        *piRet = 0;
    }
    return (sign < 0 ? 0 - ret : ret);
}

template <typename T>
void UtilString::splitString(const string &str, const string &sep, bool bStrip, vector<T> &vResult)
{
    vector<string> vStrResult;
    splitString(str, sep, bStrip, vStrResult);
    for (unsigned i = 0; i < vStrResult.size(); ++i)
    {
        T tmp = strto<T>(vStrResult[i]);
        vResult.push_back(tmp);
    }
}

template <typename T>
vector<T> UtilString::splitString(const string &str, const string &sep, bool bStrip)
{
    vector<T> vResult;
    splitString<T>(str, sep, bStrip, vResult);
    return vResult;
}

template <typename T1, typename T2>
void UtilString::splitString2(const string &str, const string &sep1, const string &sep2, vector< pair<T1, T2> > &vResult)
{
    vector<pair<string, string> > vStrResult;
    splitString2(str, sep1, sep2, vStrResult);
    for (unsigned i = 0; i < vStrResult.size(); ++i)
    {
        pair<string, string> &kv = vStrResult[i];
        vResult.push_back(make_pair(strto<T1>(kv.first), strto<T2>(kv.second)));
    }
}

template <typename T1, typename T2>
void UtilString::splitString2(const string &str, const string &sep1, const string &sep2, map<T1, T2> &mResult)
{
    vector<pair<string, string> > vStrResult;
    splitString2(str, sep1, sep2, vStrResult);
    for (unsigned i = 0; i < vStrResult.size(); ++i)
    {
        pair<string, string> &kv = vStrResult[i];
        mResult[strto<T1>(kv.first)] = strto<T2>(kv.second);
    }
}

template <typename T1, typename T2>
map<T1, T2> UtilString::splitString2(const string &str, const string &sep1, const string &sep2)
{
    map<T1, T2> mResult;
    splitString2<T1, T2>(str, sep1, sep2, mResult);
    return mResult;
}

template <typename T>
string UtilString::joinString(const vector<T> &v, const string &sep)
{
    string sResult;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (!sResult.empty())
        {
            sResult.append(sep);
        }
        sResult.append(tostr(v[i]));
    }
    return sResult;
}

template <typename T1, typename T2>
string UtilString::joinString2(const vector<pair<T1, T2> > &v, const string &sep1, const string &sep2)
{
    string sResult;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (!sResult.empty())
        {
            sResult.append(sep1);
        }

        const pair<T1, T2> &kv = v[i];
        sResult.append(tostr(kv.first));
        sResult.append(sep2);
        sResult.append(tostr(kv.second));
    }
    return sResult;
}

template <typename T1, typename T2>
string UtilString::joinString2(const map<T1, T2> &v, const string &sep1, const string &sep2)
{
    string sResult;
    for (typename map<T1, T2>::const_iterator first = v.begin(), last = v.end(); first != last; ++first)
    {
        if (!sResult.empty())
        {
            sResult.append(sep1);
        }

        sResult.append(tostr(first->first));
        sResult.append(sep2);
        sResult.append(tostr(first->second));
    }
    return sResult;
}

#endif

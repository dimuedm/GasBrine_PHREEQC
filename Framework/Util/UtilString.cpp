#include "UtilString.h"
 
string UtilString::trim(const string &str)
{
    string::size_type i = 0, j = str.size();
    while (i < j && isspace(str[i])) ++i;
    while (i < j && isspace(str[j - 1])) --j;
    return str.substr(i, j - i);
}

string UtilString::trimLeft(const string &str)
{
    string::size_type i = 0, j = str.size();
    while (i < j && isspace(str[i])) ++i;
    return str.substr(i, j - i);
}

string UtilString::trimRight(const string &str)
{
    string::size_type i = 0, j = str.size();
    while (i < j && isspace(str[j - 1])) --j;
    return str.substr(i, j - i);
}

bool UtilString::isTrimEmpty(const string &str)
{
    string sTemp = trim(str);
    return sTemp.empty();
}

void UtilString::toUpperInplace(string &str)
{
    for (string::size_type i = 0; i < str.size(); ++i)
    {
        str[i] = ::toupper(str[i]);
    }
}

void UtilString::toLowerInplace(string &s)
{
    for (string::size_type i = 0; i < s.size(); ++i)
    {
        s[i] = ::tolower(s[i]);
    }
}

string UtilString::toUpper(const string &str)
{
    string sTemp = str;
    toUpperInplace(sTemp);
    return sTemp;
}

string UtilString::toLower(const string &str)
{
    string sTemp = str;
    toLowerInplace(sTemp);
    return sTemp;
}


bool UtilString::isContains(const string &s, const string &match)
{
    return s.find(match) != string::npos;
}

bool UtilString::isBeginsWith(const string &s, const string &match)
{
    if (match.size() > s.size())
    {
        return false;
    }
    return memcmp(s.data(), match.data(), match.size()) == 0 ? true : false;
}

bool UtilString::isEndsWith(const string &s, const string &match)
{
    if (match.size() > s.size())
    {
        return false;
    }
    return memcmp(s.data() + (s.size() - match.size()), match.data(), match.size()) == 0 ? true : false;
}

bool UtilString::isdigit2(const string &s)
{
    for (string::size_type i = 0; i < s.size(); ++i)
    {
        if (!::isdigit(s[i]))
        {
            return false;
        }
    }
    return true;
}

bool UtilString::isxdigit(const string &s)
{
    for (string::size_type i = 0; i < s.size(); ++i)
    {
        if (!::isxdigit(s[i]))
        {
            return false;
        }
    }
    return true;
}

int32_t UtilString::mapToInt(const string &s)
{
    int32_t iResult = 0;
    for (string::size_type i = 0; i < s.size(); ++i)
    {
        iResult += (int32_t)s[i];
    }
    return iResult;
}

string UtilString::replace(const string &str, const string &sFrom, const string &sTo)
{
    if (sFrom.empty())
    {
        return str;
    }

    string sTemp;
    const char *begin = str.data();
    string::size_type p = 0;
    string::size_type q = str.find(sFrom);
    while (q != string::npos)
    {
        sTemp.append(begin + p, begin + q);
        sTemp.append(sTo);
        p = q + sFrom.size();
        q = str.find(sFrom, p);
    }
    sTemp.append(begin + p, begin + str.size());
    return sTemp;
}


string UtilString::format(const char *pszFormat, ...)
{
    string res;

    va_list ap;
    va_list apcopy;

    va_start(ap, pszFormat);
    va_copy(apcopy, ap);

    int iSize = vsnprintf(NULL, 0, pszFormat, ap);   // 计算需要的内存空间
    if (iSize <= 0)
    {
        return "";
    }
    ++iSize;

    int iLen = 0;
    char* pszBuf = new char[iSize];
    if (!pszBuf)
    {
        return "";
    }
    iLen = vsnprintf(pszBuf, iSize, pszFormat, apcopy);
    if (iLen < 0 || iLen >= iSize)
    {
        delete[] pszBuf;
        return "";
    }

    pszBuf[iLen] = 0;
    res.assign(pszBuf, iLen);
    delete[] pszBuf;

    va_end(apcopy);
    va_end(ap);

    return res;
}

int UtilString::remove(string& str, char ch)
{
    size_t iIndex = str.find(ch);
    int iCount = 0;
    while (iIndex != string::npos)
    {
        str.erase(iIndex, 1);
        iIndex = str.find(ch);
        ++iCount;
    }

    return iCount;
}

int UtilString::remove(string& str, const char* pszTrim)
{
    size_t iIndex = str.find_first_of(pszTrim);
    int iCount = 0;
    while (iIndex != string::npos)
    {
        str.erase(iIndex, 1);
        iIndex = str.find_first_of(pszTrim);
        ++iCount;
    }

    return iCount;
}

int UtilString::remove(string& str, const string &strTrim)
{
    return remove(str, strTrim.c_str());
}

int UtilString::removeNot(string& str, const char* pszTrim)
{
    size_t iIndex = str.find_first_not_of(pszTrim);
    int iCount = 0;
    while (iIndex != string::npos)
    {
        str.erase(iIndex, 1);
        iIndex = str.find_first_not_of(pszTrim);
        ++iCount;
    }

    return iCount;
}


int32_t UtilString::compareNoCase(const string& str1, const string& str2)
{
    return strcasecmp(str1.c_str(), str2.c_str());
}

int32_t UtilString::nCompareNoCase(const string& str1, const string& str2, size_t count)
{
    return strncasecmp(str1.c_str(), str2.c_str(), count);
}

int UtilString::nCompareCase(const string& str1, const string& str2, size_t count)
{
    return strncmp(str1.c_str(), str2.c_str(), count);
}

void UtilString::splitString(const string &str, const string &sep, bool bStrip, vector<string> &vResult)
{
    string::size_type p = 0, q = 0;
    while (true)
    {
        q = str.find_first_of(sep, p);
        if (q == string::npos)
        {
            if (!bStrip || p < str.size())
            {
                vResult.push_back(str.substr(p));
            }
            break;
        }
        if (q == p)
        {
            if (!bStrip)
            {
                vResult.push_back("");
            }
        }
        else
        {
            vResult.push_back(str.substr(p, q - p));
            p = q;
        }
        ++p;
    }
}

vector<string> UtilString::splitString(const string &str, const string &sep, bool bStrip)
{
    vector<string> vResult;
    splitString(str, sep, bStrip, vResult);
    return vResult;
}

string UtilString::joinString(const vector<string> &v, const string &sep)
{
    string sResult;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (!sResult.empty())
        {
            sResult.append(sep);
        }
        sResult.append(v[i]);
    }
    return sResult;
}

void UtilString::splitString2(const string &str, const string &sep1, const string &sep2, vector<pair<string, string> > &vResult)
{
    vector<string> vStrResult;
    splitString(str, sep1, true, vStrResult);
    
    for (unsigned i = 0; i < vStrResult.size(); ++i)
    {
        vector<string> vStrResult2;
        splitString(vStrResult[i], sep2, true, vStrResult2);
        
        if (vStrResult2.size() != 2)
        {
            continue;
        }
        vResult.push_back(make_pair(vStrResult2[0], vStrResult2[1]));
    }
}

void UtilString::splitString2(const string &str, const string &sep1, const string &sep2, map<string, string> &mResult)
{
    vector<string> vStrResult;
    splitString(str, sep1, true, vStrResult);
    
    for (unsigned i = 0; i < vStrResult.size(); ++i)
    {
        vector<string> vStrResult2;
        splitString(vStrResult[i], sep2, true, vStrResult2);
        
        if (vStrResult2.size() != 2)
        {
            continue;
        }
        mResult[vStrResult2[0]] = vStrResult2[1];
    }
}

map<string, string> UtilString::splitString2(const string &str, const string &sep1, const string &sep2)
{
    map<string, string> mResult;
    splitString2(str, sep1, sep2, mResult);
    return mResult;
}

string UtilString::joinString2(const vector<pair<string, string> > &v, const string &sep1, const string &sep2)
{
    string sResult;
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (!sResult.empty())
        {
            sResult.append(sep1);
        }

        const pair<string, string> &kv = v[i];
        sResult.append(kv.first);
        sResult.append(sep2);
        sResult.append(kv.second);
    }
    return sResult;
}

string UtilString::repeat(const string &str, uint32_t n)
{
    string r;
    for (uint32_t i = 0; i < n; ++i)
    {
        r += str;
    }
    return r;
}

string UtilString::repeat(const string &str, const string &sep, uint32_t n)
{
    string r;
    for (uint32_t i = 0; i < n; ++i)
    {
        if (i != 0)
        {
            r += sep;
        }
        r += str;
    }
    return r;
}

static const string g_sEmptyString;
const string& UtilString::getEmptyString()
{
    return g_sEmptyString;
}


uint64_t UtilString::parseHumanReadableSize(const string &s)
{
    if (s.empty())
    {
        return 0;
    }

    uint64_t iSize = UtilString::strto<uint64_t>(s);
    char c = s[s.size() - 1];
    if (c == 'k' || c == 'K')
    {
        iSize *= 1024;
    }
    else if (c == 'm' || c == 'M')
    {
        iSize *= 1024 * 1024;
    }
    else if (c == 'g' || c == 'G')
    {
        iSize *= 1024 * 1024 * 1024;
    }
    return iSize;
}


string UtilString::replaceVariable(const string &str, const map<string, string> &mVariable, REPACE_MODE mode)
{
    string sResult;
    const char *p = str.c_str(), *q = p + str.size();
    while (p < q)
    {
        const char *t1 = strstr(p, "${");
        if (t1 == NULL)
        {
            sResult.append(p, q);
            break;
        }
        sResult.append(p, t1);

        const char *t2 = strchr(t1, '}');
        if (t2 == NULL)
        {
            sResult.append(t1, q);
            break;
        }

        string sKey = string(t1 + 2, t2);
        map<string, string>::const_iterator it = mVariable.find(sKey);
        if (it != mVariable.end())
        {
            sResult.append(it->second);
        }
        else if (mode == REPLACE_MODE_ALL_MATCH)
        {
            throw std::runtime_error("replaceVariable: " + sKey);
        }
        else if (mode == REPLACE_MODE_KEEP_VARIABLE_ON_MISS)
        {
            sResult.append(t1, t2 + 1);
        }
        else if (mode == REPLACE_MODE_EMPTY_VARIABLE_ON_MISS)
        {
        }
        else if (mode == REPLACE_MODE_EMPTY_RESULT_ON_MISS)
        {
            return "";
        }

        p = t2 + 1;
    }
    return sResult;
}

string UtilString::joinURLParam(const map<string, string> &mParam, const string &sep1, const string &sep2)
{
    return joinString2<string, string>(mParam, sep1, sep2);
}

void UtilString::splitURLParam(const string &sParam, map<string, string> &mParam, const string &sep1, const string &sep2)
{
    mParam = splitString2<string, string>(sParam, sep1, sep2);
}

bool UtilString::splitIni(const string& strLine, string& strName, string& strValue, char chDel)
{
    if (strLine.size() == 0)
    {
        return false;
    }

    size_t iIndex = strLine.find(chDel);
    if (iIndex == string::npos)
    {       
        return false;
    }

    strName = trim(strLine.substr(0, iIndex));
    strValue = trim(strLine.substr(iIndex + 1, strLine.size() - iIndex - 1));

    return true;
}

void UtilString::splitPair(const string& strLine, const string& sep, string& strFirst, string& strSecond)
{
    strFirst = "";
    strSecond = "";

    if (!strLine.empty() && !sep.empty())
    {
        size_t iIndex = strLine.find(sep);
        if (iIndex != string::npos)
        {
            strFirst = trim(strLine.substr(0, iIndex));
            strSecond = trim(strLine.substr(iIndex + 1, strLine.size() - iIndex - 1));
        }
    }
}

string UtilString::str2Hex(const string &sStr)
{
    uint32_t iLen = sStr.length();
    if (iLen == 0)
    {
        return "";
    }

    static const char* const lut = "0123456789ABCDEF";

    string sOutput;
    sOutput.reserve(2 * iLen);
    for (uint32_t i = 0; i < iLen; ++i)
    {
        const unsigned char c = sStr[i];
        sOutput.push_back(lut[c >> 4]);
        sOutput.push_back(lut[c & 15]);
    }

    return sOutput;
}

string UtilString::hex2Str(const string &sHex)
{
    uint32_t iLen = sHex.length();
    if (iLen == 0)
    {
        return "";
    }

    static const char* const lut = "0123456789ABCDEF";

    string sOutput;
    sOutput.reserve(iLen);
    for (uint32_t i = 0; i < iLen; i += 2)
    {
        char a = sHex[i];
        const char* p = lower_bound(lut, lut + 16, a);
        if (*p != a)
        {
            return "";
        }

        char b = sHex[i + 1];
        const char* q = lower_bound(lut, lut + 16, b);
        if (*q != b)
        {
            return "";
        }
        sOutput.push_back(((p - lut) << 4) | (q - lut));
    }
    return sOutput;
}


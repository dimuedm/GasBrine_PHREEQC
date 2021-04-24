#include "UtilOption.h"
#include "UtilString.h"

void UtilOption::parse(int argc, char *argv[])
{
    int i = 1;
    for (; i < argc; ++i)
    {
        string s = argv[i];
        if (s == "--")
        {
            break;
        }
        if (UtilString::isBeginsWith(s, "--"))
        {
            string::size_type pos = s.find('=');
            string sKey = s.substr(2, pos - 2);
            if (pos == string::npos)
            {
                m_mOption[sKey] = "";
            }
            else
            {
                m_mOption[sKey] = s.substr(pos + 1);
            }
        }
        else
        {
            m_vParam.push_back(s);
        }
    }

    for (; i < argc; ++i)
    {
        m_vParam.push_back(argv[i]);
    }
}

bool UtilOption::hasOption(const string &sName) const
{
    return m_mOption.find(sName) != m_mOption.end();
}

string UtilOption::getOption(const string &sName) const
{
    map<string, string>::const_iterator it = m_mOption.find(sName);
    if(it != m_mOption.end())
    {
        return it->second;
    }
    return "";
}

string UtilOption::getOption(const string &sName, const string &sDefault) const
{
    map<string, string>::const_iterator it = m_mOption.find(sName);
    if(it != m_mOption.end())
    {
        return it->second;
    }
    return sDefault;
}


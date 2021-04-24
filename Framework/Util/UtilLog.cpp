#include "UtilLog.h"

// app log
INITIALIZE_EASYLOGGINGPP

class EasyLogInitialize
{
public:
    EasyLogInitialize()
    {
 /*       el::Configurations defaultConf;
        defaultConf.setToDefault();
        defaultConf.setGlobally(el::ConfigurationType::Format, "%datetime|%level|%msg");
        defaultConf.setGlobally(el::ConfigurationType::Filename, "logs//ReSoc_%datetime{%Y%M%d}.log");

        el::Loggers::reconfigureLogger("default", defaultConf);*/
 //       LOG(INFO) << "Initialize Log System OK!";
    }

    virtual ~EasyLogInitialize() { }
};

static EasyLogInitialize g_EasyLogInitialize;

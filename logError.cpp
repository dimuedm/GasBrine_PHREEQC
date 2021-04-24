#include "Framework/Util/UtilLog.h"
#include "logError.h"

logError::logError()
{
}

logError::~logError()
{
}

void logError::
LOGERROR(const std::string a_error)
{
    LOG(ERROR) << a_error << std::endl;
}
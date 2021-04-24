#include "UtilFile.h"

std::string UtilFile::getFilePath(const std::string& sFile)
{
#ifdef _DEBUG
	return "../../../Data/" + sFile;
#else
    return "../../../Data/" + sFile;
#endif
}

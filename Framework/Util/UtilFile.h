/*
* Header file for UtilFile.cpp
* Written 05/2017 by Vic
* util functions of file operations
*/

#ifndef FRAMEWORK_UTIL_UTIL_FILE_H_
#define FRAMEWORK_UTIL_UTIL_FILE_H_

#include <string>

class UtilFile
{
public:
	static std::string getFilePath(const std::string& sFile);
};

#define REAL_FILE_PATH(sFile) UtilFile::getFilePath(sFile)

#endif  // FRAMEWORK_UTIL_UTIL_FILE_H_


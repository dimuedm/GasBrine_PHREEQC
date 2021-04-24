#ifndef FRAMEWORK_COMMON_BASE_H_
#define FRAMEWORK_COMMON_BASE_H_

#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <stdint.h>
#include <stdarg.h>

#include <vector>
#include <string>
#include <list>
#include <set>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <limits>
#include <cmath>

using namespace std;

#ifdef _WIN32
// string functions
#define strcasecmp   _stricmp
#define strncasecmp  _strnicmp
#define snprintf     _snprintf
#endif

#include "Common/Macro.h"

#endif // DD_FRAMEWORK_COMMON_DD_FRAMEWORK_DEFINE__H
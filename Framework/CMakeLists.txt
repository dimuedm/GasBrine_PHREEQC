cmake_minimum_required(VERSION 2.8)

# target lib name
set(PROJECTNAME Framework)

# set the path to save the target lib file
set(LIBRARY_OUTPUT_PATH ${LIB})

# set the path of include files
include_directories(. ../ ${INC_THIRDPARTY})

# set the path of lib files
link_directories(${LIB_THIRDPARTY})

# collect the fiels target lib needed, and save to the variable : SOURCE_FILES
FILE(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h" "*.cpp" "*.cc" )

# add the target lib : framework
add_library(${PROJECTNAME} ${SOURCE_FILES})

# the libs target lib needed(Please pay attention to the order of the libs)
target_link_libraries(${PROJECTNAME})

# to create the folders for visual studio
GroupSources(${CMAKE_CURRENT_SOURCE_DIR})
SET_PROPERTY(TARGET ${PROJECTNAME} PROPERTY FOLDER ${ROOTPROJECT})

# install inlcude files
CopyFiles(${CMAKE_CURRENT_SOURCE_DIR})
CopyFiles(${CMAKE_CURRENT_SOURCE_DIR}/Common)
CopyFiles(${CMAKE_CURRENT_SOURCE_DIR}/Util)
CopyFiles(${CMAKE_CURRENT_SOURCE_DIR}/Util/Impl)





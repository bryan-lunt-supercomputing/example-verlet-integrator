cmake_minimum_required (VERSION 3.5)
project(mp3gas VERSION 1.0.0 LANGUAGES C)

# Location of additional cmake modules
#

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake
    )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_DATA_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/data)


#include_directories("${PROJECT_SOURCE_DIR}/gassim/include")

#basic library with things we will need for this project
add_library(gassim "gassim/src/gassim.c"
					"gassim/include/gassim.h"
					"gassim/src/gassim_util.c"
					"gassim/include/gassim_util.h")
include_directories(gassim "${PROJECT_SOURCE_DIR}/gassim/include")

add_library(gassim2d "gassim2d/src/gassim2d.c"
					"gassim2d/include/gassim2d.h"
					"gassim2d/src/gassim2d_util.c"
					"gassim2d/include/gassim2d_util.h")
include_directories(gassim2d "${PROJECT_SOURCE_DIR}/gassim2d/include")

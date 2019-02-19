cmake_minimum_required (VERSION 3.5)
project(gassim VERSION 1.0.0 LANGUAGES C)

# Location of additional cmake modules
#

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake
    )

set(CMAKE_C_FLAGS "-O3")
set(CMAKE_C_FLAGS_DEBUG "-g -Wall -Wextra")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_DATA_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/data)


#include_directories("${PROJECT_SOURCE_DIR}/gassim/include")

#basic library with things we will need for this project
add_library(gassim STATIC "${PROJECT_SOURCE_DIR}/src/gassim.c"
					"${PROJECT_SOURCE_DIR}/include/gassim.h"
					"${PROJECT_SOURCE_DIR}/src/gassim_util.c"
					"${PROJECT_SOURCE_DIR}/include/gassim_util.h")
include_directories(gassim INTERFACE "${PROJECT_SOURCE_DIR}/include")

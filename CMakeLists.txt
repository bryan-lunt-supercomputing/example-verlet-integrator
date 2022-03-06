cmake_minimum_required (VERSION 3.5)
project(example-verlet-integrator VERSION 1.0.0 LANGUAGES C)

# Location of additional cmake modules
#

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/cmake
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_DATA_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/data)

set(CMAKE_C_CLANG_TIDY "clang-tidy;-checks=*")

add_subdirectory("gassim")

add_subdirectory("gassim2d")

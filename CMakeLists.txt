cmake_minimum_required (VERSION 3.5)

project (openDVC)

# specify configuration file for compile-time constants
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
configure_file(config.h.in ${CMAKE_BINARY_DIR}/generated/config.h)
include_directories(${CMAKE_BINARY_DIR}/generated/)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_subdirectory (src)
#add_subdirectory (tests)


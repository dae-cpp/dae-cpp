cmake_minimum_required(VERSION 3.10)

project(dae-cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include ${BOOST_INCLUDE})

#FILE(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

#add_executable(${PROJECT_NAME} ${SOURCES})


#### Build tests ####

set(PROJECT_TEST "dae-cpp-test")

include(CTest)
enable_testing()

FILE(GLOB SOURCES_TEST ${PROJECT_SOURCE_DIR}/tests/*.cpp)

add_executable(${PROJECT_TEST} ${SOURCES_TEST})

include_directories("${PROJECT_SOURCE_DIR}/googletest/build/include" "${PROJECT_SOURCE_DIR}")

link_directories("${PROJECT_SOURCE_DIR}/googletest/build/lib")

find_library(GTEST_MAIN gtest_main HINTS ${PROJECT_SOURCE_DIR}/googletest/build/lib)
find_library(GTEST gtest HINTS ${PROJECT_SOURCE_DIR}/googletest/build/lib)

if(EXISTS "${GTEST}" AND EXISTS "${GTEST_MAIN}")
    message(STATUS "Found GoogleTest: ${GTEST}")
else(EXISTS "${GTEST}" AND EXISTS "${GTEST_MAIN}")
    message(FATAL_ERROR "Could not find GoogleTest library")
endif(EXISTS "${GTEST}" AND EXISTS "${GTEST_MAIN}")

target_link_libraries(${PROJECT_TEST} ${GTEST_MAIN} ${GTEST} pthread)

target_compile_definitions(${PROJECT_TEST} PRIVATE "TESTING")

include(GoogleTest)
gtest_discover_tests(${PROJECT_TEST})

# file(COPY ${PROJECT_SOURCE_DIR}/examples DESTINATION ${CMAKE_BINARY_DIR})


#### Install ####

# install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(TARGETS ${PROJECT_TEST} DESTINATION bin)

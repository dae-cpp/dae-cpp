cmake_minimum_required(VERSION 3.13)  # To match googletest requirement

project(dae-cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# TODO: Add examples here

#include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include ${BOOST_INCLUDE})

#FILE(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

#add_executable(${PROJECT_NAME} ${SOURCES})


#### Build tests ####

set(PROJECT_TEST "dae-cpp-test")

include(CTest)
enable_testing()

FILE(GLOB SOURCES_TEST ${PROJECT_SOURCE_DIR}/tests/*.cpp)

add_executable(${PROJECT_TEST} ${SOURCES_TEST})

set(BUILD_GMOCK OFF CACHE BOOL "Builds the googlemock subproject")
set(INSTALL_GTEST OFF CACHE BOOL "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)")

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(googletest)

include_directories("${PROJECT_SOURCE_DIR}/googletest/googletest/include" "${PROJECT_SOURCE_DIR}")

target_link_libraries(${PROJECT_TEST} gtest_main gtest)

target_compile_definitions(${PROJECT_TEST} PRIVATE "TESTING")

include(GoogleTest)
gtest_discover_tests(${PROJECT_TEST})


#### Install ####

# install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(TARGETS ${PROJECT_TEST} DESTINATION bin)

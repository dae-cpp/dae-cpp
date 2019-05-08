cmake_minimum_required(VERSION 2.8)

project(dae-cpp)

enable_language(CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(DAE_FORTRAN_STYLE "Matrices will be defined using FORTRAN style (the first index starts from 1 instead of 0)" OFF)
option(DAE_LONG_INT "Use long integer representation for huge systems (more than ~10^7 equations)" OFF)
option(DAE_SINGLE "Use single precision in the solver)" OFF)
option(DAE_BUILD_EXAMPLES "Build all examples" ON)
option(DAE_TEST "Execute all examples and test the solver" OFF)

set(DAE_MKL_DIR "/opt/intel/mkl" CACHE STRING "Path to Intel MKL library")

set(MKL_INCLUDE "${DAE_MKL_DIR}/include")
set(EXT_INCLUDE "${PROJECT_SOURCE_DIR}/src/external")

set(MKL_LIB_PATH "${DAE_MKL_DIR}/lib/intel64")

if(DAE_LONG_INT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DMKL_ILP64 -m64 -fopenmp -Wall -Wl,--no-as-needed")
else(DAE_LONG_INT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -m64 -fopenmp -Wall -Wl,--no-as-needed")
endif(DAE_LONG_INT)

include_directories(${MKL_INCLUDE} ${EXT_INCLUDE})
link_directories(${MKL_LIB_PATH})

include_directories("${PROJECT_SOURCE_DIR}/src")
configure_file(${PROJECT_SOURCE_DIR}/src/cmake_config.h.in ${PROJECT_SOURCE_DIR}/src/cmake_config.h)
add_subdirectory(src)

if(DAE_BUILD_EXAMPLES)

  set(EXAMPLE_LIST "perovskite")  # "a;b;c"
  
  #if(DAE_LONG_INT_CHANGED)
  #  set(INTEL_LP64_LIB "LIB_MKL_INTEL_LP64-NOTFOUND")
  #endif(AE_LONG_INT_CHANGED)

  if(DAE_LONG_INT)
    find_library(LIB_MKL_INTEL_LP64 mkl_intel_ilp64 HINTS ${MKL_LIB_PATH})
  else(DAE_LONG_INT)
    find_library(LIB_MKL_INTEL_LP64 mkl_intel_lp64 HINTS ${MKL_LIB_PATH})
  endif(DAE_LONG_INT)
  
  find_library(LIB_MKL_GNU_THREAD mkl_gnu_thread HINTS ${MKL_LIB_PATH})
  find_library(LIB_MKL_CORE mkl_core HINTS ${MKL_LIB_PATH})
  find_library(LIB_PTHREAD pthread)
  find_library(LIB_MATH m)
  find_library(LIB_DL dl)
  
  foreach(EXAMPLE_NAME ${EXAMPLE_LIST})
  
    FILE(GLOB SOURCES ${PROJECT_SOURCE_DIR}/examples/${EXAMPLE_NAME}/*.cpp)

    add_executable(${EXAMPLE_NAME} ${SOURCES})

    target_link_libraries(${EXAMPLE_NAME} daecpp_static)
    target_link_libraries(${EXAMPLE_NAME} ${LIB_MKL_INTEL_LP64} ${LIB_MKL_GNU_THREAD} ${LIB_MKL_CORE} ${LIB_GOMP} ${LIB_PTHREAD} ${LIB_MATH} ${LIB_DL})

    install(TARGETS ${EXAMPLE_NAME} DESTINATION bin)
    
  endforeach()
  
  if(DAE_TEST)
  
    include(CTest)
    
    foreach(EXAMPLE_NAME ${EXAMPLE_LIST})
      add_test(${EXAMPLE_NAME}_TEST ${EXAMPLE_NAME})
    endforeach()
    
  endif(DAE_TEST)
      
endif(DAE_BUILD_EXAMPLES)

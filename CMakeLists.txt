cmake_minimum_required(VERSION 2.8)

project(dae-cpp)

enable_language(CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(DAE_FORTRAN_STYLE "Matrices will be defined using FORTRAN style (one-based indexing of columns and rows)" OFF)
option(DAE_LONG_INT "Use long integer representation for huge systems (more than ~10^7 equations)" OFF)
option(DAE_SINGLE "Use single precision in the solver)" OFF)
option(DAE_TEST "Build automatic solver test" ON)

if(DAE_TEST)
  unset(DAE_BUILD_EXAMPLES CACHE)
endif(DAE_TEST)

option(DAE_BUILD_EXAMPLES "Build all examples" ON)

set(DAE_MKL_DIR "/opt/intel/mkl" CACHE STRING "Path to Intel MKL root")

unset(DAE_MKL_INCLUDE CACHE)
unset(DAE_MKL_LIB CACHE)

set(DAE_MKL_INCLUDE "${DAE_MKL_DIR}/include" CACHE STRING "Path to Intel MKL includes")
set(DAE_MKL_LIB "${DAE_MKL_DIR}/lib/intel64" CACHE STRING "Path to Intel MKL libraries")

set(EXT_INCLUDE "${PROJECT_SOURCE_DIR}/src/external")

unset(CMAKE_CXX_FLAGS CACHE)

if(DAE_LONG_INT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DMKL_ILP64 -m64 -fopenmp -Wall -Wl,--no-as-needed")
else(DAE_LONG_INT)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -m64 -fopenmp -Wall -Wl,--no-as-needed")
endif(DAE_LONG_INT)

include_directories(${DAE_MKL_INCLUDE} ${EXT_INCLUDE})
link_directories(${DAE_MKL_LIB})

include_directories("${PROJECT_SOURCE_DIR}/src")
configure_file(${PROJECT_SOURCE_DIR}/src/cmake_config.h.in ${PROJECT_SOURCE_DIR}/src/cmake_config.h)
add_subdirectory(src)

if(DAE_BUILD_EXAMPLES)

  set(EXAMPLE_LIST "perovskite")  # "a;b;c"
  
  unset(LIB_MKL_INTEL_LP64 CACHE)

  if(DAE_LONG_INT)
    find_library(LIB_MKL_INTEL_LP64 mkl_intel_ilp64 HINTS ${DAE_MKL_LIB})
  else(DAE_LONG_INT)
    find_library(LIB_MKL_INTEL_LP64 mkl_intel_lp64 HINTS ${DAE_MKL_LIB})
  endif(DAE_LONG_INT)
  
  find_library(LIB_MKL_GNU_THREAD mkl_gnu_thread HINTS ${DAE_MKL_LIB})
  find_library(LIB_MKL_CORE mkl_core HINTS ${DAE_MKL_LIB})
  find_library(LIB_PTHREAD pthread)
  find_library(LIB_MATH m)
  find_library(LIB_DL dl)
  
  foreach(EXAMPLE_NAME ${EXAMPLE_LIST})
  
    FILE(GLOB SOURCES ${PROJECT_SOURCE_DIR}/examples/${EXAMPLE_NAME}/*.cpp)

    add_executable(${EXAMPLE_NAME} ${SOURCES})

    target_link_libraries(${EXAMPLE_NAME} daecpp_static)
    target_link_libraries(${EXAMPLE_NAME} ${LIB_MKL_INTEL_LP64} ${LIB_MKL_GNU_THREAD} ${LIB_MKL_CORE} ${LIB_PTHREAD} ${LIB_MATH} ${LIB_DL})

    install(TARGETS ${EXAMPLE_NAME} DESTINATION bin)
    
  endforeach()
  
  if(DAE_TEST)
  
    include(CTest)
    
    foreach(EXAMPLE_NAME ${EXAMPLE_LIST})
      add_test(${EXAMPLE_NAME}_TEST ${EXAMPLE_NAME})
    endforeach()
    
  endif(DAE_TEST)
      
endif(DAE_BUILD_EXAMPLES)

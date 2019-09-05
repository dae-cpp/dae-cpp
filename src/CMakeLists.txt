FILE(GLOB DAE_CPP_SOURCES *.cpp)
FILE(GLOB DAE_CPP_INCLUDES *.h)
FILE(GLOB DAE_CPP_INCLUDES_PLOTTING external/matplotlib-cpp/*.h)

add_library(daecpp SHARED ${DAE_CPP_SOURCES})
add_library(daecpp_static STATIC ${DAE_CPP_SOURCES})

install(TARGETS daecpp daecpp_static DESTINATION lib)
install(FILES ${DAE_CPP_INCLUDES} DESTINATION include)
install(FILES ${DAE_CPP_INCLUDES_PLOTTING} DESTINATION include/external/matplotlib-cpp)

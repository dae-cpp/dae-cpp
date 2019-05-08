FILE(GLOB DAE_CPP_SOURCES *.cpp)

add_library(daecpp SHARED ${DAE_CPP_SOURCES})
add_library(daecpp_static STATIC ${DAE_CPP_SOURCES})

install(TARGETS daecpp daecpp_static DESTINATION lib)
install(FILES solver.h mass_matrix.h jacobian.h RHS.h DESTINATION include)

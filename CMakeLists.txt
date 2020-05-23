cmake_minimum_required(VERSION 3.9)
project(Openbps)

set(CMAKE_CXX_STANDARD 17)

#===============================================================================
# Command line options
#===============================================================================

option(optimize "Turn on O3 flag"        OFF)

#find_package(GTest REQUIRED)

include_directories(src)
include_directories(test)
include_directories(build)
include_directories(extern)

#add_library(pugixml extern/pugiData/pugixml.cpp)
add_subdirectory(extern/xtl)
add_subdirectory(extern/xtensor)
add_subdirectory(extern/xtensor-blas)
target_link_libraries(xtensor INTERFACE xtl)
target_link_libraries(xtensor-blas INTERFACE xtensor)

add_subdirectory(lib/googletest)

add_definitions(-DHAVE_CBLAS=1)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "BLAS VENDOR:    " ${BLA_VENDOR})
message(STATUS "BLAS LIBRARIES: " ${BLAS_LIBRARIES})

#Binary build
#Check compiler options
if(optimize)
  list(APPEND cxxflags -O3)
endif()
# Show flags being used
message(STATUS "C++ flags: ${cxxflags}")

list(APPEND lib_SOURCES
      "../src/chain.h"
      "../src/chain.cpp"
      "../src/parse.h"
      "../src/configure.h"
      "../src/functionals.h"
      "../src/reactions.h"
      "../src/executer.h"
      "../src/configure.cpp"
       "../src/functionals.cpp"
       "../src/reactions.cpp"
       "../src/executer.cpp"
       "../src/parse.cpp"
       "../src/form_output.cpp"
       "../src/materials.cpp"
       "../src/materials.h"
       "../extern/pugiData/pugixml.h" 
       "../extern/pugiData/pugixml.cpp")

add_library(openbps_lib STATIC ${lib_SOURCES})
# Set compile flags
target_compile_options(openbps_lib PRIVATE ${cxxflags})
# Link library
target_link_libraries(openbps_lib ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} Eigen3::Eigen xtensor xtensor-blas)


add_executable(openbps "../src/main.cpp" )
# Set compile flags
target_compile_options(openbps PRIVATE ${cxxflags})
# Link executable
target_link_libraries(openbps PUBLIC openbps_lib)

#Test build

add_executable(o_tst  "../test/test.cpp"
                      "../test/main_test.cpp"
                      "../lib/googletest/googletest/include/gtest/gtest.h"
                      ${lib_SOURCES})

add_test(NAME o_tst COMMAND o_tst)
target_link_libraries(o_tst PUBLIC openbps_lib gtest)





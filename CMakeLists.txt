cmake_minimum_required(VERSION 3.9)
project(Openbps)

set(CMAKE_CXX_STANDARD 17)

#find_package(GTest REQUIRED)

include_directories(src)
include_directories(test)
include_directories(build)
include_directories(extern)

#add_library(pugixml extern/pugiData/pugixml.cpp)
add_subdirectory(extern/xtl)
add_subdirectory(extern/xtensor)
target_link_libraries(xtensor INTERFACE xtl)

add_subdirectory(lib/googletest)

#Binary build

list(APPEND lib_SOURCES
     "../src/nuclide_class.h"
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
       "../extern/pugiData/pugixml.h" 
       "../extern/pugiData/pugixml.cpp")

add_library(openbps_lib STATIC ${lib_SOURCES})
target_link_libraries(openbps_lib xtensor)

add_executable(openbps "../src/main.cpp" )
target_link_libraries(openbps PUBLIC openbps_lib)

#Test build

add_executable(o_tst  "../test/test.cpp"
                      "../test/main_test.cpp"
                      "../lib/googletest/googletest/include/gtest/gtest.h"
                      ${lib_SOURCES})

add_test(NAME o_tst COMMAND o_tst)
target_link_libraries(o_tst PUBLIC openbps_lib gtest)






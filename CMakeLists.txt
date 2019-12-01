cmake_minimum_required(VERSION 3.9)
project(Openbps)

set(CMAKE_CXX_STANDARD 17)

#find_package(GTest REQUIRED)

include_directories(src)
include_directories(test)
include_directories(build)
include_directories(extern)

#add_library(pugixml extern/pugiData/pugixml.cpp)

add_subdirectory(lib/googletest)

#Binary build

add_executable(openbps "../src/nuclide_class.h"
                       "../src/parse.cpp"
                       "../src/form_output.cpp" 
                       "../src/main.cpp" 
                       "../extern/pugiData/pugixml.h" 
                       "../extern/pugiData/pugixml.cpp")
add_library(openbps_lib STATIC "../src/nuclide_class.h" 
                               "../src/parse.cpp"
                               "../src/form_output.cpp" 
                               "../src/main.cpp"
                               "../extern/pugiData/pugixml.h" 
                               "../extern/pugiData/pugixml.cpp")

#Test build

add_executable(o_tst  "../src/nuclide_class.h" 
                      "../src/parse.cpp"
                      "../src/form_output.cpp" 
                      "../test/test.cpp"
                      "../test/main_test.cpp"
                      "../extern/pugiData/pugixml.h" 
                      "../extern/pugiData/pugixml.cpp" 
                      "../lib/googletest/googletest/include/gtest/gtest.h")

add_test(NAME o_tst COMMAND o_tst)
target_link_libraries(o_tst PUBLIC openbps_lib gtest)






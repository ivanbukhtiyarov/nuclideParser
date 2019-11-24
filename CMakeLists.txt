cmake_minimum_required(VERSION 3.9)
project(Openbps)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)
include_directories(test)
include_directories(build)
include_directories(extern)

#add_library(pugixml extern/pugiData/pugixml.cpp)

add_subdirectory(build)
add_subdirectory(build-test)
add_subdirectory(lib/googletest)






cmake_minimum_required(VERSION 3.9)
project(ibrae2)

set(CMAKE_CXX_STANDARD 17)

add_executable(ibrae2 main.cpp nuclide_class.h pugiData/pugixml.cpp pugiData/pugixml.h pugiData/pugiconfig.h parse.cpp)
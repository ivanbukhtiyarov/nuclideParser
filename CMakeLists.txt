cmake_minimum_required(VERSION 3.9)
project(ibrae2)

set(CMAKE_CXX_STANDARD 17)

add_executable(ibrae2 src/main.cpp src/nuclide_class.h extern/pugiData/pugixml.cpp extern/pugiData/pugixml.h extern/pugiData/pugiconfig.h src/parse.cpp src/form_output.cpp)
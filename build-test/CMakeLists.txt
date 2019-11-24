set(BINARY ${CMAKE_PROJECT_NAME}_tst)
add_executable(o_tst  "../src/nuclide_class.h" "../src/parse.cpp"
               "../src/form_output.cpp" "../test/test.cpp" "../test/main_test.cpp"
               "../extern/pugiData/pugixml.h" "../extern/pugiData/pugixml.cpp" "/usr/include/gtest/gtest.h")

add_test(NAME o_tst COMMAND o_tst)
target_link_libraries(o_tst PUBLIC openbps_lib gtest)

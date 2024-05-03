add_library(f2_helper_tests OBJECT 
    f2_helper.test.cpp 
    f2_helper.cpp
)

target_link_libraries(f2_helper_tests PRIVATE Catch2::Catch2WithMain)
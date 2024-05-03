add_library(f2_helper_tests OBJECT f2_helper.test.cpp )

target_link_libraries(f2_helper_tests PRIVATE Catch2::Catch2WithMain)
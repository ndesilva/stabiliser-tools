add_library(test_sources OBJECT f2_helper.test.cpp)
target_link_libraries(test_sources PRIVATE Catch2::Catch2WithMain)
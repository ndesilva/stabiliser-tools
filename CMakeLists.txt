cmake_minimum_required(VERSION 3.10)

set(VCPKG_TARGET_TRIPLET x64-mingw-static)

project(FastStabiliser)

find_package(Catch2 3 REQUIRED)

add_subdirectory("c++ versions")

add_executable(tests)
target_link_libraries(tests test_sources)
target_link_libraries(tests Catch2::Catch2WithMain)

include(CTest)
include(Catch)

catch_discover_tests(tests)
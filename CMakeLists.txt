cmake_minimum_required(VERSION 3.13)

set(VCPKG_TARGET_TRIPLET x64-mingw-static)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow")

project(FastStabiliser)

find_package(Catch2 3 REQUIRED)

include_directories("./src")
add_subdirectory("src")
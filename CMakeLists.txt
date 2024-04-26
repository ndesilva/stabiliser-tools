cmake_minimum_required(VERSION 3.13)

set(VCPKG_TARGET_TRIPLET x64-mingw-static)

if(POLICY CMP0079 )
  cmake_policy(SET CMP0079 NEW)
endif()

project(FastStabiliser)

find_package(Catch2 3 REQUIRED)

include_directories("./src")
add_subdirectory("src")
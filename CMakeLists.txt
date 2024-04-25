cmake_minimum_required(VERSION 3.10)

set(VCPKG_TARGET_TRIPLET x64-mingw-static)

project(HelloWorld)

find_package(fmt CONFIG REQUIRED)

add_executable(HelloWorld helloworld.cpp)

target_link_libraries(HelloWorld PRIVATE fmt::fmt)
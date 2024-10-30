cmake_minimum_required( VERSION 3.28 )

set ( CMAKE_COMPILE_WARNING_AS_ERROR ON )
set ( CMAKE_CXX_STANDARD 23 )

project( FastStabiliser CXX )

# This currently doesnt work, because by the time this if-else statement is reached
# vpckg --install has already been run by the project() call, so the vcpkg triplet is already
# wrong. Unfortunately, it seems that the MSVC and MINGW flags are set by the project() call
# so that moving project() beneath this if-else statement means both branches are not run.
# As a workaround for MINGW, one can paste the MINGW branch above the project() command, and then
# compile once, causing vcpkg to install with the correct triplet. As long as the
# vcpkg_installed/x64-mingw-static folder is not subsequently deleted, or the vcpkg dependancies change,
# the additional set() commands can then be deleted - future builds will use the wrong triplet for vcpkg
# --install, but the correct installed folder already exists.
if ( MSVC )
	add_compile_options( /W4 /Zc:__cplusplus /Zc:enumTypes /Zc:inline /Zc:templateScope /Oi )
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
	message("Compiler detected as MSVC")
elseif( MINGW ) 
	add_compile_options( -Wall -Wextra )
	set( VCPKG_TARGET_TRIPLET x64-mingw-static )
	set( DVCPKG_TARGET_TRIPLET x64-mingw-static ) 
	message("Compiler detected as MINGW")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    message( "64 bit architecture" )
elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
	message( "32 bit architecture" )
endif()

find_package( Catch2 3 REQUIRED )
find_package( Python3 COMPONENTS Interpreter Development REQUIRED )
find_package( pybind11 REQUIRED CONFIG )

include( CTest )
add_subdirectory( "cpp" )
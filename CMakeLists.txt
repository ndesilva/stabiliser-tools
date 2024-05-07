cmake_minimum_required( VERSION 3.28 )
project( FastStabiliser CXX )

set( CMAKE_COMPILE_WARNING_AS_ERROR ON )
set( CMAKE_CXX_STANDARD 23 )

if ( MSVC )
	add_compile_options( /W4 /Zc:__cplusplus /Zc:enumTypes /Zc:inline /Zc:templateScope /Oi )
else()
	add_compile_options( -Wall -Wextra )
	set( VCPKG_TARGET_TRIPLET x64-mingw-static )
endif()

include( CTest )
add_subdirectory( "src" )

cmake_minimum_required( VERSION 3.28 )
project( FastStabiliser CXX )

set( CMAKE_COMPILE_WARNING_AS_ERROR ON )
set( CMAKE_CXX_STANDARD 20 )

if ( MSVC )
	add_compile_options( /W4 /Zc:__cplusplus /Zc:enumTypes /Zc:inline /Zc:templateScope /Oi )
else()
	add_compile_options( -Wall -Wextra -Wshadow )
	set( VCPKG_TARGET_TRIPLET x64-mingw-static )
endif()

find_package( Catch2 3 REQUIRED )

include_directories( "./src" )
add_subdirectory( "src" )

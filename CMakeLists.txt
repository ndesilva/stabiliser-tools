cmake_minimum_required( VERSION 3.28 )

set( CMAKE_COMPILE_WARNING_AS_ERROR ON )
set( CMAKE_CXX_STANDARD 23 )

if ( MSVC )
	add_compile_options( /W4 /Zc:__cplusplus /Zc:enumTypes /Zc:inline /Zc:templateScope /Oi )
else()
	add_compile_options( -Wall -Wextra )
	set( VCPKG_TARGET_TRIPLET x64-mingw-static )
	set( DVCPKG_TARGET_TRIPLET x64-mingw-static )
	message("here")
endif()

project( FastStabiliser CXX )

message("now here")
find_package( Catch2 3 REQUIRED )

include_directories( "./src" )
add_subdirectory( "src" )

#pragma once

#define STRINGIFY_X( X ) #X
#define STRINGIFY( X ) STRINGIFY_X( X )

#ifdef _MSC_VER
#define MSVC_PUSH_AND_DISABLE_WARNINGS( ... ) \
	_Pragma( "warning( push )" ) \
	_Pragma( STRINGIFY( warning( disable : ##__VA_ARGS__ ## ) ) )

#define MSVC_POP_WARNINGS _Pragma( "warning( pop )" )

#else
#define MSVC_PUSH_AND_DISABLE_WARNINGS( ... )
#define MSVC_POP_WARNINGS
#endif

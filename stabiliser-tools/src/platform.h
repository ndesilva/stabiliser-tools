#pragma once

#include <version>

#define STRINGIFY_X(X) #X
#define STRINGIFY(X) STRINGIFY_X(X)

#ifdef _MSC_VER
#define MSVC_PUSH_AND_DISABLE_WARNINGS(...) \
	_Pragma("warning( push )")              \
		_Pragma(STRINGIFY(warning(disable :##__VA_ARGS__##)))

#define MSVC_POP_WARNINGS _Pragma("warning( pop )")

#else
#define MSVC_PUSH_AND_DISABLE_WARNINGS(...)
#define MSVC_POP_WARNINGS
#endif

#ifdef __GNUC__
#define GCC_PUSH_AND_DISABLE_WARNINGS(...) \
	_Pragma("GCC diagnostic push")         \
		_Pragma(STRINGIFY(GCC warning ignored##__VA_ARGS__))

#define GCC_POP_WARNINGS _Pragma("GCC warning pop")
#else
#define GCC_PUSH_AND_DISABLE_WARNINGS(...)
#define GCC_POP_WARNINGS
#endif

#ifndef __cpp_size_t_suffix
// Added in C++23, we can mimic it ourselves
// Disable warning for user defined literal not starting with _, we know this one is safe since we compile with C++20
MSVC_PUSH_AND_DISABLE_WARNINGS(4455)
MSVC_PUSH_AND_DISABLE_WARNINGS(4244)
GCC_PUSH_AND_DISABLE_WARNINGS("-Wno-c++23-extensions")
constexpr std::size_t operator"" uz(unsigned long long n)
{
	return n;
}
MSVC_POP_WARNINGS
GCC_POP_WARNINGS
#endif

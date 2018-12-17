#ifndef TARGET_H
#define TARGET_H

// Android target-specific headers and macros

#ifdef ANDROID_NDK

#include <android/log.h>

#define  LOG_TAG    "Satellite Safari"
#define  LOG_INFO(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)
#define  LOG_ERROR(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)

#endif // ANDROID_NDK

// Windows target-specific headers and macros

#ifdef WIN32
#include <float.h>
#define isinf(x)	(!_finite(x))
#define isnan(x)	_isnan(x)
#define snprintf	sprintf_s
#define strcasecmp	_strcmpi
#define strncasecmp	_strnicmp
#define strlcpy		strncpy
#ifdef ASTROLIBDLL_EXPORTS
#define ASTROLIBDLL_API __declspec(dllexport)
#else
#define ASTROLIBDLL_API __declspec(dllimport)
#endif
#else
#define ASTROLIBDLL_API
#endif	// WIN32

#ifdef LINUX
#define strlcpy strncpy
#define strlcat strncat
#endif

#endif	// TARGET_H

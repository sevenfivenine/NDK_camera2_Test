#include <jni.h>
#include <string>
#include <math.h>
#include <iostream>
#include <android/log.h>
#include "AstroLib.h"
#include "Astrometry.h"

#define  LOG_TAG    "NDKCameraTest"
#define  LOG_INFO(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)
#define  LOG_ERROR(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)

#define STARRADIUS		0.003			// radius to use for star flux measurement (fraction of image width)
#define BACKRADIUS		0.01			// radius to use for background measurement (fraction of image width)
#define PEAKSIGMA		5.0				// peak pixel must be this many standard deviations above background.
#define MEANSIGMA		0.1				// mean flux must be this many standard deviations above background.
#define MAX_MAJ_AXIS	INFINITY		// maximum allowed object major axis in pixels
#define MAX_AXIS_RATIO	1.5				// maximum allowed ratio between object major and minor axes.
#define MAX_OBJECTS		50				// maximum number of objects to indentify (brightest N in image).
#define MIN_REF_ALT		DEG_TO_RAD(10)	// minimun altitude above horizon for references.
#define MIN_OBJ_ALT		DEG_TO_RAD(10)	// minimum altitude above horizon for objects.
#define TIMEOUT			60				// timeout for plate-solving routine in seconds; 0 = no timeout.
#define REF_MAG_LIMIT	4.5				// visual magnitude limit for reference objecrs.

extern "C" JNIEXPORT jstring

JNICALL
Java_com_darkfuturestudios_ndktest_MainActivity_stringFromJNI(
        JNIEnv *env,
        jobject /* this */) {
    std::string hello = "Hello from C++";
    return env->NewStringUTF(hello.c_str());
}

extern "C" JNIEXPORT jstring

JNICALL
Java_com_darkfuturestudios_ndktest_MainActivity_computeOrbitParams(
        JNIEnv *env,
        jobject,
        jfloat a,
        jfloat e) {
    float rp = (1 - e)*a;
    float ra = (1 + e)*a;
    double mu = 3.986*pow(10, 14);
    double period = 2*M_PI*sqrt(pow(a, 3)/mu);
    std::string rpString = std::to_string(rp);
    std::string raString = std::to_string(ra);
    std::string periodString = std::to_string(period);

    std::string outputString = "Perigee: " + rpString + " meters\nApogee: " + raString
                               + " meters\nOrbital Period: " + periodString + " seconds";

    return env->NewStringUTF(outputString.c_str());

}

extern "C" JNIEXPORT jint

JNICALL
Java_com_darkfuturestudios_ndktest_MainActivity_processImageBuffer(
        JNIEnv *env,
        jobject,
        jintArray pixelArray,
        jint width,
        jint height,
        jfloat widthAngle,
        jfloat heightAngle)
{
    jsize len = env->GetArrayLength ( pixelArray );
    uint32_t *pixels = (uint32_t *) env->GetIntArrayElements ( pixelArray, 0 );

    double sumR = 0.0, sumG = 0.0, sumB = 0.0, sumA = 0;
    uint8_t r = 0, g = 0, b = 0 , a = 0;

    for ( int i = 0; i < len; i++ )
    {
        r = pixels[i] & 0x000000ff;
        g = ( pixels[i] >> 8 ) & 0x000000ff;
        b = ( pixels[i] >> 16 ) & 0x000000ff;
        a = ( pixels[i] >> 24 ) & 0x000000ff;

        sumR += r;
        sumG += g;
        sumB += b;
        sumA += a;

        pixels[i] = r | ( (int) g << 8 ) | ( (int) b << 16 ) | ( (int) a << 24 );
    }

    // Initialize sky image struct with data from the camera

    A3Image skyImage = { 0 };
    double jd = AACurrentUTC();
    double lon = DEG_TO_RAD ( -122.0 ); // hard-coded approximate longitude of San Francisco
    double lat = DEG_TO_RAD ( 37.0 );   // hard-coded approximate latitude of San Francisco
    double fov = DEG_TO_RAD ( MAX ( widthAngle, heightAngle ) );

    A3ImageInit ( &skyImage, lon, lat, jd, HUGE_VAL, HUGE_VAL, fov, width, height, A3IMAGE_RGBA8, pixels );

    // Compute image overall statistics.
    // If image is too bright or has insufficient contrast, refuse to find any objects.

    A3ImageStats	stats = { 0 };

    A3ImageStatsInRectangle ( &skyImage, &stats, 0, 0, width, height );
    if ( stats.fMedian > 192 || stats.fMax - stats.fMin < 2 )
    {
        A3ImageFree ( &skyImage );
        LOG_INFO ( "Image too bright (median=%.0f) or has insufficient contrast (%.0f); exiting.\n", stats.fMedian, stats.fMax - stats.fMin );
        return ( -1 );
    }

    int backRadius = MAX ( width, height ) * BACKRADIUS;
    int starRadius = MAX ( width, height ) * STARRADIUS;

    A3FindObjectParams	findParams = { 0, 0, width, height, MIN_OBJ_ALT, starRadius, backRadius, PEAKSIGMA, MEANSIGMA, MAX_MAJ_AXIS, MAX_AXIS_RATIO };

    int nFound = A3ImageFindObjects ( &skyImage, &findParams );

    LOG_INFO ( "Current JD=%.6f", jd );
    LOG_INFO ( "sumR = %.0f, sumG = %0.f, sumB = %.0f, sumA = %.0f", sumR, sumG, sumB, sumA );

    // Release memory for sky image data

    A3ImageFree ( &skyImage );
    return ( 0 );
}
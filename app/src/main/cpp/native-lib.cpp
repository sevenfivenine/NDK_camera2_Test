#include <jni.h>
#include <string>
#include <math.h>
#include <iostream>

#include <android/log.h>

#define  LOG_TAG    "NDKCameraTest"
#define  LOG_INFO(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)
#define  LOG_ERROR(...)  __android_log_print(ANDROID_LOG_ERROR,LOG_TAG,__VA_ARGS__)

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
    int *pixels = env->GetIntArrayElements ( pixelArray, 0 );

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
    }

    LOG_INFO ( "sumR = %.0f, sumG = %0.f, sumB = %.0f, sumA = %.0f", sumR, sumG, sumB, sumA );
    return ( 0 );
}
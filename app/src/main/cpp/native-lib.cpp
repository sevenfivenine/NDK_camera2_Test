#include <jni.h>
#include <string>
#include <math.h>
#include <iostream>

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

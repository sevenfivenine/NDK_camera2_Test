#include <jni.h>
#include <string>
#include <math.h>
#include <iostream>
#include <android/log.h>
#include "AstroLib.h"
#include "Astrometry.h"
#include "BrightStars.h"

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

int AddReferenceSolarSystemObjects ( A3Image *pImage );

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

    LOG_INFO ( "sumR = %.0f, sumG = %0.f, sumB = %.0f, sumA = %.0f", sumR, sumG, sumB, sumA );

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

    // Now find objects in the image, but use only the brightest MAX_OBJECTS.

    int backRadius = MAX ( width, height ) * BACKRADIUS;
    int starRadius = MAX ( width, height ) * STARRADIUS;

    A3FindObjectParams	findParams = { 0, 0, width, height, MIN_OBJ_ALT, starRadius, backRadius, PEAKSIGMA, MEANSIGMA, MAX_MAJ_AXIS, MAX_AXIS_RATIO };

    int nFound = A3ImageFindObjects ( &skyImage, &findParams );
    if ( nFound > MAX_OBJECTS )
        skyImage.nObjects = MAX_OBJECTS;

    LOG_INFO ( "Total %d objects found.\n", skyImage.nObjects );

    // If image center is unknown, use reference stars to mag 3.5; otherwise use stars to 4.5.

    float magLimit = REF_MAG_LIMIT;
    if ( isinf ( skyImage.ra0 ) || isinf ( skyImage.dec0 ) )
        magLimit = REF_MAG_LIMIT - 1.0;

    for ( int i = 0; i < gNumBrightStars; i++ )
    {
        if ( gBrightStars[i].mag > magLimit )
            continue;

        A3Reference reference = { 0 };

        reference.ra = HOUR_TO_RAD ( gBrightStars[i].ra );
        reference.dec = DEG_TO_RAD ( gBrightStars[i].dec );
        reference.mag = gBrightStars[i].mag;
        reference.data = NULL;

        strncpy ( reference.name, gBrightStars[i].bayer, sizeof ( reference.name ) );
        A3ImageAddReference ( &skyImage, &reference, MIN_REF_ALT );
    }

    int nRefs = AddReferenceSolarSystemObjects ( &skyImage );

    // Release memory for sky image data

    A3ImageFree ( &skyImage );
    return ( 0 );
}

// Add Sun, Moon and planets Mercury thru Saturn to the global array of refernce stars.

int AddReferenceSolarSystemObjects ( A3Image *pImage )
{
    int 	i = 0;
    double	jd = pImage->jd, lst = pImage->lst;
    double	eclmat[3][3] = { 0.0 }, xyz[3] = { 0 }, geoXYZ[3] = { 0 }, helXYZ[3] = { 0 }, relXYZ[7][3] = { 0 };
    double	lon[7] = { 0.0 }, lat[7] = { 0.0 }, rad[7] = { 0 }, helRad = 0.0;
    double	ra[7] = { 0.0 }, dec[7] = { 0.0 }, dst[7] = { 0.0 }, pha[7] = { 0.0 }, mag[7] = { 0.0 };
    static const char *names[7] = { "Sun", "Mercury", "Venus", "Moon", "Mars", "Jupiter", "Saturn" };

    // Compute matrix for precessing current-equinox ecliptic coords to equatorial.
    // Compute Earth's heliocentric XYZ position in current-epoch equatorial frame.
    // Compute viewer's geocentric XYZ in current-epoch equatorial frame.

    AASetEclipticRotationMatrix ( eclmat, AAObliquity ( jd ), -1 );
    VFPEarth ( jd, &lon[3], &lat[3], &helRad );
    AASphericalToXYZVector ( lon[3], lat[3], helRad, helXYZ );
    AATransformVector ( eclmat, helXYZ );
    AAGeodeticToGeocentricXYZ ( lst, pImage->lat, 0.0, 1.0, EARTH_FLATTENING, &geoXYZ[0], &geoXYZ[1], &geoXYZ[2] );

    // Compute current-equinox ecliptic longitude, latitude, and heliocentric distance in AU
    // of Earth and major planets Mercury thru Saturn.

    VFPMercury ( jd, &lon[1], &lat[1], &rad[1] );
    VFPVenus ( jd, &lon[2], &lat[2], &rad[2] );
    VFPMoon ( jd, &lon[3], &lat[3], &rad[3] );
    VFPMars ( jd, &lon[4], &lat[4], &rad[4] );
    VFPJupiter ( jd, &lon[5], &lat[5], &rad[5] );
    VFPSaturn ( jd, &lon[6], &lat[6], &rad[6] );

    // For each solar system object, compute J2000 equatorial RA, Dec, phase angle, and distance.

    for ( i = 0; i < 7; i++ )
    {
        AASphericalToXYZVector ( lon[i], lat[i], rad[i], xyz );
        AATransformVector ( eclmat, xyz );
        AAVectorDifference ( xyz, i == 3 ? geoXYZ : helXYZ, relXYZ[i] );
        pha[i] = AAPhaseAngle ( i == 3 ? helXYZ : xyz, relXYZ[i] );
        AAUnTransformVector ( pImage->premat, relXYZ[i] );
        AAXYZVectorToSpherical ( relXYZ[i], &ra[i], &dec[i], &dst[i] );
    }

    // Compute visual magnitudes from distances to sun, Earth, and phase angle

    mag[0] = AASunMagnitude ( helRad );
    mag[1] = AAMercuryMagnitude ( rad[1], dst[1], pha[1] );
    mag[2] = AAVenusMagnitude ( rad[2], dst[2], pha[2] );
    mag[3] = AAMoonMagnitude ( helRad, dst[3] / EARTH_RADII_PER_AU, pha[3] );
    mag[4] = AAMarsMagnitude ( rad[4], dst[4], pha[4] );
    mag[5] = AAJupiterMagnitude ( rad[5], dst[5], pha[5] );
    mag[6] = AASaturnMagnitude ( rad[6], dst[6], pha[6], AASaturnRingPlaneInclination ( relXYZ[6], dst[6] ) );

    A3Reference reference = { 0 };

    reference.ra = ra[i];
    reference.dec = dec[i];
    reference.mag = mag[i];
    reference.data = NULL;
    strncpy ( reference.name, names[i], sizeof ( reference.name ) );

    for ( i = 0; i < 7; i++ )
        A3ImageAddReference ( pImage, &reference, MIN_REF_ALT );

    return ( pImage->nReferences );
}

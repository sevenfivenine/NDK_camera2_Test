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
int CompareSkyReferenceMagnitudes ( const void *p1, const void *p2 );

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
Java_com_darkfuturestudios_ndktest_MainActivity_stackImageBuffers(
        JNIEnv *env,
        jobject,
        jintArray pixelArray,
        jint width,
        jint height,
        jintArray stackedPixelArray)
{
    jsize len = env->GetArrayLength(pixelArray);
    jsize stackedLen = env->GetArrayLength(stackedPixelArray);

    if (len != stackedLen || len != width * height) {
        LOG_ERROR ("stackImageBuffers: inconsistent array lengths (%d vs %d x %d vs %d",
                   len, width, height, stackedLen);
        return -1;
    }

    uint32_t *pixels = (uint32_t *) env->GetIntArrayElements(pixelArray, 0);
    uint32_t *stackedPixels = (uint32_t *) env->GetIntArrayElements(stackedPixelArray, 0);

    uint32_t r = 0, g = 0, b = 0;
    uint32_t stackedR = 0, stackedG = 0, stackedB = 0;

    for (int i = 0; i < len; i++)
    {
        r = pixels[i] & 0x000000FF;
        g = (pixels[i] >> 8) & 0x000000FF;
        b = (pixels[i] >> 16) & 0x000000FF;

        stackedR = stackedPixels[i] & 0x000000FF;
        stackedG = (stackedPixels[i] >> 8) & 0x000000FF;
        stackedB = (stackedPixels[i] >> 16) & 0x000000FF;

        stackedR = MIN ( stackedR + r, 255 );
        stackedG = MIN ( stackedG + g, 255 );
        stackedB = MIN ( stackedB + b, 255 );

        stackedPixels[i] = 0xFF000000 | (stackedB << 16) | (stackedG << 8) | stackedR;
    }

    return 0;
}

extern "C" JNIEXPORT jint
JNICALL
Java_com_darkfuturestudios_ndktest_MainActivity_processStackedImage(
        JNIEnv *env,
        jobject,
        jintArray pixelArray,
        jint width,
        jint height,
        jfloat widthAngle,
        jfloat heightAngle,
        jfloat exposureSecs,
        jint iso,
        jstring imagePath,
        jstring logDir)
{
    const char *logDirString = env->GetStringUTFChars ( logDir, 0 );
    char logPath[256] = { 0 };

    // Open log file

    strlcpy ( logPath, logDirString, sizeof ( logPath ) );
    strlcat ( logPath, "/AstrometryLog.txt", sizeof ( logPath ) );
    FILE *logFile = fopen ( logPath, "a" );
    if ( logFile == NULL )
        LOG_ERROR ( "Can't open log file %s!", logPath );

    env->ReleaseStringUTFChars ( logDir, logDirString );

    // Print image file path to log file, if we have one.

    if ( logFile )
    {
        const char *imagePathString = env->GetStringUTFChars ( imagePath, 0 );
        fprintf ( logFile, "Saved image %s\n", imagePathString );
        env->ReleaseStringUTFChars ( imagePath, imagePathString );
    }

    // Get length of, and pointer to, 32-bit ARGB image pixel data from the camera

    jsize len = env->GetArrayLength ( pixelArray );
    uint32_t *pixels = (uint32_t *) env->GetIntArrayElements ( pixelArray, 0 );

    // Initialize sky image struct with data from the camera

    A3Image skyImage = { 0 };
    double jd = AACurrentUTC();
    double lon = DEG_TO_RAD ( -122.0 ); // hard-coded approximate longitude of San Francisco
    double lat = DEG_TO_RAD ( 37.0 );   // hard-coded approximate latitude of San Francisco
    double fov = DEG_TO_RAD ( MAX ( widthAngle, heightAngle ) );

    A3ImageInit ( &skyImage, lon, lat, jd, HUGE_VAL, HUGE_VAL, fov, width, height, A3IMAGE_RGBA8, pixels );
    skyImage.exposure = exposureSecs;
    skyImage.iso = iso;
    skyImage.pLogFile = logFile;

    A3ImagePrint ( &skyImage, logFile );

    // Compute image overall statistics.
    // If image is too bright or has insufficient contrast, refuse to find any objects.

    A3ImageStats	stats = { 0 };

    A3ImageStatsInRectangle ( &skyImage, &stats, 0, 0, width, height );
    if ( stats.fMedian > 192 || stats.fMax - stats.fMin < 2 )
    {
        A3ImageFree ( &skyImage );
        if ( logFile )
        {
            fprintf(logFile,
                    "Image too bright (median=%.0f) or has insufficient contrast (%.0f); exiting.\n\n",
                    stats.fMedian, stats.fMax - stats.fMin);
            fclose(logFile);
            return(-1);
        }
    }

    // Now find objects in the image, but use only the brightest MAX_OBJECTS.

    int backRadius = MAX ( width, height ) * BACKRADIUS;
    int starRadius = MAX ( width, height ) * STARRADIUS;

    A3FindObjectParams	findParams = { 0, 0, width, height, MIN_OBJ_ALT, starRadius, backRadius, PEAKSIGMA, MEANSIGMA, MAX_MAJ_AXIS, MAX_AXIS_RATIO };

    int nFound = A3ImageFindObjects ( &skyImage, &findParams );
    if ( nFound > MAX_OBJECTS )
        skyImage.nObjects = MAX_OBJECTS;

    if ( logFile )
        fprintf ( logFile, "Total %d objects found.\n", skyImage.nObjects );

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
    qsort ( skyImage.pReferences, skyImage.nReferences, sizeof ( A3Reference ), CompareSkyReferenceMagnitudes );
    if ( logFile )
    {
        fprintf(logFile, "Total %d references.\n", skyImage.nReferences);
    }

    // Now attempt to identify objects in the image.  If we only attempt to identify the 6 brightest
    // objects, we'll get an answer much faster, and this usually works.  But we may fail to correctly
    // identify the objects if there are many bright "false stars" present.  So for now, we do an
    // exhaustive attempt to identify objects using all possible triplets.

    A3MatchParams	params = { 0 };

    params.width = width;
    params.height = height;
    params.tolerance = starRadius;
    params.timeout = TIMEOUT;

    int success = A3ImageIdentifyObjects ( &skyImage, &params );
    if ( success )
    {
        // after successful plate-solve, print out solved image center, scale, etc.

        if ( logFile )
        {
            fprintf(logFile, "%d objects matched, mean error is %.2f°.\n",
                    params.nObjsMatched,
                    RAD_TO_DEG (params.meanError));
            A3ImagePrint(&skyImage, logFile);
        }
    }
    else
    {
        if (logFile)
            fprintf(logFile, "Failed to identify enough objects to determine field of view.\n");
    }

    if ( logFile )
    {
        A3ImagePrintIdentifiedObjects(&skyImage, FALSE, logFile);
        fprintf(logFile, "\n" );
    }

    // Release memory for sky image data

    A3ImageFree ( &skyImage );
    fclose ( logFile );
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

    for ( i = 0; i < 7; i++ )
    {
        A3Reference reference = { 0 };

        reference.ra = ra[i];
        reference.dec = dec[i];
        reference.mag = mag[i];
        reference.data = NULL;
        strncpy ( reference.name, names[i], sizeof ( reference.name ) );

        A3ImageAddReference(pImage, &reference, MIN_REF_ALT);
    }

    return ( pImage->nReferences );
}

// Comparison function for sorting reference stars by magnitude

int CompareSkyReferenceMagnitudes ( const void *p1, const void *p2 )
{
    A3Reference *pRef1 = (A3Reference *) p1;
    A3Reference *pRef2 = (A3Reference *) p2;

    if ( pRef1->mag < pRef2->mag )
        return -1;
    else if ( pRef1->mag > pRef2->mag )
        return 1;
    else
        return 0;
}

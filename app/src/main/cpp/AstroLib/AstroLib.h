/*** COPYRIGHT NOTICE *********************************************************

Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.

******************************************************************************/

/*** Prevent multiple #includes! ***/

#ifndef ASTROLIB_H
#define ASTROLIB_H

#ifdef __cplusplus
extern "C" {
#endif

/*** ANSI headers ***/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>

/*** Place target platform-specific headers and macros into a separate
     "target.h" file - feel free to comment this line out if you don't
     have any, or if your build environment gives you a better way. ***/

#include "Target.h"

/*** Infinities and Undefineds ***/

#ifndef INFINITY
#define INFINITY	HUGE_VAL
#endif

/*** Numerical routines headers ***/

#include "NMatrix.h"

/*** Define all versions of true and false ***/

#ifndef false
#define false 0
#endif
#ifndef true
#define true (!false)
#endif

#ifndef FALSE
#define FALSE false
#endif
#ifndef TRUE
#define TRUE true
#endif

/************************ constants and macros  ******************************/

#define PI							3.141592653589						// 180 degrees = 3.141592653589 radians
#define HALF_PI						(0.5*PI)							//  90 degrees = 1.570796326795 radians
#define THREE_HALVES_PI				(1.5*PI)							// 270 degrees = 4.712388980385 radians
#define TWO_PI						(2.0*PI)							// 360 degrees = 6.28318530718 radians
#define THREE_PI					(3.0*PI)							// 540 degrees = 9.42477796077 radians
#define FOUR_PI						(4.0*PI)							// 720 degrees = 12.56637061436 radians
#define RAD_PER_DEG					(PI/180.0)
#define DEG_PER_RAD					(180.0/PI)							// Degrees per Radian = 57.295779513082
#define HOUR_PER_RAD				(12.0/PI)							// Hours per Radian = 3.819718634205
#define RAD_PER_HOUR				(PI/12.0)
#define ARCMIN_PER_RAD				(10800.0/PI)
#define RAD_PER_ARCMIN				(PI/10800.0)
#define ARCSEC_PER_RAD				(648000.0/PI)						// Arcseconds per Radian = 206264.8062471
#define RAD_PER_ARCSEC				(PI/648000.0)						
#define DEG_TO_RAD(x)				((x)*RAD_PER_DEG)					// converts degrees to radians
#define RAD_TO_DEG(x)				((x)*DEG_PER_RAD)					// converts radians to degrees
#define RAD_TO_HOUR(x)				((x)*HOUR_PER_RAD)					// converts radians to hours
#define HOUR_TO_RAD(x)				((x)*RAD_PER_HOUR)					// converts hours to radians
#define RAD_TO_ARCSEC(x)			((x)*ARCSEC_PER_RAD)				// converts radians to arcseconds
#define ARCSEC_TO_RAD(x)			((x)*RAD_PER_ARCSEC)				// converts arcseconds to radians
#define RAD_TO_ARCMIN(x)			((x)*ARCMIN_PER_RAD)				// converts radians to arcminutes
#define ARCMIN_TO_RAD(x)			((x)*RAD_PER_ARCMIN)				// converts arcminutes to radians

#ifndef MIN
#define MIN(x,y)					((x)<(y)?(x):(y))					// minimum of x or y
#endif

#ifndef MAX
#define MAX(x,y)					((x)>(y)?(x):(y))					// maximum of x or y
#endif

#ifndef CLIP
#define CLIP(x,min,max)				MIN(MAX(x,min),max)					// clips x to minimum and maximum
#endif

#ifndef SWAP
#define SWAP(type,x,y)	{ type t; t = x; x = y; y = t; }				// swaps values of x and y
#endif
	
#define B1900                		2415020.31352						// JD of standard Besseilian epoch B1900
#define B1950                		2433282.423							// JD of standard Besselian epoch B1950
#define J2000                		2451545.0							// JD of standard Julian epoch J2000
#define J1970						2440587.5							// JD of standard UNIX time base 1.0 January 1970 UTC
#define SEC_PER_DAY					86400								// Seconds per day
#define MIN_PER_DAY					1440								// Minutes per day
#define HRS_PER_DAY					24									// Hours per day
#define DAY_PER_JYEAR				365.25								// Days per Julian year
#define DAY_PER_BYEAR				365.242198781						// Days per Besselian year
#define SEC_PER_JYEAR				31557600							// Seconds per Julian year
#define SEC_PER_BYEAR				31556925.974678						// Seconds per Besselian year
#define DAY_PER_SIDEREAL_DAY		0.99726957							// Days per Sidereal Day
#define DAY_PER_SIDEREAL_MONTH		27.321661							// Days per Sidereal Month
#define DAY_PER_SYNODIC_MONTH		29.531								// Days per Synodic Month
#define DAY_PER_LUNAR_YEAR			(12 * DAY_PER_SYNODIC_MONTH)		// Days per Lunar Year
#define DAY_PER_SIDEREAL_YEAR		365.256363004						// Days per Sidereal Year
#define DAY_PER_SAROS				6585.322							// Days per Saros
#define SIDEREAL_PER_SOLAR			1.00273790934						// Sidereal per Solar days 

#define KM_PER_LIGHT_SEC			299792.458							// Speed of Light in Kilometers per Second
#define KM_PER_AU					149597870.0							// Kilometers per Astronomical Unit
#define AU_PER_KM					(1.0/KM_PER_AU)						// Astronomical Units per Kilometer = 0.00000000668459
#define AU_PER_PC					ARCSEC_PER_RAD						// Astronomical Units per Parsec
#define PC_PER_AU					RAD_PER_ARCSEC						// Parsecs per Astronomical Unit
#define KM_PER_LIGHT_DAY			(KM_PER_LIGHT_SEC*SEC_PER_DAY)		// Kilometers per light day = 25902068371.2
#define AU_PER_LIGHT_DAY			(KM_PER_LIGHT_DAY/KM_PER_AU)		// Astronomical Units per Light-Day = 173.14463348442
#define LIGHT_DAY_PER_AU			(KM_PER_AU/KM_PER_LIGHT_DAY)		// Light-Days per Astronomical Unit = 0.005775518
#define KM_PER_LY					(SEC_PER_JYEAR*KM_PER_LIGHT_SEC)	// Kilometers per Light-Year = 9460730472580.8
#define KM_PER_PC					(KM_PER_AU*AU_PER_PC)				// Kilometers per Parsec = 30856775670528.86
#define LY_PER_PC					(KM_PER_PC/KM_PER_LY)				// Light-Years per Parsec = 3.261563761906
#define PC_PER_LY					(KM_PER_LY/KM_PER_PC)				// Parsecs per Light-Year = 0.30660139522
#define AU_PER_LY					(KM_PER_LY/KM_PER_AU)				// Astronomical Units per Light-Year = 63241.0773801846
#define LY_PER_AU					(KM_PER_AU/KM_PER_LY)				// Light-Years per Astronomical Unit = 0.00001581250734

#define KM_PER_SOLAR_RADII			696342.0							// from SOHO spacecraft timings of Mercury transits.  Previously 695500.  See http://en.wikipedia.org/wiki/Solar_radius
#define AU_PER_SOLAR_RADII			0.004652							// Solar radius in AU

#define KM_PER_LUNAR_RADII			1738.1								// Moon's equatorial radius in km
#define AU_PER_LUNAR_RADII			0.000011618480932					// Moon's equatorial radius in AU

#define KM_PER_JUPITER_RADII		71492.0								// IAU nominal equatorial Jovian radius, 2015 value
#define AU_PER_JUPITER_RADII		0.000477894504781					// Jovian equatorial radius in AU

#define HELIO_GAUSS_CONST			0.01720209895						// days
#define GEO_GAUSS_CONST				0.0743669161						// minutes
#define EARTHS_PER_SOLAR_MASS		332946.0
#define KM_PER_EARTH_RADII			6378.140							// IAU 1976 reference ellipsoid
#define EARTH_RADII_PER_AU			(KM_PER_AU/KM_PER_EARTH_RADII)
#define AU_PER_EARTH_RADII			(KM_PER_EARTH_RADII/KM_PER_AU)
#define EARTH_FLATTENING			(1.0/298.257)						// IAU 1976 reference ellipsoid
#define EARTH_ROTATION_RATE			6.300387487008						// sidereal diurnal angular rotation rate, radians per day
#define EARTH_J2					0.00108263
#define EARTH_J3					-0.254E-5
#define EARTH_J4					-0.161E-5
#define MARS_FLATTENING				0.00647630
#define JUPITER_FLATTENING			0.0648744
#define SATURN_FLATTENING			0.0979624
#define URANUS_FLATTENING			0.0229273
#define NEPTUNE_FLATTENING			0.0171

#define PHOBOS_FLATTENING_X			0.0									// 26.8 km long (Mars-pointing) axis
#define PHOBOS_FLATTENING_Y			(1.0-22.4/26.8)						// 22.4 km axis pointing in direction of orbital motion
#define PHOBOS_FLATTENING_Z			(1.0-18.4/26.8)						// 18.4 km polar axis

#define DEIMOS_FLATTENING_X			0.0									// 15 km long (Mars-pointing) axis
#define DEIMOS_FLATTENING_Y			(1.0-12.2/15.0)						// 12.2 km axis pointing im direction of orbital motion
#define DEIMOS_FLATTENING_Z			(1.0-10.4/15.0)						// 10.4 km polar axis

#define AMALTHEA_FLATTENING_X		0.0									// 250 km long (Jupiter-pointing) axis
#define AMALTHEA_FLATTENING_Y		(1.0-146.0/250.0)					// 146 km axis pointing im direction of orbital motion
#define AMALTHEA_FLATTENING_Z		(1.0-128.0/250.0)					// 128 km polar axis

#define HYPERION_FLATTENING_X		0.0									// 360.2 km long axis
#define HYPERION_FLATTENING_Y		(1.0-266/360.2)						// 266 km
#define HYPERION_FLATTENING_Z		(1.0-205.4/360.2)					// 205.4 km polar axis
	
#define PHOEBE_FLATTENING_X			0.0									// 218.8 km long axis
#define PHOEBE_FLATTENING_Y			(1.0-217/218.8)						// 217 km
#define PHOEBE_FLATTENING_Z			(1.0-203.6/218.8)					// 203.6 km polar axis

#define CERES_FLATTENING_X			0.0									// 965.2 km long axis
#define CERES_FLATTENING_Y			(1.0-961.2/965.2)					// 961.2 km
#define CERES_FLATTENING_Z			(1.0-891.2/965.2)					// 891.2 km polar axis

#define VESTA_FLATTENING_X			0.0									// 572.6 km long axis
#define VESTA_FLATTENING_Y			(1.0-557.2/572.6)					// 557.2 km
#define VESTA_FLATTENING_Z			(1.0-446.4/572.6)					// 446.4 km polar axis

#define KG_PER_EARTH_MASS			5.9742e+24							// Earth's mass in kilograms
#define KG_PER_SOLAR_MASS			1.9891e+30							// Sun's mass in kilograms
#define KG_PER_LUNAR_MASS			7.3477e+22							// Moon's mass in kilograms
	
#define MASS_SUN					332946.050894783					// Sun's mass relative to Earth's
#define MASS_MERCURY				0.05527359899309					// Mercury's mass relative to Earth's
#define MASS_VENUS					0.81499810842015					// Venus's mass relative to Earth's
#define MASS_EARTH					1.00000000000000					// Earth's mass relative to Earth's
#define MASS_MARS					0.10744688443521					// Mars's mass relative to Earth's
#define MASS_JUPITER				317.828391263478					// Jupiter's mass relative to Earth's
#define MASS_SATURN					95.1609799938971					// Saturn's mass relative to Earth's
#define MASS_URANUS					14.5357348649260					// Uranus's mass relative to Earth's
#define MASS_NEPTUNE				17.1477661508462					// Neptune's mass relative to Earth's
#define MASS_PLUTO					0.0025611							// Pluto's mass relative to Earth's
#define MASS_MOON					0.01230003827772					// Moon's mass relative to Earth's
#define MASS_CHARON					0.000254							// Charon's mass relative to Earth's
#define MASS_EARTH_SYSTEM			1.01230003827772					// Earth + Moon system mass relative to Earth (alone)
#define MASS_JUPITER_SYSTEM			317.894205324555					// Jupiter + moon system mass relative to Earth (alone)
#define MASS_SATURN_SYSTEM			95.1846082689706					// Saturn + moon system mass relative to Earth (alone)
#define MASS_URANUS_SYSTEM			14.5372371147677					// Uranus + moon system mass relative to Earth (alone)
#define MASS_NEPTUNE_SYSTEM			17.1513463098694					// Neptune + moon system mass relative to Earth (alone)

#define GALACTIC_N_POLE_RA			192.25								// in degrees, B1950
#define GALACTIC_N_POLE_DEC			27.4								// in degrees, B1950
#define GALACTIC_LON_ASC_NODE		33.0								// in degrees

/**********************  structs and typedefs  ***************************/

typedef struct VFPTerm				VFPTerm;
typedef struct GuideStar			GuideStar, *GuideStarPtr;
typedef struct GuideStarList		GuideStarList, *GuideStarListPtr;

typedef double	XYZVector[3];
typedef double	XYZMatrix[3][3];
	
typedef double (*ImageSurfaceFuncPtr)	( long, long, double[] );
typedef double (*ImageModelFuncPtr)		( long, long, double[], double[] );

/*** Event codes for Jupiter and moon events ***/
	
#define JUPITER_SPOT_TRANSIT	0x00	// Jupiter Great Red Spot transits
#define JUPITER_MOON_ECLIPSE	0x01	// moon eclipse by Jupiter begins/ends
#define JUPITER_MOON_OCCULT		0x02	// moon occultation behind Jupiter begins/ends
#define JUPITER_MOON_TRANSIT	0x03	// moon transits in front of Jupiter begins/ends
#define JUPITER_MOON_SHADOW		0x04	// moon shadow transits across Jupiter begins/ends
#define JUPITER_EVENT_MASK		0x0F	// bit mask for extracting event type from event code
#define JUPITER_EVENT_END		0x80	// if set, this is the end of an event

typedef struct JupiterMoonEvent
{
	double	jd;		// Julian date on which event occurs
	int		moon;	// 0 = Jupiter (GRS), 1 = Io, 2 = Europa, 3 = Ganymede, 4 = Callisto, -1 = multiple moon shadow event
	int		type;	// JUPITER_MOON_ECLIPSE ... JUPITER_MOON_SHADOW, possibly bitwise ORed with JUPITER_MOON_EVENT_END
}
JupiterMoonEvent;
	
/**********************  Functions in Angle.c *****************************/

/******************************  atan2pi  *********************************

  DESCRIPTION: Returns the arctangent of (y/x) in the range zero to TWO_PI.

     SYNOPSIS: double atan2pi ( double x, double y )

          (y): Y coordinate
		  (x): X coordinate

      RETURNS: arctangent of (y/x) in the range 0 to TWO_PI.  If the angle
	           is undefined, returns zero.

****************************************************************************/

double ASTROLIBDLL_API atan2pi ( double, double );

/*********************************  Mod2Pi  ********************************

	Reduces an angle to a specific range.

	double ModPi ( double x )
	double Mod2Pi ( double x )
	double Mod12h ( double x )
	double Mod24h ( double x )
	double Mod180 ( double x )
	double Mod360 ( double x )

	(x): Value of the angle, in radians.

	Returns value of the angle, reduced to the following ranges:
	
	-PI to PI for ModPi();
	0 to TWO_PI for Mod2Pi();
	-12 to 12 for Mod12h();
	0 to 24 for Mod24h();
	-180 to 180 for Mod180();
	0 to 360 for Mod360();
	
*****************************************************************************/

double ASTROLIBDLL_API ModPi ( double );
double ASTROLIBDLL_API Mod2Pi ( double );
double ASTROLIBDLL_API Mod24h ( double );
double ASTROLIBDLL_API Mod360 ( double );

/*** AADegMinSecToDecimal **************************************************

	Converts an angle in degrees, minutes, and seconds to its equivalent
	value in decimal degrees.

	double AADegMinSecToDecimal ( short deg, short min, double sec, char sign )

	(deg): Degree part of the angle.
	(min): Minutes part of the angle.
	(sec): Seconds part of the angle.
	(sign): Sign of the angle.

	Returns the decimal-degree value of the angle.  The angle is assumed
	to be positive regardless of the sign on (deg), (min), or (sec), unless
	the value passed in (sign) is '-'.

****************************************************************************/

double ASTROLIBDLL_API AADegMinSecToDecimal ( short, short, double, char );

/*** AADecimalToDegMinSec **************************************************

	Converts an angle in decimal degrees to its equivalent value
	in degrees, minutes, and seconds.

	void AADecimalToDegMinSec ( double dec,
	     short *deg, short *min, double *sec, char *sign )

	(dec): Decimal-degree value of the angle.
	(deg): Receives degree part of the angle.
	(min): Receives minutes part of the angle.
	(sec): Receives seconds part of the angle.
	(sgn): Receives sign of the angle, as '+' is the angle is zero
	       or positive, or as '-' if the angle is negative.

	This function returns nothing.
	
****************************************************************************/

void ASTROLIBDLL_API AADecimalToDegMinSec ( double, short *, short *, double *, char * );
void ASTROLIBDLL_API AADecimalToDegMin ( double decimal, short *deg, double *min, char *sign );
	
/*** AASeparation *******************************************************

	Computes the angular separation between two points on a sphere.

	double AASeparation ( double a1, double d1, double a2, double d2 )

	(a1): Longitude coordinate of the first point, in radians.
	(d1): Latitude coordinate of the first point, in radians.
	(a2): Longitude coordinate of the second point, in radians.
	(d2): Latitude coordinate of the second point, in radians.

	Returns angular separation between points (a1,d1) and (a2,d2) in radians.

***************************************************************************/

double ASTROLIBDLL_API AASeparation ( double, double, double, double );

/*** AAPositionAngle *****************************************************

	Computes the direction from one point on a sphere to another.

	double AAPositionAngle ( double a1, double d1, double a2, double d2 )

	(a1): Longitude coordinate of the first point, in radians.
	(d1): Latitude coordinate of the first point, in radians.
	(a2): Longitude coordinate of the second point, in radians.
	(d2): Latitude coordinate of the second point, in radians.

	Returns the direction from point (a1,d1) to (a2,d2), in radians.

	Position angle is measured eastward (counterclockwise) from north:
	north = 0, east = 90, south = 180, west = 270. (Note that these values
	would be returned in radians!)

**************************************************************************/

double ASTROLIBDLL_API AAPositionAngle ( double, double, double, double );

/*** AAVectorSeparation ************************************************

	Computes the angular separation between two points on a sphere,
	where both are expressed as rectangular unit vectors.
	
	double AAVectorPositionAngle ( double u[3], double v[3] )

	(u): Elements of the first vector.
	(v): Elements of the second vector.
	     
	Returns angular separation between vectors (u) and (v), in radians.
	
***********************************************************************/

double ASTROLIBDLL_API AAVectorSeparation ( double [3], double [3] );

/*** AAVectorPositionAngle  ********************************************

	Computes the position angle (direction) from one point to another
	on a sphere, where both are expressed as rectangular unit vectors.
	
	double AAVectorPositionAngle ( double u[3], double v[3] )

	(u): Elements of the first vector.
	(v): Elements of the second vector.
	     
	Returns position angle from vector (u) to vector (v), in radians.

	Position angle is measured eastward (counterclockwise)
	from north: north = 0, east = 90, south = 180, west = 270.
	
	(Note that these values would be returned in radians!)
	
***********************************************************************/

double ASTROLIBDLL_API AAVectorPositionAngle ( double [3], double [3] );

/*********************  Functions in CoordSys.c ***********************/

/***  AASphericalToXYZ  **************************************************

	Converts spherical coordinates to rectangular.
  
	void AASphericalToXYZ ( double l, double b, double r,
	double *x, double *y, double *z )

	void AASphericalToXYZVector ( double l, double b, double r, double v[3] )

	(l): Longitude coordinate, in radians.
	(B): Latitude coordinate, in radians.
	(r): Radial coordinate.
	(x): Receives x coordinate.
	(y): Receives y coordinate.
	(z): Receives z coordinate.
	(v): Receives array of (x,y,z) coordinates
	          
	Return nothing.

***************************************************************************/

void ASTROLIBDLL_API AASphericalToXYZ ( double, double, double, double *, double *, double * );
void ASTROLIBDLL_API AASphericalToXYZVector ( double, double, double, double[3] );

/***  AAXYZToSpherical  ****************************************************

	Converts rectangular coordinates to spherical.
  
	void AAXYZToSpherical ( double x, double y, double z,
	double *l, double *b, double *r )

	void AAXYZVectorToSpherical ( double x, double y, double z,
	double *l, double *b, double *r )

	(x): X coordinate.
	(y): Y coordinate.
	(z): Z coordinate.
	(v): Vector containing X-Y-Z coordinates.
	(l): Receives longitude coordinate, in radians from 0 to TWO_PI.
	(b): Receives latitude coordinate, in radians from -HALF_PI to HALF_PI.
	(r): Receives radial coordinate.  If you do not wish to compute
	     the radial coordinate, pass NULL instead of a pointer.
          
	Returns nothing.  If the (r) contains zero, (a) and (d) will contain
	zero.

**************************************************************************/

void ASTROLIBDLL_API AAXYZToSpherical ( double, double, double, double *, double *, double * );
void ASTROLIBDLL_API AAXYZVectorToSpherical ( double[3], double *, double *, double * );

/***  AASphericalToXYZMotion  ***********************************************

	Converts position and motion in spherical coordinates to position
	and motion in rectangular coordinates.

	void AASphericalToXYZMotion ( double l, double b, double r,
	double dl, double db, double dr, double *x, double *y, double *z,
	double *dx, double *dy, double *dz )

	 (l): Longitude coordinate, in radians.
	 (b): Latitude coordinate, in radians.
	 (r): Radial coordinate.
	(dl): Change in longitude coordinate per unit time.
	(db): Change in longitude coordinate per unit time.
	(dr): Change in radial coordinate per unit time.
	 (x): Receives x coordinate.
	 (y): Receives y coordinate.
	 (z): Receives z coordinate.
	(dx): Receives x velocity.
	(dy): Receives y velocity.
	(dz): Receives z velocity.

	Returns nothing.

	The longitude velocity (dl) should be expressed as the change
	in the longitude coordinate (l) per unit time.  To convert this
	value to radians per unit time, you should multiply it by the
	cosine of the latitude coordinate (b).  Note that latitude velocity
	is already expressed in radians per unit time.
	
*************************************************************************/
	
void ASTROLIBDLL_API AASphericalToXYZMotion ( double, double, double, double, double, double,
     double *, double *, double *, double *, double *, double *);

/***  AAXYZToSphericalMotion  ***********************************************

	Converts position and motion in rectangular coordinates to position
	and motion in spherical coordinates.
               
	void AAXYZToSphericalMotion ( double x, double y, double z,
	double dx, double dy, double dz, double *l, double *b, double *r,
	double *dl, double *db, double *dr )

	 (x): X coordinate.
	 (y): Y coordinate.
	 (z): Z coordinate.
	(dx): X velocity.
	(dy): Y velocity.
	(dz): Z velocity.
	 (l): Receives longitude coordinate, in radians.
	 (b): Receives latitude coordinate, in radians.
	 (r): Receives radial coordinate.
	(dl): Receives change in longitude coordinate per unit time.
	(db): Receives change in latitude coordinate per unit time.
	(dr): Receives change in radial coordinate per unit time.

	Returns nothing.  For singular cases, which occur if (r) or both
	(x) and (y) are zero, the function will return zero for all values.
	
	The longitude velocity (dl) is expressed as the change in the
	longitude coordinate (l) per unit time.  To convert this value
	to radians per unit time, you should multiply it by the cosine
	of the latitude coordinate (b).  Note that latitude velocity is 
	already expressed in radians per unit time.
	
***************************************************************************/

void ASTROLIBDLL_API AAXYZToSphericalMotion ( double, double, double, double, double, double,
double *, double *, double *, double *, double *, double * );

void ASTROLIBDLL_API AAXYZVectorToSphericalMotion ( double v[3], double dv[3],
double *a, double *d, double *r, double *da, double *dd, double *dr );

/*** AASetVector **********************************************************

	Sets elements of a vector.

	void AASetVector ( double v[3], double x, double y, double z )
               
	(v): Contains the 3 elements of the vector to set.
    (x): contains X coordinate to copy into vector.
    (y): contains Y coordinate to copy into vector.
    (z): contains Z coordinate to copy into vector.
    
	Returns nothing.                   
  
***************************************************************************/

#define AASetVector(v,x,y,z)	{ (v)[0] = (x); (v)[1] = (y); (v)[2] = (z); }

/*** AACopyVector **********************************************************

	Copies elements of one vector into another.

	void AACopyVector ( double v2[3], double v1[3] )
               
	(v1): Source vector to copy elements from.
    (v2): Destination vector to copy elements into.
    
	Returns pointer to destination vector (v2).                   
  
***************************************************************************/

#define AACopyVector3(v2,v1)		{ (v2)[0] = (v1)[0]; (v2)[1] = (v1)[1]; (v2)[2] = (v1)[2]; }

double ASTROLIBDLL_API *AACopyVector ( double v2[3], double v1[3] );
	
/*** AAEqualVector **********************************************************

	Tests whether or not two vectors are equal.

	void AAEqualVector ( double v2[3], double v1[3] )
               
	(v1): First vector.
    (v2): Second vector.
    
	Returns nothing.                   
  
*****************************************************************************/

#define AAEqualVector(v2,v1)	( (v1)[0] == (v2)[0] && (v1)[1] == (v2)[1] && (v1)[2] == (v2)[2] )

/*** AAVectorMagnitude ****************************************************

	Computes the magnitude of a vector.

	double AAVectorMagnitude ( double v[3] );
               
	(v): Contains the 3 elements of the vector.
          
	Returns magnitude of the vector.                   
  
***************************************************************************/

#define AAVectorMagnitude2(v)	(((v)[0]*(v)[0])+((v)[1]*(v)[1])+((v)[2]*(v)[2]))
#define AAVectorMagnitude3(v)	(sqrt(AAVectorMagnitude2(v)))
	
double ASTROLIBDLL_API AAVectorMagnitude ( double [3] );

/*** AAScaleVector ********************************************************

	Multiplies a vector by a constant.

	double *AAScaleVector ( double v[3], double k )
               
	(v): Contains the 3 elements of the vector.
    (k): Constant by which to multiply vector.
          
	The function returns a pointer to the given vector (v).                   
  
***************************************************************************/

#define AAScaleVector3(v,k)	{ (v)[0] *= (k); (v)[1] *= (k); (v)[2] *= (k); }

double ASTROLIBDLL_API *AAScaleVector ( double [3], double );

/*** AANormalizeVector ****************************************************

	Scales a vector to unit length and returns its original magnitude.
	
	double AANormalizeVector ( double v[3] );
	
	(v): contains the 3 elements of the vector.
	
	The function returns the vector's original magnitude.
	
***************************************************************************/

#define AANormalizeVector3(v)	{ double m = AAVectorMagnitude3(v); if ( m ) { m = 1.0 / m; AAScaleVector3 ( v, 1.0 / m ); } }

double ASTROLIBDLL_API AANormalizeVector ( double [3] );

/*** AANegateVector ********************************************************

	Multiplies the elements of a vector by -1.

	double AANegateVector ( double v[3] )
	
	Returns nothing.
	
****************************************************************************/

#define AANegateVector(v)	{ (v)[0] = -(v)[0]; (v)[1] = -(v)[1]; (v)[2] = -(v)[2]; }

/*** AAVectorSum **********************************************************

	Computes the sum of two vectors.
	
	double *AAVectorSum ( double u[3], double v[3], double w[3] )

	(u): Contains the 3 elements of the first vector.
	(v): Contains the 3 elements of the second vector.
	(w): Recieves the 3 elements of the vector sum u + v.
	
	The function returns a pointer to the vector that receives the
	sum (w).  Note that (w) can be the same as either (u) or (v).
          
***************************************************************************/

double ASTROLIBDLL_API *AAVectorSum ( double [3], double [3], double [3] );
double ASTROLIBDLL_API *AAScale1VectorSum ( double ka, double a[3], double b[3], double c[3] );
double ASTROLIBDLL_API *AAScale2VectorSum ( double ka, double a[3], double kb, double b[3], double c[3] );

#define AAVectorSum3(u,v,w)		{ (w)[0] = (u)[0] + (v)[0]; w[1] = (u)[1] + (v)[1]; (w)[2] = (u)[2] + (v)[2]; }

/*** AAVectorDifference *****************************************************

	Computes the difference between two vectors.

	double *AAVectorDifference ( double u[3], double v[3], double w[3] );
               
	(u): Contains the 3 elements of the first vector.
	(v): Contains the 3 elements of the second vector.
	(w): Recieves the 3 elements of the vector difference u - v.
          
	The function returns a pointer to the vector that receives the
	difference (w).  Note that (w) can be the same as either (u) or (v).
  
***************************************************************************/

double ASTROLIBDLL_API *AAVectorDifference ( double [3], double [3], double [3] );

#define AAVectorDifference3(u,v,w)	{ (w)[0] = (u)[0] - (v)[0]; w[1] = (u)[1] - (v)[1]; (w)[2] = (u)[2] - (v)[2]; }

/*** AAVectorDistance ****************************************************
	 
	Computes the distance from one vector to another, i.e. the magnitude
	of the difference between the vectors.
	 
	double AAVectorDistance ( double u[3], double v[3] );
	 
	(u): Contains the 3 elements of the first vector.
	(v): Contains the 3 elements of the second vector.
	 
	Returns magnitude of the vector.                   
	 
***************************************************************************/
	
double ASTROLIBDLL_API AAVectorDistance ( double [3], double [3] );
	
/*** AADotProduct **********************************************************

	Computes the dot product of two vectors.
               
	double AADotProduct ( double u[3], double v[3]);

	(u): Contains the 3 elements of the first vector.
	(v): Contains the 3 elements of the second vector.
          
	Returns dot product of the two vectors.                   
  
***************************************************************************/

double ASTROLIBDLL_API AADotProduct ( double [3], double [3] );
	
#define AADotProduct3(u,v)	( (u)[0] * (v)[0] + (u)[1] * (v)[1] + (u)[2] * (v)[2] )

/***  AACrossProduct ********************************************************

	Computes the cross product of two vectors.
               
	double *AACrossProduct ( double u[3], double v[3], double w[3] );

	(u): Contains the 3 elements of the first vector.
	(v): Contains the 3 elements of the second vector.
	(w): Recieves the 3 elements of the cross product u x v.
          
	This function returns a pointer to the vector which receives the
	cross product (w).  Note that the vector cross product is not
	commutative, i.e. u x v does not equal v x u; rather, u x v = - v x u.                
  
***************************************************************************/

double ASTROLIBDLL_API *AACrossProduct ( double[3], double[3], double[3] );

/*** AANormalVector ********************************************************

	Given an input vector, generates another vector which is normal
	(perpendicular) to that vector in an arbitrary direction.
               
	double *AANormalVector ( double v[3], double p[3] );

	(v): Contains the 3 elements of the input vector.
	(p): Recieves the 3 elements of the perpendicular output vector.
         
	This function returns a pointer to the vector which receives the normal/
	perpendicular vector (p).  The magnitude of the output normal vector will
	vary depending on the magnitude of the input vector, but will be > 0
	as long as the magnitude of the input vector is > 0.
	
***************************************************************************/

double ASTROLIBDLL_API *AANormalVector ( double v[3], double p[3] );

/***  AATransformVector **********************************************************

	Multiplies an XYZ coordinate vector by a 3x3 rotation matrix
	to transform it into a new coordinate frame.

	double *AATransformVector ( double m[3][3], double v[3] )
	
	(m): Contains the 3 rows of the rotation matrix sequentially;
	     must consist of 9 consecutive doubles.
	(v): Contains the 3 elements of the vector.
	
	The function returns a pointer to the vector (v).  Upon return,
	the elements of (v) will be multiplied by (m).                   
  
********************************************************************************/

#define AATransformVector3(w,m,v)	{ \
	w[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2]; \
	w[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2]; \
	w[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]; \
};

double ASTROLIBDLL_API *AATransformVector ( double [3][3], double[3] );

/***  AAUnTransformVector  *******************************************************

	Multiplies an XYZ coordinate vector by the transpose of a 3x3 rotation
	matrix.  Performs the opposite transformation as AATransformVector().

	void AAUnTransformVector ( double m[3][3], double v[3] )
               
	(m): Contains the 3 rows of the rotation matrix sequentially;
	     must consist of 9 consecutive doubles.
	(v): Contains the 3 elements of the vector.
          
	The function returns a pointer to the vector (v).  Upon return, the elements
	of (v) will be multiplied by the transpose of (m).                   
  
*********************************************************************************/

#define AAUnTransformVector3(v,m,w)	{ \
	v[0] = m[0][0] * w[0] + m[1][0] * w[1] + m[2][0] * w[2]; \
	v[1] = m[0][1] * w[0] + m[1][1] * w[1] + m[2][1] * w[2]; \
	v[2] = m[0][2] * w[0] + m[1][2] * w[1] + m[2][2] * w[2]; \
};
	
double ASTROLIBDLL_API *AAUnTransformVector ( double [3][3], double [3] );

/***  AATransformRotationMatrix  ***************************************************

	Multiplies one 3x3 rotation matrix by another to make a single matrix
	which is the equivalent of the first rotation followed by the second.

	void AATransformRotationMatrix ( double b[3][3], double a[3][3] )
               
	(a): Contains the first rotation matrix.
	(b): Contains the second rotation matrix.
          
	Returns nothing.  Upon return, the columns of (a) will be rotated by (b).
	Note that the product b * a is not the same as a * b!!!
	
**********************************************************************************/
	
void ASTROLIBDLL_API AATransformRotationMatrix ( double [3][3], double[3][3] );

/*** AAUnTransformRotationMatrix *************************************************

	Multiplies one 3x3 rotation matrix by the transpose of another.
               
	void AAUnTransformRotationMatrix ( double a[3][3], double b[3][3] );

	(a): Contains the first rotation matrix.
	(b): Contains the second rotation matrix.
	
	Returns nothing.  Upon return, the columns of (a) will be multipled
	by the columns of (b), i.e. rotated by the matrix which is the
	inverse of b.
               
**********************************************************************************/

void ASTROLIBDLL_API AAUnTransformRotationMatrix ( double [3][3], double [3][3] );

/*** AACopyRotationMatrix ************************************************

	Copies elements from one rotation matrix into another.
  
	void CopyRotationMatrix ( double b[3][3], double a[3] )

	(b): Pointer to destination matrix.
	(a): Pointer to source matrix.
          
	Both matrices must contain 9 doubles (the matrix rows).
               
	Returns nothing.

**************************************************************************/

#define AACopyRotationMatrix3(m2,m1)	{ \
	AACopyVector(m2[0],m1[0]); \
	AACopyVector(m2[1],m1[1]); \
	AACopyVector(m2[2],m1[2]); \
}
	
void ASTROLIBDLL_API AACopyRotationMatrix ( double [3][3], double [3][3] );

/*** AATransposeRotationMatrix *******************************************

	Transposes a 3x3 rotation matrix in place.
  
	void AATransposeRotationMatrix ( double m[3][3], double t[3][3] )
	
	(m): Address of the matrix, which must contain 9 doubles
	     (the 3 rows of the matrix sequentially).
    (t): Address of the matrix to recieve transpose.
               
	Returns nothing.  Upon return, the matrix in (t) is replaced
	with the transpose of (m).  If both (m) and (t) point to the same
 	matrix, then the matrix is transposed in place.

**************************************************************************/

void ASTROLIBDLL_API AATransposeRotationMatrix ( double m[3][3], double t[3][3] );

/*** AASetRotationMatrix **************************************************

	Calculates the elements for a generalized rotation matrix representing
	an arbitrary number of rotations about the principal XYZ axes.
               
	void AASetRotationMatrix ( double m[3][3], int n, ... )

	  (m): An array of 9 doubles; will receive the elements of the matrix.
	  (n): The total number of rotations the rotation matrix represents.
	(...): Describes the actual rotations.  There must be (n) pairs of
	       arguments following (n).  The first argument in each pair
	       is an integer, representing the axis about which that rotation
	       is supposed to take place: 0 = X axis, 1 = Y axis, and 2 = Z axis.
	       The second argument in each pair is a double, representing the
	       angular size of the rotation in radians.
      
	Returns nothing.
               
	For example,
	
	AASetRotationMatrix ( m, 3, 0, 0.5, 1, -0.6, 2, 0.7 );

	would return in (m) the elements of the rotation matrix that
	performs the following 3 transformations:
  
	1) rotation around the X axis by 0.5 rad;
	2) rotation around the Y axis by -0.6 rad;
	3) rotation around the Z axis by 0.7 rad.

****************************************************************************/

void ASTROLIBDLL_API AASetRotationMatrix ( double [3][3], int, ... );

/***  AASetIdentityRotationMatrix  *****************************************

	Creates an identity rotation matrix.

	void AASetIdentityRotationMatrix ( double m[3][3] )
  
	(m): Address of the matrix, which must contain 9 doubles (the 3 rows
	     of the matrix sequentially).
               
	Returns nothing.  On return, (m) contains a 3x3 identity matrix.

**************************************************************************/

#define AASetIdentityRotationMatrix3(m)	{ \
	AASetVector(m[0],1,0,0); \
	AASetVector(m[1],0,1,0); \
	AASetVector(m[2],0,0,1); \
}

void ASTROLIBDLL_API AASetIdentityRotationMatrix ( double [3][3] );

/*** AARotationMatrixToOpenGLMatrix ******************************************
 
	Converts a 3x3 rotation matrix to an OpenGL 4x4 transformation matrix.
 
	double *AARotationMatrixToOpenGLMatrix ( double m[3][3], double p[16] )

	m: pointer to 3x3 rotation matrix
	p: pointer to buffer to store OpenGL 4x4 transformation matrix.
 
	The function returns a pointer to the output matrix (p).
	If p is NULL, return a pointer to an internal static buffer which
	holds the converted 3x3 matrix.
 
	The converted 4x4 matrix is stored in column-major order,
	conventional for OpenGL.
 
**************************************************************************/

double *AARotationMatrixToOpenGLMatrix ( double m[3][3], double p[16] );

/*** AARotationMatrixFromGLMatrix ******************************************
 
	Converts an OpenGL 4x4 transformation matrix to a 3x3 rotation matrix.
 
	double *AARotationMatrixFromOpenGLMatrix ( double m[3][3], double p[16] )
 
	m: recieves output 3x3 rotation matrix
	p: input OpenGL 4x4 transformation matrix.
 
	The function returns nothing.  It assumes the input 4x4 matrix is
	stored in column-major order, as conventional for OpenGL.
 
**************************************************************************/

void AARotationMatrixFromOpenGLMatrix ( double m[3][3], double p[16] );

/**********************  Functions in Time.c *****************************/

/*** AALocalJD ***********************************************************

	Returns the current value of the computer's system clock, expressed
	as a Julian Date.
	
	double AALocalJD ( void )
	
	This function takes no arguments, and returns the current value of
	the computer's system clock expressed as a Julian Day number.  No
	attempt is made to convert the time from the computer's local time
	zone to UTC, so the value returned is in local time.  To convert
	this value to UTC, subtract the local time zone, expressed as a
	fraction of a day (i.e. the number of hours east of Greenwich
	divided by 24). 
	
	This function makes use of the standard C library timing functions
	provided in <time.h>; this allows the computer's system clock to be
	read to a precision of one second (although the clock may measure
	time more to a higher precision internally.)
	
****************************************************************************/

double ASTROLIBDLL_API AALocalJD ( void );

/*** AACurrentJD ***********************************************************
	 
	Returns the current value of the computer's system clock, expressed
	as a Julian Date.
	 
	double AACurrentJD ( void )
	 
	This function takes no arguments, and returns the current value of
	the computer's system clock expressed as a Julian Day number.  The
	computer's internal time zone setting is used to convert local time
	to UTC; ensure the computer's time zone is set correctly.
	 
	This function makes use of the standard C library timing functions
	provided in <time.h>; this allows the computer's system clock to be
	read to a precision of one second (although the clock may measure
	time more to a higher precision internally.)
	 
****************************************************************************/
	
double ASTROLIBDLL_API AACurrentJD ( void );

/*** AACurrentUTC *************************************************************

	Returns the current Universal Coordinated Time, expressed
	as a Julian Date.
	
	double AACurrentUTC ( void )
	
	This function takes no arguments, and returns the current Universal
	Coordinated Time, expressed as a Julian Day number.  The UTC is taken
	from the computer's system clock, and converted to UTC using the time
	zone value stored in the operating system.
	
	This function makes use of OS-specific timing functions which allows
 	the computer's system clock to be read to a precision of one millisecond.

*******************************************************************************/

double ASTROLIBDLL_API AACurrentUTC ( void );

/*** AACurrentTimeZone ********************************************************

	Returns current difference between local time and UTC, and optionally
	whether or nut daylight saving time is in effect.
 
	double AACurrentTimeZone ( int *isDST )
 
	(isDST): if not NULL, recieves TRUE or FALSE depending on whether daylight
	         saving time is in effect.
 
	The function returns the difference between local time and UTC in days
 	east of greenwhich (west is negative) determined from the computer system
 	clock at the current time and location.
 
*******************************************************************************/
	
double ASTROLIBDLL_API AACurrentTimeZone ( int *isDST );

/*** AADateTimeToJD ********************************************************

	Converts a local date and time in the Gregorian or Julian calendar
	to a Julian Day number.
               
	double AADateTimeToJD ( int year, short month, double day, short hour,
	short min, double sec, double zone, short calendar )

	    (year): Year.
	   (month): Month.
	     (day): Day, optionally including fractional part.
	    (hour): Hour.
	     (min): Minute.
	     (sec): Second.
	    (zone): Time zone (local time - UTC), as fraction of a day.
	(calendar): Calendar system, AA_CALENDAR_JULIAN ... AA_CALENDAR_MAYAN.
               
	Returns Julian Day number corresponding to the calendar date.
	
	For the Mayan calendar, the Long Count should be passed as the year,
	and the Haab date should be passed as the month and day.  Only the
	fractional part of the date is used, since the Long Count completely
	specifies a day.  (The Haab "month" value is actually ignored.)
	
	If the calendar system provided (calendar) is not recognized, the date
	will be assumed to be in the Gregorian calendar if after 4 October 1582,
	and in the Julian calendar if before 15 October 1482.
      
	References:
	
	Jean Meeus, "Astronomical Algorithms", p. 61, 71-76.
	The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	http://www.fourmilab.ch/documents/calendar/	
	
*****************************************************************************/

#define AA_CALENDAR_JULIAN				0
#define AA_CALENDAR_GREGORIAN			1
#define AA_CALENDAR_GREGORIAN_JULIAN	2
#define AA_CALENDAR_HEBREW				3
#define AA_CALENDAR_ISLAMIC				4
#define AA_CALENDAR_PERSIAN				5
#define AA_CALENDAR_INDIAN				6
#define AA_CALENDAR_CHINESE				7
#define AA_CALENDAR_MAYAN				8

double ASTROLIBDLL_API AADateTimeToJD ( int, short, double, short, short, double, double, short );

/*** AAJDToDateTime ***************************************************************

	Converts a Julian date to a local date and time in the Gregorian or Julian
	calendar.
               
	AAJDToDateTime ( double jd, double zone, int *year, short *month, double *day,
	short *hour, short *min, double *sec, short calendar )

	      (jd): Julian Day number.
	    (zone): Time zone (local time - UTC), as a fraction of a day.
	    (year): Receives year.
	   (month): Receives month.
	     (day): Receives day, including fractional part.
	    (hour): Receives hour.
	     (min): Receives minute.
	     (sec): Receives second.
	(calendar): Calendar system, AA_CALENDAR_JULIAN ... AA_CALENDAR_MAYAN.
               
	Returns nothing.  The algorithm presented in the reference has been modified
	to provide correct results for negative Julian dates as well as positive ones,
	unlike the original algorithm.
	
	For the Mayan calendar, the Long Count date corresponding to the Julian Date
	is returned as the year.  The Haab date is returned as the month and date.
	
	If the calendar system provided (calendar) is not recognized, the date
	will be returned in the Gregorian calendar if after JD 2299160.5, and in
	the Julian calendar if before then.

	References:
	
	Jean Meeus, "Astronomical Algorithms", p. 63, 71-76.
	The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	http://www.fourmilab.ch/documents/calendar/
	
**********************************************************************************/

void ASTROLIBDLL_API AAJDToDateTime ( double, double, int *, short *, double *, short *, short *, double *, short );

/*** AAGregorianToJD ********************************************************

	Converts a Gregorian calendar date to a Julian Day number.
               
	double AAGregorianToJD ( int y, int m, double d )

	(y): Year.
	(m): Month.
	(d): Day, optionally including fractional part.
               
	Returns the Julian Day number.  The algorithm presented in the source
	reference is valid for all Gregorian calendar dates corresponding to
	JD >= 0, i.e. dates after -4713 November 23, but has been verified to
	produce valid results for all year values.
      
	Reference:
	
	The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	
******************************************************************************/

double ASTROLIBDLL_API AAGregorianToJD ( int, int, double );

/*** AAJDToGregorian ************************************************************

	Converts a Julian Day number to a date in the Gregorian calendar.
               
	void AAJDToGregorian ( double jd, int *y, short *m, double *d )

	(jd): Julian Day number.
	 (y): Receives year.
	 (m): Receives month.
	 (d): Receives day, including fractional part.
               
	Returns nothing.  The algorithm presented in the source reference is valid
	for all Gregorian calendar dates corresponding to JD >= 0, i.e. dates after
	-4713 November 23, but has been modified to produce valid results for all
	Julian dates, both positive and negative.
      
	Reference:
	
	The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	
**********************************************************************************/

void ASTROLIBDLL_API AAJDToGregorian ( double, int *, short *, double * );

/*** AAJulianToJD ***********************************************************

	Converts a Julian calendar date to a Julian Day number.
               
	double AAJulianToJD ( int y, int m, double d )

	(y): Year.
	(m): Month.
	(d): Day, optionally including fractional part.
               
	Returns the Julian Day number.  The algorithm presented in the source
	reference is valid for all Julian calendar dates corresponding to
	JD >= 0, i.e. years after -4712, but has been modified to be valid
	for all Julian dates (both positive and negative).
      
	References:
	
	The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	
******************************************************************************/

double ASTROLIBDLL_API AAJulianToJD ( int, int, double );

/*** AAJDToJulian ****************************************************************

	Converts a Julian date to a date in the Julian calendar.
               
	void AAJDToJulian ( double jd, int *y, short *m, double *d )

	(jd): Julian Day number.
	 (y): Receives year.
	 (m): Receives month.
	 (d): Receives day, including fractional part.
               
	Returns nothing.  The algorithm presented in the source reference is valid
	for all Julian calendar dates corresponding to JD >= 0, i.e. years after
	-4712.
      
	References:
	
	The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	
**********************************************************************************/

void ASTROLIBDLL_API AAJDToJulian ( double, int *, short *, double * );

/*** AAHebrewToJD ***********************************************************

	Converts a Hebrew calendar date to a Julian Day number.
               
	double AAHebrewToJD ( int y, short m, double d )

	(y): Year.
	(m): Month.
	(d): Day, optionally including fractional part.
               
	Returns the Julian Day number.  The algorithm is based on JavaScript code
	from www.fourmilab.ch, but has been modified to include the fractional part
	of the day properly, and has been verified to produce valid results for
	all Hebrew year values (both positive and negative).
	
	References: 
	
	http://www.fourmilab.ch/documents/calendar/
	Jean Meeus, "Astronomical Algorithms", pp. 71-76.
	The Explanatory Supplement to the Astronomical Almanac, pp. 584-589.
	
******************************************************************************/

double ASTROLIBDLL_API AAHebrewToJD ( int, short, double );

/*** AAJDToHebrew ***********************************************************

	Converts a Julian date to a date in the Hebrew calendar.
               
	void AAJDToHebrew ( double jd, int *y, short *m, double *d )

	(jd): Julian Day number.
	 (y): Receives year.
	 (m): Receives month.
	 (d): Receives day, including fractional part.
               
	Returns nothing.  The algorithm is based on JavaScript code from
	www.fourmilab.ch, but has been modified to include the fractional part
	of the day properly, and has been verified to produce valid results for
	all Julian Dates (both positive and negative).
	      
	References:
	
	http://www.fourmilab.ch/documents/calendar/
	Jean Meeus, "Astronomical Algorithms", pp. 71-87.
	The Explanatory Supplement to the Astronomical Almanac, pp. 584-589.
	
******************************************************************************/

void ASTROLIBDLL_API AAJDToHebrew ( double, int *, short *, double * );

/*** AAIslamicToJD ***********************************************************

	Converts an Islamic tabular calendar date to a Julian Day number.
               
	double AAIslamicToJD ( int y, int m, double d )

	(y): Year.
	(m): Month.
	(d): Day, optionally including fractional part.
               
	Returns the Julian Day number.  The algorithm presented in the source
	reference is valid for all years >= 1, corresponding to JD >= 1948440
	but has been modified to be valid for all year values (both positive
	and negative).
      
	References:
	
	The Explanatory Supplement to the Astronomical Almanac, pp. 589-591.
	Jean Meeus, "Astronomical Algorithms", pp. 71-76.
	
******************************************************************************/

double ASTROLIBDLL_API AAIslamicToJD ( int, int, double );

/*** AAJDToIslamic ************************************************************

	Converts a Julian date to a date in the Islamic tabular calendar.
               
	void AAJDToIslamic ( double jd, int *y, short *m, double *d )

	(jd): Julian Day number.
	 (y): Receives year.
	 (m): Receives month.
	 (d): Receives day, including fractional part.
               
	Returns nothing.  The algorithm presented in the source reference is valid
	for all years >= 1, corresponding to JD >= 1948440, but has been modified
	to produce valid results for all year values.
	
	References:
	
	Jean Meeus, "Astronomical Algorithms", pp. 71-76.
	The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	
**********************************************************************************/

void ASTROLIBDLL_API AAJDToIslamic ( double, int *, short *, double * );

/*** AAPersianToJD ***********************************************************

	Converts a Persian calendar date to a Julian Day number.
               
	double AAPersianToJD ( int y, short m, double d )

	(y): Year.
	(m): Month.
	(d): Day, optionally including fractional part.
               
	Returns the Julian Day number.  The algorithm is based on JavaScript code
	from www.fourmilab.ch, but has been modified to include the fractional part
	of the day properly, and to treat years before 1 as 0, -1, -2, ... instead
	of 1, -1, -2 ...  It has been verified to produce valid results for
	all Persian year values (both positive and negative).
	
	References: 
	
	http://www.fourmilab.ch/documents/calendar/
	
******************************************************************************/

double ASTROLIBDLL_API AAPersianToJD ( int, short, double );

/*** AAJDToPersian ***********************************************************

	Converts a Julian date to a date in the Persian calendar.
               
	void AAJDToPersian ( double jd, int *y, short *m, double *d )

	(jd): Julian Day number.
	 (y): Receives year.
	 (m): Receives month.
	 (d): Receives day, including fractional part.
               
	Returns nothing.  The algorithm is based on JavaScript code from
	www.fourmilab.ch, but has been modified to include the fractional part
	of the day properly, and to treat years before 1 as 0, -1, -2, ...
	instead of 1, -1, -2 ...

	It has been verified to produce valid results for all Julian dates
	(both positive and negative).
	      
	References:
	
	http://www.fourmilab.ch/documents/calendar/
	
******************************************************************************/

void ASTROLIBDLL_API AAJDToPersian ( double, int *, short *, double * );

/*** AAIndianToJD ***********************************************************

	Converts an Indian civil calendar date to a Julian Day number.
               
	double AAIndianToJD ( int y, int m, double d )

	(y): Year.
	(m): Month.
	(d): Day, optionally including fractional part.
               
	Returns the Julian Day number.  The algorithm presented in the source
	reference is valid for all years >= 1, corresponding to JD >= 17499995,
	but has been modified to be valid for both positive and negative years.
      
	Reference: The Explanatory Supplement to the Astronomical Almanac,
	           pp. 603-606.
	
******************************************************************************/

double ASTROLIBDLL_API AAIndianToJD ( int, int, double );

/*** AAJDToIndian ************************************************************

	Converts a Julian date to a date in the Indian civil calendar.
               
	void AAJDToJulian ( double jd, int *y, short *m, double *d )

	(jd): Julian Day number.
	 (y): Receives year.
	 (m): Receives month.
	 (d): Receives day, including fractional part.
               
	Returns nothing.  The algorithm presented in the source reference is valid
	for all years >= 1, or jd >= 1749995.
      
	Reference: The Explanatory Supplement to the Astronomical Almanac, pp. 603-606.
	
**********************************************************************************/

void ASTROLIBDLL_API AAJDToIndian ( double, int *, short *, double * );

/*** AAMayanToJD ***********************************************************

	Converts a Mayan Long Count date to a Julian Day number.
               
	double AAMayanLongCountToJD ( int lcount );
	double AAMayanToJD ( short baktun, short katun, short tun, short uinal, short kin )

	(lcount): Long Count sum.
	(baktun): Baktun component of long count.
	(katun):  Katun component of long count.
	(tun):    Tun component of long count.
	(uinal):  Uinal component of long count.
	(kin):    Kin component of long count.
               
	Returns the Julian Day number.  The algorithm is based on JavaScript code
	from www.fourmilab.ch, and has been verified to produce valid results for
	all Mayan long count values (both positive and negative).  Note that these
	functions discard the fractional part of the Julian Date and return the
	long count as an integer.
	
	References: 
	
	http://www.fourmilab.ch/documents/calendar/
	
******************************************************************************/

double ASTROLIBDLL_API AAMayanLongCountToJD ( int );
double ASTROLIBDLL_API AAMayanToJD ( short, short, short, short, short );

/*** AAJDToMayan ***********************************************************

	Converts a Julian date to a Mayan Long Count date.
               
	long AAJDToMayanLongCount ( double jd );
	void AAJDToMayan ( double jd, short *baktun, short *katun, short *tun, short *uinal, short *kin )
	void AAJDToMayanHaab ( double jd, short *haab, double *day );

	(jd):     Julian Date.
	(baktun): Receives baktun part of long count.
	(katun):  Receives katun part of long count.
	(tun):    Receives un part of long count.
	(uinal):  Receives uinal part of long count.
	(kin):    Receives kin part of long count.
	(haab):   Receives Haab "month".
	(day):    Receives day within Haab, including fractional part.
	
	The first function returns the long count sum.  The second returns
	nothing, but computes the long count parsed into its component parts
	(baktun...kin).  The third function computes the 20-day Haab "month"
	and day within the Haab that the specified Julian date corresponds to.
	Note that the day is counted from 0-19, not 1-20.
		
	The algorithm used by both functions is based on JavaScript code
	from www.fourmilab.ch.  It has been modified to produce valid results
	for all Julian dates (both positive and negative).
	      
	References:
	
	http://www.fourmilab.ch/documents/calendar/
	
******************************************************************************/

int	 ASTROLIBDLL_API AAJDToMayanLongCount ( double );
void ASTROLIBDLL_API AAJDToMayan ( double, short *, short *, short *, short *, short * );
void ASTROLIBDLL_API AAJDToMayanHaab ( double, short *, double * );

/*** AJulianYearToJD  ***************************************************

	Converts a Julian or Besselian epoch to a Julian date.
  
	double AAJulianYearToJD ( double year )
	double AABesselianYearToJD ( double year )
	
	(year): The Julian or Besselian year.

	Returns the Julian date corresponding to the specified Julian
	or Besselian year.
      
	Reference: the Astronomical Almanac for the Year 1990, p. B4.
   
*************************************************************************/

double ASTROLIBDLL_API AAJulianYearToJD ( double );
double ASTROLIBDLL_API AABesselianYearToJD ( double );

/***  AAJDToJulianYear  *****************************************************

	Converts a Julian date to a Julian or Besselian year.
  
	double AAJDToJulianYear ( double jd )
	double AAJDToBesselianYear ( double jd )
	
	(jd): The Julian date.
	
	Returns the Julian or Besselian year, as desired.
      
	Reference: the Astronomical Almanac for the Year 1990, p. B4.

***************************************************************************/

double ASTROLIBDLL_API AAJDToJulianYear ( double );
double ASTROLIBDLL_API AAJDToBesselianYear ( double );

/*** ALocalWeekDay *******************************************************

	Given a Julian date, computes the day of the week in a particular
	local time zone.

	int AALocalWeekDay ( double jd, double zone )
	
	(jd):   Julian Date, expressed in Universal time.
	(zone): Local time zone, expressed as a fraction of a day east of
	        Greenwich; east time zones are positive.
	
	The function returns an integer from 0-7, indicating the day of
	the week: 0 = sunday, 1 = monday, 2 = tuesday ... 6 = saturday.
	
	Reference: Jean Meeus, "Astronomical Algorithms", p. 65
		
***************************************************************************/

int ASTROLIBDLL_API AALocalWeekDay ( double, double );

/*** AADaylightSavingsTime **************************************************

	Determines whether a Julian date falls within Daylight Saving Time
	in a particular local time zone.

	int AADaylightSavingsTime ( double jd, double zone, int usa,
	    double *start, double *end )
	
	(jd):    Julian Date, expressed in Universal time.
	(zone):  Local time zone, expressed as a fraction of a day east of
	         Greenwich; east time zones are positive.
	(rule):  specifies the daylist savings time rules to use.
	(start): receives JD of start of daylight saving time.
	(end):   receives JD of end of daylight saving time.
	
	The function returns TRUE if the given JD lies within daylight saving
	time, and FALSE otherwise.  The variables (start) and (end) receive the
	Julian dates (expressed in Universal time) corresponding to the
	beginning and end of Daylight Saving Time in the specified local time
	zone, for the same year as the given JD.
	
	The daylight saving time rules to use are specified in the (rule)
	parameter, which can be one of the values #defined below (AA_DLST_USA_CANADA
	... AA_DLST_AUSTRALIA).  For example, in the United States and Canada
	(AA_DLST_USA_CANADA),  Daylight Saving Time begins at 02:00:00 local time
	on the first Sunday in April, and ends at 02:00:00 local time on the last
	Sunday in October.  In 2007, the USA and Canada switched to using the
	second Sunday in March and the first Sunday in November.  Mexico continues
	to use the same rules as for the pre-2007 USA, so if (rule) is AA_DLST_MEXICO,
	the function will not change behavior after 2006.  Similar logic for the
	other regional daylight saving time rules are implemented in the function.
	
	Note that for countries in the Southern Hemisphere, DLST rules may imply
	a start date for DLST that is after the end date in the same calendar year.
	For example, in Australia, DLST starts on the last Sunday in October and
	ends on the last Sunday in March.
	
	References: 
	
	Wikipedia, "Daylight Saving Time Around the World", 2007
	(http://en.wikipedia.org/wiki/Daylight_saving_time_around_the_world)
	
	WebExhibits, "When is Daylight Saving Time Worldwide?", 2007
	(http://webexhibits.org/daylightsaving/g.html)
			
***************************************************************************/

#define AA_DLST_ALWAYS			-1
#define AA_DLST_NONE			0
#define AA_DLST_USA_CANADA		1
#define AA_DLST_MEXICO			2
#define AA_DLST_BRAZIL			10
#define AA_DLST_CHILE			11
#define AA_DLST_PARAGUAY		12
#define AA_DLST_EUROPE			20
#define AA_DLST_RUSSIA			21
#define AA_DLST_IRAN			22
#define AA_DLST_SYRIA			23
#define AA_DLST_AUSTRALIA		30
#define AA_DLST_NEW_ZEALAND		31
#define AA_DLST_FIJI			32
#define AA_DLST_MOROCCO			40
#define AA_DLST_NAMIBIA			41

int ASTROLIBDLL_API AADaylightSavingsTime ( double, double, int, double *, double * );

/***  AADeltaT  *************************************************************

	Computes the value of Delta T, the difference between Universal Coordinated
	Time (UTC) and Terrestrial Dynamic Time (TDT) or Ephemeris Time (ET) for any
	given date between the years -1999 and +3000.
 
	double AADeltaT ( double jd )
	
	(jd): Julian date at which to compute Delta T.
 
	Returns the value of TDT - UTC in seconds at the specified Julian date
	(jd).  This should be added to UTC to obtain TDT.
	
	Prior to 1984, TDT is referred to as ET.
	
	Delta T is computed using a series of polynomial expressions derived
	from the historical record and from direct observations.  The values
	returned by this function differ from published observations of Delta
	T by less than 1 second over the interval from 1700 to 2010, and when
	extrapolated forward, should differ by less than 10 seconds through
	the end of the 21st century.
 
	References:
 
	Espenak, F. and Meeus, J. "Five Millennium Canon of Solar Eclipses:
	-1999 to +3000" (NASA Tecnical Publication 2006-214141); see also
	http://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
 
	The Astronomical Almanac for the Years 1984-2004, p. K8
 
***************************************************************************/

double AADeltaT ( double );

/*** AALocalSiderealTime **************************************************
	
	Calculates Greenwich mean sidereal time and an observer's local
	mean sidereal time.
               
	double AAGreenwichMeanSiderealTime ( double jd );
	double AALocalMeanSiderealTime ( double jd, double lon )
	
	(jd):  Julian date expressed in universal time (UT) not dynamic time.
    (lon): observer's longitude, in radians; positive to the east.
        
	Returns the observer's local mean sidereal time at the time specified
	by the julian date (jd).  The result will be returned in radians in
	the range 0 to TWO_PI.

	To obtain Greenwich mean sidereal time, use the first function or
 	pass a value of zero for the longitude parameter.
	
	To convert mean sidereal time to true sidereal time, you should add
	the value of the nutation in longitude computed for the same Julian
	date; see Nutation().
	
	Reference: the Astronomical Almanac for the Year 1984.
	
***************************************************************************/

double ASTROLIBDLL_API AAGreenwichMeanSiderealTime ( double );
double ASTROLIBDLL_API AALocalMeanSiderealTime ( double, double );

/*** AASemiDiurnalArc *****************************************************

	Computes the angular distance from rising or setting to transit for an
	object with any given declination as seen from any latitude.
        
	double SemiDiurnalArc ( double a, double d, double f);

	(a): The object's altitude when rising/setting, in radians.
	(d): The object's declination, in radians.
	(f): The observer's latitude, in radians.       
               
	Returns the object's hour angle when it is at altitude (a) as viewed from
	latitude (f), in radians.  If the object's altitude is always greater than
	(a), SemiDiurnalArc() returns PI.  If the altitude is always less than (a),
	it returns zero.
  
****************************************************************************/
   
double ASTROLIBDLL_API AASemiDiurnalArc ( double, double, double );

/*** AARiseSetTime *****************************************************

	Given an object's celestial coordinates at a particular instant,
	computes the object's time of an rising, transit, or setting that
	is closest to that instant, as seen from a particular location.

	double AARiseSetTime ( double ra, double dec, double jd, int sign,
	double lon, double lat, double alt )

	(ra):   object's right ascension, in radians.
	(dec):  object's declination, in radians.
	(jd):   Julian date for which object's coordinates are given.
	(sign): event to compute: -1 = rising, 0 = transit, +1 = setting.
	(lon):  observer's east longitude, in radians.
	(lat):  observer's north latitude, in radians.
	(alt):  observer's local horizon altitude, in radians.
	
	The function returns a Julian date which represents the object's rise/
	transit/set time that is closest to the time specified for the object's
	coordinates (jd).  This will always be within half a day of the value
	passed in the (jd) parameter, with two exceptions: if the object does
	not set below the horizon as seen from the specified location, the
	function returns INFINTY; and if the object does rise above the horizon
	at the specified location, it returns -INFINITY.
	
	The object's equatorial coordinates, (ra) and (dec), should be given
	for their current precessional epoch (jd).
	
	The (alt) parameter specifies the altitude of the local horizon, in
	radians; the rise and set times returned by this function correspond
	to the instant that the object's coordinates lie at this altitude.
	For point objects, a horizon altitude of -0.5 degrees is typically
	used to account for atmospheric refraction.  For the sun and moon,
	which are extended objects, a horizon altitude of -50 arcminutes
	is conventionally used to compute rise/set times.  The times of
	civil, nautical, and astronomical twilight can be computed by using
	a horizon altitude of -6, -12, and -18 degrees, respectively.
	
	Note that this function does not take into account the object's motion
	over the course of the day.  For objects, like stars, whose equatorial
	coordinates remain relatively fixed over the course of a day, this should
	pose no problem.  However, for fast-moving objects like the sun and moon,
	you should call this function iteratively, recomputing the object's
	coordinates at the rise/set time predicted by the function on each
	iteration.
	
*******************************************************************************/

double ASTROLIBDLL_API AARiseSetTime ( double, double, double, int, double, double, double );

/*** AARiseSetTimeSearch ****************************************************

	Computes the time of an object's rise, transit, or set that is closest
	to an initial starting date.
	
	double AARiseSetTimeSearch ( RiseSetProcPtr proc, void *data, double jd,
	int sign, double lon, double lat, double alt, double precis, int imax )

	(proc):   pointer to a function which compute the object's position.
	(data):   pointer to data structure containing parameters to pass to proc.
	(jd):     starting date, expressed as a julian date.
	(sign):   event to find: -1 = rise, 0 = transit, +1 = set.
	(lon):    observer's east longitude, in radians.
	(lat):    observer's north latitude, in radians.
	(alt):    observer's local horizon altitude, in radians.
	(precis): desired precision of result, as a fraction of a day.
	(imax):   maximum number of iterations to perform.
	
	The function returns the  Julian date of the object's rising, transit,
	or setting that is closest to the initial starting date (jd), to the
	specified precision.  It works by iterately computing the object's
	position and using it to predict increasingly accurate estimates
	of the object's rise/transit/set time.  The maximum number of
	iterations to perform is specified by the (imax) parameter.  If the
	object reaches a position from which it does not rise or does not
	set on the specified date, the function returns -INFINITY or
	INFINITY, respectively.
	
	The (sign) parameter specifies the phenomenon to search for: -1
	specifies a rise time, 0 specifies a transit time, and +1 specifies
	a set time.
	
	The (proc) parameter is a pointer to a function which computes the
	object's right ascension and declination on a given Julian date.
	This function has the following prototype:
	
	void RiseSetSearchProc ( double jd, double *ra, double *dec, void *data )
	
	The Julian date at which the object's position should be computed is
	given in the (jd) parameter; its right ascension and declination
	should be returned in the (dec) parameter.  For highest accuracy,
	the right ascension and declination should be computed for the same
	precessional epoch as (jd), and should be corrected for nutation,
	aberration, light time, and geocentric parallax at the observer's
	position.  The (data) parameter is a pointer to a caller-defined
	structure which contains any additional information needed to compute
	the object's apparent position (for example, the object's identifier,
	the observer's location, etc.)

******************************************************************************/

typedef void (*AARiseSetSearchProcPtr) ( double, double *, double *, void * );

double ASTROLIBDLL_API AARiseSetTimeSearch ( AARiseSetSearchProcPtr, void *, double, int,
double, double, double, double, int );

/*** AADailyRiseSetTimeSearch ************************************************

	Computes the time an object's rises, transits, or sets on a particular
	local day.
	
	double AADailyRiseSetTimeSearch ( RiseSetProcPtr proc, void *data, double jd,
	int sign, double lon, double lat, double alt, double precis, int imax )

	(proc):   pointer to a function which compute the object's position.
	(data):   pointer to data structure containing parameters to pass to proc.
	(jd):     starting date, expressed as a julian date.
	(sign):   event to find: -1 = rise, 0 = transit, +1 = set.
	(zone):   observer's local time zone, days east of Greenwich.
	(lon):    observer's east longitude, in radians.
	(lat):    observer's north latitude, in radians.
	(alt):    observer's local horizon altitude, in radians.
	(precis): desired precision of result, as a fraction of a day.
	(imax):   maximum number of iterations to perform.
	
	The function returns the  Julian date of the object's rising, transit,
	or setting that is closest to the starting time (jd) on a particular
	local day.  It works by iterately computing the object's position and
	using it to predict increasingly accurate estimates of the object's
	rise/transit/set time, until the estimate of the time converges to the
	specified precision.  The maximum number of iterations to perform is
	specified by the (imax) parameter.  If the object does not rise or set
	on the specified local day, the function returns -INFINITY or INFINITY,
	respectively.
	
	This function is guaranteed to work only for objects which rise
	and set at most once per day.  It will not work for objects which
	rise and set multiple times, e.g. artifical earth satellites.
	
	The (sign) parameter specifies the phenomenon to search for: -1
	specifies a rise time, 0 specifies a transit time, and +1 specifies
	a set time.
	
	The (proc) parameter is a pointer to a function which computes the
	object's right ascension and declination on a given Julian date.
	See RiseSetTimeSearch() for a complete description of this function.

******************************************************************************/

double ASTROLIBDLL_API AADailyRiseSetTimeSearch ( AARiseSetSearchProcPtr, void *, double, int,
double, double, double, double, double, int );

/*** Functions in Reduce.c ***********************************************/

/*** AAPrecession *********************************************************

	Computes the precessional constants.
               
	void AAPrecession ( double jd0, double jd1, double *zeta, double *z,
	double *theta, double *eta, double *pi, double *p )

	  (jd0): Julian date of the initial epoch.
	  (jd1): Julian date of the final epoch.
	 (zeta): Receives Newcomb's precessional constant zeta, in radians.
	    (z): Receives Newcomb's precessional constant z, in radians.
	(theta): Receives Newcomb's precessional constant theta, in radians.
	  (eta): Receives angle of inclination between the two ecliptics.
	   (pi): Receives angle from initial equinox to intersection of ecliptics.
	    (p): Receives general precession in ecliptic longitude.

	Returns nothing.  To supress computation of any of the above constants,
	pass NULL instead of a pointer.
      
	References: Explanatory Supplement to the Astronomical Almanac, p. 104.
	            Jean Meeus, "Astronomical Algorithms", p. 126-128.
               
***************************************************************************/
                
void ASTROLIBDLL_API AAPrecession ( double, double, double *, double *, double *, double *, double *, double * );

/*** AALongTermPrecession **************************************************

	Computes the precessional constants and obliquity of the ecliptic, valid
	for a timespan of +/-500,000 years from J2000.
               
	void AALongTermPrecession ( double jd, double *chi, double *omega,
	                            double *psi, double *epsilon )

	     (jd): Julian date at which to compute angles.
	    (chi): Receives accumulated planetary precession, in radians.
	  (omega): Receives inclination of equator of date to ecliptic of J2000, in radians.
	    (psi): Receives accumulated lunisolar precession, in radians.
	(epsilon): Receives obliquity of ecliptic, in radians.

	Returns nothing.  To supress computation of any of the above constants,
	pass NULL instead of a pointer.
      
	Reference: William M. Owen, Jr. (wmo@jpl.nasa.gov), "A Theory of the
	           Earth's Precession Relative to the Invariable Plane of the
	           Solar System", Ph.D. dissertation, University of Florida, 1990.
               
***************************************************************************/

void ASTROLIBDLL_API AALongTermPrecession ( double, double *, double *, double *, double * );

/*** AASetLongTermPrecessionMatrix ********************************************

	Computes a rotation matrix for precessing equatorial coordinates using
	long-term expessions for the precessional constants valid for a span of
	+/- 500,000 years from J2000.
  
	void AASetLongTermPrecessionMatrix ( double m[3][3], double jd0,
	     double jd1, int n )

	  (m): array of 9 doubles to recieve rows of the matrix.
	(jd0): Julian date of the initial epoch.
	(jd1): Julian date of the final epoch.
	  (n): if non-zero, the matrix will also contain the effects of nutation.
               
	Returns nothing.  If nutation is included, the matrix returned will
	transform a vector from the the mean equatorial frame of epoch JD0
	to the true equatorial frame of epoch JD1.  If not, the matrix will
	only transform to the mean frame of epoch JD1.      
      
	Reference: William M. Owen, Jr. (wmo@jpl.nasa.gov), "A Theory of the
	           Earth's Precession Relative to the Invariable Plane of the
	           Solar System", Ph.D. dissertation, University of Florida, 1990.
               
***************************************************************************/

void ASTROLIBDLL_API AASetLongTermPrecessionMatrix ( double [3][3], double, double, int );

/*** AAObliquity ********************************************************

	Computes the mean obliquity of the ecliptic for a given epoch.

	double AAObliquity ( double jd )
               
	(jd): Julian Day number of the epoch.
         
	Returns the mean obliquity of the ecliptic for epoch (jd), in radians.
               
	References: The Astronomical Almanac for the Year 1984.

*************************************************************************/

double ASTROLIBDLL_API AAObliquity ( double );

/*** AANutation **********************************************************

	Computes the nutation in longitude and in obliquity for a given epoch.

	void AANutation ( double jd, double *dl, double *de )
               
	(jd): Julian Day number of the epoch.
	(dl): Receives nutation in longitude, in radians.
	(de): Receives nutation in obliquity, in radians.
         
	Returns nothing.  Values are computed using the 1980 IAU theory of
	nutation.  Only linear terms in the fundamental arguments, and the
	four largest periodic terms, have been retained, providing an accuracy
	for (dl) and (de) of 0".5 and 0".1.
               
	References: Explanatory Supplement to the Astronomical Almanac, pp. 109-116.
                Jean Meeus, "Astronomical Algorithms", pp. 131-134.

***************************************************************************/

void ASTROLIBDLL_API AANutation ( double, double *, double * );

/*** AAEarthVelocity *****************************************************

	Computes the earth's heliocentric velocity.

	AAEarthVelocity ( double jd, double *vx, double *vy, double *vz )
  
	      (jd): Julian date on which to compute velocity.
	(vx,vy,vz): Recieves components of the earth's velocity vector
                with respect to the solar system barycenter in the
                J2000 equatorial frame, in units of AU per day.
               
	Returns nothing.
               
	References: Explanatory Supplement to the Astronomical Almanac,
	            pp. 130, 484.
   
**************************************************************************/

void ASTROLIBDLL_API AAEarthVelocity ( double, double *, double *, double * );

/*** AARelativisticAberration ********************************************

	Applies aberration of light to an object's geometric position vector.
  
	double AARelativisticAberration ( double u[3], double v[3], double w[3] )

	(u): Object's geometric position vector relative to viewer.
	(v): Observer's velocity vector, in units of the speed of light.
    (w): Recieves apparent direction vector.
    
    The function returns the object's doppler shift.  To obtain the
    observed frequency of electromagnetic radiation emitted from the
    object, divide the true frequency emitted at the object by the
    doppler shift; similarly, to obtain the observed wavelenth,
    multiply the true wavelength by this value.
    
    The input vector (u) represents the object's geometric position
    relative to the viewer.  The vector (v) represents the viewer's
    velocity relative to the origin of the coordinate system in which
    both the viewer's and object's positions are measured.  While
    vector (u) can take any units and have any magnitude, vector (v)
    must be given in units of the speed of light; its magnitude must
    therefore be less than 1.
    
    The vector returned in (w) will have the same magnitude as the
    original vector (u), but its direction will be corrected for
    aberration of light.  If (w) points to the same array as (u),
    the original vector will be overwritten.
    
    This function is adapted from the relativistic expression for
    aberration given in the Explanatory Supplement.
              
	References:
	
	The Explanatory Supplement to the Astronomical Almanac, p. 129.
	
***************************************************************************/

double ASTROLIBDLL_API AARelativisticAberration ( double *, double *, double * );

/*** AAGeocentricVelocity  ************************************************
	 
	Computes velocity of a point on Earth's surface relative to Earth's
 	center due to Earth's diurnal rotation.
	 
	void AAGeocentricVelocity ( double l, double p, double h, double a, double f, double w,
	double *x, double *y, double *z )
	 
	(l): Geodetic longitude, in radians.
	(b): Geodetic latitude, in radians.
	(h): Geodetic height.
	(a): Equatorial radius of reference ellipsoid, in same units as (h).
	(f): Flattening of reference ellipsoid: f = ( a - b ) / a,
	     where a and b are ellipsoid's equatorial and polar radii.
	(w): Earth's sidereal angular rotation rate, radians per day.
	(x): Receives geocentric X velocity component, in same units as (h) per day.
	(y): Receives geocentric Y velocity component, in same units as (h) per day.
	(z): Receives geocentric Z velocity component, in same units as (h) per day.
	 
	Returns nothing.
	The values for a, f, and w adopted in the 1976 IAU reference ellipsoid
	are a = 6378140 m, f = 1 / 298.257, and w = 6.300387487008 rad/day
	 
	Reference: Explanatory Supplement to the Astronomical Almanac, pp. 132-133.
     
*********************************************************************************/

void ASTROLIBDLL_API AAGeocentricVelocity ( double l, double p, double h, double a, double f, double w,
						   double *x, double *y, double *z );

/*** AAGeodeticToGeocentricXYZ  ************************************************

	Converts geodetic latitude, longitude, and height to geocentric rectangular
	coordinates.
               
	void AAGeodeticToGeocentricXYZ ( double l, double b, double h, double a,
	double f, double *x, double *y, double *z )

	(l): Geodetic longitude, in radians.
	(b): Geodetic latitude, in radians.
	(h): Geodetic height.
	(a): Equatorial radius of reference ellipsoid, in same units as (h).
	(f): Flattening of reference ellipsoid: f = ( a - b ) / a,
	     where a and b are ellipsoid's equatorial and polar radii. 
	(x): Receives geocentric X-coordinate, in same units as (h).
	(y): Receives geocentric Y-coordinate, in same units as (h).
	(z): Receives geocentric Z-coordinate, in same units as (h).

	Returns nothing.  Geocentric longitude, latitude, and radius can be
	computed by converting (x), (y), and (z) to spherical coordinates.
	The values for a and f adopted in the 1976 IAU reference ellipsoid
	are a = 6378140 m and f = 1 / 298.257.
               
    Reference: The Astronomical Almanac for the Year 1990, pp. K11-K13.
     
*********************************************************************************/

void ASTROLIBDLL_API AAGeodeticToGeocentricXYZ (double, double, double, double, double,
double *, double *, double * );

/*** AAGeocentricXYZToGeodetic **************************************************

	Converts geocentric rectangular coordinates to geodectic longitude, latitude,
	and height.
               
	void AAGeocentricXYZToGeodetic ( double x, double y, double z,
	double a, double f, double *l, double *b, double *h )

	(x): Geocentric X-coordinate.
	(y): Geocentric Y-coordinate.
	(z): Geocentric Z-coordinate.
	(a): Equatorial radius of reference ellipsoid.
	(f): Flattening of reference ellipsoid: f = ( a - b ) / a,
	     where a and b are ellipsoid's equatorial and polar radii. 
	(l): Recieves geodetic longitude, in radians from -PI to PI.
	(b): Recieves geodetic latitude, in radians from -PI/2 to PI/2.
	(h): Receives geodetic height, in same units as (a).
          
	Returns nothing.  The values for a and f adopted in the 1976 IAU reference
	ellipsoid are a = 6378140 m and f = 1 / 298.257.
     
	Reference: The Astronomical Almanac for the Year 1990, pp. K11-K13.
     
********************************************************************************/

void ASTROLIBDLL_API AAGeocentricXYZToGeodetic ( double, double, double, double, double,
double *, double *, double * );

/*** AARefractionAngle ***************************************************************

	Computes altitude correction due to atmospheric refraction.
  
	double AARefractionAngle ( double h, double p, double t, int a )

	(h): The object's altitude in radians.
	(p): Barometric pressure in millibars (1010 is standard).
	(t): Termperature in degrees centigrade (10 is standard).
	(a): If true, h is true altitude; otherwise h is apparent altitude.
 
	Returns refraction correction, in radians.  If a is true (i.e. h is true
	altitude) this should be added to h to obtain the apparent altitude;
	otherwise the result should be subtracted from the apparent altitude h
	to obtain the true altitude.
 
	Reference: Jean Meeus, "Astronomical Algorithms", pp. 105 - 108.
	
*********************************************************************************/

double ASTROLIBDLL_API AARefractionAngle ( double, double, double, int );

/*** AAHorizonToEquatorial ***************************************************
 
	Converts horizon coordinates to equatorial.
 
	void AAHorizonToEquatorial ( double azm, double alt, double lon, double lat,
	double jd, double *ra, double *dec )
 
	(azm): Azimuth in radians, measured from North = 0, East = pi / 2.
	(alt): Altitude in radians, from horizon 0 to zenith = pi / 2.
	(lon): Geographic longitude in radians, east is positive.
	(lat): Gepographc latitude in radians, north is positive.
	(jd):  Julian Date
	(ra):  Receives right ascension in radians from 0 to 2 * pi
	(dec): Receives declination in radians from -pi / 2 to +pi / 2
 
	The input altitude is assumed to be geometric, not apparent/refracted.
 
	The output right ascension and declination are for the current equinox,
	not corrected for precession to J2000 or nutation.
 
******************************************************************************/

void ASTROLIBDLL_API AAHorizonToEquatorial ( double azm, double alt,
	double lon, double lat, double jd, double *ra, double *dec );

/*** AAEquatorialToHorizon ***************************************************
 
	Converts equatorial coordinates to horizon.
 
	void AAEquatorialToHorizon ( double ra, double dec, double lon, double lat,
	double jd, double *azm, double *alt )
 
	(ra):  Right Ascension in radians from 0 to 2 * pi.
	(dec): Declination in radians from -pi / 2 to +pi / 2.
	(lon): Geographic longitude in radians, east is positive.
	(lat): Gepographc latitude in radians, north is positive.
	(jd):  Julian Date
	(azm): Receives azimuth in radians, measured from North = 0, East = pi / 2.
	(alt): Receives altitude in radians, from horizon = 0 to zenith = pi / 2.
 
	The input right ascension and declination are assumed to be for the current
	equinox, not J2000, and not corrected for nutation.
 
	The returned altitude is geometric, not apparent; refraction is NOT taken
	into account.
 
******************************************************************************/

void ASTROLIBDLL_API AAEquatorialToHorizon ( double ra, double dec,
	double lon, double lat, double jd, double *azm, double *alt );

/*** AASetEclipticRotationMatrix ********************************************

	Computes rotation matrix for transforming XYZ vectors between the
	equatorial and ecliptic frames.

	void AASetEclipticRotationMatrix ( double m[3][3], double e, int i )
  
	(e): Obliquity of the ecliptic, in radians.
	(m): Array of 9 doubles to recieve rows of the matrix.
	(i): Direction of the operation desired for the matrix:
	     if i > 0, the matrix will transform equatorial -> ecliptic;
	     if i < 0, the matrix will transform ecliptic -> equatorial.
               
	This function returns nothing.
  
***************************************************************************/

void ASTROLIBDLL_API AASetEclipticRotationMatrix ( double [3][3], double, int );

/*** AASetHorizonRotationMatrix ********************************************

	Computes rotation matrix for transforming XYZ vectors between the
	equatorial and local horizon frames.

	void AASetHorizonRotationMatrix ( double m[3][3], double l, double b, int i )
          
	(m): Array of 9 doubles to recieve rows of the matrix.     
	(l): The local sidereal time, in radians.
	(b): The local latitude, in radians.
	(i): Direction of the operation desired for the matrix:
	     if i > 0, the matrix will transform equatorial -> horizon;
         if i < 0, the matrix will transform horizon -> equatorial.

	Returns nothing.
      
	NOTE: Azimuth increases eastward from north, i.e. N = 0, E = 90,
	S = 180, W = 270, i.e. clockwise as viewed from above.  Therefore,
	a vector transformed by this matrix will be mirror imaged
	with respect to the sky, which we view from below. 
               
******************************************************************************/

void ASTROLIBDLL_API AASetHorizonRotationMatrix ( double [3][3], double, double, int );

/*** AASetGalacticRotationMatrix ************************************************

	Computes matrix for transforming between equatorial and galactic coordinates.
               
	void AASetGalacticRotationMatrix ( double m[3][3], double jd, int i );

	 (m): Array of 9 doubles to recieve rows of matrix.
	(jd): Julian date of precessional epoch to which equatorial coordinates
	      are referred.
	 (i): Direction of the operation desired for the matrix:
	      if i > 0, the matrix will transform equatorial -> galactic;
	      if i < 0, the matrix will transform galactic -> equatorial.

	This function returns nothing.
      
	NOTE: The galactic north pole is fixed by convention at B1950 coordinates
	R.A. 12h 49m, Dec +27 24'.  It does not move with precession, but lies fixed
	with the stars.  The matrix returned by this function will correctly precess
	an equatorial position vector to the B1950 frame before transforming it to
	the galactic frame, and vice-versa.
               
	References: Jean Meeus, "Astronomical Algorithms", pp. 89-91.
   
********************************************************************************/

void ASTROLIBDLL_API AASetGalacticRotationMatrix ( double [3][3], double, int );

/*** AASetPrecessionRotationMatrix ************************************************

	Computes a rotation matrix for precessing equatorial coordinates.
  
	void AASetPrecessionRotationMatrix ( double m[3][3], double jd0, double jd1, int n )

	  (m): array of 9 doubles to recieve rows of the matrix.
	(jd0): Julian date of the initial epoch.
	(jd1): Julian date of the final epoch.
	  (n): if non-zero, the matrix will also contain the effects of nutation.
               
	Returns nothing.  If nutation is included, the matrix returned will
	transform a vector from the the mean equatorial frame of epoch JD0
	to the true equatorial frame of epoch JD1.  If not, the matrix will
	only transform to the mean frame of epoch JD1.      

	References: Explanatory Supplement to the Astronomical Almanac, pp. 101-105.
                Jean Meeus, "Astronomical Algorithms", pp. 126-128.
  
	[Note: the signs of the rotations in the matrix given by the Supplement
	appear to be inconsistent.  The correct matrix is given below.]
                 
***********************************************************************************/

void ASTROLIBDLL_API AASetPrecessionRotationMatrix ( double [3][3], double, double, int );

/*** AASetEclipticPrecessionRotationMatrix ******************************************

	Computes a rotation matrix for precessing ecliptic coordinates.

	void AASetEclipticPrecessionRotationMatrix ( double m[3][3], double jd0, double jd1 )
  
	  (m): array of 9 doubles to recieve rows of the matrix.
	(jd0): Julian date of the initial epoch.
	(jd1): Julian date of the final epoch.
  
	This function returns nothing.         

	References: Explanatory Supplement to the Astronomical Almanac, pp. 99-101.
                Jean Meeus, "Astronomical Algorithms", pp. 128-129.
                 
*************************************************************************************/

void ASTROLIBDLL_API AASetEclipticPrecessionRotationMatrix ( double [3][3], double, double );

/*** AASetNutationRotationMatrix ******************************************************

	Calculates the matrix which applies nutation to transform between the mean
	and true equatorial frame of any epoch.
               
	void AASetNutationRotationMatrix ( double m[3][3], double e, double dl, double de, int i )

	 (m): Array of 9 doubles to recieve elements of the matrix.
	 (e): Mean obliquity of the ecliptic at epoch, in radians.
	(dl): Nutation in longitude, in radians.
	(de): Nutation in obliquity, in radians.
	 (i): Sign of the operation desired for the matrix:
	      if i > 0, the matrix will transform mean->true;
          if i < 0, the matrix will transform true->mean.
	
	This function returns nothing.
	  
	References: Explanatory Supplement to the Astronomical Almanac, p. 114.

**************************************************************************************/

void ASTROLIBDLL_API AASetNutationRotationMatrix ( double m[3][3], double e, double dl, double de, int i );

/**********************  Functions in Orbit.c **********************************/

#ifdef WIN32

/*** asinh ***************************************************************

	Computes the inverse hyperbolic sine, cosine, and tangent of its argument.

	double asinh ( double x )
	double acosh ( double x )
	double atanh ( double x )

	(x): argument, in radians.
          
	These function return the inverse hyperbolic sine, cosine, and tangent
	of their arguments.

	These functions are only compiled for Win32.  They seem to be part of
	the standard library on other platforms.

****************************************************************************/

double asinh ( double );
double acosh ( double );
double atanh ( double );

#endif

/*** AAMeanMotion *************************************************************

	Computes the mean motion of an object in a 2-body orbit, or the periapse
	distance from the mean motion.
  
	double AAMeanMotion ( double mu, double q, double e )
	double AAPeriapseDistance ( double mu, double n, double e )

	(mu): Reduced mass of orbiting bodies.
	 (q): Periapse distance of the orbit.
	 (e): Eccentriciy of the orbit.
     (n): Mean motion, in radians per unit time.
          
	Returns mean motion of the orbiting body; units depend on (k).
	(mu) is related to the masses of the orbiting bodies by
               
	mu = G * ( M + m ) = k * k * ( M + m )
               
	where G is the Newtonian constant of gravitation, and th masses of the
	bodies are M and m.  G is related to the Gaussian gravitational constant,
	k, by G = k * k.
               
	In the 1976 IAU heliocentric system of units, the unit of length is the AU,
	the unit of mass is the solar mass, the unit of time is the SI day, and k
	takes the value 0.01720209895.
               
	In the geocentric system, the unit of length is the Earth's equatorial
	radius, the unit of mass is the earth's mass, and the unit of time is the
	SI minute, and k takes the value 0.07436680.  
               
*******************************************************************************/

double ASTROLIBDLL_API AAMeanMotion ( double, double, double );
double ASTROLIBDLL_API AAPeriapseDistance ( double, double, double );

/*** J2MeanMotion *************************************************************

	Computes J2-perturbed mean motion and secular variations in angular orbital
	elements due to J2.
  
	void J2MeanMotion ( double ae, double j2, double a, double e, double i,
	double dm, double *dw, double *dn, double *dmm )

	 (ae): Radius of primary, in arbitrary units.
	 (j2): J2 harmonic coefficient, dimensionless.
	  (a): Semimajor axis of the orbit, in same units as (ae).
	  (e): Eccentricity of the orbit.
	  (i): Inclination of the orbit plane to the equatorial plane of primary.
	 (dm): Unperturbed mean motion.
	 (dw): Recieves secular variation in argument of periapse.
	 (dn): Recieves secular variation in longitude of ascending node.
	(dmm): Recieves perturbation to mean motion due to J2.
        
	Returns nothing.  J2 does not cause any secular variations in the elements
	a, e, and i.
               
	References:
	
	Kelso, T.S. "SPACETRACK REPORT No. 3: Models for Propogation of
	Norad Element Sets", 31 December 1988, pp. 3-4. National Technical
	Information Service, Springfield, VA. 
      
	Danby, J.M.A., "Fundamentals of Celestial Mechanics", p. 347.
	Willmann-Bell, Richmond VA, 1992.
               
**********************************************************************************/

void ASTROLIBDLL_API J2MeanMotion ( double, double, double, double, double, double,
double *, double *, double * );

/*** AASolveKeplersEqn ***********************************************************

	Solves Kepler's equation for elliptical, parabolic, and hyperbolic orbits.
               
	void AASolveKeplersEqn ( double m, double e, double q, double *nu, double *r )

	 (m): Mean anomaly, in radians.
	 (e): Eccentricity.
	 (q): Periapse distance.
	(nu): Receives true anomaly, in radians.
	 (r): Receives true distance, in same units as (q).
           
	Returns nothing.  If convergence is not achieved with the maximum number of
	iterations (defined below), erroneous results may be returned.
	
	For elliptical orbits, true anomaly (nu) is always returned in the range 0 to
	TWO_PI radians.  For parabolic and hyperbolic orbits, nu may have any value.
	
	Note: The method used to solve Kepler's equation depends on the eccentricity
	of the orbit.  For closed (elliptical) orbits, Newton's method usually provides
	the fastest solution.  However, if the orbit is highly eccentric (e > 0.975),
	Newton's method may not converge, so we substitute a slower but foolproof binary
	search algorithm.  For parabolic orbits (e = 1.0), we use the the standard
	Newton method.

************************************************************************************/

void ASTROLIBDLL_API AASolveKeplersEqn ( double, double, double, double *, double * );

/*** AAInverseKeplersEqn ***********************************************************

	Solves the inverse of Kepler's equation for elliptical, parabolic, and hyperbolic
	orbits.
               
	void AAInverseKeplersEqn ( double nu, double e, double r, double *m, double *q )

	(nu): True anomaly, in radians.
	 (e): Eccentricity.
	 (r): True distance.
	 (m): Receives mean anomaly, in radians.
	 (q): Receives periapse distance, in same units as (r).
           
	Returns nothing.

	For elliptical orbits, mean anomaly (m) is always returned in the range 0 to
	TWO_PI radians.  For parabolic and hyperbolic orbits, m may have any value.

************************************************************************************/

void ASTROLIBDLL_API AAInverseKeplersEqn ( double, double, double, double *, double * );

/*** AAOrbitToSpherical ****************************************************

	Computes spherical coordinates for an object in a 2-body orbit
	from classical Keplerian elements.
               
	void AAOrbitToSpherical ( double q, double e, double w, double i, double n,
	double m, double *l, double *b, double *r )

	(q): Periapse distance.
	(e): Eccentricity of the orbit.
	(w): Argument of the perihelion, in radians.
	(i): Inclination of the orbit, in radians.
	(n): Longitude of the ascending node, in radians.
	(m): Mean anomaly, in radians.
	(l): Receives the object's longitude coordinate, in radians.
	(b): Receives object's latitude coordinate, in radians.
	(r): Receives object's radial coordinate, in same units as (q).
          
********************************************************************************/

void ASTROLIBDLL_API AAOrbitToSpherical ( double, double, double, double, double, double,
double *, double *, double * );

/*** AAOrbitToXYZ ************************************************************

	Computes rectangular coordinates for an object in a 2-body orbit from
	classical elements.
               
	void OrbitToXYZ ( double q, double e, double i, double w, double n, double m,
	     double *x, double *y, double *z )

	void OrbitToXYZVector ( double q, double e, double i, double w, double n, double m,
	     double v[3] )

	(q): Periapse distance.
	(e): Eccentricity of the orbit.
	(i): Inclination of the orbit, in radians.
	(w): Argument of the perihelion, in radians.
	(n): Longitude of the ascending node, in radians.
	(m): Mean anomaly, in radians.
	(x): Receives the object's X-coordinate, in same units as (q).
	(y): Receives object's Y-coordinate, in same units as (q).
	(z): Receives object's Z-coordinate, in same units as (q).
	(z): Receives object's XYZ coordinates, in same units as (q).
          
	Returns nothing.
                 
********************************************************************************/

void ASTROLIBDLL_API AAOrbitToXYZ ( double, double, double, double, double, double,
	double *, double *, double * );

void ASTROLIBDLL_API AAOrbitToXYZVector ( double, double, double, double, double, double,
	double[3] );

/*** AAOrbitToXYZMotion ******************************************************

	Computes rectangular coordinates and motion for an object in a 2-body
	orbit from classical elements.
               
	void AAOrbitToXYZMotion ( double q, double e, double i, double w, double n,
	double m, double dm, double *x, double *y, double *z, double *dx,
	double *dy, double *dz, double *mu )

	 (q): Periapse distance.
	 (e): Eccentricity of the orbit.
	 (i): Inclination of the orbit, in radians.
	 (w): Argument of the perihelion, in radians.
	 (n): Longitude of the ascending node, in radians.
	 (m): Mean anomaly, in radians.
	(dm): Mean motion, in radians per unit time.
	 (x): Receives X coordinate.
	 (y): Receives Y coordinate.
	 (z): Receives Z coordinate.
	(dx): Receives X component of velocity.
	(dy): Receives Y component of velocity.
	(dz): Receives Z component of velocity.
	(mu): Receives combined mass of the bodies multiplied by the
	      Newtonian constant of gravitation: mu = G * ( M + m )
          
********************************************************************************/

void ASTROLIBDLL_API AAOrbitToXYZMotion ( double, double, double, double, double, double, double,
double *, double *, double *, double *, double *, double *, double * );

void ASTROLIBDLL_API AAOrbitToXYZMotionVector ( double, double, double, double, double, double, double,
double [3], double [3], double * );

/*** AAXYZMotionToOrbit ***********************************************************

	Computes classical 2-body orbital elements from rectangular position
	and velocity vectors.

	void AAXYZMotionToOrbit ( double, double, double, double, double, double, double,
	double *x, double *y, double *z, double *dx, double *dy, double *dz, double *mu )

	 (x): X component of position vector.
	 (y): Y component of position vector.
	 (z): Z component of position vector.
	(dx): X component of velocity vector.
	(dy): Y component of velocity vector.
	(dz): Z component of velocity vector.
	(mu): Combined masses of the two bodies, m and M, multiplied by the
	      Newtonian constant of gravitation, G; i.e. mu = G * ( M + m ).              
	 (q): Receives periapse distance, in same units as (X,Y,Z).
	 (e): Receives orbital eccentricity.
	 (i): Receives orbital inclination, in radians.
	 (w): Receives argument of periapse, in radians.
	 (n): Receives longitude of ascending node, in radians.
	 (m): Receives mean anomaly, in radians.
	(dm): Receives mean motion, in radians per unit time.
      
	Returns nothing.
      
	For combinations of position and velocity that imply an eccentricity
	within 1.0E-8 of 1.0, the orbit is assumed to be parabolic.
               
	Note that (w) is undetermined when (e) = 0.0 and (n) is undetermined
	when (i) = 0.0.
            
*******************************************************************************/

void ASTROLIBDLL_API AAXYZMotionToOrbit ( double, double, double, double, double, double, double,
double *, double *, double *, double *, double *, double *, double * );

void ASTROLIBDLL_API AAXYZMotionVectorToOrbit ( double xyz[3], double dxyz[3], double mu,
double *q, double *e, double *i, double *w, double *n, double *m, double *dm );

/*** AATransformOrbit **********************************************************

	Transforms angular orbital elements from one coordinate frame to another.
  
	void AATransformOrbit ( double m[3][3], double *i, double *w, double *n )

	(m): rotation matrix describing transformation from initial
	     to final frame.
	(i): inclination angle, in radians.
	(w): argument of periapse, in radians.
	(n): longitude of ascending node, in radians.
          
	Returns nothing.  On return, the values contained in (i,w,n) will be replaced
	by their values in the new frame.
      
*********************************************************************************/

void ASTROLIBDLL_API AATransformOrbit ( double [3][3], double *, double *, double * );

/***  AASetOrbitRotationMatrix  *************************************************

	Sets matrix for transforming rectangular coordinates from orbit plane to
	reference plane.
  
	void AASetOrbitRotationMatrix ( double m[3][3], double i, double w,
	     double n, int s )

	(m): receives rotation matrix.
	(i): inclination angle, in radians.
	(w): argument of periapse, in radians.
	(n): longitude of ascending node, in radians.
	(e): obliquity of ecliptic, in radians.
    (s): sign of transformation: if > 0, resulting matring transforms
         from orbit plane to reference plane; otherwise from reference plane
         to orbit plane.
    
	Returns nothing.  If s is positive, the matrix returned in (m) will transform
	from orbit plane coordinates (where the x,y plane is the orbit plane, and the
	periapse lies along the positive x axis) to reference plane coordinates (e.g.
	the J2000 mean ecliptic).  If s is negative, the returned matrix will transform
	in the opposite direction (i.e. from reference plane coordinates to orbit
	plane coordinates.)
	
	If (e) is non-zero, an additional "twist" will be included in the resulting
	rotation matrix, which transforms from the ecliptic to the equatorial reference
	frame (if s is positive) or vice-versa (if s is negative).  This is useful
	if the orbital elements (i,w,n) are referred to the ecliptic (as is typical),
	but you want orbit positions computed relative to the equatorial frame.
	
*********************************************************************************/

void ASTROLIBDLL_API AASetOrbitRotationMatrix ( double [3][3], double, double, double, double, int );

/*** AABinaryStarOrbitPosition ***************************************************

	Computes separation and position angle of binary star system components from
	their orbital elements.
 
	void AABinaryStarOrbitPosition ( double m, double a, double e, double i, double w, double n,
		 double *r, double *rp, double *pa )

	(m): current mean anomaly of binary star component, in radians.
	(a): semimajor axis of true binary orbit, in arcseconds or radians.
	(e): eccentricity of true binary orbit, dimensionless.
	(i): inclination of binary orbit to the plane of the sky, in radians.
	(w): argument of periastron, in radians.
	(n): position angle of the ascending node of the binary orbit on the sky plane, in radians.
	(r): recieves true separation of binary star system components.
   (rp): recieves apparent separation of binary star system components.
   (pa): recieves position angle of binary star system components, in radians.
 
	The plane of the sky is the plane perpendicular to the line-of-sight to the
	binary system.  Note that the binary system's true orbit differs from its
	apparent orbit because it is inclined to the plane of the sky.  Similarly,
	the true and apparent separation returned by this function differ for the same
	reason: the apparent separation is the true separation projected onto the
	plane of the sky.  For example, if one component appears directly in front
	of the other, their apparent separation may be zero, but their true separation
	quite large.  The true and apparent separations will be returned in the same
	units as the semimajor axis; e.g. if the semimajor axis was provided in arc-
	seconds, then the separations will be returned in arcseconds.
 
	REFERENCES:
 
	Jean Meeus, "Astronomical Algorithms", pp. 397-400.
 
**********************************************************************************/

void ASTROLIBDLL_API AABinaryStarOrbitPosition ( double m, double a, double e, double i, double w, double n,
double *r, double *rp, double *pa );

/*** VFPlanet ***************************************************************

	Calculates the ecliptic longitude, latitude, and radius vector for the
    Sun, Moon, and major planets on a given Julian date.
               
	void VFPSun ( double jd, double *l, double *b, double *r )
	void VFPMoon ( double jd, double *l, double *b, double *r )

	void VFPMercury ( double jd, double *l, double *b, double *r )
	void VFPVenus ( double jd, double *l, double *b, double *r )
	void VFPEarth ( double jd, double *l, double *b, double *r )
	void VFPMars ( double jd, double *l, double *b, double *r )
	void VFPJupiter ( double jd, double *l, double *b, double *r )
	void VFPSaturn ( double jd, double *l, double *b, double *r )
	void VFPUranus ( double jd, double *l, double *b, double *r )
	void VFPNeptune ( double jd, double *l, double *b, double *r )
	void VFPPluto ( double jd, double *l, double *b, double *r )

	(jd): Julian date.
	 (l): Receives ecliptic longitude, in radians.
	 (b): Receives ecliptic latitude, in radians.
	 (r): Receives radius vector, in AU.
         
	Returns nothing.  (l) and (b) are referred to the mean ecliptic and equinox
	of date, and should be correct to 1' for any date within 300 years of the
    present.

    For the sun and moon, the longitude and latitude returned are geocentric.
    For the major planets Mercury - Pluto, the longitude and latitude returned
    are heliocentric.

	Note that VFPEarth simply computes the sun's geocentric position and
	reverses its coordinates by negating the latitude and adding 180 degrees
    to the longitude.

	Reference:
	
	T. C. Van Flandern and F. K. Pulkkinnen, "Low Precision Formulae for Planetary
	Positions", Astrophysical Journal Supplement Series 41: 391-411, 1979 November.

************************************************************************************/
    
void ASTROLIBDLL_API VFPSun ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPMoon ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPMercury ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPVenus ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPEarth ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPMars ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPJupiter ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPSaturn ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPUranus ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPNeptune ( double, double *, double *, double * );
void ASTROLIBDLL_API VFPPluto ( double, double *, double *, double * );

/*** VSOP87 *****************************************************************

	Calculates the heliocentric ecliptic longitude, latitude, and radius vector
    for the major planets Mercury - Neptune on a given Julian date.
               
	void VSOP87Mercury ( double jd, double *l, double *b, double *r )
	void VSOP87Venus ( double jd, double *l, double *b, double *r )
	void VSOP87Earth ( double jd, double *l, double *b, double *r )
	void VSOP87Mars ( double jd, double *l, double *b, double *r )
	void VSOP87Jupiter ( double jd, double *l, double *b, double *r )
	void VSOP87Saturn ( double jd, double *l, double *b, double *r )
	void VSOP87Uranus ( double jd, double *l, double *b, double *r )
	void VSOP87Neptune ( double jd, double *l, double *b, double *r )

	(jd): Julian date.
	 (l): Receives heliocentric ecliptic longitude, in radians.
	 (b): Receives heliocentric ecliptic latitude, in radians.
	 (r): Receives radius vector, in AU.
         
	Returns nothing.  (l) and (b) are referred to the mean ecliptic and equinox
	of date.

    For Mercury, Venus, Earth-Moon barycenter, and Mars, positions should be
    computed with a precision of 1" for 4000 years before and after J2000.
	The same precision is ensured for Jupiter and Saturn over 2000 years and
    for Uranus and Neptune over 6000 years before and after J2000.

	References:
	
	Bretagnon P., Francou G., : 1988, Astron. Astrophys., 202, 309.
	Meeus, J., "Astronomical Algorithms", pp. 205-208.

************************************************************************************/

void ASTROLIBDLL_API VSOP87Mercury ( double, double *, double *, double * );
void ASTROLIBDLL_API VSOP87Venus ( double, double *, double *, double * );
void ASTROLIBDLL_API VSOP87Earth ( double, double *, double *, double * );
void ASTROLIBDLL_API VSOP87Mars ( double, double *, double *, double * );
void ASTROLIBDLL_API VSOP87Jupiter ( double, double *, double *, double * );
void ASTROLIBDLL_API VSOP87Saturn ( double, double *, double *, double * );
void ASTROLIBDLL_API VSOP87Uranus ( double, double *, double *, double * );
void ASTROLIBDLL_API VSOP87Neptune ( double, double *, double *, double * );

/**** ELP2000Moon ************************************************************

	Calculates the geocentric ecliptic longitude, latitude, and radius vector
    of the Moon on a given Julian date.
               
	void ELP2000Moon ( double jd, double *l, double *b, double *r )

	(jd): Julian date.
	 (l): Receives geocentric ecliptic longitude, in radians.
	 (b): Receives geocentric ecliptic latitude, in radians.
	 (r): Receives radius vector, in AU.
         
	Returns nothing.  (l) and (b) are referred to the mean ecliptic and
	equinox of date.

	References:
	
	Meeus, J. "Astronomical Algorithms". pp. 307-314.

	Chapront-Touze M., Chapront J. "ELP 2000-85: a semi-analytical lunar
	ephemeris adequate for historical times."  Astronomy & Astrophysics,
	vol. 190, p. 342 (1988).

	Chapront-Touze M., Chapront J. "The Lunar Ephemeris ELP 2000".
    Astronomy & Astrophysics, vol. 124, pp. 50-62 (1983).

************************************************************************************/

void ASTROLIBDLL_API ELP2000Moon ( double, double *, double *, double * );

/*** ELP2000SphericalToJ2000XYZ ****************************************************

	Converts ELP2000 spherical coordinates referred to the ecliptic of date
	to rectangular coordinates referred to the J2000 ecliptic.
               
	void ELP2000SphericalToJ2000XYZ ( double jd, double l, double b, double r,
	double *x, double *y, double *z )

	(jd): Julian date.
	 (l): Geocentric ecliptic longitude, in radians.
	 (b): Geocentric ecliptic latitude, in radians.
	 (r): Radius vector, in AU.
     (x): Recieves J2000 ecliptic X coordinate, in A.U.
     (y): Recieves J2000 ecliptic Y coordinate, in A.U.
     (z): Recieves J2000 ecliptic Z coordinate, in A.U.
         
	Returns nothing.  The function implements Laskar's precession matrix as
	presented in the ELP 2000 papers described below, specifically for converting
	ELP2000 lunar coordinates from the ecltipic and equinox of date to the 
	ecltipic of J2000 in order to compare results against published values.

	References:
	
	Meeus, J. "Astronomical Algorithms". pp. 307-314.

	Chapront-Touze M., Chapront J. "ELP 2000-85: a semi-analytical lunar
	ephemeris adequate for historical times."  Astronomy & Astrophysics,
	vol. 190, p. 342 (1988).

	Chapront-Touze M., Chapront J. "The Lunar Ephemeris ELP 2000".
    Astronomy & Astrophysics, vol. 124, pp. 50-62 (1983).

************************************************************************************/

void ASTROLIBDLL_API ELP2000SphericalToJ2000XYZ ( double, double, double, double t,
double *, double *, double * );

/*** PLUTO95Pluto ******************************************************************

	Calculates Pluto's heliocentric rectangular coordinates on a given date.
               
	void PLUTO95Pluto ( double jd, double *x, double *y, double *z )

	(jd):    Julian date.
	(x,y,z): Receives Pluto's J2000 rectangular equatorial coordinates
	         on the given date.

	Returns nothing.  The rectangular coordinates (x,y,z) are returned in
	Astronomical Units, referred to the Earth's equator and equinox of J2000.
	
	The expressions used to compute Pluto's position are an analytic fit
	to the JPL numerical DE200/LE200 ephemerides.  The interval of validity
	of the fit is JD 2341972.5 - 2488092.5 (1 Jan 1700 - 24 Jan 2100).
	
	References:
	
	J. Chapront, G. Francou, "Representation of planetary ephemerides by
	frequency analysis."  Astronomy & Astrophysics Supplemental Series,
	109, 191 (1995).

************************************************************************************/

void ASTROLIBDLL_API PLUTO95Pluto ( double, double *, double *, double * );

/*** AAPluto *********************************************************************

  	Calculates Pluto's heliocentric rectangular coordinates on a given date.
               
	void AAPluto ( double jd, double *lon, double *lat, double *rad )

	(jd):  Julian date.
	(lon): Receives Pluto's J2000 heliocentric ecliptic longitude, in radians
	(lat): Receives Pluto's J2000 heliocentric ecliptic latitude, in radians
	(rad): Receives Pluto's J2000 heliocentric distance, in AU.

	Returns nothing.
	
	The expressions used to compute Pluto's position are an analytic fit
	to the JPL numerical DE405 ephemeris.  The interval of validity
	of the fit is 1885 - 2099.
	
	References:
	
	Jean Meeus, "Astronomical Algorithms", pp. 263 - 267.

**********************************************************************************/

void ASTROLIBDLL_API AAPluto ( double, double *, double *, double * );

/*** AAMercury *******************************************************************

	Calculates the heliocentric ecliptic longitude, latitude, and radius vector
    for the major planets Mercury - Neptune on a given Julian date.
               
	void AAMercury ( double jd, double *l, double *b, double *r )
	void AAVenus ( double jd, double *l, double *b, double *r )
	void AAEarth ( double jd, double *l, double *b, double *r )
	void AAMars ( double jd, double *l, double *b, double *r )
	void AAJupiter ( double jd, double *l, double *b, double *r )
	void AASaturn ( double jd, double *l, double *b, double *r )
	void AAUranus ( double jd, double *l, double *b, double *r )
	void AANeptune ( double jd, double *l, double *b, double *r )

	(jd): Julian date.
	 (l): Receives heliocentric ecliptic longitude, in radians.
	 (b): Receives heliocentric ecliptic latitude, in radians.
	 (r): Receives radius vector, in AU.
         
	Returns nothing.  (l) and (b) are referred to the mean ecliptic and equinox
	of date.

	For Julian Dates more than 10,000 years from J2000, the higher-order
	polynomial terms are truncted at T = +/- 10,000 years from the present,
	with the exception of the linear terms for longitude.  This provides
	results which are realistic (if not scientifically accurate) over an
	unlimited time interval from J2000.
	
	References:
	
	Bretagnon P., Francou G., : 1988, Astron. Astrophys., 202, 309.
	Meeus, J., "Astronomical Algorithms", pp. 205-208.

**********************************************************************************/

void ASTROLIBDLL_API AAMercury ( double, double *, double *, double * );
void ASTROLIBDLL_API AAVenus ( double, double *, double *, double * );
void ASTROLIBDLL_API AAEarth ( double, double *, double *, double * );
void ASTROLIBDLL_API AAMars ( double, double *, double *, double * );
void ASTROLIBDLL_API AAJupiter ( double, double *, double *, double * );
void ASTROLIBDLL_API AASaturn ( double, double *, double *, double * );
void ASTROLIBDLL_API AAUranus ( double, double *, double *, double * );
void ASTROLIBDLL_API AANeptune ( double, double *, double *, double * );

/*** AASunEclipticCoords/AAMoonEclipticCoords *************************************
	 
	Computes the Sun and Moon's approximate geocentric longitude, latitude,
	and distance for a given Julian Date, referred to the ecliptic of that date.
	 
	void AASunEclipticCoords ( double jd, double *lon, double *lat, double *rad )
	void AAMoonEclipticCoords ( double jd, double *lat, double *lat, double *rad )

	(jd): Julian Date at which to compute the Sun or Moon's coordinates.
	(lon): Receives ecliptic longitude, in radians.
	(lat): Receives ecliptic latitude, in radians (zero for Sun).
	(rad): Receives geocentric distance, in AU (Sun) or Earth-radii (Moon).
	 
	Returns nothing. To avoid computing lon, lat, or rad, pass NULL for any of
    those parameters.

	These are fast, approximate formulae.  They compute the Sun's ecliptic longitude
	to a precision of 0.01 degrees between 1950 and 2050.  For the Moon, errors will
	rarely exceed 0.3 degrees in longitude, 0.2 degrees in latitude, abd 0.2 Earth-
	radii in distance.
 
	Reference:
	 
	The Astronomical Almanac for the Year 1984, page C24 and D46.
	 
***********************************************************************************/
	
void ASTROLIBDLL_API AAMoonEclipticCoords ( double jd, double *lon, double *lat, double *rad );
void ASTROLIBDLL_API AASunEclipticCoords ( double jd, double *lon, double *lat, double *rad );

/*** AAEclipticCoordsToEquatorialXYZVector *************************************
	 
	Converts ecliptic spherical coordinates (longitude, latitude, and distance)
    to an equatorial rectangular XYZ vector for a particular Julian Date.
	 
	void AAEclipticCoordsToEqatorialXYZVector ( double jd,
 	     double lon, double lat, double rad, double xyz[3] );
 
	(jd): Julian Date of coordinate system.
	(lon): Ecliptic longitude, in radians.
	(lat): Ecliptic latitude, in radians.
	(rad): Distance from origin of XYZ coordinate system.
	 
	Returns nothing.  The Julian date (jd) determines the obliquity of the ecliptic;
 	if high precision is not needed, passing J2000 will give results accurate
	to 0.001 degrees for the years 1950 to 2050.
 
	Reference:
	 
	The Astronomical Almanac for the Year 1984, page C24 and D46.
	 
***********************************************************************************/

void ASTROLIBDLL_API AAEclipticCoordsToEquatorialXYZVector ( double jd, double lon, double lat, double rad, double xyz[3] );

/*** AAMoon ***********************************************************************

	Computes the Moon's geocentric longitude, latitude, and distance for a given
	Julian Date, referred to the ecliptic and equinox of that date.
	
	void AAMoon ( double jd, double *l, double *b, double *r )
	
	(jd): Julian Date at which to compute the Moon's position.
	(l):  Receives Moon's ecliptic longitude, in radians.
	(b):  Receives Moon's ecliptic latitude, in radians.
	(r):  Receives Moon's geocentric distance, in kilometers.

	Returns nothing.

	The algorithm for computing the Moon's position is based on a simplification
	of the complete ELP2000-82 theory, and should provide results accurate to
	approximately 10" in ecliptic longitude, and 4" in latitude.  However, it
	does not contain any high-order polynomial terms in time, and thus provides
	realistic (if not scientifically accurate) results for an unlimited time
	interval.
	
	References:
	
	Jean Meeus, "Astronomical Algorithms", pp. 337 - 344.
	
***********************************************************************************/

void ASTROLIBDLL_API AAMoon ( double, double *, double *, double * );

/*** AANextMoonPhase ******************************************************************

	Computes the Julian Date of the next lunar phase event after a particular JD.
 
	void AANextMoonPhase ( double jd, int event )
 
	(jd): Julian Date from which to compute the next lunar phase event.
	(event): Event identifier: AA_NEW_MOON thru AA_LAST_QUARTER.

	Returns the Julian Ephemeris Date of the next event following the specified
	Julian Date.

	References:
 
	Jean Meeus, "Astronomical Algorithms", 2nd ed., pp. 349 - 354.
 
***********************************************************************************/

#define AA_NEW_MOON			0
#define AA_FIRST_QUARTER	1
#define AA_FULL_MOON		2
#define AA_LAST_QUARTER		3
	
double ASTROLIBDLL_API AANextMoonPhase ( double, int );

/*** AAMoonOrbit ***************************************************************
 
	Calculates the mean orbital elements of the Moon for a given Julian Date,
 	referred to the ecliptic and equinox of that date.
 
	void AAMoonOrbit ( double jd, double *a, double *e, double *i, double *w,
	double *n, double *m, double *dm )

	(jd): Julian date for which mean orbital elements are to be computed.
	(a):  Receives semimajor axis in AU.
	(e):  Receives eccentricity.
	(i):  Receives inclination to J2000 ecliptic, in radians.
	(p):  Receives longitude of perihelion, in radians.
	(n):  Receives longitude of ascending node, in radians.
	(m):  Receives mean anomaly at specified epoch, in radians.
	(dm): Receives mean motion, in radians per day.
 
	Returns nothing.
 
	The expressions used to compute the mean elements of the Moon's orbit
	include polynomial representations of their long-term secular variations,
	but do not include any short-term perturbations.  The polynomials are
	truncated at time t = +/- 10,000 years from J2000, except for the first-order
	terms for (w,n,m), which are allowed to increase linarly with time beyond
	that interval.  This should provide results which are realistic (if not
	scientifically accurate) over an unlimited timespan.
 
	References:
 
	Jean Meeus, "Astronomical Algorithms", pp. 338 - 344.
 	Explanatory Supplement to the Astronomical Alamanac, p. 701.
 
***********************************************************************************/

void ASTROLIBDLL_API AAMoonOrbit ( double jd, double *a, double *e, double *i, double *w,
double *n, double *m, double *dm );

/*** AAPlanetOrbit ***************************************************************
 
	Calculates the approximate orbital elements of the planets Mercury - Pluto
	for a given Julian Ephemeris Date, referred to the J2000 ecliptic.
 
	void AAPlanetOrbit ( double jd, double *a, double *e double *i, double *w,
	double *n, double *m, double *dm );
 
	(jd): Julian ephemeris date for which mean orbital elements are to be computed.
	(a):  Receives semimajor axis in AU.
	(e):  Receives eccentricity.
	(i):  Receives inclination to J2000 ecliptic, in radians.
	(w):  Receives argument of perihelion, in radians.
	(n):  Receives longitude of ascending node, in radians.
	(m):  Receives mean anomly, in radians.
	(dm): Receives mean motion, in radians per day.
 
	Returns nothing.
 
	To convert semimajor axis (a) to perihelion distance (q): q = a * ( 1 - e ).
	To convert longitude of perihelion (l) to argument of perihelion (w): w = p - n.
	To convert mean longitude (l) to mean anomaly (m): m = l - p.
 
	The expressions used to compute the approximate elements of the planetary orbits
	are the result of a best fit, and are not valid outside the time interval -3000
	to 3000.  For the interval 1850 - 2000, positions predicted using these orbits
	are accurate to about 1 arcminute; outside that interval, accuracy is about 10
	arcminutes.
 
	Reference:
 
	E. M. Standish, "Keplerian Elements for Approximate Positions of the Major Planets",
	Solar System Dynamics Group, JPL/Caltech.
	
***********************************************************************************/

void ASTROLIBDLL_API AAMercuryOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AAVenusOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AAEarthOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AAMarsOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AAJupiterOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AASaturnOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AAUranusOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AANeptuneOrbit ( double, double *, double *, double *, double *, double *, double *, double * );
void ASTROLIBDLL_API AAPlutoOrbit ( double, double *, double *, double *, double *, double *, double *, double * );

/* Physical.c */

/*** AAPlanetographicCoordinates **************************************************

  DESCRIPTION: Determines the planetographic coordinatates of the center
               of a planet's disk and position angle of its north pole.
           
          (a): The planet's apparent right ascension, in radians.
          (d): The planet's apparent declination, in radians.
         (a1): Right ascension of the planet's north pole, in radians.
         (d1): Declination of the planet's north pole, in radians.
         (w1): Argument of the planet's prime meridian, in radians:
               positive for planets where longitude increases clockwise,
               negative for planets where longitude decreases clockwise.
               See note below.
          (l): Receieves longitude of the central meridian of the
               planet's disk, in radians.
          (b): Receievs planetocentric declination of the central
               point on the planet's disk, in radians.
          (p): Receives position angle of the north pole of the planet's
               disk, in radians.
               
      RETURNS: Nothing.  One should give (w1) for the time of observation
               antedated by the light time from the planet.
               
         NOTE: Currently, the IAU has embraced several different, conflicting
               conventions for measuring planetary longitudes.  Under the
               current definitions,

               1) A planet's north pole is the pole which lies north of the
                  invariable plane of the solar system (roughly the ecliptic).

               2) For the Earth, Sun, and Moon, longitude increases in the direction
                  of planetary rotation.  For all other planets, longitude increases
                  in the direction -opposite- planetary rotation.

               These definitions lead to several problems.  First, because of
               definition (1), some planets (Venus, Uranus, Pluto) rotate
               clockwise when viewed from above their north poles, while others
               rotate anticlockwise.  Worse, depending on whether the planet itself
               rotates clockwise or anticlockwise, convention (2) can make longitudes
               increase either clockwise -or- anticlockwise when viewed from above.
               For Uranus, it's even more confusing because the north pole lies
               south of the Earth's equator, which reverses things once again.
               
               In summary, for the Earth, Sun, Moon, Venus, and Pluto, longitude
               increases anticlockwise when viewed from above the north pole, so
               for these bodies the value given for (w1) should be negative.
               For the other planets, it should be positive.
               
               Also, note that (b) differs slightly from the latitude of the
               central point on the planet's disk, due to the planet's oblateness.
               
   REFERENCES: The Astronomical Almanac for the Year 1990.
               The Explanatory Supplement to the Astronomical Almanac, p. 387.
               
************************************************************************************/

void ASTROLIBDLL_API AAPlanetographicCoordinates(double, double, double, double, double, double *, double *, double *);

/*** AASetPlanetographicRotationMatrix **********************************************

	Computes a rotation matrix for transforming from J2000 equatorial coordinates
	to planetographic coordiantes, or vice-versa.
	
	void AASetPlanetographicRotationMatrix ( double m[3][3], double a1, double d1,
	     double w1, int i )
	     
	(m):  receives pointer to rotation matrix.
	(a1): Right ascension of the planet's north pole, in radians.
	(d1): Declination of the planet's north pole, in radians.
	(w1): Argument of the planet's prime meridian, in radians.
	(i):  Direction of the operation desired for the matrix:
	      if i > 0, the matrix will transform equatorial -> planetographic;
	      if i < 0, the matrix will transform planetographic -> equatorial.

	RETURNS: Nothing.  
	
	NOTE: The sign of the argument of the prime meridian (w1) is ignored by this
	      function (i.e. the function takes the abolute value of w1 before computing
	      the matrix.)
	
************************************************************************************/

void ASTROLIBDLL_API AASetPlanetographicRotationMatrix ( double [3][3], double, double, double, int );

/*** AAVectorPlanetographicCoordinates  ***********************************************

  DESCRIPTION: Determines the planetographic coordinatates of the center
               of a planet's disk and position angle of its north pole.
           
          (v): The planet's apparent direction, expressed as a unit vector
               in the J2000 equatorial frame.
         (a1): Right ascension of the planet's north pole, in radians.
         (d1): Declination of the planet's north pole, in radians.
         (w1): Argument of the planet's prime meridian, in radians;
               positive for planets where longitude increases clockwise,
               negative for planets where longitude decreases clockwise.
               See note above.
          (l): Receieves longitude of the central meridian of the
               planet's disk, in radians.
          (b): Receievs planetocentric declination of the central
               point on the planet's disk, in radians.
          (p): Receives position angle of the north pole of the planet's
               disk, in radians.
               
      RETURNS: Nothing.  One should give (w1) for the time of observation
               antedated by the light time from the planet.
               
         NOTE: See PlanetographicCoordinates() for more information on
               sign conventions regarding planetary longitudes.  This
               function should give identical results, but significant
               computation can be saved by using a vector-matrix method
               rather than the spherical-trigonometric method given in the
               Astronomical Almanace.  The method is one of my own devising
               and is detailed in the comments below.

************************************************************************************/

void AAVectorPlanetographicCoordinates ( double[3], double, double, double, double *, double *, double * );

/*** AAXXXRotation *************************************************************

  DESCRIPTION: Computes the rotational elements of the Sun, Moon, and major
               planets at a given time.
 
               void AASunRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AAMercuryRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AAVenusRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AAEarthRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AAMarsRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AAJupiterRotation ( double jd, double *a1, double *d1, double *w1, double *wd, char sys );
               void AASaturnRotation ( double jd, double *a1, double *d1, double *w1, double *wd, char sys );
               void AAUranusRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AANeptuneRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AAPlutoRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
               void AAMoonRotation ( double jd, double *a1, double *d1, double *w1, double *wd );
 
         (jd): Desired time, expressed as a Julian Date.
         (a1): Receives J2000 right ascension of object's north pole, in radians.
         (d1): Receives J2000 declination of object's north pole, in radians.
         (w1): Receives the argument of object's prime meridian, in radians
               from 0 to TWO_PI.  Negative indicates anticlockwise rotation.
		 (wd): Receives the planet's rotation rate, in radians per day.
		       Negative indicates anticlockwise rotation.
        (sys): For Jupiter, the longitude system that (w1) should be computed for:
               1 = System I, the mean equatorial atmospheric rotation
               2 = System II, the mean atmospheric rotation north of the
                   south component of the north equatorial belt and south
                   of the north component of the south equatorial belt.
               3 = System III, the magnetic field rotation, which corresponds
                   to the rotation of the deep interior.
               All other values will produce meaningless results.
               For Saturn, the longitude system that (w1) should be computed for:
               1 = System I, the mean atmospheric rotation.  This system has
                   been found to be of little use since the Voyager encounters,
                   and is no longer supported by the IAU.
               3 = System III, the magnetic field rotation, which corresponds
                   to the rotation of the deep interior.
               All other values will produce meaningless results.
          
               For the Sun, Venus, Earth, Moon, and Pluto, the value returned in (w1)
               will be negative to indicate that longitude (by convention) increases
               anticlockwise when viewed from above the planet's north pole.  For
               all other planets, the value returned in (w1) will be positive, to
               indicate that longitude increases clockwise when viewed from above
               the planet's north pole.
               
      RETURNS: Nothing.

        NOTES: See PlanetographicCoordinates().
               Notes and references for individual objects follow:
               
               Sun:
               
    REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 398.

               Mercury:
               
        NOTES: Uses the formulae of the 1985 IAU working group on cartographic
               and rotational elements, identical to the 1991 formulae.
   REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 402.

               Venus:
               
         NOTE: See PlanetographicCoordinates().
               Uses the method of the 1985 IAU working group on cartographic
               and rotational elements.  The 1991 forumla for (w) differs.
   REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 402.

               Earth:
               
         NOTE: See PlanetographicCoordinates().  Additionally, the constant
               term in the expression for (w1) had to be changed from 100.21
               to 190.21 degrees to obtain results consistent with current
               terrestrial longitude convention.
   REFERENCES: The Astronomical Almanac for the Year 1984, p. S30.

               Moon:
               
         NOTE: See PlanetographicCoordinates().
               Uses the formulae of the 1985 IAU working group on cartographic
               and rotational elements.
   REFERENCES: The Astronomical Almanac for the Year 1984, p. S31.

               Mars:
               
        NOTES: Uses the formulae of the 1991 IAU working group on cartographic
               and rotational elements, identical to the 1985 formulae.
   REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 403.
               The Astronomical Almanac for the Year 1984, p. S30.

               Jupiter:
               
        NOTES: Uses the formulae of the 1985 IAU working group on cartographic
               and rotational elements, identical to the 1991 formulae.
   REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 403.
               The Astronomical Almanac for the Year 1984, p. S30.

               Saturn:
               
        NOTES: Uses the method of the 1985 IAU working group on cartographic
               and rotational elements, identical to the 1991 formulae.
   REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 404.
               The Astronomical Almanac for the Year 1984, p. S30.

               Uranus:
               
         NOTE: Referred to the ecliptic, the invariable plane of the solar system,
               or Uranus's orbit plane, Uranus's rotation rate is negative.
               Referred to the Earth's equator, Uranus's rotation rate is positive,
               since its north pole lies south of the Earth's equator.
               Uses the method of the 1991 IAU working group on cartographic
               and rotational elements.  The 1985 formula for (w1) differs.
   REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 405.
               The Astronomical Almanac for the Year 1984, p. S30.

               Neptune:
               
        NOTES: Uses the formulae of the 1991 IAU working group on cartographic
               and rotational elements.  The formulae of the 1985 working group
               are included but commented out.
               As with Saturn, only the rotation rate for the magnetic
               field (System III) is defined for Neptune.
   REFERENCES: Explanatory Supplement to the Astronomical Almanac, p. 406.
               The Astronomical Almanac for the Year 1984, p. S30.

               Pluto:
               
   REFERENCES: Uses the formulae of the 2009 IAU working group on cartographic
               and rotational elements.

*******************************************************************************/

void ASTROLIBDLL_API AASunRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAMercuryRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAVenusRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAEarthRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAMarsRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAJupiterRotation(double, double *, double *, double *, double *, char);
void ASTROLIBDLL_API AASaturnRotation(double, double *, double *, double *, double *, char);
void ASTROLIBDLL_API AAUranusRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AANeptuneRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAPlutoRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAMoonRotation(double, double *, double *, double *, double *);

/*** AAAsteroidRotation *************************************************************
	 
	DESCRIPTION: Computes rotational elements for numbered asteroids at a given time.
	 
	void AAAsteroidRotation ( int number, double jd, double *a1, double *d1, double *w1, double *wd );

 	(number): asteroid number (1=Ceres, 2=Pallas, 4=Vesta, ...)
		(jd): Desired time, expressed as a Julian Date.
	    (a1): Receives J2000 right ascension of object's north pole, in radians.
	    (d1): Receives J2000 declination of object's north pole, in radians.
	    (w1): Receives the argument of object's prime meridian, in radians.
	    (wd): Receives the object's rotation rate, in radians per day.

	The function returns TRUE if rotational elements for the specified asteroid
	number are known, and FALSE otherwise.  If the function returns FALSE, output
    variables a1, d1, w1, wd will not be modified.
 
	Negative values for w1 and wd indicate anticlockwise rotation.

 REFERENCES: IAU Working Group Report on Cartographic and Rotational Elements, 2009.
             IAU Working Group Reccommended Coordinate System for (4) Vesta, 2011.

************************************************************************************/
	
int ASTROLIBDLL_API AAAsteroidRotation(int, double, double *, double *, double *, double *);
	
/*** AAAngularDiameter *******************************************************

	Computes an object's angular radius or diameter as seen from a given
	distance.
	
	double AAAngularDiameter ( double diameter, double distance )
	double AAAngularRadius ( double radius, double distance )
	
	(radius): object's physical radius. 
	(diameter): object's physical diameter. 
	(distance): object's distance, in same units as diameter.
	
	The function returns the object's angular diameter at the distance
	given in the (distance) parameter.  Note that the distance refers
	to the distance from the center of the object, not its surface.
 
	If seen from "inside" the object (i.e. distance < radius), the angular
	radius returned will be PI and returned diameter will be TWO_PI, e.g.
	the object fills the entire sky.
 
**************************************************************************/

double ASTROLIBDLL_API AAAngularDiameter ( double, double );
double ASTROLIBDLL_API AAAngularRadius ( double, double );

/**************************  AAPhaseAngle  **********************************

	Computes the angle between the sun and the viewer as seen from the
	object the viewer is observing.
	
	double AAPhaseAngle ( double sun[3], double viewer[3] )
	
	(sun):    sun's position vector relative to the object
	(viewer): viewer's position vector relative to the object
	
	The function returns the phase angle, in radians.  This is the
	angle between the sun and the viewer as seen from the object,
	and it determines the amount of the object's apparent disk that
	appears illuminated.  An object is fully illminated at a phase
	angle of zero degrees; it is completely dark at a phase angle
	of 180 degrees.
	
	Note that the input vectors may be expressed in any coordinate system,
	as long as they are the same for both vectors.  The units of the
	vectors need not be consistent; this function is only concerned
	with the directions of these vectors, not their magnitudes.

**************************************************************************/

double ASTROLIBDLL_API AAPhaseAngle ( double *, double * );

/*** AAIlluminatedFraction ***********************************************

	Computes the illuminated fraction of a planet's disk from its phase
	angle.
	
	double AAIlluminatedFraction ( double phase )
	
	(phase): phase angle, in radians from zero to PI.
	
	The function returns the illuminated fraction of the planet's
	apparent disk as a floating-point value from zero to one.

**************************************************************************/

double ASTROLIBDLL_API AAIlluminatedFraction ( double );

/*** AAPlanetMagnitude *********************************************************

	Computes the apparent magnitude of the major planets, given their distances
	from the sun and viewer, and phase angle.
	
	double AAMercuryMagnitude ( double r, double d, double phase )
	double AAVenusMagnitude ( double r, double d, double phase )
	double AAMarsMagnitude ( double r, double d, double phase )
	double AAJupiterMagnitude ( double r, double d, double phase  )
	double AASaturnMagnitude ( double r, double d, double phase, double b0 )
	double AAUranusMagnitude ( double r, double d, double phase )
	double AANeptuneMagnitude ( double r, double d, double phase )
	double AAPlutoMagnitude ( double r, double d, double phase )
	
	(r):     object's distance from sun, in AU.
	(d):     object's distance from viewer, in AU.
	(phase): object's phase angle, in radians.
	(b0):    Saturn's central latitude as seen by viewer.
	
	The functions return the planet's apparent magnitude using the formulae
	given in the Explanatory Supplement to the Astronomical Almanac, with
	the following modification: the absolute magnitudes used for the planets
	are those given in the Astronomical Almanac for the Year 1984, p. E88,
	and in Meeus, p. 270.  These values provide results which agree more
	closely with those actually tabulated in Astronomical Almanacs published
	after 1984 than the equivalent values listed in the Explanatory Supplement.
	
	For Mercury, the formula is only designed to produce realistic results
	for phase angles between 3 and 123 degrees.  Outside that range, Mercury
	is too close to the sun to be observed photometrically.  (See the Explanatory
	Supplement to the Astronomical Almanac, page 402).
	
	For Earth, the forumla uses an absolue magitude of -3.86 for a fully-illuminated
	Earth at a distance of 1 AU from both Sun and observer.  It takes the Earth's
	distance from both Sun and observer into account, but does not correct for phase
	angle (usually not an issue, since one rarely observes the Earth from outside).
 
	Note that for Saturn, the apparent magnitude is strongly dependent on
	the amount of light reflected by the rings, which in turn depends on the
	apparent inclination of the ring plane to the viewer.  Since the rings
	orbit in Saturn's equatorial plane, their inclination is equivalent to the
	planetographic latitude at the center of Saturn's apparent disk.  You
	can compute this value accurately enough for magnitude estimates with
	the function AASaturnRingPlaneInclination().
	    
	References:
   
	The Explanatory Supplement to the Astronomical Almanac, pp. 401-407.
	The Astronomical Almanac for the Year 1984, p. E88.
	Jean Meeus, "Astronomical Algorithms", pp. 269-270.
	    
*******************************************************************************/

double ASTROLIBDLL_API AAMercuryMagnitude ( double, double, double );
double ASTROLIBDLL_API AAVenusMagnitude ( double, double, double );
double ASTROLIBDLL_API AAEarthMagnitude ( double, double, double );
double ASTROLIBDLL_API AAMarsMagnitude ( double, double, double );
double ASTROLIBDLL_API AAJupiterMagnitude ( double, double, double );
double ASTROLIBDLL_API AASaturnMagnitude ( double, double, double, double );
double ASTROLIBDLL_API AAUranusMagnitude ( double, double, double );
double ASTROLIBDLL_API AANeptuneMagnitude ( double, double, double );
double ASTROLIBDLL_API AAPlutoMagnitude ( double, double, double );

/**** AASaturnRingPlaneInclination *******************************************

	Computes the apparent inclination of Saturn's ring plane.
	
	double AASaturnRingPlaneInclination ( double p[3], double d )
	
	(p): Saturn's position vector relative to viewer.
	(d): magnitude of vector p.
	
	The function returns the inclination of the ring plane to the viewer,
	in radians.  This is equivalent to the Saturnographic latitude at the
	center of Saturn's disk as seen from the viewer's position.  A return
	value of zero indicates that the rings are being viewed edge-on; a value
	of PI / 2 radians (i.e. 90 degrees) would indicate that the rings are
	being viewed from above Saturn's north pole.
	
	Saturn's position vector relative to the viewer (p) should be expressed
	in the J2000 equatorial coordinate frame.  The parameter (d) contains
	the magnitude of the vector, i.e. Saturn's distance from the viewer in
	the same units as (p).
	
	The function computes the ring plane inclination assuming fixed J2000
	coordinates of 40.66 degrees and 83.52 degrees for the right ascension
	and declination of Saturn's north pole.
	
	References:
	
	The Astronomical Almanac for the Year 1984, p. S30.
	
******************************************************************************/

double ASTROLIBDLL_API AASaturnRingPlaneInclination ( double *, double );

/*** AASunMagnitude ***********************************************************

	Computes the apparent magnitude of the sun.
	
	double AASunMagnitude ( double d )

	(d): sun's distance from viewer, in AU.

	The function returns the Sun's apparent magnitude.  It is computed assuming
	an apparent magnitude of -26.72 at a distance of 1 AU.

*******************************************************************************/

double ASTROLIBDLL_API AASunMagnitude ( double );

/*** AAMoonMagnitude *********************************************************

	Computes the apparent magnitude of the moon, given its distance
	from the sun and viewer, and phase angle.

	double AAMoonMagnitude ( double r, double d, double phase )
	
	(d):     moon's distance from viewer, in AU.
	(r):     moon's distance from sun, in AU.
	(phase): moon's phase angle, in radians.
	
	The function returns the Moon's apparent magnitude.  The Explanatory
	Supplement to the Astronomical Almanac does not give an expression
	for the Moon's apparent magnitude.  This function computes the Moon's
	apparent magnitude as if the Moon were an asteroid with an absolute
	magnitude (H) of +0.21 and a phase parameter (G) of +0.25.  The absolute
	magnitude is taken from the Astronomical Almanac for the Year 1984,
	p. E88, and G = +0.25 is a typical value assumed for most asteroids.
	
	References:
	
	The Astronomical Almanac for the Year 1984, p. E88.
	The Explanatory Supplement to the Astronomical Almanac, p. 311.
	Jean Meeus, "Astronomical Algorithms", p. 217.
	
*******************************************************************************/
	    	
double ASTROLIBDLL_API AAMoonMagnitude ( double, double, double );

/*** AAAsteroidMagnitude  ******************************************************

	Computes the apparent magnitude of an asteroid, given its distance from
	the sun and viewer, phase angle, absolute magnitude, and phase parameter.
	
	double AAAsteroidMagnitude ( double r, double d, double phase, double h, double g )
	
	(d):     asteroid's distance from viewer, in AU.
	(r):     asteroid's distance from sun, in AU.
	(phase): asteroid's phase angle, in radians.
 	(h):     asteroid's absolute magnitude parameter.
	(g):     asteroid's phase parameter.
	
	The functions returns the asteroid's apparent magnitude.  This formula
	is only valid for phase angles less than 120 degrees.  It may return
	INFINITY if the input parameters result in an invalid magnitude.
		    
	The asteroid's absolute magnitude (h) is defined as its apparent magnitude
	at full illumination, at a distance of 1 AU from both the viewer and sun.
	The phase parameter (g) has no physical meaning, but describes how the
	magnitude varies with the phase angle.
	
	References:
   
	The Explanatory Supplement to the Astronomical Almanac, p. 311.
	Jean Meeus, "Astronomical Algorithms", p. 217.
	    
*******************************************************************************/

double ASTROLIBDLL_API AAAsteroidMagnitude ( double, double, double, double, double );

/*** AACometMagnitude ***********************************************************

	Computes the apparent magnitude of a comet, given its distance from
	the sun and viewer, phase angle, absolute magnitude, and heliocentric
	distance parameter.
	
	double AACometMagnitude ( double r, double d, double h, double k )
	
	(d):     comet's distance from viewer, in AU.
	(r):     comet's distance from sun, in AU.
 	(h):     comet's absolute magnitude parameter.
	(k):     comet's heliocentric distance parameter.
	
	The functions returns the comet's apparent magnitude.  Note that comet
	magnitudes are extremely difficult to predict accurately, since they are
	heavily dependent on the amount of material released from the comet's
	surface.  Values computed by this function should only be considered
	rough estimates.
		    
	The comet's absolute magnitude (h) is defined as its apparent magnitude
	at full illumination, at a distance of 1 AU from both the viewer and sun.
	The phase parameter (k) describes how the absolute magnitude varies with
	the comet's distance from the sun.
	
	References:
   
	Jean Meeus, "Astronomical Algorithms", p. 216.
	    
*******************************************************************************/

double ASTROLIBDLL_API AACometMagnitude ( double, double, double, double );

/*** AASatelliteMagnitude *****************************************************

	Computes the apparent magnitude of an artificial Earth-orbiting satellite.

	double AASatelliteMagnitude ( double d, double b, double k )

	(d): satellite's distance from observer, in km.
	(b): Sun-satellite-observer phase angle, in radians.
	(h): satellite's absolute magnitude.

	The function returns the satellite's apparent magnitude at the specified
	distance and phase angle.  Note that this function does not take into account
	whether or not the satellite is in the Earth's shadow, in which case the
	satellite would appear infinitely faint.
 
	The satellite's absolute magnitude (h) is defined as its apparent magnitude
	at 50% illumination (phase angle PI/2 radians or 90 degrees), at a distance
	of 1000 km from the observer.

*******************************************************************************/
	
double ASTROLIBDLL_API AASatelliteMagnitude ( double, double, double );
	
/*** AAAbsoluteMagnitude ******************************************************

	Returns an object's absolute magnitude, given its apparent magnitude and
	distance in parsecs.
	
	double AAAbsoluteMagnitude ( double m, double d )

	(m): apparent magnitude
	(d): distance in parsecs
	
	The function returns the absolute magnitude.  If the distance is zero or
	infinite, the function returns infinity.
	
*******************************************************************************/

double ASTROLIBDLL_API AAAbsoluteMagnitude ( double m, double d );

/*** AAApparentMagnitude ******************************************************

	Returns an object's apparent magnitude, given its absolute magnitude and
	distance in parsecs.
	
	double AAApparentMagnitude ( double m, double d )

	(m): absolute magnitude
	(d): distance in parsecs
	
	The function returns the apparent magnitude.  If the distance is zero or
	infinite, the function returns infinity.
	
*******************************************************************************/

double ASTROLIBDLL_API AAApparentMagnitude ( double m, double d );

/*** AAMagnitudeDistance ******************************************************
 
	Returns an object's distance in parsecs from the difference between its
	apparent and absolute magnitudes.
 
	double AAMagnitudeDistance ( double m, double M )
 
	(m): apparent magnitude
	(M): absolute magnitude
 
	The function returns the distance in parsecs, which may be infinite.
 
*******************************************************************************/

double ASTROLIBDLL_API AAMagnitudeDistance ( double m, double M );
	
/*** AAMagnitudeRatio **********************************************************

	Returns the brightness ratio that corresponds to the difference between
	two different magnitudes.

	double AAMagnitudeRatio ( double m1, double m2 )
	
	(m1): first magnitude
	(m2): second magnitude
	
	The function returns the brightness ratio that corresponds to the
	magnitude difference (m2 - m1).  If m2 > m1, the ratio will be > 1;
	if m2 < m1, the ratio will be < 1.  If m1 is infinitely faint, the
	ratio will be zero; if m2 is infinitely faint, the ratio will be
	infinite.
	
********************************************************************************/

double ASTROLIBDLL_API AAMagnitudeRatio ( double m1, double m2 );

/*** AACombinedMagnitude *******************************************************

	Returns the combined magnitude of two objects.

	double AACombinedMagnitude ( double m1, double m2 )
	
	(m1): magnitude of first object
	(m2): magnitude of second object
	
	The function returns the magnitude corresponds to the combined light
	of (m1 + m2).  Note that this function is commutative, i.e.:
	
	AACombinedMagnitude ( m1, m2 ) == AACombinedMagnitude ( m2, m1 )

	If either magnitude is infinitely faint, the function returns the
	other magnitude.
	
********************************************************************************/

double ASTROLIBDLL_API AACombinedMagnitude ( double m1, double m2 );

/*** AALineEllipsoidIntersect ****************************************************

	Determines if a line from an external point intersects an oblate ellipsoid.

	int AALineEllipsoidIntersect ( double pxyz[3], double uxyz[3],
		 double exyz[3], double re, double f, double *d, double dxyz[3] )

	(pxyz): coordinates of external point
 	(uxyz): unit vector describing direction in which line extends from point
	(exyz): coordinates of ellipsoid center
 	(re):   equatorial radius of ellipsoid
 	(f):	ellipsoid flattening factor
 	(d):	receives distance from pxyz to intersection point (if any)
 	(dxyz): receives coordinates of intersection point (if any)
 
 	Returns TRUE if ray intersects the ellipsoid, or FALSE otherwise.
	If vector extends away from ellipsoid center, returns FALSE.
 
	Adapted from https://gis.stackexchange.com/questions/20780/point-of-intersection-for-a-ray-and-earths-surface

**********************************************************************************/

int ASTROLIBDLL_API AALineEllipsoidIntersect ( double pxyz[3], double uxyz[3],
												double exyz[3], double re, double f,
												double *d, double dxyz[3] );

/*** AARedshiftToRadialVelocity *************************************************
 
	Given a positive or negative redshift (z), returns the equivalent radial
	velocity as a fraction of the speed of light (v), and vice-versa.
 
	double AARedshiftToRadialVelocity ( double z )
	double AARadialVelocityToRedshift ( double v )
 
	Both functions use relativistic formulae, so redshifts greater than 1.0
	correspond to radial velocities less than 1.0.
 
********************************************************************************/

double ASTROLIBDLL_API AARedshiftToRadialVelocity ( double z );
double ASTROLIBDLL_API AARadialVelocityToRedshift ( double v );

/*** functions in Spectrum.c ***/

#define SPECTRAL_TYPE_W0	  0
#define SPECTRAL_TYPE_O0	 10
#define SPECTRAL_TYPE_B0 	 20
#define SPECTRAL_TYPE_A0	 30
#define SPECTRAL_TYPE_F0	 40
#define SPECTRAL_TYPE_G0	 50
#define SPECTRAL_TYPE_K0	 60
#define SPECTRAL_TYPE_M0	 70
#define SPECTRAL_TYPE_R0	 80
#define SPECTRAL_TYPE_N0	 90
#define SPECTRAL_TYPE_S0	100
#define SPECTRAL_TYPE_C0	110
	
#define LUMINOSITY_CLASS_IA0	1
#define LUMINOSITY_CLASS_IA		2
#define LUMINOSITY_CLASS_IAB	3
#define LUMINOSITY_CLASS_IB		4
#define LUMINOSITY_CLASS_II		5
#define LUMINOSITY_CLASS_III	6
#define LUMINOSITY_CLASS_IV		7
#define LUMINOSITY_CLASS_V		8
	
int ASTROLIBDLL_API ParseSpectrumString ( char *, int *, int * );
int ASTROLIBDLL_API GetSpectralType ( char * );
int ASTROLIBDLL_API GetLuminosityClass ( char * );
int ASTROLIBDLL_API GetSpectralTypeData ( int spec, int lum, float *mv, float *bmv, float *temp, float *bc, float *mass );
int ASTROLIBDLL_API SpectralTypeToAbsoluteMagnitude ( int, int, float * );
int ASTROLIBDLL_API SpectralTypeToColorIndex ( int, int, float * );
int ASTROLIBDLL_API SpectralTypeToTemperature ( int, int, float * );
int ASTROLIBDLL_API SpectroscopicParallaxDistance ( int, int, float, float, float * );

/* MarsMoon.c */

void ASTROLIBDLL_API AASetPhobosMatrix ( double[3][3], double);
void ASTROLIBDLL_API AASetDeimosMatrix ( double[3][3], double);
void ASTROLIBDLL_API AAPhobosOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AADeimosOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAPhobosRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AADeimosRotation(double, double *, double *, double *, double *);

/* JupiMoon.c */

void ASTROLIBDLL_API AASetJupiterMoonMatrix ( double[3][3] );
void ASTROLIBDLL_API AAIoXYZ(double, double *, double *, double *);
void ASTROLIBDLL_API AAEuropaXYZ(double, double *, double *, double *);
void ASTROLIBDLL_API AAGanymedeXYZ(double, double *, double *, double *);
void ASTROLIBDLL_API AACallistoXYZ(double, double *, double *, double *);
void ASTROLIBDLL_API AAIoRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAEuropaRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAGanymedeRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AACallistoRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAAmaltheaRotation(double, double *, double *, double *, double *);
int ASTROLIBDLL_API AAFindJupiterMoonEvents ( double jd0, double jd1, JupiterMoonEvent *events, int maxEvents );
double ASTROLIBDLL_API AAJupiterGRSLongitude ( double jd );

/* SatuMoon.c */

void ASTROLIBDLL_API AASetSaturnMoonMatrix ( double[3][3] );
void ASTROLIBDLL_API AAMimasOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAEnceladusOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AATethysOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AADioneOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AARheaOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AATitanOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAHyperionOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAIapetusOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAPhoebeOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAMimasRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAEnceladusRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AATethysRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AADioneRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AARheaRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AATitanRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAHyperionRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAIapetusRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAPhoebeRotation(double, double *, double *, double *, double *);

/* UranMoon.c */

void ASTROLIBDLL_API AASetUranusMoonMatrix ( double[3][3] );
void ASTROLIBDLL_API AAMirandaOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAArielOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAUmbrielOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AATitaniaOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAOberonOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AAMirandaRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAArielRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAUmbrielRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AATitaniaRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AAOberonRotation(double, double *, double *, double *, double *);

/* NeptMoon.c */

void ASTROLIBDLL_API AASetTritonMatrix ( double [3][3], double);
void ASTROLIBDLL_API AASetNereidMatrix ( double [3][3] );
void ASTROLIBDLL_API AATritonOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AANereidOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AACharonOrbit(double, double *, double *, double *, double *, double *, double *, double *);
void ASTROLIBDLL_API AATritonRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AANereidRotation(double, double *, double *, double *, double *);
void ASTROLIBDLL_API AACharonRotation(double, double *, double *, double *, double *);

#ifdef __cplusplus
}
#endif

#endif	/* ASTROLIB_H */

/*** COPYRIGHT NOTICE *********************************************************
 
	Copyright (c) 1992-2013 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

#include "AstroLib.h"

/*** AAPlanetographicCoordinates ***************************************************/

void AAPlanetographicCoordinates (double a, double d, double a1, double d1, double w1, double *l, double *b, double *p)
{
  double sa, ca, sd, cd, sd1, cd1, k;
	
  sa = sin ( a1 - a );
  ca = cos ( a1 - a );
  sd = sin ( d );
  cd = cos ( d );
  sd1 = sin ( d1 );
  cd1 = cos ( d1 );
	
  *b = asin ( - sd1 * sd - cd1 * cd * ca );
  *p = atan2pi ( cd1 * sa, sd1 * cd - cd1 * sd * ca );

  k = atan2pi ( - cd1 * sd + sd1 * cd * ca, cd * sa );
  if ( w1 >= 0.0 )
    *l = Mod2Pi ( w1 - k );
  else
    *l = Mod2Pi ( k + w1 );
}

/*** AASetPlanetographicRotationMatrix ***********************************************/

void AASetPlanetographicRotationMatrix ( double m[3][3], double a1, double d1, double w1, int i )
{
	double j1, n1;
	
	/*** Compute inclination of planet's equator to J2000 equator, and
	     longitude of ascending node of planet's equator on J2000 equator. ***/

	j1 = HALF_PI - d1;
	n1 = a1 + HALF_PI;
	
	/*** Take absolute value of argument of prime meridian (negative values
	 indicate anticlockwise longitude convention; this is irrelevant here
	 since we need to use the actual, physical orientation of the prime
	 meridian in the J2000 equatorial coordinate system. ***/
	
	if ( w1 < 0.0 )
		w1 = -w1;
	
	/*** Compute matrix which transforms coordinates from J2000 equatorial frame
	     to planetographic frame, or vice-versa. ***/
		
	if ( i > 0 )	
		AASetRotationMatrix ( m, 3, 2, -n1, 0, -j1, 2, -w1 );
	else
		AASetRotationMatrix ( m, 3, 2, w1, 0, j1, 2, n1 );
}

/*** AAVectorPlanetographicCoordinates  ***********************************************/

void AAVectorPlanetographicCoordinates ( double v[3], double a1, double d1, double w1,
double *l0, double *b0, double *pN )
{
	double m[3][3], v0[3];
	
	/*** Compute viewer's direction vector as seen from the planet. ***/
	
	v0[0] = -v[0];
	v0[1] = -v[1];
	v0[2] = -v[2];

	/*** Compute matrix which transforms coordinates from J2000 equatorial frame
	     to planetographic frame.  Then transform viewer's direction vector to
	     planetographic frame. ***/

	AASetPlanetographicRotationMatrix ( m, a1, d1, w1, 1 );
	AATransformVector ( m, v0 );
	AAXYZToSpherical ( v0[0], v0[1], v0[2], l0, b0, NULL );
	
	/*** For planets where longitude increases clockwise when viewed from above
	     their north pole, negate the longitude (the transformation matrix was
	     computed assuming a right-handed coordinate system; clockwise-increasing
	     longitude implies a left-handed coordinate system). ***/
	
	if ( w1 > 0.0 )
		*l0 = TWO_PI - *l0;

	/*** Compute position angle from planet's direction vector to planet's north
	     pole vector (which corresponds to bottom row, i.e. Z-axis vector, of
	     J2000-to-planetographic matrix computed above. ***/
	     
	*pN = AAVectorPositionAngle ( v, m[2] );
}

/*** AASunRotation ****************************************************************/
         
void AASunRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, t;
	
	d = jd - 2451545.0;
	t = d / 36525.0;
	
	*a1 = 286.13 * RAD_PER_DEG;
	*d1 =  63.87 * RAD_PER_DEG;
	*w1 = -Mod2Pi ( ( 84.10 + 14.1844000 * d ) * RAD_PER_DEG );
	*wd = 14.1844000 * RAD_PER_DEG;
}

/*** AAMercuryRotation *****************************************************/

void AAMercuryRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, t;

	d = jd - 2451545.0;
	t = d / 36525.0;

	*a1 = ( 281.01 - 0.033 * t ) * RAD_PER_DEG;
	*d1 = (  61.45 - 0.005 * t ) * RAD_PER_DEG;
	*w1 = Mod2Pi ( ( 329.71 + 6.1385025 * d ) * RAD_PER_DEG );
	*wd = 6.1385025 * RAD_PER_DEG;
}

/*** AAVenusRotation *******************************************************/

void AAVenusRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, t;

	d = jd - 2451545.0;
	t = d / 36525.0;

	*a1 = 272.78 * RAD_PER_DEG;
	*d1 =  67.21 * RAD_PER_DEG;	
	*w1 = -Mod2Pi ( ( 159.91 - 1.4814205 * d ) * RAD_PER_DEG );
	*wd = -1.4814205 * RAD_PER_DEG;
	
  /*** *w1 = -Mod2Pi ( ( 160.26 - 1.4813596 * d ) * RAD_PER_DEG ); 1991 value ***/
}

/*** AAEarthRotation *******************************************************/

void AAEarthRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, t;

	d = jd - 2451545.0;
	t = d / 36525.0;

	*a1 = (  0.00 - 0.641 * t ) * RAD_PER_DEG;
	*d1 = ( 90.00 - 0.557 * t ) * RAD_PER_DEG;
	*w1 = -Mod2Pi ( ( 190.21 + 360.9856123 * d ) * RAD_PER_DEG );
	*wd = 360.9856123 * RAD_PER_DEG;
}

/*** AAMarsRotation  **********************************************************/

void AAMarsRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, t;
	
	d = jd - 2451545.0;
	t = d / 36525.0;
	
	*a1 = ( 317.681 - 0.108 * t ) * RAD_PER_DEG;
	*d1 = (  52.886 - 0.061 * t ) * RAD_PER_DEG;
	*w1 = Mod2Pi ( ( 176.868 + 350.8919830 * d ) * RAD_PER_DEG );
	*wd = 350.8919830 * RAD_PER_DEG;
}

/*** AAJupiterRotation ******************************************************/

void AAJupiterRotation ( double jd, double *a1, double *d1, double *w1, double *wd, char system )
{
	double d, t;

	d = jd - 2451545.0;
	t = d / 36525.0;

	*a1 = ( 268.05 - 0.009 * t ) * RAD_PER_DEG;
	*d1 = (  64.49 + 0.003 * t ) * RAD_PER_DEG;
	if ( system == 1 )
	{
		*w1 = Mod2Pi ( (  67.10 + 877.9000000 * d ) * RAD_PER_DEG );
		*wd = 877.9000000 * RAD_PER_DEG;
	}
	else if ( system == 2 )
	{
		*w1 = Mod2Pi ( (  43.30 + 870.2700000 * d ) * RAD_PER_DEG );
		*wd = 870.2700000 * RAD_PER_DEG;
	}
	else
	{
		*w1 = Mod2Pi ( ( 284.95 + 870.5360000 * d ) * RAD_PER_DEG );
		*wd = 870.5360000 * RAD_PER_DEG;
	}
}

/*** AASaturnRotation **********************************************************/

void AASaturnRotation ( double jd, double *a1, double *d1, double *w1, double *wd, char system )
{
	double d, t;

	d = jd - 2451545.0;
	t = d / 36525.0;

	*a1 = ( 40.66 - 0.036 * t ) * RAD_PER_DEG;
	*d1 = ( 83.52 - 0.004 * t ) * RAD_PER_DEG;
	if ( system == 1 )
	{
		*w1 = Mod2Pi ( ( 227.2037 + 844.3000000 * d ) * RAD_PER_DEG );
		*wd = 844.3000000 * RAD_PER_DEG;
	}
	else
	{
		*w1 = Mod2Pi ( (  38.9000 + 810.7939024 * d ) * RAD_PER_DEG );
		*wd = 810.7939024 * RAD_PER_DEG;
	}
}

/*** AAUranusRotation **********************************************************/

void AAUranusRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, t;

	d = jd - 2451545.0;
	t = d / 36525.0;

	*a1 = 257.43 * RAD_PER_DEG;
	*d1 = -15.10 * RAD_PER_DEG;
	*w1 = Mod2Pi ( ( 203.81 - 501.1600928 * d ) * RAD_PER_DEG );
	*wd = -501.1600928 * RAD_PER_DEG;
	
	/**************************** 1985 value **********************************
  
	*w1 = Mod2Pi ( ( 261.62 - 554.9130000 * d ) * RAD_PER_DEG );
  
	***************************************************************************/
}

/*** AANeptuneRotation ***************************************************************/

void AANeptuneRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, t, n;
	
	d = jd - 2451545.0;
	t = d / 36525.0;

	n = Mod2Pi ( 357.85 + 52.316 * t ) * RAD_PER_DEG;
	*a1 = ( 299.36 + 0.07 * sin ( n ) ) * RAD_PER_DEG;
	*d1 = (  43.46 - 0.51 * cos ( n ) ) * RAD_PER_DEG;
	*w1 = Mod2Pi ( ( 253.18 + 536.3128492 * d - 0.48 * sin ( n ) ) * RAD_PER_DEG );
	*wd = 536.3128492 * RAD_PER_DEG;
	
	/******************************  1985 formulae  **********************************

	*a1 = 295.33 * RAD_PER_DEG;
	*d1 =  40.65 * RAD_PER_DEG;
	*w1 = Mod2Pi ( ( 107.21 + 48.75000000 * d ) * RAD_PER_DEG );

	**********************************************************************************/
}

/*** AAPlutoRotation ******************************************************************
 
	Revised formulae from 2009 IAU COSPAR working group.
 
***************************************************************************************/

void AAPlutoRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d;
	
	d = jd - 2451545.0;

	*a1 = 132.993 * RAD_PER_DEG;
	*d1 =  -6.163 * RAD_PER_DEG;
	*w1 = Mod2Pi ( ( 237.305 + 56.3625225 * d ) * RAD_PER_DEG );
	*wd = 56.3625225 * RAD_PER_DEG;
}

/*** AAMoonRotation ***************************************************************/

void AAMoonRotation ( double jd, double *a1, double *d1, double *w1, double *wd )
{
	double d, e1, e2, e3, e4, e5;

	d = jd - 2451545.0;
		
	e1 = Mod2Pi ( ( 125.045 -  0.052992 * d ) * RAD_PER_DEG );
	e2 = Mod2Pi ( ( 249.390 -  0.105948 * d ) * RAD_PER_DEG );
	e3 = Mod2Pi ( ( 196.694 - 13.012000 * d ) * RAD_PER_DEG );
	e4 = Mod2Pi ( ( 176.630 + 13.340716 * d ) * RAD_PER_DEG );
	e5 = Mod2Pi ( ( 358.219 -  0.985600 * d ) * RAD_PER_DEG );
	
	*a1 = ( 270.000
		- 3.878 * sin ( e1 )
		- 0.120 * sin ( e2 )
		+ 0.070 * sin ( e3 )
		- 0.017 * sin ( e4 ) ) * RAD_PER_DEG;
	*d1 = ( 66.534
		+ 1.543 * cos ( e1 )
		+ 0.024 * cos ( e2 )
		- 0.028 * cos ( e3 )
		+ 0.007 * cos ( e4 ) ) * RAD_PER_DEG;
	*w1 = -Mod2Pi ( ( 38.314 + 13.1763581 * d
		+ 3.558 * sin ( e1 )
		+ 0.121 * sin ( e2 )
		- 0.064 * sin ( e3 )
		+ 0.016 * sin ( e4 )
		+ 0.025 * sin ( e5 ) ) * RAD_PER_DEG );
	
	*wd = 13.1763581 * RAD_PER_DEG;
}

/*** AAAsteroidRotation ***************************************************************/

int AAAsteroidRotation ( int number, double jd, double *a1, double *d1, double *w1, double *wd )
{
	double	d = jd - J2000;
	
	if ( number == 1 )	// Ceres
	{
		*a1 = 291 * RAD_PER_DEG;
		*d1 =  59 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 170.90 + 952.1532 * d ) * RAD_PER_DEG );
		*wd = 952.1532 * RAD_PER_DEG;
	}
	else if ( number == 2 )	// Pallas
	{
		*a1 = 33 * RAD_PER_DEG;
		*d1 = -3 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 38 + 1105.8036 * d ) * RAD_PER_DEG );
		*wd = 1105.8036 * RAD_PER_DEG;
	}
	else if ( number == 4 )	// Vesta, reccommended coordinate system 2011
	{
		*a1 = 309.031 * RAD_PER_DEG;
		*d1 =  42.235 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 285.39 + 1617.3329428 * d ) * RAD_PER_DEG );
		*wd = 1617.3329428 * RAD_PER_DEG;
	}
	else if ( number == 21 ) // Lutetia
	{
		*a1 = 52 * RAD_PER_DEG;
		*d1 = 12 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 94 + 1057.7515 * d ) * RAD_PER_DEG );
		*wd = 1057.7515 * RAD_PER_DEG;
	}
	else if ( number == 243 ) // Ida
	{
		*a1 = 168.76 * RAD_PER_DEG;
		*d1 =  -2.88 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 265.95 + 1864.6280070 * d ) * RAD_PER_DEG );
		*wd = 1864.6280070 * RAD_PER_DEG;
	}
	else if ( number == 433 ) // Eros
	{
		*a1 = 11.35 * RAD_PER_DEG;
		*d1 = 17.22 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 326.07 + 1639.38864745 * d ) * RAD_PER_DEG );
		*wd = 1639.38864745 * RAD_PER_DEG;
	}
	else if ( number == 511 ) // Davida
	{
		*a1 = 297 * RAD_PER_DEG;
		*d1 =   5 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 268.1 + 1684.4193549 * d ) * RAD_PER_DEG );
		*wd = 1684.4193549 * RAD_PER_DEG;
	}
	else if ( number == 951 ) // Gaspra
	{
		*a1 =  9.47 * RAD_PER_DEG;
		*d1 = 26.70 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 83.67 + 1226.9114850 * d ) * RAD_PER_DEG );
		*wd = 1226.9114850 * RAD_PER_DEG;
	}
	else if ( number == 25143 ) // Itokawa
	{
		*a1 =  90.53 * RAD_PER_DEG;
		*d1 = -66.30 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 000 + 712.143 * d ) * RAD_PER_DEG );
		*wd = 712.143 * RAD_PER_DEG;
	}
	else if ( number == 134340 ) // Pluto
	{
		*a1 = 132.993 * RAD_PER_DEG;
		*d1 =  -6.163 * RAD_PER_DEG;
		*w1 = Mod2Pi ( ( 237.305 + 56.3625225 * d ) * RAD_PER_DEG );
		*wd = 56.3625225 * RAD_PER_DEG;
	}
	else
		return ( FALSE );
	
	return ( TRUE );
}

/*** AAAngularDiameter *************************************************/

double AAAngularDiameter ( double diameter, double distance )
{
	return ( 2.0 * AAAngularRadius ( diameter / 2.0, distance ) );
}

/*** AAAngularRadius ***************************************************/

double AAAngularRadius ( double radius, double distance )
{
	if ( distance < radius )
		return ( PI );
	else if ( radius == 0.0 )
		return ( 0.0 );
	else
		return ( asin ( radius / distance ) );
}

/*** AAPhaseAngle *****************************************************/

double AAPhaseAngle ( double sun[3], double viewer[3] )
{
	double scale, phase = 0.0;
	
	scale = AAVectorMagnitude ( sun ) * AAVectorMagnitude ( viewer );
	if ( scale > 0.0 )
	{
		phase = AADotProduct ( sun, viewer ) / scale;
		
		if ( phase < -1.0 )
			phase = PI;
		else if ( phase > 1.0 )
			phase = 0.0;
		else
			phase = acos ( phase );
	}

	return ( phase );
}

/*** AAIlluminatedFraction ********************************************/

double AAIlluminatedFraction ( double phase )
{
	return ( ( 1.0 + cos ( phase ) ) / 2.0 );
}

/*** AASunMagnitude ****************************************************/

double AASunMagnitude ( double d )
{
	return ( -26.72 + 5.0 * log10 ( d ) );
}

/*** AAMoonMagnitude ***************************************************/

double AAMoonMagnitude ( double r, double d, double b )
{
	return ( AAAsteroidMagnitude ( r, d, b, 0.21, 0.25 ) );
}

/*** AAMercuryMagnitude *************************************************/
   
double AAMercuryMagnitude ( double r, double d, double b )
{
	double i, i2, i3;
	double m;
	
	i = b * DEG_PER_RAD;
	i2 = i * i;
	i3 = i2 * i;
	m = -0.42 + 5.0 * log10 ( r * d ) + 0.0380 * i - 0.000273 * i2 + 0.000002 * i3;
	
	return ( m );
}

/*** AAVenusMagnitude **************************************************/

double AAVenusMagnitude ( double r, double d, double b )
{
	double i, i2, i3, m;
	
	i = b * DEG_PER_RAD;
	i2 = i * i;
	i3 = i2 * i;
	m = -4.40 + 5.0 * log10 ( r * d ) + 0.0009 * i + 0.000239 * i2 - 0.00000065 * i3;
	
	return ( m );
}

/*** AAEarthMagnitude **************************************************/

double AAEarthMagnitude ( double r, double d, double b )
{
	double m;
	
	m = -3.86 + 5.0 * log10 ( r * d );
	
	return ( m );
}

/*** AAMarsMagnitude **************************************************/

double AAMarsMagnitude ( double r, double d, double b )
{
	double m, i;
	
	i = b * DEG_PER_RAD;
	m = -1.52 + 5.0 * log10 ( r * d ) + 0.016 * i;
	
	return ( m );
}

/*** AAJupiterMagnitude *************************************************/

double AAJupiterMagnitude ( double r, double d, double b )
{
	double m, i;
	
	i = b * DEG_PER_RAD;
	m = -9.40 + 5.0 * log10 ( r * d ) + 0.005 * i;
	return ( m );
}

/*** AASaturnMagnitude *********************************************/

double AASaturnMagnitude ( double r, double d, double beta, double b0 )
{
	double	m, i, sb;
	
	sb = sin ( b0 );
	i = beta * DEG_PER_RAD;
	m = -8.88 + 5.0 * log10 ( r * d ) + 0.044 * i - 2.60 * fabs ( sb ) + 1.25 * sb * sb;
	
	return ( m );
}

/*** AASaturnRingPlaneInclination **********************************/

double AASaturnRingPlaneInclination ( double p[3], double d )
{
	double	ps[3] = { -0.0856116575103, -0.0735337181506, -0.99361131052 };
	
	/*** The inclination of the ring plane is equivalent to the
	     Saturnocentric latitude at the center of the planet's disk.
	     This value, in turn, is equivalent to 90 degrees minus the
	     angle between the Saturn-to-viewer vector and Saturn's north
	     pole vector, or the angle between the viewer-to-Saturn vector
	     (which is what we have!) and Saturn's south pole vector.  We
	     can find this angle from the dot product of the two vectors.
	     
	     The elements of Saturn's south pole vector in the J2000 equatorial
	     frame were computed above using constant values of 40.66 degrees
	     and 83.52 degrees for the right ascensision and declination of
	     Saturn's north pole, respectively (see the Astronomical Almanac
	     for the Year 1984, p. S30).  They can also be obtained by
	     un-commenting the following line:

	     SphericalToXYZ ( 40.66 * RAD_PER_DEG, 83.52 * RAD_PER_DEG, -1.0,
	     &ps[0], &ps[1], &ps[2] ); ***/
	
	return ( HALF_PI - acos ( AADotProduct ( p, ps ) / d ) );
}

/*** AAUranusMagnitude **********************************************/

double AAUranusMagnitude ( double r, double d, double b )
{
	double m, i;
	
	i = b * DEG_PER_RAD;
	m = -7.19 + 5.0 * log10 ( r * d ) + 0.0028 * i;
	
	return ( m );
}

/*** AANeptuneMagnitude *********************************************/

double AANeptuneMagnitude ( double r, double d, double b )
{
	double m, i;
	
	i = b * DEG_PER_RAD;
	m = -6.87 + 5.0 * log10 ( r * d );
	
	return ( m );
}

/*** AAPlutoMagnitude ************************************************/

double AAPlutoMagnitude ( double r, double d, double b )
{
	double m, i;
	
	i = b * DEG_PER_RAD;
	m = -1.01 + 5.0 * log10 ( r * d ) + 0.041 * i;
	
	return ( m );
}

/*** AAAsteroidMagnitude ****************************************************/

double AAAsteroidMagnitude ( double r, double d, double b, double h, double g )
{
	double phi1, phi2, m;
	
	phi1 = exp ( -3.33 * pow ( tan ( b / 2.0 ), 0.63 ) );
	phi2 = exp ( -1.87 * pow ( tan ( b / 2.0 ), 1.22 ) );
	
	m = ( 1.0 - g ) * phi1 + g * phi2;
	
	if ( m > 0.0 )
		m = h + 5.0 * log10 ( r * d ) - 2.5 * log10 ( m );
	else
		m = INFINITY;
		
	return ( m );
}

/*** AACometMagnitude **********************************************/

double AACometMagnitude (double r, double d, double h, double k)
{
	double m;
	
	m = h + 5.0 * log10 ( d ) + 2.5 * k * log10 ( r );
	
	return ( m );
}

/*** AASatelliteMagnitude ****************************************************/

double AASatelliteMagnitude ( double dist, double phase, double stdmag )
{
	double mag;
	
	// Formula from http://www.prismnet.com/~mmccants/tles/mccdesc.html
	// Standard magnitude is at 1000 km range, and 50% illuminated
	
	if ( phase < PI )
		mag = stdmag - 15.75 + 2.5 * log10 ( dist * dist / ( ( 1.0 + cos ( phase ) ) / 2.0 ) );
	else
		mag = INFINITY;
	
	return ( mag );
}

/*** AAAbsoluteMagnitude *****************************************************/

double AAAbsoluteMagnitude ( double m, double d )
{
	if ( d > 0.0 && d < INFINITY )
		return ( m - 5.0 * ( log10 ( d ) - 1.0 ) );
	else
		return ( -INFINITY );
}

/*** AAApparentMagnitude ******************************************************/

double AAApparentMagnitude ( double m, double d )
{
	if ( d > 0.0 && d < INFINITY )
		return ( m + 5.0 * ( log10 ( d ) - 1.0 ) );
	else if ( d <= 0.0 )
		return ( -INFINITY );
	else
		return ( INFINITY );
}

/*** AAMagnitudeDistance *******************************************************/

double AAMagnitudeDistance ( double m, double M )
{
	return ( pow ( 10.0, 1.0 + ( m - M ) / 5.0 ) );
}

/*** AAMagnitudeRatio **********************************************************/

double AAMagnitudeRatio ( double m1, double m2 )
{
	if ( isinf ( m1 ) )
		return ( 0.0 );
	else if ( isinf ( m2 ) )
		return ( INFINITY );
	else
		return ( pow ( 10.0, ( m2 - m1 ) / 2.5 ) );
}

/*** AACombinedMagnitude ******************************************************/

double AACombinedMagnitude ( double m1, double m2 )
{
	double r = AAMagnitudeRatio ( m1, m2 );
	
	if ( isinf ( m1 ) )
		return ( m2 );
	else if ( isinf ( m2 ) )
		return ( m1 );
	else
		return ( m1 - 2.5 * log10 ( r + 1.0 ) );
}

/*** AALineEllipsoidIntersect ****************************************************/

int AALineEllipsoidIntersect ( double pxyz[3], double uxyz[3], double exyz[3], double re, double f, double *d, double dxyz[3] )
{
	double x = pxyz[0] - exyz[0];
	double y = pxyz[1] - exyz[1];
	double z = pxyz[2] - exyz[2];
	double x2 = x * x, y2 = y * y, z2 = z * z;
	double u = uxyz[0], v = uxyz[1], w = uxyz[2];
	double u2 = u * u, v2 = v * v, w2 = w * w;
	double a = re, a2 = a * a;
	double b = re * ( 1.0 - f ), b2 = b * b;
	double t = b2 * ( u * x + v * y ) + a2 * w * z;
	
    // Compute vector from exterior point to center of ellipsoid.
    // If unit vector points away from center of ellipsoid, line doesn't intersect.
    
    AAVectorDifference ( exyz, pxyz, dxyz );
    *d = AADotProduct ( dxyz, uxyz );
    if ( *d < 0.0 )
        return ( FALSE );

    // Any point (x,y,z) along a line from the external point (x0,y0,z0)
	// satisfies the equation (x,y,z) = (x0,y0,z0) + t * (u,v,w).
	// The ellipsoid's equatorial and polar radii are a and b, and the ellipsoid
	// is described by the equation (x^2/a^2) + (y^2/a^2) + (z^2/b^2) = 1.
	// Plug the first equation into the second, and we get a quadratic equation,
	// which we solve for t. There is no solution if the line does not intersect
	// the ellipsoid at all, one solution if the line is exactly tangent to the
	// ellipsoid, and two solutions for most intersection cases (where the line
	// enters the ellipsoid on one side, and exits out the opposite side).
	
	t = t * t - ( b2 * ( u2 + v2 ) + a2 * w2 ) * ( b2 * ( -a2 + x2 + y2 ) + a2 * z2 );
	if ( t < 0 )
		return ( FALSE );
	
	// We choose the solution which is closest to the external point (x0,y0,z0).
	
	t = ( -1.0 / ( b2 * ( u2 + v2 ) + a2 * w2 ) )
	  * ( b2 * ( u * x + v * y ) + a2 * w * z + sqrt ( t ) );
	
	*d = t;
	
	// Plug this value into the first equation to obtain the intersection point.
	
	dxyz[0] = pxyz[0] + t * u;
	dxyz[1] = pxyz[1] + t * v;
	dxyz[2] = pxyz[2] + t * w;
	
	return ( TRUE );
}

/*** AARedshiftToRadialVelocity ************************************************/

double AARedshiftToRadialVelocity ( double z )
{
	double z12 = ( z + 1.0 ) * ( z + 1.0 );
	double rv = ( z12 - 1.0 ) / ( z12 + 1.0 );
	
	if ( z < 0.0 )
		z = z;
	
	return ( rv );
}

/*** AARadialVelocityToRedshift ************************************************/

double AARadialVelocityToRedshift ( double rv )
{
	double z = sqrt ( ( 1.0 + rv ) / ( 1.0 - rv ) ) - 1.0;
	
	return ( z );
}

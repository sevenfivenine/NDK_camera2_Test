/*** COPYRIGHT NOTICE *********************************************************
 
	Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

#include "AstroLib.h"

/******************************  atan2pi  *********************************/

double atan2pi ( double y, double x )
{
  if ( x == 0.0 )
  {
	 if ( y > 0.0 )
	  return ( HALF_PI );
	if ( y < 0.0 )
	  return ( -HALF_PI );
  }
  else
  {
	 if ( y < 0.0 )
	  return ( atan2 ( y, x ) + TWO_PI );
	else
	  return ( atan2 ( y, x ) );
  }

  return ( 0.0 );
}

/*********************************  Mod2Pi  ************************************/

double Mod2Pi ( double x )
{
	return ( x - TWO_PI * floor ( x / TWO_PI ) );
}

double Mod24h ( double x )
{
	return ( x - 24.0 * floor ( x / 24.0 ) );
}

double Mod360 ( double x )
{
	return ( x - 360.0 * floor ( x / 360.0 ) );
}

double ModPi ( double x )
{
	x = Mod2Pi ( x );
	if ( x > PI )
		x -= TWO_PI;
	
	return ( x );
}

double Mod12h ( double x )
{
	x = Mod24h ( x );
	if ( x > 12.0 )
		x -= 24.0;
		
	return ( x );
}

double Mod180 ( double x )
{
	x = Mod360 ( x );
	if ( x > 180.0 )
		x -= 360.0;
		
	return ( x );
}

/*** AADegMinSecToDecimal *************************************************/

double AADegMinSecToDecimal ( short deg, short min, double sec, char sign )
{
	double decimal;

	decimal = abs ( deg ) + min / 60.0 + sec / 3600.0;
	if ( sign == '-' )
	  decimal = - decimal;

	return ( decimal );
}

/*** AADecimalToDegMinSec **************************************************/

void AADecimalToDegMinSec ( double decimal, short *deg, short *min, double *sec, char *sign )
{
	if ( decimal < 0.0 )
		*sign = '-';
	else
		*sign = '+';

	decimal = fabs ( decimal );
	*deg = floor ( decimal );
	*min = floor ( 60.0 * ( decimal - *deg ) );
	*sec = 3600.0 * ( decimal - *deg - *min / 60.0 );
	
	/*** Occasionally, when the seconds should be exactly zero, they can
	     instead be slightly negative due to floating point round-off error.
	     Correct these cases before returning their value. ***/
	
	if ( *sec <= 0.0 )
		*sec = 0.0;
}

/**** AADecimalToDegMin ******************************************************/

void AADecimalToDegMin ( double decimal, short *deg, double *min, char *sign )
{
	if ( decimal < 0.0 )
		*sign = '-';
	else
		*sign = '+';
	
	decimal = fabs ( decimal );
	*deg = floor ( decimal );
	*min = 60.0 * ( decimal - *deg );
	
	/*** Occasionally, when the seconds should be exactly zero, they can
	 instead be slightly negative due to floating point round-off error.
	 Correct these cases before returning their value. ***/
	
	if ( *min < 0.0 )
		*min = 0.0;
}

/*** AASeparation ********************************************************/

double haversin ( double x )
{
	double s = sin ( x * 0.5 );
	return ( s * s );
}

double AASeparation ( double a1, double d1, double a2, double d2 )
{
	double s = 0.0;

//	This spherical triangle formula is numerically inaccurate for small angles.
//	s = CLIP ( sin ( d1 ) * sin ( d2 ) + cos ( d1 ) * cos ( d2 ) * cos ( a2 - a1 ), -1.0, 1.0 );
//	s = acos ( s );
	
	s = haversin ( d2 - d1 ) + cos ( d1 ) * cos ( d2 ) * haversin ( a2 - a1 );
	s = CLIP ( s, 0.0, 1.0 );
	s = 2.0 * asin ( sqrt ( s ) );
	
	return ( s );
}

/*** AAPositionAngle ******************************************************/

double AAPositionAngle ( double a1, double d1, double a2, double d2 )
{
  double eta, xi, pa;

  eta = cos ( d2 ) * sin ( a2 - a1 );
  xi = cos ( d1 ) * sin ( d2 ) - sin ( d1 ) * cos ( d2 ) * cos ( a2 - a1 );
  pa = atan2pi ( eta, xi );

  return ( pa );
}	

/*** AAVectorSeparation ************************************************/

double AAVectorSeparation ( double u[3], double v[3] )
{
	double	w[3], d;
	
	w[0] = u[0] - v[0];
	w[1] = u[1] - v[1];
	w[2] = u[2] - v[2];
	
	d = sqrt ( w[0] * w[0] + w[1] * w[1] + w[2] * w[2] );
	d = 2.0 * asin ( d / 2.0 );
	
	return ( d );
}

/*** AAVectorPositionAngle*********************************************/

double AAVectorPositionAngle ( double u[3], double v[3] )
{
	double n[3], e[3], edotv, ndotv, pa;

	n[2] = sqrt ( 1.0 - u[2] * u[2] );
	if ( n[2] == 0.0 )
		return ( 0.0 );

	n[0] = -u[0] * u[2] / n[2];
	n[1] = -u[1] * u[2] / n[2];

	e[0] = -u[1] / n[2];
	e[1] =  u[0] / n[2];
	e[2] = 0.0;
	
	edotv = e[0] * v[0] + e[1] * v[1];
	ndotv = n[0] * v[0] + n[1] * v[1] + n[2] * v[2];
	
	if ( edotv == 0.0 && ndotv == 0.0 )
		pa = 0.0;
	else
		pa = atan2pi ( edotv, ndotv );

	return ( pa );
}

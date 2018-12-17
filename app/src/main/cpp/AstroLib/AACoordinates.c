/*** COPYRIGHT NOTICE *********************************************************
 
	Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

#include "AstroLib.h"

/*** AASphericalToXYZ **************************************************/

void AASphericalToXYZ ( double a, double d, double r, double *x, double *y, double *z )
{
	*x = r * cos ( a ) * cos ( d );
	*y = r * sin ( a ) * cos ( d );
	*z = r * sin ( d );
}

/*** AASphericalToXYZVector ********************************************/

void AASphericalToXYZVector ( double a, double d, double r, double v[3] )
{
	AASphericalToXYZ ( a, d, r, &v[0], &v[1], &v[2] );
}

/*** AAXYZToSpherical  **************************************************/

void AAXYZToSpherical ( double x, double y, double z, double *ap, double *dp, double *rp )
{
	double r;
  
	r = sqrt ( x * x + y * y + z * z );
	if ( rp )
		*rp = r;

	if ( dp )
		*dp = ( r == 0.0 ) ? 0.0 : asin ( z / r );

	if ( ap )
		*ap = ( x == 0.0 && y == 0.0 ) ? 0.0 : atan2pi ( y, x );
}

/*** AAXYZVectorToSpherical  **********************************************/

void AAXYZVectorToSpherical ( double v[3], double *ap, double *dp, double *rp )
{
	AAXYZToSpherical ( v[0], v[1], v[2], ap, dp, rp );
}

/*** AASphericalToXYZMotion ***********************************************/
	
void AASphericalToXYZMotion ( double a, double d, double r, double da, double dd,
double dr, double *x, double *y, double *z, double *dx, double *dy, double *dz )
{
  double ca, sa, cd, sd;
	
  ca = cos ( a );
  sa = sin ( a );
  cd = cos ( d );
  sd = sin ( d );
	
  *x = r * ca * cd;
  *y = r * sa * cd;
  *z = r * sd;
	
  *dx = r * ( - cd * sa * da - ca * sd * dd ) + ca * cd * dr;
  *dy = r * ( ca * cd * da - sa * sd * dd ) + cd * sa * dr;
  *dz = r * cd * dd + sd * dr; 
}

/*** AAXYZToSphericalMotion ***************************************************/

void AAXYZToSphericalMotion ( double x, double y, double z, double dx, double dy,
double dz, double *a, double *d, double *r, double *da, double *dd, double *dr )
{
  double cd;
  
  *r = sqrt ( x * x + y * y + z * z );
  if ( *r == 0.0 || ( x == 0.0 && y == 0.0 ) )
  {
    *a = *d = *dr = *dd = *da = 0.0;
  }
  else
  {
    *d = asin ( z / *r );
	*a = atan2pi ( y, x );

    cd = cos ( *d );
  
    *dr = ( x * dx + y * dy + z * dz ) / *r;
    *dd = ( *r * dz - z * *dr ) / ( cd * *r * *r );
    *da = ( x * dy - y * dx ) / ( x * x + y * y );
  }
}

/*** AAXYZVectorToSphericalMotion ***************************************************/

void AAXYZVectorToSphericalMotion ( double v[3], double dv[3],
double *a, double *d, double *r, double *da, double *dd, double *dr )
{
	return ( AAXYZToSphericalMotion ( v[0], v[1], v[2],
		dv[0], dv[1], dv[2], a, d, r, da, dd, dr ) );
}

/*** AAVectorMagnitude ***************************************************/

double AAVectorMagnitude ( double v[3] )
{
	return ( sqrt ( v[0] * v[0] + v[1] * v[1] + v[2] * v[2] ) );
}

/*** AAVectorDistance ***************************************************/

double AAVectorDistance ( double u[3], double v[3] )
{
	double w[3];
	
	AAVectorDifference ( u, v, w );
	
	return ( AAVectorMagnitude ( w ) );
}

/*** AACopyVector *******************************************************/

double *AACopyVector ( double v2[3], double v1[3] )
{
	v2[0] = v1[0];
	v2[1] = v1[1];
	v2[2] = v1[2];
	
	return ( v2 );
}

/*** AAVectorSum **********************************************************/

double *AAVectorSum ( double a[3], double b[3], double c[3] )
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
	
	return ( c );
}

/*** AAScale1VectorSum **********************************************************/

double *AAScale1VectorSum ( double ka, double a[3], double b[3], double c[3] )
{
	c[0] = ka * a[0] + b[0];
	c[1] = ka * a[1] + b[1];
	c[2] = ka * a[2] + b[2];
	
	return ( c );
}

/*** AAScale2VectorSum **********************************************************/

double *AAScale2VectorSum ( double ka, double a[3], double kb, double b[3], double c[3] )
{
	c[0] = ka * a[0] + kb * b[0];
	c[1] = ka * a[1] + kb * b[1];
	c[2] = ka * a[2] + kb * b[2];
	
	return ( c );
}

/*** AAVectorDifference **************************************************/

double *AAVectorDifference ( double a[3], double b[3], double c[3] )
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
	
	return ( c );
}

/*** AAScaleVector *******************************************************/

double *AAScaleVector ( double v[3], double k )
{
	v[0] *= k;
	v[1] *= k;
	v[2] *= k;
	
	return ( v );
}

/*** AANormalizeVector ****************************************************/

double AANormalizeVector ( double v[3] )
{
	double	m = AAVectorMagnitude ( v );
	
	if ( m > 0.0 )
		AAScaleVector ( v, 1.0 / m );

	return ( m );
}

/*** AANormalVector **************************************************/

double *AANormalVector ( double v[3], double p[3] )
{
	double	u[3] = { v[1], v[2], v[0] };
	
	AACrossProduct ( u, v, p );
	
	return ( p );
}

/*** AADotProduct *******************************************************/

double AADotProduct ( double a[3], double b[3] )
{
	return ( a[0] * b[0] + a[1] * b[1] + a[2] * b[2] );
}

/*** AACrossProduct *****************************************************/

double *AACrossProduct ( double a[3], double b[3], double c[3] )
{
	double	u[3] = { a[0], a[1], a[2] };
	double	v[3] = { b[0], b[1], b[2] };
	
	c[0] = u[1] * v[2] - u[2] * v[1];
	c[1] = u[2] * v[0] - u[0] * v[2];
	c[2] = u[0] * v[1] - u[1] * v[0];
	
	return ( c );
}

/*** AATransformVector ******************************************************/

double *AATransformVector ( double m[3][3], double v[3] )
{
	double w[3];
	short i;
  
	w[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
	w[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
	w[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];

	for ( i = 0; i < 3; i++ )
		v[i] = w[i];
		
	return ( v );
}

/*** AAUnTransformVector ***************************************************/

double *AAUnTransformVector ( double m[3][3], double w[3] )
{
	double		v[3];
	short		i;
  
	v[0] = m[0][0] * w[0] + m[1][0] * w[1] + m[2][0] * w[2];
	v[1] = m[0][1] * w[0] + m[1][1] * w[1] + m[2][1] * w[2];
	v[2] = m[0][2] * w[0] + m[1][2] * w[1] + m[2][2] * w[2];

	for ( i = 0; i < 3; i++ )
		w[i] = v[i];
		
	return ( w );
}

/*** AATransformRotationMatrix *********************************************/

void AATransformRotationMatrix ( double b[3][3], double a[3][3] )
{
	double	c[3][3];
	int		i, j, k;

	for ( i = 0; i < 3; i++ )
		for ( k = 0; k < 3; k++ )
			for ( c[i][k] = 0.0, j = 0; j < 3; j++ )
        		c[i][k] += b[i][j] * a[j][k];

	for ( i = 0; i < 3; i++ )
		for ( k = 0; k < 3; k++ )
    		a[i][k] = c[i][k];
}

/*** AAUnTransformRotationMatrix ******************************************/

void AAUnTransformRotationMatrix ( double b[3][3], double a[3][3] )
{
	double	c[3][3];
	int		i, j, k;
  
	for ( i = 0; i < 3; i++ )
		for ( k = 0; k < 3; k++ )
      		for ( c[i][k] = 0.0, j = 0; j < 3; j++ )
        		c[i][k] += b[j][i] * a[j][k];

	for ( i = 0; i < 3; i++ )
		for ( k = 0; k < 3; k++ )
			a[i][k] = c[i][k];
}

/*** AACopyRotationMatrix **************************************************/
   
void AACopyRotationMatrix ( double b[3][3], double a[3][3] )
{
	int i, j;
  
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			b[i][j] = a[i][j];
}

/*** AATransposeRotationMatrix **********************************************/
   
void AATransposeRotationMatrix ( double matrix[3][3], double transpose[3][3] )
{
	double	temp[3][3];
	int		i, j;
	
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			temp[j][i] = matrix[i][j];

	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			transpose[i][j] = temp[i][j];
}

/*** AASetRotationMatrix *****************************************************/

void AASetRotationMatrix ( double m[3][3], int n, ... )
{
	int i, k;
	double a, c, s, r[3][3];
	va_list ap;
  
	va_start ( ap, n );

	AASetIdentityRotationMatrix ( m );
  
	for ( k = 0; k < n; k++ )
	{
		i = va_arg ( ap, int );
		a = va_arg ( ap, double );
		c = cos ( a );
		s = sin ( a );
    
		switch ( i )
		{
			case 0:
				r[0][0] = 1.0;  r[0][1] = 0.0;  r[0][2] = 0.0;
				r[1][0] = 0.0;  r[1][1] =   c;  r[1][2] =  -s;
				r[2][0] = 0.0;  r[2][1] =   s;  r[2][2] =   c;
				break;
				
			case 1:
				r[0][0] =   c;  r[0][1] = 0.0;  r[0][2] =  -s;
				r[1][0] = 0.0;  r[1][1] = 1.0;  r[1][2] = 0.0;
				r[2][0] =   s;  r[2][1] = 0.0;  r[2][2] =   c;
				break;
				
			case 2:
				r[0][0] =   c;  r[0][1] =  -s;  r[0][2] = 0.0;
				r[1][0] =   s;  r[1][1] =   c;  r[1][2] = 0.0;
				r[2][0] = 0.0;  r[2][1] = 0.0;  r[2][2] = 1.0;
				break;
		}
     
		AATransformRotationMatrix ( r, m ); 
	}
   
	va_end ( ap ); 
}

/*** AASetIdentityRotationMatrix ************************************************/

void AASetIdentityRotationMatrix ( double m[3][3] )
{
	int i, j;
	
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			if ( i == j )
				m[i][j] = 1.0;
			else
				m[i][j] = 0.0;
}

/*** AARotationMatrixToOpenGLMatrix **************************************************/

double *AARotationMatrixToOpenGLMatrix ( double m[3][3], double p[16] )
{
	static double temp[16] = { 0 };
	
	if ( p == NULL )
		p = temp;
	
	p[0] = m[0][0]; p[4] = m[0][1]; p[8]  = m[0][2]; p[12] = 0.0;
	p[1] = m[1][0]; p[5] = m[1][1]; p[9]  = m[1][2]; p[13] = 0.0;
	p[2] = m[2][0]; p[6] = m[2][1]; p[10] = m[2][2]; p[14] = 0.0;
	p[3] = 0.0;     p[7] = 0.0;     p[11] = 0.0;     p[15] = 1.0;
	
	return p;
}

/*** AARotationMatrixFromOpenGLMatrix **************************************************/

void AARotationMatrixFromOpenGLMatrix ( double m[3][3], double p[16] )
{
	m[0][0] = p[0]; m[0][1] = p[4]; m[0][2] = p[8];
	m[1][0] = p[1]; m[1][1] = p[5]; m[1][2] = p[9];
	m[2][0] = p[2]; m[2][1] = p[6]; m[2][2] = p[10];
}


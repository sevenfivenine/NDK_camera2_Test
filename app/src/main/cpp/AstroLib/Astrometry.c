/*** COPYRIGHT NOTICE *********************************************************
 
	Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

#include "Astrometry.h"

/*** local utility function declarations and data types ***/

static short	A2ComputeNumTerms ( short );
static void		A2ComputeXYPolynomial ( double, double, double[10] );

/*** A2RADecToXiEta ************************************************************/

void A2RADecToXiEta ( double ra, double dec, double ra0, double dec0, double *xi, double *eta )
{
	double h;
	
	h = sin ( dec ) * sin ( dec0 ) + cos ( dec ) * cos ( dec0 ) * cos ( ra - ra0 );
	
	*xi = cos ( dec ) * sin ( ra - ra0 ) / h;
	*eta = ( sin ( dec ) * cos ( dec0 ) - cos ( dec ) * sin ( dec0 ) * cos ( ra - ra0 ) ) / h;
}

/*** A2XiEtaToRADec ************************************************************/

void A2XiEtaToRADec ( double xi, double eta, double ra0, double dec0, double *ra, double *dec )
{
	double delta, gamma;
	
	delta = cos ( dec0 ) - eta * sin ( dec0 );
	gamma = sqrt ( xi * xi + delta * delta );
	*ra = atan2 ( xi, delta ) + ra0;
	*ra = *ra - TWO_PI * floor ( *ra / TWO_PI );
	*dec = atan ( ( sin ( dec0 ) + eta * cos ( dec0 ) ) / gamma );
}

/*** A2XiEtaToXY *******************************************************************/

void A2XiEtaToXY ( double xi, double eta, double **c, double *x, double *y )
{
	int		i;
	double	p[10];
	
	A2ComputeXYPolynomial ( xi, eta, p );
	
	for ( *x = *y = 0.0, i = 0; i < 10; i++ )
	{
		*x += c[i][0] * p[i];
		*y += c[i][1] * p[i];
	}	 
}

/*** A2XiEtaToXY *****************************************************************/

void A2XYToXiEta ( double x, double y, double **c, double *xi, double *eta )
{
	int		i;
	double	p[10];
	
	A2ComputeXYPolynomial ( x, y, p );
	
	for ( *xi = *eta = 0.0, i = 0; i < 10; i++ )
	{
		*xi  += c[i][0] * p[i];
		*eta += c[i][1] * p[i];
	}	 
}

/*** A2XYToRADec ****************************************************************/

void A2XYToRADec ( double x, double y, double **c, double ra0, double dec0,
double *ra, double *dec )
{
	double xi, eta;
	
	A2XYToXiEta ( x, y, c, &xi, &eta );
	A2XiEtaToRADec ( xi, eta, ra0, dec0, ra, dec );
}

/*** A2RADecToXY ***************************************************************/

void A2RADecToXY ( double ra, double dec, double ra0, double dec0, double **c,
double *x, double *y )
{
	double xi, eta;
	
	A2RADecToXiEta ( ra, dec, ra0, dec0, &xi, &eta );
	A2XiEtaToXY ( xi, eta, c, x, y );
}

/*** A2NewAstrometricSolution ***************************************************/

int A2NewAstrometricSolution ( double ***a, double ***b, double ***c, double ***d )
{
	if ( a )
	{
		*a = NMatrix ( double, 10, 10 );
		if ( *a == NULL )
		{
			return ( FALSE );
		}
	}
	
	if ( b )
	{
		*b = NMatrix ( double, 10, 2 );
		if ( *b == NULL )
		{
			if ( a )
				NDestroyMatrix ( (void **) *a );
				
			return ( FALSE );
		}
	}
	
	if ( c )
	{
		*c = NMatrix ( double, 10, 10 );
		if ( *c == NULL )
		{
			if ( a )
				NDestroyMatrix ( (void **) *a );
				
			if ( b )
				NDestroyMatrix ( (void **) *b );
				
			return ( FALSE );
		}
	}
	
	if ( d )
	{
		*d = NMatrix ( double, 10, 2 );
		if ( *d == NULL )
		{
			if ( a )
				NDestroyMatrix ( (void **) *a );
			
			if ( b )	
				NDestroyMatrix ( (void **) *b );
			
			if ( c )
				NDestroyMatrix ( (void **) *c );
				
			return ( FALSE );
		}
	}
		
	return ( TRUE );
}

/*** A2DeleteAstrometricSolution ****************************************************/

void A2DeleteAstrometricSolution ( double **a, double **b, double **c, double **d )
{
	if ( a )
		NDestroyMatrix ( (void **) a );
	
	if ( b )	
		NDestroyMatrix ( (void **) b );
	
	if ( c )	
		NDestroyMatrix ( (void **) c );
		
	if ( d )
		NDestroyMatrix ( (void **) d );
}

/*** A2InitializeAstrometricSolution ************************************************/

int A2InitializeAstrometricSolution ( double **a, double **b, double **c, double **d )
{
	int i, j;
	
	if ( a )
		for ( i = 0; i < 10; i++ )
			for ( j = 0; j < 10; j++ )
				a[i][j] = 0.0;
		
	if ( b )
		for ( i = 0; i < 10; i++ )
			for ( j = 0; j < 2; j++ )
				b[i][j] = 0.0;

	if ( c )
		for ( i = 0; i < 10; i++ )
			for ( j = 0; j < 10; j++ )
				c[i][j] = 0.0;

	if ( d )
		for ( i = 0; i < 10; i++ )
			for ( j = 0; j < 2; j++ )
				d[i][j] = 0.0;

	return ( TRUE );
}

/*** A2CopyAstrometricSolution **********************************************************/

void A2CopyAstrometricSolution ( double **a, double **b, double **c, double **d )
{
	int i, j;
	
	if ( a && c )
		for ( i = 0; i < 10; i++ )
			for ( j = 0; j < 10; j++ )
				a[i][j] = c[i][j];
				
	if ( b && d )
		for ( i = 0; i < 10; i++ )
			for ( j = 0; j < 2; j++ )
				b[i][j] = d[i][j];
}

/*** A2SetAstrometricSolution ************************************************************/

void A2SetAstrometricSolution ( double x, double y, double a, double s, double **b, double **d )
{
	int i;
	
	if ( b )
	{
		b[0][0] =  x;
		b[0][1] =  y;
		b[1][0] =  cos ( a ) / s;
		b[1][1] =  sin ( a ) / s;
		b[2][0] = -sin ( a ) / s;
		b[2][1] =  cos ( a ) / s;
		
		for ( i = 3; i < 10; i++ )
			b[i][0] = b[i][1] = 0.0;
	}
	
	if ( d )
	{	
		d[0][0] = -s * cos ( a ) * x - s * sin ( a ) * y;
		d[0][1] =  s * sin ( a ) * x - s * cos ( a ) * y;
		d[1][0] =  s * cos ( a );
		d[1][1] =  s * sin ( a );
		d[2][0] = -s * sin ( a );
		d[2][1] =  s * cos ( a );

		for ( i = 3; i < 10; i++ )
			d[i][0] = d[i][1] = 0.0;
	}
}

/*** A2GetAstrometricSolution ****************************************************************/

void A2GetAstrometricSolution ( double **b, double **d, double *x, double *y, double *a, double *s )
{
	if ( b )
	{
		*s = b[1][0] * b[2][1] - b[2][0] * b[1][1];
		if ( *s > 0.0 )
		{
			*s = 1.0 / sqrt ( *s );
			*a = atan2pi ( b[1][1] - b[2][0], b[1][0] + b[2][1] );
			*x = b[0][0];
			*y = b[0][1];
		}
	}

	if ( d )
	{
		*s = d[1][0] * d[2][1] - d[2][0] * d[1][1];
		if ( *s > 0.0 )
		{
			*s = sqrt ( *s );
			*a = atan2pi ( d[2][0] - d[1][1], d[1][0] + d[2][1] );
			*x = ( -cos ( *a ) * d[0][0] + sin ( *a ) * d[0][1] ) / *s;
			*y = ( -sin ( *a ) * d[0][0] - cos ( *a ) * d[0][1] ) / *s;
		}
	}
}

/*** A2AugmentAstrometricSolution *******************************************************/

void A2AugmentAstrometricSolution ( double **a, double **b, double **c, double **d,
short order, double ra0, double dec0, double ra, double dec, double x, double y )
{
	int n = A2ComputeNumTerms ( order );
	double xi, eta, p[10], q[2];
	
	A2RADecToXiEta ( ra, dec, ra0, dec0, &xi, &eta );

	q[0] = x;
	q[1] = y;
	
	A2ComputeXYPolynomial ( xi, eta, p );
	NAugmentNormalEqns ( n, 2, p, q, a, b );

	q[0] = xi;
	q[1] = eta;
	
	A2ComputeXYPolynomial ( x, y, p );
	NAugmentNormalEqns ( n, 2, p, q, c, d );
}

/*** A2FitAstrometricSolution ************************************************************/

int A2FitAstrometricSolution ( double **a, double **b, double **c, double **d, short order )
{
	int result1, result2, n = A2ComputeNumTerms ( order );
	
	result1 = NGaussJordanSolveMatrixEqn ( a, n, b, 2 );
	result2 = NGaussJordanSolveMatrixEqn ( c, n, d, 2 );
	
	return ( result1 && result2 );
}

/*** A2ComputeXYPolynomial ***************************************************************

	Computes polynomial terms up to 3rd order for use in astrometric solutions.
	
	(x,y): coordinate variables.
	(p): recieves polynomial terms.
	
*****************************************************************************************/

void A2ComputeXYPolynomial ( double x, double y, double p[10] )
{
	p[0] = 1.0;
	p[1] = x;
	p[2] = y;
	p[3] = x * x;
	p[4] = y * y;
	p[5] = x * y;
	p[6] = x * x * x;
	p[7] = y * y * y;
	p[8] = x * x * y;
	p[9] = y * y * x;
}

/*** A2ComputeNumTerms ***************************************************************/

short A2ComputeNumTerms ( short order )
{
	short	n = 0;
	
	if ( order == 1 )
		n = 3;
		
	if ( order == 2 )
		n = 6;
		
	if ( order == 3 )
		n = 10;

	return ( n );		
}

/*** A3AngleToScale ********************************************************************/

double A3AngleToScale ( double angle, double width )
{
	double scale = 2.0 * tan ( angle / 2.0 ) / width;
	
	return ( scale );
}

/*** A3ScaleToAngle ********************************************************************/

double A3ScaleToAngle ( double scale, double width )
{
	double angle = 2.0 * atan ( width * scale / 2.0 );
	
	return ( angle );
}


/*** A3ObjectDistance2 ****************************************************************/

double A3ObjectDistance2 ( A3Object *pObject1, A3Object *pObject2 )
{
	double	dx = pObject2->x - pObject1->x;
	double	dy = pObject2->y - pObject1->y;
	
	return ( dx * dx + dy * dy );
}

/*** A3ObjectPixelDistance ****************************************************************/

double A3ObjectPixelDistance ( A3Object *pObject1, A3Object *pObject2 )
{
	return ( sqrt ( A3ObjectDistance2 ( pObject1, pObject2 ) ) );
}

// Comparison function for sorting objects by amplitude.

int A3ObjectCompareAmplitudes ( const void *p1, const void *p2 )
{
	A3Object *pObj1 = (A3Object *) p1;
	A3Object *pObj2 = (A3Object *) p2;
	
	if ( pObj1->peak < pObj2->peak )
		return 1;
	else if ( pObj1->peak > pObj2->peak )
		return -1;
	else
		return 0;
}

// Comparison function for sorting object by total flux.

int A3ObjectCompareFluxes ( const void *p1, const void *p2 )
{
	A3Object *pObj1 = (A3Object *) p1;
	A3Object *pObj2 = (A3Object *) p2;
	
	if ( pObj1->flux < pObj2->flux )
		return 1;
	else if ( pObj1->flux > pObj2->flux )
		return -1;
	else
		return 0;
}

// Given two image objects, the image center position (x0,y0) in pixels,
// and the image scale in radians per pixel at the image center, this function computes the
// angular separation between the two objects, in radians.

double A3ObjectAngularDistance ( A3Object *pObject1, A3Object *pObject2, double x0, double y0, double scale )
{
	double uvw1[3] = { 0 }, uvw2[3] = { 0 };
	
	uvw1[0] = 1.0;
	uvw1[1] = ( pObject1->x - x0 ) * scale;
	uvw1[2] = ( y0 - pObject1->y ) * scale;
	
	uvw2[0] = 1.0;
	uvw2[1] = ( pObject2->x - x0 ) * scale;
	uvw2[2] = ( y0 - pObject2->y ) * scale;
	
	AANormalizeVector ( uvw1 );
	AANormalizeVector ( uvw2 );
	
	return ( AAVectorSeparation ( uvw1, uvw2 ) );
}

/*** A3ImageInit ******************************************************************/

void A3ImageInit ( A3Image *pImage, double lon, double lat, double jd, double ra, double dec, double angle, int width, int height, int depth, void *data )
{
	int row, col;
	
	A3ImageFree ( pImage );
	memset ( pImage, 0, sizeof ( A3Image ) );
	
	// Set time and location. If they are known, calculate
	// J2000 -> current equatorial (precession) matrix and
	// current equatorial -> horizon transformation matrix.
	// Set image center and image transformation matrix.
	
	A3ImageSetLonLatJD ( pImage, lon, lat, jd );
	A3ImageSetCenterRADec ( pImage, ra, dec, 0.0 );

	// Calculate image scale if field-of-view angle is known.
	
	if ( ! isinf ( angle ) && angle > 0.0 )
		pImage->scale = 2.0 * tan ( angle / 2.0 ) / MAX ( width, height );
	else
		pImage->scale = 0.0;
	
	// Allocate image data matrix.  Return on failure.
	
	pImage->pData = (unsigned char **) NCreateMatrix ( sizeof ( unsigned char ), height, width );
	if ( pImage->pData == NULL )
		return;
	
	// Initialize exposure/ISO to "unknown" values
	
	pImage->exposure = pImage->iso = 0;
	
	// Store image data
	
	pImage->width = width;
	pImage->height = height;
	pImage->x0 = pImage->width / 2;
	pImage->y0 = pImage->height / 2;
	
	// Copy image data, converting RGB to monochrome.
	
	for ( row = 0; row < height; row++ )
	{
		unsigned char	*pPixel = (unsigned char *) data + row * width * depth;
		unsigned short	*pPix16 = (unsigned short *) pPixel;
		
		for ( col = 0; col < width; col++ )
		{
			if ( depth == A3IMAGE_MONO8 )
				pImage->pData[row][col] = pPixel[0];
			else if ( depth == A3IMAGE_MONO16 )
				pImage->pData[row][col] = MIN ( pPix16[0], 255 );
			else if ( depth == A3IMAGE_RGB8 || depth == A3IMAGE_RGBA8 )
				pImage->pData[row][col] = pPixel[0] * 0.3 + pPixel[1] * 0.6 + pPixel[2] * 0.1;
			else if ( depth == A3IMAGE_RGB16 || depth == A3IMAGE_RGBA16 )
				pImage->pData[row][col] = MIN ( pPix16[0] * 0.3 + pPix16[1] * 0.6 + pPix16[2] * 0.1, 255 );
			
			pPixel += depth;
			pPix16 += depth;
		}
	}
}

/*** A3ImageFree ******************************************************************/

void A3ImageFree ( A3Image *pImage )
{
    // Delete any previously existing image data, any previously-found objects,
    // and any reference objects, then zero out the whole structure.

    free ( pImage->pObjects );
    pImage->pObjects = NULL;

    free ( pImage->pReferences );
    pImage->pReferences = NULL;

    NDestroyMatrix ( (void **) pImage->pData );
    pImage->pData = NULL;
}

/*** A3ImageSetLonLatJD **********************************************************/

void A3ImageSetLonLatJD ( A3Image *pImage, double lon, double lat, double jd )
{
	pImage->jd = jd;
	pImage->lon = lon;
	pImage->lat = lat;
	pImage->lst = AALocalMeanSiderealTime ( jd, lon );
	
	if ( ! isinf ( lon ) && ! isinf ( lat ) && ! isinf ( jd ) )
	{
		AASetPrecessionRotationMatrix ( pImage->premat, J2000, pImage->jd, false );
		AASetHorizonRotationMatrix ( pImage->hormat, pImage->lst, pImage->lat, 1 );
	}
}

// Initializes image statistics structure prior to statistical analysis

void A3ImageStatsInit ( A3ImageStats *pStats )
{
	memset ( pStats, 0, sizeof ( A3ImageStats ) );
	
	pStats->fMin = HUGE_VAL;
	pStats->fMax = -HUGE_VAL;
}

// Updates image statistics with a pixel value (fVal) located at (iCol,iRow).
// Call this once for every pixel in the pat of the image you want to analyze.

void A3ImageStatsUpdate ( A3ImageStats *pStats, int iCol, int iRow, float fVal )
{
	float	fMean = pStats->fMean;

	// Mean and standard deviation calculation is based on
	// Welford's method, described in Knuth's "Algorithms":
	// http://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
	
	pStats->nPixels++;
	pStats->fMean += ( fVal - fMean ) / pStats->nPixels;
	pStats->fSum2 += ( fVal - pStats->fMean ) * ( fVal - fMean );
	pStats->fSum  += fVal;
	
	if ( fVal < pStats->fMin )
	{
		pStats->fMin = fVal;
		pStats->iMinRow = iRow;
		pStats->iMinCol = iCol;
	}
	
	if ( fVal > pStats->fMax )
	{
		pStats->fMax = fVal;
		pStats->iMaxRow = iRow;
		pStats->iMaxCol = iCol;
	}
	
	if ( fVal >= 0.0 && fVal < 256.0 )
		pStats->nHistogram[ (int) fVal ]++;
}

// Finalizes image statistic computation (mean, sigma, median, mode).
// Call this after updating image statistics on each pixel in a region.

void A3ImageStatsFinal ( A3ImageStats *pStats )
{
	if ( pStats->nPixels > 1 )
		pStats->fStdDev = sqrt ( pStats->fSum2 / ( pStats->nPixels - 1 ) );
	
	int iSum = 0, iMax = 0;
	pStats->fMode = 0;
	
	for ( int iVal = 0; iVal < 256; iVal++ )
	{
		iSum += pStats->nHistogram[iVal];
		if ( iSum < pStats->nPixels / 2 )
			pStats->fMedian = iVal;
		
		if ( pStats->nHistogram[iVal] > iMax )
		{
			iMax = pStats->nHistogram[iVal];
			pStats->fMode = iVal;
		}
	}
}

/*** A3ImageStatsInRectangle *******************************************************/

void A3ImageStatsInRectangle ( A3Image *pSkyImage, A3ImageStats *pStats, int iLeft, int iTop, int iWidth, int iHeight )
{
	int				iCol, iRow;
	int 			iRight = iLeft + iWidth;
	int				iBottom = iTop + iHeight;
	
	if ( iTop < 0 )
		iTop = 0;
	
	if ( iBottom > pSkyImage->height )
		iBottom = pSkyImage->height;
	
	if ( iLeft < 0 )
		iLeft = 0;
	
	if ( iRight > pSkyImage->width )
		iRight = pSkyImage->width;

	A3ImageStatsInit ( pStats );
	
	for ( iRow = iTop; iRow < iBottom; iRow++ )
		for ( iCol = iLeft; iCol < iRight; iCol++ )
			A3ImageStatsUpdate ( pStats, iCol, iRow, pSkyImage->pData[iRow][iCol] );
	
	A3ImageStatsFinal ( pStats );
}

/*** A3ImageStatsInEllipse ***********************************************************/

void A3ImageStatsInEllipse ( A3Image *pSkyImage, A3ImageStats *pStats, int iCol0, int iRow0, int iColRad, int iRowRad )
{
	int		iRow, iCol;
	int 	iLeft = iCol0 - iColRad;
	int 	iTop = iRow0 - iRowRad;
	int 	iRight = iCol0 + iColRad;
	int		iBottom = iRow0 + iRowRad;
	
	if ( iTop < 0 )
		iTop = 0;
	
	if ( iBottom >= pSkyImage->height )
		iBottom = pSkyImage->height - 1;
	
	if ( iLeft < 0 )
		iLeft = 0;
	
	if ( iRight >= pSkyImage->width )
		iRight = pSkyImage->width - 1;
	
	A3ImageStatsInit ( pStats );
	
	for ( iRow = iTop; iRow <= iBottom; iRow++ )
	{
		for ( iCol = iLeft; iCol <= iRight; iCol++ )
		{
			float fDelRow = (float) ( iRow - iRow0 ) / iRowRad;
			float fDelCol = (float) ( iCol - iCol0 ) / iColRad;
			
			if ( fDelRow * fDelRow + fDelCol * fDelCol > 1.0 )
				continue;
			
			A3ImageStatsUpdate ( pStats, iCol, iRow, pSkyImage->pData[iRow][iCol] );
		}
	}
	
	A3ImageStatsFinal ( pStats );
}

/*** A3ImageStatsInPolygon *************************************************************/

void A3ImageStatsInPolygon ( A3Image *pSkyImage, A3ImageStats *pStats, int nPoints, int pPointX[], int pPointY[] )
{
	int		iPixY, iPixX;
	int 	iLeft = 0;
	int 	iTop = pSkyImage->height;
	int		iRight = pSkyImage->width;
	int		iBottom = 0;
	int  	nNodes, *pNodeX, i, j, temp;
	
	pNodeX = (int *) calloc ( nPoints, sizeof ( int ) );
	if ( pNodeX == NULL )
		return;
	
	A3ImageStatsInit ( pStats );
	
	// from http://alienryderflex.com/polygon_fill/
	
	for ( i = 0; i < nPoints; i++ )
	{
		if ( pPointY[i] < iTop )
			iTop = pPointY[i];
		
		if ( pPointY[i] > iBottom )
			iBottom = pPointY[i];
	}
	
	if ( iBottom > pSkyImage->height )
		iBottom = pSkyImage->height;
	
	if ( iTop < 0 )
		iTop = 0;
	
	for ( iPixY = iTop; iPixY < iBottom; iPixY++ )
	{
		//  Build a list of nodes.
		
		nNodes = 0;
		j = nPoints - 1;
		
		for ( i = 0; i < nPoints; i++ )
		{
			if ( ( pPointY[i] < iPixY && pPointY[j] >= iPixY ) || ( pPointY[j] < iPixY && pPointY[i] >= iPixY ) )
				pNodeX[nNodes++] = pPointX[i] + (float) ( iPixY - pPointY[i] ) / ( pPointY[j] - pPointY[i] ) * ( pPointX[j] - pPointX[i] );
			
			j = i;
		}
		
		//  Sort the nodes, via a simple bubble sort.
		
		i = 0;
		while ( i < nNodes - 1 )
		{
			if ( pNodeX[i] > pNodeX[i+1] )
			{
				temp = pNodeX[i];
				pNodeX[i] = pNodeX[i+1];
				pNodeX[i+1] = temp;
				if ( i )
					i--;
			}
			else
			{
				i++;
			}
		}

		// Use the pixels between node pairs.
		
		for ( i = 0; i < nNodes; i += 2 )
		{
			if ( pNodeX[i] >= iRight )
				break;
			
			if ( pNodeX[i+1] >= iLeft )
			{
				if ( pNodeX[i] < iLeft )
					pNodeX[i] = iLeft;
				
		  		if ( pNodeX[i+1] >= iRight )
					pNodeX[i+1] = iRight - 1;
				
		  		for ( iPixX = pNodeX[i]; iPixX < pNodeX[i+1]; iPixX++ )
					A3ImageStatsUpdate ( pStats, iPixX, iPixY, pSkyImage->pData[iPixY][iPixX] );
			}
		}
	}

	free ( pNodeX );
	A3ImageStatsFinal ( pStats );
}

/*** A3ImageXYToUVW ******************************************************************/

void A3ImageXYToUVW ( A3Image *pImage, float x, float y, double uvw[3] )
{
	uvw[0] = 1.0;
	uvw[1] = ( pImage->x0 - x ) * pImage->scale;
	uvw[2] = ( pImage->y0 - y ) * pImage->scale;
	
	AANormalizeVector ( uvw );
}

/*** A3ImageUVWToXY ******************************************************************/

void A3ImageUVWToXY ( A3Image *pImage, double uvw[3], float *x, float *y )
{
	*x = pImage->x0 - ( uvw[1] / uvw[0] ) / pImage->scale;
	*y = pImage->y0 - ( uvw[2] / uvw[0] ) / pImage->scale;
}

/*** A3ImageUVWToRADec ******************************************************************/

void A3ImageUVWToRADec ( A3Image *pImage, double uvw[3], double *ra, double *dec )
{
	AAUnTransformVector ( pImage->imgmat, uvw );
	AAXYZVectorToSpherical ( uvw, ra, dec, NULL );
}

/*** A3ImageRADecToUVW ******************************************************************/

void A3ImageRADecToUVW ( A3Image *pImage, double ra, double dec, double uvw[3] )
{
	AASphericalToXYZVector ( ra, dec, 1.0, uvw );
	AATransformVector ( pImage->imgmat, uvw );
}

/*** A3ImageXYToRADec ******************************************************************/

void A3ImageXYToRADec ( A3Image *pImage, float x, float y, double *ra, double *dec )
{
	double uvw[3] = { 0 };
	
	A3ImageXYToUVW ( pImage, x, y, uvw );
	A3ImageUVWToRADec ( pImage, uvw, ra, dec );
}

/*** A3ImageRADecToXY ********************************************************************/

void A3ImageRADecToXY ( A3Image *pImage, double ra, double dec, float *x, float *y )
{
	double uvw[3] = { 0 };
	
	A3ImageRADecToUVW ( pImage, ra, dec, uvw );
	A3ImageUVWToXY ( pImage, uvw, x, y );
}

/*** A3ImageXYToAzAlt ******************************************************************/

void A3ImageXYToAzAlt ( A3Image *pImage, float x, float y, double *az, double *alt )
{
	double uvw[3] = { 0 };
	
	A3ImageXYToUVW ( pImage, x, y, uvw );
	AAUnTransformVector ( pImage->imgmat, uvw );
	AATransformVector ( pImage->premat, uvw );
	AATransformVector ( pImage->hormat, uvw );
	AAXYZVectorToSpherical ( uvw, az, alt, NULL );
}

/*** A3ImageAzAltToXY ********************************************************************/

void A3ImageAzAltToXY ( A3Image *pImage, double az, double alt, float *x, float *y )
{
	double uvw[3] = { 0 };
	
	AASphericalToXYZVector ( az, alt, 1.0, uvw );
	AAUnTransformVector ( pImage->hormat, uvw );
	AAUnTransformVector ( pImage->premat, uvw );
	AATransformVector ( pImage->imgmat, uvw );
	A3ImageUVWToXY ( pImage, uvw, x, y );
}

/*** A3ImageRADecToAzAlt *****************************************************************/

void A3ImageRADecToAzAlt ( A3Image *pImage, double ra, double dec, double *az, double *alt )
{
	double xyz[3] = { 0 };
	
	AASphericalToXYZVector ( ra, dec, 1.0, xyz );
	AATransformVector ( pImage->premat, xyz );
	AATransformVector ( pImage->hormat, xyz );
	AAXYZVectorToSpherical ( xyz, az, alt, NULL );
}

/*** A3ImageAzAltToRADec *****************************************************************/

void A3ImageAzAltToRADec ( A3Image *pImage, double az, double alt, double *ra, double *dec )
{
	double xyz[3] = { 0 };
	
	AASphericalToXYZVector ( az, alt, 1.0, xyz );
	AAUnTransformVector ( pImage->hormat, xyz );
	AAUnTransformVector ( pImage->premat, xyz );
	AAXYZVectorToSpherical ( xyz, ra, dec, NULL );
}

/*** A3ImageSetCenterRADec ****************************************************************/

void A3ImageSetCenterRADec ( A3Image *pImage, double ra, double dec, double pa )
{
	pImage->ra0 = ra;
	pImage->dec0 = dec;
	pImage->pa = pa;
	
	AASetRotationMatrix ( pImage->imgmat, 3, 2, -pImage->ra0, 1, -pImage->dec0, 0, pImage->pa );
}

/*** A3ImageSetCenterAzAlt ****************************************************************/

void A3ImageSetCenterAzAlt ( A3Image *pImage, double az, double alt, double pa )
{
	double	matrix[3][3] = { 0 };
	
	// Compute matrix for transforming J2000 equatorial coordinates
	// to current local horizon coordinates
	
	AACopyRotationMatrix ( pImage->imgmat, pImage->premat );
	AATransformRotationMatrix ( pImage->hormat, pImage->imgmat );
	
	// note that alt-azimuth is a left-handed coordinate system
	// so negate the middle row of the horizon matrix and use the positive azimuth
	// when computing the view matrix

	AASetRotationMatrix ( matrix, 3, 2, az, 1, -alt, 0, pa );
	AANegateVector ( pImage->imgmat[1] );
	AATransformRotationMatrix ( matrix, pImage->imgmat );

	A3ImageGetCenterFromMatrix ( pImage->imgmat, pImage );
}

/*** A3ImageSetMatrixFromCenter ***********************************************************/

void A3ImageSetMatrixFromCenter ( A3Image *pImage, double matrix[3][3] )
{
	AASetRotationMatrix ( matrix, 3, 2, -pImage->ra0, 1, -pImage->dec0, 0, pImage->pa );
}

/*** A3ImageGetCenterFromMatrix ***********************************************************/

void A3ImageGetCenterFromMatrix ( double matrix[3][3], A3Image *pImage )
{
	AAXYZVectorToSpherical ( matrix[0], &pImage->ra0, &pImage->dec0, NULL );
	pImage->pa = AAVectorPositionAngle ( matrix[0], matrix[2] );
}

// Given an image, this function determines if the pixel at (iCol0,iRow0) is
// >= other pixels within a circle of (iRadius) pixels around (iCol0,iRow).
// Returns true if no pixels within radius are brighter tnan (iRow0,iCol0),
// and (iRow0,iCol0) is brigher than at least one pixel in the circle.
// Returns false if any pixels in circle are brighter than (iRow0,iCol0)
// or if all pixels in circle are equal to (iRow0,iCol0).

int A3ImageLocalMaximum ( A3Image *pImage, int iRow0, int iCol0, int iRadius )
{
	int				iCol = 0, iRow = 0;
	int 			iLeft = iCol0 - iRadius;
	int 			iTop = iRow0 - iRadius;
	int 			iRight = iCol0 + iRadius;
	int				iBottom = iRow0 + iRadius;
	int				iRadius2 = iRadius * iRadius;
	unsigned char	**pData = pImage->pData, iVal0 = pData[iRow0][iCol0];
	int				isMaximum = FALSE;
	
	for ( iRow = iTop; iRow <= iBottom; iRow++ )
	{
		for ( iCol = iLeft; iCol <= iRight; iCol++ )
		{
			int iRowDelta = iRow - iRow0;
			int iColDelta = iCol - iCol0;
			int iDelta2 = iRowDelta * iRowDelta + iColDelta * iColDelta;
			
			if ( iDelta2 > iRadius2 )
				continue;
			
			if ( pData[iRow][iCol] > iVal0 )
				return ( FALSE );
			else if ( pData[iRow][iCol] < iVal0 )
				isMaximum = TRUE;
		}
	}
	
	return ( isMaximum );
}

// returns true if measurement indicates a valid star image, or false otherwise.

void A3ImageMeasureBackground ( A3Image *pImage, int iCol0, int iRow0, A3FindObjectParams *pParams, A3Object *pObject )
{
	int 	iRow, iCol, nBack = 0;
	int		iObjRad2 = pParams->objectRadius * pParams->objectRadius;
	int		iBackRad2 = pParams->backgroundRadius * pParams->backgroundRadius;
	float	fValue = 0, fSumBack2 = 0, fBack = 0, fStdDev = 0;
	unsigned char **pData = pImage->pData;
	
	// First, compute background level (mean) and noise (std dev)
	// in an annular ragion between the object to background radii.
	
	for ( iRow = iRow0 - pParams->backgroundRadius; iRow <= iRow0 + pParams->backgroundRadius; iRow++ )
	{
		for ( iCol = iCol0 - pParams->backgroundRadius; iCol <= iCol0 + pParams->backgroundRadius; iCol++ )
		{
			int	iDeltaRow = iRow - iRow0;
			int iDeltaCol = iCol - iCol0;
			int iDelta2 = iDeltaRow * iDeltaRow + iDeltaCol * iDeltaCol;
			
			// Ignore pixels inside the object radius or outside the background radius.
			
			if ( iDelta2 <= iObjRad2 || iDelta2 > iBackRad2 )
				continue;
			
			fValue = pData[iRow][iCol];
			
			// Background mean and standard deviation calculation is based on
			// Welford's method, described in Knuth's "Algorithms"; see here:
			// http://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
			
			float fPrevBack = fBack;
			
			nBack++;
			fBack += ( fValue - fBack ) / nBack;
			fSumBack2 += ( fValue - fBack ) * ( fValue - fPrevBack );
		}
	}
	
	fStdDev = sqrt ( fSumBack2 / ( nBack - 1 ) );

	pObject->back = fBack;
	pObject->noise = fStdDev;
	pObject->nback = nBack;
}

// Sets all pixel values below a specific altitude to zero.

int A3ImageMaskBelowAltitude ( A3Image *pImage, double dMaskAlt )
{
	int 	iRow, iCol;
	double	dAzm, dAlt;
	
	for ( iRow = 0; iRow < pImage->height; iRow++ )
	{
		for ( iCol = 0; iCol < pImage->width; iCol++ )
		{
			A3ImageXYToAzAlt ( pImage, iCol, iRow, &dAzm, &dAlt );
			if ( dAlt < dMaskAlt )
				pImage->pData[iRow][iCol] = 0;
		}
	}
	
	return ( 0 );
}

// returns true if measurement indicates a valid star image, or false otherwise.

int A3ImageMeasureObject ( A3Image *pImage, int iCol0, int iRow0, A3FindObjectParams *pParams, A3Object *pObject )
{
	int 	iRow, iCol, nBack;
	int		iObjRad2 = pParams->objectRadius * pParams->objectRadius;
	float	fValue = 0, fBack = 0, fNoise = 0;
	unsigned char **pData = pImage->pData;
	
	// First determine local background and noise levels.
	
	A3ImageMeasureBackground ( pImage, iCol0, iRow0, pParams, pObject );
	
	fBack = pObject->back;
	fNoise = pObject->noise;
	nBack = pObject->nback;

	// Now determine object peak, flux, and centroid
	// within a circular region inside the object radius.

	int 	nStar = 0;
	float	fSumStar = 0, fSumStarX = 0, fSumStarY = 0, fMaxStar = 0;

	for ( iRow = iRow0 - pParams->objectRadius; iRow <= iRow0 + pParams->objectRadius; iRow++ )
	{
		for ( iCol = iCol0 - pParams->objectRadius; iCol <= iCol0 + pParams->objectRadius; iCol++ )
		{
			int	iDeltaRow = iRow - iRow0;
			int iDeltaCol = iCol - iCol0;
			int iDelta2 = iDeltaRow * iDeltaRow + iDeltaCol * iDeltaCol;
			
			// Ignore pixels outside object radius

			if ( iDelta2 > iObjRad2 )
				continue;
			
			// Ignore pixels fainter than background
			
			fValue = pData[iRow][iCol] - fBack;
			if ( fValue < 0 )
				continue;
			
			if ( fValue > fMaxStar )
				fMaxStar = fValue;
			
			nStar++;
			fSumStar += fValue;
			fSumStarX += fValue * iCol;
			fSumStarY += fValue * iRow;
		}
	}

	if ( nStar < 1 || fSumStar < 1 )
		return FALSE;
	
	pObject->peak = fMaxStar;
	pObject->flux = fSumStar;
	pObject->nflux = nStar;
	pObject->x = fSumStarX / fSumStar;
	pObject->y = fSumStarY / fSumStar;
	
	// Reject object if its peak or mean flux are below the required
	// number of standard deviations above the background level.

	if ( pObject->peak < pParams->peakSigma * pObject->noise )
		return ( FALSE );
	
	if ( pObject->flux / pObject->nflux < pParams->meanSigma * pObject->noise )
		return ( FALSE );
	
	// Now determine object major and minor axes and orientation angle.
	
	float fSum = 0, fSigXX = 0, fSigYY = 0, fSigXY = 0, fSigR2 = 0;
	
	for ( iRow = iRow0 - pParams->objectRadius; iRow <= iRow0 + pParams->objectRadius; iRow++ )
	{
		for ( iCol = iCol0 - pParams->objectRadius; iCol <= iCol0 + pParams->objectRadius; iCol++ )
		{
			int	iDeltaRow = iRow - iRow0;
			int iDeltaCol = iCol - iCol0;
			int iDelta2 = iDeltaRow * iDeltaRow + iDeltaCol * iDeltaCol;
			
			// Ignore pixels outside object radius
			
			if ( iDelta2 > iObjRad2 )
				continue;
			
			// Ignore pixels fainter than background

			fValue = pData[iRow][iCol] - fBack;
			if ( fValue < 0 )
				continue;
			
			float fDelX = iCol - pObject->x;
			float fDelY = iRow - pObject->y;
			
			fValue *= fValue;
			
			fSum += fValue;
			fSigXX += fValue * fDelX * fDelX;
			fSigYY += fValue * fDelY * fDelY;
			fSigXY += fValue * fDelX * fDelY;
			fSigR2 += fValue * ( fDelX * fDelX + fDelY * fDelY );
		}
	}
	
	fSigXX /= fSum;
	fSigYY /= fSum;
	fSigXY /= fSum;
	fSigR2 /= fSum;
	
	pObject->angle = 0.5 * atan2 ( 2.0 * fSigXY, ( fSigXX - fSigYY ) );
	
	double fCosAng = cos ( pObject->angle );
	double fCosAng2 = fCosAng * fCosAng;
	double fSinAng2 = 1.0 - fCosAng2;
	double fSin2Ang = sin ( 2.0 * pObject->angle );
	
	pObject->maj = 2.0 * sqrt ( fSigXX * fCosAng2 + fSigYY * fSinAng2 + fSigXY * fSin2Ang );
	pObject->min = 2.0 * sqrt ( fSigXX * fCosAng2 + fSigYY * fSinAng2 - fSigXY * fSin2Ang );

	// Reject objects larger than the maximum size
	
	if ( pParams->maxMajorAxis > 0.0 )
		if ( pObject->maj > pParams->maxMajorAxis )
			return ( FALSE );
	
	// Reject objects which are non-circular by more than the maximum allowed ellipticity
	
	if ( pParams->maxAxisRatio > 0.0 )
		if ( pObject->maj / pObject->min > pParams->maxAxisRatio )
			return ( FALSE );

	return ( TRUE );
}

/*** A3ImageFindObjects *********************************************************************/

int A3ImageFindObjects ( A3Image *pImage, A3FindObjectParams *pParams )
{
	int iLeft = MAX ( pParams->left, pParams->backgroundRadius );
	int iTop = MAX ( pParams->top, pParams->backgroundRadius );
	int iRight = MIN ( pParams->left + pParams->width, pImage->width - pParams->backgroundRadius - 1 );
	int iBottom = MIN ( pParams->top + pParams->height, pImage->height - pParams->backgroundRadius - 1 );
	int	iCol = 0, iRow = 0;
	unsigned char **pData = pImage->pData;
	
	for ( iRow = iTop; iRow < iBottom; iRow++ )
	{
		for ( iCol = iLeft; iCol < iRight; iCol++ )
		{
			if ( pData[iRow][iCol] == 0 )
				continue;
			
			if ( A3ImageLocalMaximum ( pImage, iRow, iCol, pParams->objectRadius ) )
			{
				A3Object	object = { 0 };
				
				if ( ! A3ImageMeasureObject ( pImage, iCol, iRow, pParams, &object ) )
					continue;
				
				if ( ! isinf ( pParams->minAlt ) && ! isinf ( pImage->ra0 ) && ! isinf ( pImage->dec0 ) && ! isinf ( pImage->pa ) && ! isinf ( pImage->scale ) )
				{
					double	az = 0, alt = 0;
					
					A3ImageXYToAzAlt ( pImage, object.x, object.y, &az, &alt );
					if ( alt < pParams->minAlt )
						continue;
				}
					
				int index = A3ImageFindObjectIndex ( pImage, iCol, iRow, pParams->objectRadius );
				if ( index >= 0 && index < pImage->nObjects )
				{
					if ( object.flux > pImage->pObjects[index].flux )
						A3ImageAddObject ( pImage, &object, index );
				}
				else
				{
					A3ImageAddObject ( pImage, &object, pImage->nObjects );
				}
			}
		}
	}
	
	qsort ( pImage->pObjects, pImage->nObjects, sizeof ( A3Object ), A3ObjectCompareFluxes );
	return ( pImage->nObjects );
}

/*** A3ImageAddObject ***********************************************************/

int A3ImageAddObject ( A3Image *pImage, A3Object *pObject, int index )
{
	if ( index < 0 || index >= pImage->nObjects )
	{
		A3Object *pNewObjects = (A3Object *) realloc ( pImage->pObjects, sizeof ( A3Object ) * ( pImage->nObjects + 1 ) );
		if ( pNewObjects == NULL )
			return 0;
		
		pNewObjects[ pImage->nObjects ] = *pObject;
		pImage->pObjects = pNewObjects;
		pImage->nObjects++;
	}
	else
	{
		pImage->pObjects[index] = *pObject;
	}
	
	A3ImageXYToUVW ( pImage, pObject->x, pObject->y, pObject->uvw );
	return ( pImage->nObjects );
}

/*** A3ImageAddReference ***********************************************************/

int A3ImageAddReference ( A3Image *pImage, A3Reference *pReference, double minAlt )
{
	double				prexyz[3] = { 0 };
	double				horxyz[3] = { 0 };
	
	AASphericalToXYZVector ( pReference->ra, pReference->dec, 1.0, pReference->xyz );
	
	// If the image's location and time are known, compute the reference's azimuth and altitude
	// from where/when the image was taken.  Use this to eliminate references below the horizon.

	if ( ! isinf ( pImage->lon ) && ! isinf ( pImage->lat ) && ! isinf ( pImage->jd ) && ! isinf ( minAlt ) )
	{
		double azm, alt;
		
		AACopyVector ( prexyz, pReference->xyz );
		AATransformVector ( pImage->premat, prexyz );
		AACopyVector ( horxyz, prexyz );
		AATransformVector ( pImage->hormat, horxyz );
		AAXYZVectorToSpherical ( horxyz, &azm, &alt, NULL );
		
		// Eliminate references below the nminumum altitude when the image was taken.
		
		if ( alt < minAlt )
			return ( -1 );
	}

	// Eliminate references significantly outside the image's field of view.
		
	if ( ! isinf ( pImage->ra0 ) && ! isinf ( pImage->dec0 ) && pImage->scale > 0.0 )
	{
		double sep = AASeparation ( pReference->ra, pReference->dec, pImage->ra0, pImage->dec0 );
		double rad = hypot ( pImage->width, pImage->height ) / 2.0;
		
		rad = atan ( pImage->scale * rad );
		if ( sep > rad )
			return ( -2 );
	}
	
	// Now add the reference to the end of the existing array of references.
	
	A3Reference *pNewRefs = (A3Reference *) realloc ( pImage->pReferences, sizeof ( A3Reference ) * ( pImage->nReferences + 1 ) );
	if ( pNewRefs == NULL )
		return 0;
	
	pNewRefs[ pImage->nReferences ] = *pReference;
	pImage->pReferences = pNewRefs;
	pImage->nReferences++;
	
	return ( pImage->nReferences );
}

/*** A3ImageFindObjectIndex *******************************************************/

int A3ImageFindObjectIndex ( A3Image *pImage, float x, float y, float radius )
{
	for ( int i = 0; i < pImage->nObjects; i++ )
	{
		float dx = pImage->pObjects[i].x - x;
		float dy = pImage->pObjects[i].y - y;
		
		if ( dx * dx + dy * dy < radius * radius )
			return ( i );
	}
	
	return ( -1 );
}

/*** A3ImageFreeObjects *************************************************************/

void A3ImageFreeObjects ( A3Image *pImage )
{
	free ( pImage->pObjects );
	pImage->pObjects = NULL;
	pImage->nObjects = 0;
}

/*** A3ImageFreeReferences *************************************************************/

void A3ImageFreeReferences ( A3Image *pImage )
{
	free ( pImage->pReferences );
	pImage->pReferences = NULL;
	pImage->nReferences = 0;
}

/*** A3ImageObjectXYToUVW ***********************************************************/

void A3ImageObjectXYToUVW ( A3Image *pImage )
{
	int 		i, nObjects = pImage->nObjects;
	A3Object	*pObjects = pImage->pObjects;
	
	for ( i = 0; i < nObjects; i++ )
		A3ImageXYToUVW ( pImage, pObjects[i].x, pObjects[i].y, pObjects[i].uvw );
}

/*** A3VectorDistance2 *******************************************************************/

double A3VectorDistance2 ( double xyz1[3], double xyz2[3] )
{
	double dx = xyz2[0] - xyz1[0];
	double dy = xyz2[1] - xyz1[1];
	double dz = xyz2[2] - xyz1[2];
	
	return ( dx * dx + dy * dy + dz * dz );
}

// Given two star centroid positions in pixels, the image center position (x0,y0) in pixels,
// and the known angular separation between the two stars in radians, this function computes
// the image scale in radians per pixel at the image center.
// Returns HUGE_VAL (infinity) if computation fails to converge on a reasonable value.

double A3FindImageScale ( A3Object *pObject1, A3Object *pObject2, double x0, double y0, double sep )
{
	double sep1, scale = sep / A3ObjectPixelDistance ( pObject1, pObject2 ), ratio = 1.0;
	
	do
	{
		sep1 = A3ObjectAngularDistance ( pObject1, pObject2, x0, y0, scale );
		ratio = sep / sep1;
		scale *= ratio;
	}
	while ( fabs ( ratio - 1.0 ) > 0.000001 && scale < 1.0 );
	
	return ( scale < 1.0 ? scale : HUGE_VAL );
}

/*** A3SolveXYZToUVWMatrix *****************************************************************/

int A3SolveXYZToUVWMatrix ( double xyz1[3], double xyz2[3], double xyz3[3], double uvw1[3], double uvw2[3], double uvw3[3], double matrix[3][3] )
{
	int		result = false;
	double	**b = (double **) NCreateMatrix ( sizeof ( double ), 3, 3 );
	double	**a = (double **) NCreateMatrix ( sizeof ( double ), 3, 3 );
	
	b[0][0] = uvw1[0];	b[0][1] = uvw1[1];	b[0][2] = uvw1[2];
	b[1][0] = uvw2[0];	b[1][1] = uvw2[1];	b[1][2] = uvw2[2];
	b[2][0] = uvw3[0];	b[2][1] = uvw3[1];	b[2][2] = uvw3[2];
	
	a[0][0] = xyz1[0];	a[0][1] = xyz1[1];	a[0][2] = xyz1[2];
	a[1][0] = xyz2[0];	a[1][1] = xyz2[1];	a[1][2] = xyz2[2];
	a[2][0] = xyz3[0];	a[2][1] = xyz3[1];	a[2][2] = xyz3[2];
	
	result = NGaussJordanSolveMatrixEqn ( a, 3, b, 3 );
	if ( result )
	{
		matrix[0][0] = b[0][0];
		matrix[0][1] = b[1][0];
		matrix[0][2] = b[2][0];
		
		matrix[1][0] = b[0][1];
		matrix[1][1] = b[1][1];
		matrix[1][2] = b[2][1];
		
		matrix[2][0] = b[0][2];
		matrix[2][1] = b[1][2];
		matrix[2][2] = b[2][2];
	}
	
	NDestroyMatrix ( (void **) a );
	NDestroyMatrix ( (void **) b );
	
	return ( result );
}

#define DELTA(x,y)	fabs((x-y)/MAX(x,y))
  
/*** A3ImageTripletsMatch *****************************************************************/

int A3ImageTripletsMatch ( A3Image *pImage, int i1, int i2, int i3, int j1, int j2, int j3, double tolerance )
{
	double	di12, di13, di23;
	double	dj12, dj13, dj23;
	double	s12, s13, s23;
	A3Object *pObjects = pImage->pObjects;
	A3Reference	*pRefs = pImage->pReferences;
	
	// If the image has a known scale, we can use that scale to save some computation
	// and prevent false matches.
	
	if ( ! isinf ( pImage->scale ) && pImage->scale > 0 )
	{
		double	tol2 = 0.5 * 0.5;

		di12 = A3VectorDistance2 ( pObjects[i1].uvw, pObjects[i2].uvw );
		dj12 = A3VectorDistance2 (  pRefs[j1].xyz,  pRefs[j2].xyz );
		if ( DELTA ( di12, dj12 ) > tol2 )
			return ( FALSE );
		
		di13 = A3VectorDistance2 ( pObjects[i1].uvw, pObjects[i3].uvw );
		dj13 = A3VectorDistance2 (  pRefs[j1].xyz,  pRefs[j3].xyz );
		if ( DELTA ( di13, dj13 ) > tol2 )
			return ( FALSE );
		
		di23 = A3VectorDistance2 ( pObjects[i2].uvw, pObjects[i3].uvw );
		dj23 = A3VectorDistance2 (  pRefs[j2].xyz,  pRefs[j3].xyz );
		if ( DELTA ( di23, dj23 ) > tol2 )
			return ( FALSE );
	}
	else
	{
		// Angular distance from reference stars 1 to 2 and 1 to 3 must be
		// less than 90 degrees of each other, so the distance between their
		// XYZ vectors must be < sqrt(2).
		
		dj12 = A3VectorDistance2 ( pRefs[j1].xyz,  pRefs[j2].xyz );
		dj13 = A3VectorDistance2 ( pRefs[j1].xyz,  pRefs[j3].xyz );
		dj23 = A3VectorDistance2 ( pRefs[j2].xyz,  pRefs[j3].xyz );
		if ( dj12 > 2.0 || dj13 > 2.0 || dj23 > 2.0 )
			return ( FALSE );
		
		// Compute square of distance from star 1 to 2, 1 to 3, and 2 to 3.
		// Compute approximate image scale (radians/pixel) for each pair.
		// Image scale must be within 10% of identical for each pair.
		
		di12 = A3ObjectDistance2 ( &pObjects[i1], &pObjects[i2] );
		di13 = A3ObjectDistance2 ( &pObjects[i1], &pObjects[i3] );
		di23 = A3ObjectDistance2 ( &pObjects[i2], &pObjects[i3] );
		
		s12 = sqrt ( dj12 / di12 );
		s13 = sqrt ( dj13 / di13 );
		s23 = sqrt ( dj23 / di23 );
		
		if ( DELTA ( s12, s13 ) > 0.1 || DELTA ( s12, s23 ) > 0.1 || DELTA ( s13, s23 ) > 0.1 )
			return ( FALSE );
	}
	
	// Compute precise image scale.
	
	dj12 = 2.0 * asin ( sqrt ( dj12 ) / 2.0 );
	dj13 = 2.0 * asin ( sqrt ( dj13 ) / 2.0 );
	dj23 = 2.0 * asin ( sqrt ( dj23 ) / 2.0 );
	
	s12 = A3FindImageScale ( &pObjects[i1], &pObjects[i2], pImage->x0, pImage->y0, dj12 );
	s13 = A3FindImageScale ( &pObjects[i1], &pObjects[i3], pImage->x0, pImage->y0, dj13 );
	s23 = A3FindImageScale ( &pObjects[i2], &pObjects[i3], pImage->x0, pImage->y0, dj23 );
	
	if ( isinf ( s12 ) || isinf ( s13 ) || isinf ( s23 ) )
		return ( FALSE );
	
	// Precise image scale must be within 1% of identical for each pair.
	
	if ( DELTA ( s12, s13 ) > 0.01 || DELTA ( s12, s23 ) > 0.01 || DELTA ( s13, s23 ) > 0.01 )
		return ( FALSE );
	
	pImage->scale = ( s12 + s13 + s23 ) / 3.0;
	return ( TRUE );
}

/*** A3ImageSolveMatrixFromTriplet *****************************************************************/

int A3ImageSolveMatrixFromTriplet ( A3Image *pImage, int i1, int i2, int i3, int j1, int j2, int j3 )
{
	A3Object	*pObjects = pImage->pObjects;
	A3Reference	*pRefs = pImage->pReferences;
	
	// If the triplets match, use the derived image scale to recompute UVW coordinates of all objects in the image.
	// Then use the triplet's XYZ and UVW coordinates to solve the XYZ -> UVW transformation matrix.
	
	if ( A3ImageTripletsMatch ( pImage, i1, i2, i3, j1, j2, j3, 0.01 ) )
	{
		A3ImageObjectXYToUVW ( pImage );
		A3ImageUnmatchObjects ( pImage );
		
		pObjects[i1].pRefStar = &pRefs[j1];
		pObjects[i2].pRefStar = &pRefs[j2];
		pObjects[i3].pRefStar = &pRefs[j3];
		
		A3SolveXYZToUVWMatrix ( pRefs[j1].xyz, pRefs[j2].xyz, pRefs[j3].xyz, pObjects[i1].uvw, pObjects[i2].uvw, pObjects[i3].uvw, pImage->imgmat );
		return ( TRUE );
	}

	return ( FALSE );
}

/*** A3ImageSolveMatrix ***************************************************/

int A3ImageSolveMatrix ( A3Image *pImage, double matrix[3][3] )
{
	int		success = FALSE;
	double	**a = NULL, **b = NULL;

	if ( NCreateNormalEqns ( 3, 3, &a, &b ) == FALSE )
		return ( FALSE );
	
	for ( int i = 0; i < pImage->nObjects; i++ )
	{
		A3Object	*pObj = &pImage->pObjects[i];
		A3Reference	*pRef = pObj->pRefStar;
		
		if ( pRef == NULL )
			continue;
		
		NAugmentNormalEqns ( 3, 3, pRef->xyz, pObj->uvw, a, b );
	}
	
	success = NGaussJordanSolveMatrixEqn ( a, 3, b, 3 );
	if ( success )
	{
		matrix[0][0] = b[0][0];
		matrix[0][1] = b[1][0];
		matrix[0][2] = b[2][0];
		
		matrix[1][0] = b[0][1];
		matrix[1][1] = b[1][1];
		matrix[1][2] = b[2][1];
		
		matrix[2][0] = b[0][2];
		matrix[2][1] = b[1][2];
		matrix[2][2] = b[2][2];
	}
	
	// delete normal equations
	
	NDestroyMatrix ( (void **) a );
	NDestroyMatrix ( (void **) b );
	
	return ( success );
}

/*** A3ImageUnmatchObjects ************************************************************************/

void A3ImageUnmatchObjects ( A3Image *pImage )
{
	for ( int j = 0; j < pImage->nObjects; j++ )
	{
		pImage->pObjects[j].pRefStar = NULL;
		pImage->pObjects[j].refError = 0.0;
	}
}

/*** A3ImageMatchObjects **************************************************************************/

void A3ImageMatchObjects ( A3Image *pImage, A3MatchParams *pParams, int maxObjects, int maxRefs )
{
	double			uvw[3] = { 0 }, d = 0;
	float			x = 0, y = 0, dx = 0, dy = 0, d2 = 0, tolerance2 = pParams->tolerance * pParams->tolerance;
	float			fSumRefError = 0;
	int				nObjectsMatched = 0;
	int 			nObjects = MIN ( maxObjects, pImage->nObjects ), nRefs = MIN ( maxObjects, pImage->nReferences );
	A3Object		*pObjects = pImage->pObjects;
	A3Reference		*pRefs = pImage->pReferences;
	
	A3ImageUnmatchObjects ( pImage );
	
	// For each object found in the image...
	
	for ( int j = 0; j < nObjects; j++ )
	{
		// Find the reference which is closest to the object.
		
		int min_i = -1;
		double min_d = PI;
		
		for ( int i = 0; i < nRefs; i++ )
		{
			// Figure out where that reference object should be located in the sky image.
			// Reject reference objects that fall within margin pixels of the image boundaries.
			
			AACopyVector ( uvw, pRefs[i].xyz );
			AATransformVector ( pImage->imgmat, uvw );
			A3ImageUVWToXY ( pImage, uvw, &x, &y );
			if ( x < pParams->left || x >= pParams->left + pParams->width || y < pParams->top || y >= pParams->top + pParams->height )
				continue;
			
			dx = pObjects[j].x - x;
			dy = pObjects[j].y - y;
			d2 = dx * dx + dy * dy;
			
			if ( d2 < tolerance2 )
			{
				d = AAVectorSeparation ( uvw, pObjects[j].uvw );
				if ( d < min_d )
				{
					min_i = i;
					min_d = d;
				}
			}
		}

		// if we found a matching reference, link it and store angular error.
		// Increment match count and accumulate total angular error.
		
		if ( min_i >= 0 )
		{
			pObjects[j].pRefStar = &pRefs[min_i];
			pObjects[j].refError = min_d;
			
			fSumRefError += min_d;
			nObjectsMatched++;
		}
	}
	
	// return results
	
	pParams->nObjsMatched = nObjectsMatched;
	pParams->meanError = nObjectsMatched ? fSumRefError / nObjectsMatched : 0.0;
}

#if 1

/*** A3ImageIdentifyObjects **************************************************************************/

int A3ImageIdentifyObjects ( A3Image *pImage, A3MatchParams *pParams )
{
	int				blind = isinf ( pImage->ra0 ) || isinf ( pImage->dec0 );	// if we have no approximate image center we must do a "blind" solution.
	double			scale = 0, bestScale = 0;
	int				i1, i2, i3;
	int				j1, j2, j3;
	double			jd0 = AACurrentUTC();
	int 			nObjects = MIN ( pImage->nObjects, 20 ), nRefs = MIN ( pImage->nReferences, 20 );
	A3Object		*pBestObjects = NULL;
	int				maxObjMatched = 0;
	float			minRefErr = HUGE_VAL;
	
	pBestObjects = (A3Object *) calloc ( sizeof ( A3Object ), pImage->nObjects );

	// If image scale is known beforehand, use it to compute UVW coordinates
	// of all objects found in the image.
	
	if ( pImage->scale > 0.0 )
		A3ImageObjectXYToUVW ( pImage );
	
	for ( i1 = 0; i1 < nObjects; i1++ )
		for ( i2 = i1 + 1; i2 < nObjects; i2++ )
			for ( i3 = i2 + 1; i3 < nObjects; i3++ )
			{
				for ( j1 = 0; j1 < nRefs; j1++ )
					for ( j2 = 0; j2 < nRefs; j2++ )
						for ( j3 = 0; j3 < nRefs; j3++ )
							if ( j1 != j2 && j1 != j3 && j2 != j3 )
							{
								// Save image scale. If triplet matches, scale will be recomputed from the triplet.
								// Then solve image matrix using triplet, and use it to match all objects with references.
								
								scale = pImage->scale;

								if ( A3ImageSolveMatrixFromTriplet ( pImage, i1, i2, i3, j1, j2, j3 ) )
								{
//									printf ( "\n\ntriplet matched: implied image scale = %.1f\"/pix!\n", pImage->scale * ARCSEC_PER_RAD );
//									A3ImagePrintIdentifiedObjects ( pImage, TRUE, stdout );
									
									// If at least 5 objects match (i.e. the triplet plus two others), find a least-squares best-fit image matrix
									// using all matched objects, not just the triplet.  Then use that matrix to make a better match.

									A3ImageMatchObjects ( pImage, pParams, blind ? pImage->nObjects : nObjects, blind ? pImage->nReferences : nRefs );
									if ( pParams->nObjsMatched >= 5 )
									{
//										printf ( "from triplet, %d of brightest %d objects matched, implied image scale = %.1f\"/pix!\n", pParams->nObjsMatched, nObjects, pImage->scale * ARCSEC_PER_RAD );
//										A3ImagePrintIdentifiedObjects ( pImage, TRUE, stdout );

										A3ImageSolveMatrix ( pImage, pImage->imgmat );
										A3ImageMatchObjects ( pImage, pParams, pImage->nObjects, pImage->nReferences );
										
//										printf ( "from above, %d of brightest %d objects matched, implied image scale = %.1f\"/pix!\n", pParams->nObjsMatched, pImage->nObjects, pImage->scale * ARCSEC_PER_RAD );
//										A3ImagePrintIdentifiedObjects ( pImage, TRUE, stdout );

										// If more objects match than the previous best match, or the average error decreases,
										// save the image scale, object array, number of matching objects, and average error.
										
										if ( pParams->nObjsMatched > maxObjMatched || ( pParams->nObjsMatched == maxObjMatched && pParams->meanError < minRefErr ) )
										{
											bestScale = pImage->scale;
											memcpy ( pBestObjects, pImage->pObjects, sizeof ( A3Object ) * pImage->nObjects );
											maxObjMatched = pParams->nObjsMatched;
											minRefErr = pParams->meanError;
										}
									}
									
									// restore previous image scale; if non-zero,
									// restore previous UVW coordinates of all image objects.
									
									pImage->scale = scale;
									if ( pImage->scale > 0.0 )
										A3ImageObjectXYToUVW ( pImage );
								}
							}

				// check for timeout, if we have a non-zero time limit
				
//				if ( pParams->timeout > 0.0 )
//					if ( AACurrentUTC() - jd0 > pParams->timeout / SEC_PER_DAY )
//						goto failure;
			}
	
	// success!
	// extract the image center coordinates from the matrix and declare success.

	if ( bestScale > 0.0 )
	{
		pImage->scale = bestScale;
		
		memcpy ( pImage->pObjects, pBestObjects, sizeof ( A3Object ) * pImage->nObjects );
		free ( pBestObjects );
		
		A3ImageSolveMatrix ( pImage, pImage->imgmat );
		A3ImageGetCenterFromMatrix ( pImage->imgmat, pImage );

		pParams->nObjsMatched = maxObjMatched;
		pParams->meanError = minRefErr;
		return ( TRUE );
	}
	
	// failure: unmatch any image objects before returning

failure:
	free ( pBestObjects );
	A3ImageUnmatchObjects ( pImage );
	return ( FALSE );
}

#else	// EXPERIMENTAL - DOES NOT WORK CORRECTLY 2018-06-12

/*** A3ImagePairMatch *****************************************************************/

int A3ImagePairMatch ( A3Image *pImage, int i1, int i2, int j1, int j2, double tolerance )
{
	double	di12, dj12, s12;
	A3Object *pObjects = pImage->pObjects;
	A3Reference	*pRefs = pImage->pReferences;
	
	// If the image scale is unknown, all pairs match - by definition!
	
	if ( isinf ( pImage->scale ) || pImage->scale == 0 )
	{
		// Angular distance from reference stars 1 to 2 must be less than 90 degrees,
		// so the distance between their XYZ vectors must be < sqrt(2).
		
		dj12 = A3VectorDistance2 ( pRefs[j1].xyz,  pRefs[j2].xyz );
		if ( dj12 < 2.0 )
			dj12 = 2.0 * asin ( sqrt ( dj12 ) / 2.0 );
		else
			return ( FALSE );
		
		// Compute square of distance from star 1 to 2.
		// Compute approximate image scale (radians/pixel) for each pair.
		
		di12 = A3ObjectDistance2 ( &pObjects[i1], &pObjects[i2] );
		s12 = A3FindImageScale ( &pObjects[i1], &pObjects[i2], pImage->x0, pImage->y0, dj12 );
		if ( isinf ( s12 ) )
			return ( FALSE );
	}
	else
	{
		double	tol2 = 0.5 * 0.5;
		
		di12 = A3VectorDistance2 ( pObjects[i1].uvw, pObjects[i2].uvw );
		dj12 = A3VectorDistance2 (  pRefs[j1].xyz,  pRefs[j2].xyz );
		if ( DELTA ( di12, dj12 ) > tol2 )
			return ( FALSE );
		
		// Compute precise image scale.
		
		dj12 = 2.0 * asin ( sqrt ( dj12 ) / 2.0 );
		s12 = A3FindImageScale ( &pObjects[i1], &pObjects[i2], pImage->x0, pImage->y0, dj12 );
		if ( isinf ( s12 ) )
			return ( FALSE );
		
		// Precise image scale must be within 10% of assumed scale.
		
		if ( DELTA ( s12, pImage->scale ) > 0.1 )
			return ( FALSE );
	}
	
	pImage->scale = s12;
	return ( TRUE );
}

/*** A3ImageSolveMatrixFromPair *****************************************************************/

int A3ImageSolveMatrixFromPair ( A3Image *pImage, int i1, int i2, int j1, int j2 )
{
	A3Object	*pObjects = pImage->pObjects;
	A3Reference	*pRefs = pImage->pReferences;
	double		xyz3[3] = { 0 }, uvw3[3] = { 0 };
	
	// If the pair matches, use the derived image scale to recompute UVW coordinates of all objects in the image.
	// Then use the triplet's XYZ and UVW coordinates to solve the XYZ -> UVW transformation matrix.
	
	if ( A3ImagePairMatch ( pImage, i1, i2, j1, j2, 0.01 ) )
	{
		A3ImageObjectXYToUVW ( pImage );
		A3ImageUnmatchObjects ( pImage );
		
		pObjects[i1].pRefStar = &pRefs[j1];
		pObjects[i2].pRefStar = &pRefs[j2];
		
		AACrossProduct ( pObjects[i1].uvw, pObjects[i2].uvw, uvw3 );
		AANormalizeVector ( uvw3 );
		
		AACrossProduct ( pRefs[j2].xyz, pRefs[j1].xyz, xyz3 );
		AANormalizeVector ( xyz3 );
		
		A3SolveXYZToUVWMatrix ( pRefs[j1].xyz, pRefs[j2].xyz, xyz3, pObjects[i1].uvw, pObjects[i2].uvw, uvw3, pImage->imgmat );
		
		double m0 = AAVectorMagnitude ( pImage->imgmat[0] );
		double m1 = AAVectorMagnitude ( pImage->imgmat[1] );
		double m2 = AAVectorMagnitude ( pImage->imgmat[2] );
		
		double d00 = AADotProduct ( pImage->imgmat[0], pImage->imgmat[1] );
		double d02 = AADotProduct ( pImage->imgmat[0], pImage->imgmat[2] );
		double d12 = AADotProduct ( pImage->imgmat[1], pImage->imgmat[2] );
		
		return ( TRUE );
	}
	
	return ( FALSE );
}

/*** A3ImageIdentifyObjects **************************************************************************/

int A3ImageIdentifyObjects ( A3Image *pImage, A3MatchParams *pParams )
{
	double			scale = 0, bestScale = 0;
	int				i1, i2, j1, j2;
	double			jd0 = AACurrentUTC();
	int 			nObjects = MIN ( pImage->nObjects, 20 ), nRefs = pImage->nReferences; // MIN ( pImage->nReferences, 20 );
	A3Object		*pBestObjects = NULL;
	int				maxObjMatched = 0;
	float			minRefErr = HUGE_VAL;
	
	pBestObjects = (A3Object *) calloc ( sizeof ( A3Object ), pImage->nObjects );
	
	// If image scale is known beforehand, use it to compute UVW coordinates
	// of all objects found in the image.
	
	if ( pImage->scale > 0.0 )
		A3ImageObjectXYToUVW ( pImage );
	
	for ( i1 = 0; i1 < nObjects; i1++ )
		for ( i2 = i1 + 1; i2 < nObjects; i2++ )
		{
			for ( j1 = 0; j1 < nRefs; j1++ )
				for ( j2 = 0; j2 < nRefs; j2++ )
					if ( j1 != j2 )
					{
						// Save image scale. If triplet matches, scale will be recomputed from the triplet.
						// Then solve image matrix using triplet, and use it to match all objects with references.
						
						scale = pImage->scale;
						
						if ( A3ImageSolveMatrixFromPair ( pImage, i1, i2, j1, j2 ) )
						{
							printf ( "\n\npair matched: implied image scale = %.1f\"/pix!\n", pImage->scale * ARCSEC_PER_RAD );
							A3ImagePrintIdentifiedObjects ( pImage, TRUE, stdout );
							
							// If at least 5 objects match (i.e. the triplet plus two others), find a least-squares best-fit image matrix
							// using all matched objects, not just the triplet.  Then use that matrix to make a better match.
							
							A3ImageMatchObjects ( pImage, pParams, pImage->nObjects /*20*/, pImage->nReferences /*20*/ );
							if ( pParams->nObjsMatched >= 5 )
							{
//								printf ( "from triplet, %d of brightest 20 objects matched, implied image scale = %.1f\"/pix!\n", pParams->nObjsMatched, pImage->scale * ARCSEC_PER_RAD );
//								A3ImagePrintIdentifiedObjects ( pImage, TRUE, stdout );
								
								A3ImageSolveMatrix ( pImage, pImage->imgmat );
								A3ImageMatchObjects ( pImage, pParams, 50, 50 );
								
//								printf ( "from above, %d of brightest %d objects matched, implied image scale = %.1f\"/pix!\n", pParams->nObjsMatched, pImage->nObjects, pImage->scale * ARCSEC_PER_RAD );
//								A3ImagePrintIdentifiedObjects ( pImage, TRUE, stdout );
								
								// If more objects match than the previous best match, or the average error decreases,
								// save the image scale, object array, number of matching objects, and average error.
								
								if ( pParams->nObjsMatched > maxObjMatched || ( pParams->nObjsMatched == maxObjMatched && pParams->meanError < minRefErr ) )
								{
									bestScale = pImage->scale;
									memcpy ( pBestObjects, pImage->pObjects, sizeof ( A3Object ) * pImage->nObjects );
									maxObjMatched = pParams->nObjsMatched;
									minRefErr = pParams->meanError;
//									printf ( "%d objects matched, %.2f mean err - A BEST MATCH!\n", pParams->nObjsMatched, RAD_TO_DEG ( pParams->meanError ) );
								}
								else
								{
//									printf ( "%d objects matched, %.2f mean err - not best match\n", pParams->nObjsMatched, RAD_TO_DEG ( pParams->meanError ) );
								}
							}
							
							// restore previous image scale; if non-zero,
							// restore previous UVW coordinates of all image objects.
							
							pImage->scale = scale;
							if ( pImage->scale > 0.0 )
								A3ImageObjectXYToUVW ( pImage );
						}
					}
			
			// check for timeout, if we have a non-zero time limit
			
			//				if ( pParams->timeout > 0.0 )
			//					if ( AACurrentUTC() - jd0 > pParams->timeout / SEC_PER_DAY )
			//						goto failure;
		}
	
	// success!
	// extract the image center coordinates from the matrix and declare success.
	
	if ( bestScale > 0.0 )
	{
		pImage->scale = bestScale;
		
		memcpy ( pImage->pObjects, pBestObjects, sizeof ( A3Object ) * pImage->nObjects );
		free ( pBestObjects );
		
		A3ImageSolveMatrix ( pImage, pImage->imgmat );
		A3ImageGetCenterFromMatrix ( pImage->imgmat, pImage );
		
		pParams->nObjsMatched = maxObjMatched;
		pParams->meanError = minRefErr;
		return ( TRUE );
	}
	
	// failure: unmatch any image objects before returning
	
failure:
	free ( pBestObjects );
	A3ImageUnmatchObjects ( pImage );
	return ( FALSE );
}

#endif

void A3ImagePrint ( A3Image *pImage, FILE *pFile )
{
	int		year;
	short	month, hour, min;
	double	day, sec;
	short	rah, ram, decd, decm;
	double	ras, decs;
	char	sign;
	double	azm, alt, xyz0[3] = { 0 };
	
	AADecimalToDegMinSec ( RAD_TO_HOUR ( pImage->ra0 ), &rah, &ram, &ras, &sign );
	AADecimalToDegMinSec ( RAD_TO_DEG ( pImage->dec0 ), &decd, &decm, &decs, &sign );
	
	fprintf ( pFile, "Width:     %d\n", pImage->width );
	fprintf ( pFile, "Height:    %d\n", pImage->height );
	
	if ( pImage->iso > 0 )
		fprintf ( pFile, "ISO:       %d\n", pImage->iso );
	
	if ( ! isinf ( pImage->exposure ) && pImage->exposure > 0.0 )
		fprintf ( pFile, "Exposure:  %.3f\n", pImage->exposure );
	
	if ( ! isinf ( pImage->jd ) )
	{
		AAJDToDateTime ( pImage->jd, 0.0, &year, &month, &day, &hour, &min, &sec, AA_CALENDAR_GREGORIAN );
		if ( sec > 59.5 )
			AAJDToDateTime ( pImage->jd + 0.5 / SEC_PER_DAY, 0.0, &year, &month, &day, &hour, &min, &sec, AA_CALENDAR_GREGORIAN );
		
		fprintf ( pFile, "JD:        %.6f\n", pImage->jd );
		fprintf ( pFile, "Date:      %d/%02hd/%02.0f\n", year, month, floor ( day ) );
		fprintf ( pFile, "Time:      %02hd:%02hd:%02.0f UTC\n", hour, min, sec );
	}
	
	if ( ! isinf ( pImage->lon ) )
		fprintf ( pFile, "Longitude: %+.4f°\n", RAD_TO_DEG ( pImage->lon ) );
	
	if ( ! isinf ( pImage->lat ) )
		fprintf ( pFile, "Latitude:  %+.4f°\n", RAD_TO_DEG ( pImage->lat ) );
	
	if ( ! isinf ( pImage->ra0 ) )
		fprintf ( pFile, "RA:        %02hdh %02hdm %04.1fs\n", rah, ram, ras );
	
	if ( ! isinf ( pImage->dec0 ) )
		fprintf ( pFile, "Dec:       %c%02hd° %02hd' %02.0f\"\n", sign, decd, decm, decs );
	
	if ( ! isinf ( pImage->pa ) )
		fprintf ( pFile, "PA:        %.1f°\n", RAD_TO_DEG ( pImage->pa ) );

	if ( ! isinf ( pImage->scale ) && pImage->scale > 0.0 )
	{
		fprintf ( pFile, "Scale:     %.1f\"/pix\n", RAD_TO_ARCSEC ( pImage->scale ) );
		fprintf ( pFile, "Field:     %.1f° x %.1f°\n",
			RAD_TO_DEG ( A3ScaleToAngle ( pImage->scale, pImage->width ) ),
			RAD_TO_DEG ( A3ScaleToAngle ( pImage->scale, pImage->height ) ) );
	}

	if ( ! isinf ( pImage->dec0 ) && ! isinf ( pImage->jd ) && ! isinf ( pImage->lon ) && ! isinf ( pImage->lat ) )
	{
		AASphericalToXYZVector ( pImage->ra0, pImage->dec0, 1.0, xyz0 );
		AATransformVector ( pImage->premat, xyz0 );
		AATransformVector ( pImage->hormat, xyz0 );
		AAXYZVectorToSpherical ( xyz0, &azm, &alt, NULL );
		
		if ( ! isinf ( azm ) )
			fprintf ( pFile, "Azimuth:   %.4f°\n", RAD_TO_DEG ( azm ) );
		
		if ( ! isinf ( alt ) )
			fprintf ( pFile, "Altitude:  %+.4f°\n", RAD_TO_DEG ( alt ) );
	}
}

void A3ImagePrintObjects ( A3Image *pImage, FILE *pFile )
{
	A3Object *pObjects = pImage->pObjects;
	int nObjects = pImage->nObjects;
	
	fprintf ( pFile, "%3s  %6s  %6s  %4s  %4s  %5s  %5s  %7s  %5s  %5s\n",
			 "Obj", "X", "Y", "Maj", "Min", "Ang", "Peak", "Flux", "Back", "Noise" );
	
	for ( int i = 0; i < nObjects; i++ )
		fprintf ( pFile, "%3d  %6.1f  %6.1f  %4.1f  %4.1f  %5.1f  %5.1f  %7.1f  %5.1f  %5.1f\n",
				 i + 1, pObjects[i].x, pObjects[i].y, pObjects[i].maj, pObjects[i].min, RAD_TO_DEG ( pObjects[i].angle ), pObjects[i].peak, pObjects[i].flux, pObjects[i].back, pObjects[i].noise );
}

void A3ImagePrintReferences ( A3Image *pImage, FILE *pFile )
{
	A3Reference *pRefs = pImage->pReferences;
	int nRefs = pImage->nReferences;
	
	fprintf ( pFile, "%3s  %10s  %9s  ", "Ref", "RA", "Dec" );
	
	if ( ! isinf ( pImage->lon ) && ! isinf ( pImage->lat ) && ! isinf ( pImage->jd ) )
		fprintf ( pFile, "%8s  %8s  ", "Azm", "Alt" );
	
	fprintf ( pFile, "%5s  %s\n", "Mag", "Name" );

	for ( int i = 0; i < nRefs; i++ )
	{
		short	rah, ram, decd, decm;
		double	ras, decs;
		char	sign;
		
		AADecimalToDegMinSec ( RAD_TO_HOUR ( pRefs[i].ra ), &rah, &ram, &ras, &sign );
		AADecimalToDegMinSec ( RAD_TO_DEG ( pRefs[i].dec ), &decd, &decm, &decs, &sign );
		
		fprintf ( pFile, "%3d  %02hd %02hd %04.1f  %c%02hd %02hd %02.0f  ",
			i + 1, rah, ram, ras, sign, decd, decm, decs );

		if ( ! isinf ( pImage->lon ) && ! isinf ( pImage->lat ) && ! isinf ( pImage->jd ) )
		{
			double	xyz[3] = { 0 }, azm = 0, alt = 0;
			
			AACopyVector ( xyz, pRefs[i].xyz );
			AATransformVector ( pImage->premat, xyz );
			AATransformVector ( pImage->hormat, xyz );
			AAXYZVectorToSpherical ( xyz, &azm, &alt, NULL );
			
			fprintf ( pFile, "%8.4f  %+8.4f  ", RAD_TO_DEG ( azm ), RAD_TO_DEG ( alt ) );
		}
		
		fprintf ( pFile, "%5.2f  %s\n", pRefs[i].mag, pRefs[i].name );
	}
}

void A3ImagePrintIdentifiedObjects ( A3Image *pImage, int onlyIdentified, FILE *pFile )
{
	A3Reference *pRef = NULL;
	A3Object *pObjects = pImage->pObjects;
	int nObjects = pImage->nObjects;
	
	fprintf ( pFile, "%3s  %6s  %6s  %4s  %4s  %5s  %7s  %5s  %5s  %3s  %5s  %10s  %9s  %5s  %s\n",
			 "Obj", "X", "Y", "Maj", "Min", "Peak", "Flux", "Back", "Noise", "Ref", "Err", "RA", "Dec", "Mag", "Name" );
	
	for ( int i = 0; i < nObjects; i++ )
	{
		short	rah, ram, decd, decm;
		double	ras, decs;
		char	sign;
		
		pRef = pObjects[i].pRefStar;
		if ( pRef == NULL && onlyIdentified )
			continue;
		
		fprintf ( pFile, "%3d  %6.1f  %6.1f  %4.1f  %4.1f  %5.1f  %7.1f  %5.1f  %5.1f  ",
				 i + 1, pObjects[i].x, pObjects[i].y, pObjects[i].maj, pObjects[i].min, pObjects[i].peak, pObjects[i].flux, pObjects[i].back, pObjects[i].noise );
		
		if ( pRef == NULL )
		{
			fprintf ( pFile, "\n" );
			continue;
		}
		
		AADecimalToDegMinSec ( RAD_TO_HOUR ( pRef->ra ), &rah, &ram, &ras, &sign );
		AADecimalToDegMinSec ( RAD_TO_DEG ( pRef->dec ), &decd, &decm, &decs, &sign );
		
		fprintf ( pFile, "%3ld  %5.2f  %02hd %02hd %04.1f  %c%02hd %02hd %02.0f  %5.2f  %s\n",
				 (long) (pRef - pImage->pReferences + 1), RAD_TO_DEG ( pObjects[i].refError ), rah, ram, ras, sign, decd, decm, decs, pRef->mag, pRef->name );
	}
}

/*** COPYRIGHT NOTICE ********************************************************
 
	Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.
 
	This file contains functions for converting 2D (x,y) positions in an image
 	to celestial coordinates (ra,dec) using known reference objects: a problem
 	traditionally called "plate solving" even though astronomers haven't used
	glass plates for decades!
 
	The "A2" functions implement traditional plate solving and work well at
 	small image scales, when the image only covers a small part of the sky
 	and the image center coordinates are well known.
 
	The "A3" functions and data types are my original implementation.  They
	are intended for large sky images with significant spherical distortion,
	where the field of view is unknown. These functions include methods to
 	match triplets of light sources on the image with known reference objects,
 	test the match, and determine the image scale and center coordinates.
 
******************************************************************************/

#ifndef ASTROMETRY_H
#define ASTROMETRY_H

#ifdef __cplusplus
extern "C" {
#endif

#include "AstroLib.h"

// Represents a star, planet, or other known reference object in the sky

typedef struct A3Reference
{
	double		ra, dec;	// Right Ascension and declination of reference star, radians
	float		mag;		// Magnitude of reference star
	double		xyz[3];		// reference star RA/Dec converted to rectangular XYZ unit vector
	char		name[32];	// ASCII name string, zero-terminated.
	void		*data;		// Pointer to client-specified data (could be pointer to struct w/ additional params, etc.)
}
A3Reference;

// Represents a pointlike object (a light source) observed in a 2-dimensional image

typedef struct A3Object
{
	float		x, y;			// position of object centroid on image
	float		maj, min;		// major and minor axes of object ellipsoid, in pixels
	float		angle;			// orientation angle of object major axis in radians; 0 = parallel to X; HALF_PI = parallel to Y
	float		peak;			// peak value of object, minus background.
	float		back;			// background level
	float		flux;			// total light flux from object, minus background.
	float		noise;			// background noise level
	int			nback, nflux;	// number of pixels used in background and flux measurement.
	double		uvw[3];			// apparent direction to object centroid in image coordinate system.
	A3Reference	*pRefStar;		// pointer to reference object matched with this image object, or NULL if no reference object has been matched with this image object.
	float		refError;		// angular distance in radians from this image object's computed position to reference object's known position; zero if pRefStar is NULL.
}
A3Object;

// Contains information about a 2-dimensional raster image of part of the sky.

typedef struct A3Image
{
	int				width;			// horizontal pixel count, i.e. number of colums
	int				height;			// vertical pixel count, i.e. number of rows
	int				iso;			// ISO rating (dimensionless, but assume 100 = gain of 1.0); 0 if unknown
	float			exposure;		// exposure length in seconds; 0 if unknown
	float			x0;				// horizontal position of image center, pixels
	float			y0;				// vertical position of image center, pixels
	double			scale;			// image scale in radians per pixel at center of image
	double			ra0;			// right ascension at image center, radians
	double			dec0;			// declination at image center, radians
	double			pa;				// position angle from north in sky to vertical on image at image center, radians
	double			jd;				// Julian Date of mid-exposure (i.e. when image was captured); INFINITY if unknown
	double			lon;			// Longitude of location where image was taken, radians (east positive); INFINITY if unknown
	double			lat;			// Latitude of location where image was taken, radians (north positive); INFINITY if unknown
	double			lst;			// Local sidereal time when image was taken, radians; INFINITY if unknown
	double			premat[3][3];	// Rotation matrix for transforming from J2000 equatorial to current equatorial coordinates.
	double			hormat[3][3];	// Rotation matrix for transforming from current equatorial to local horizon coordinates.
	double			imgmat[3][3];	// Rotation matrix for transforming from J2000 equatorial coordinates to image coordinates.
	A3Object		*pObjects;		// pointer to array of A3Object structs describing pointlike objects observed in the image
	int				nObjects;		// number of A3Object structs contained in above array.
	A3Reference		*pReferences;	// pointer to array of A3Reference structs describing known sky objects that may appear in the image.
	int				nReferences;	// number of A3Reference structs contained in above array.
	unsigned char	**pData;		// array of pointers to rows of raster image data
	FILE			*pLogFile;		// pointer to file for logging debug information; can be stdout, stderr, a file, or NULL (for no logging)
}
A3Image;

// Contains statistical information about a rectangular or circular region of an image.

typedef struct A3ImageStats
{
	int		nPixels;				// total number of pixels in region
	float	fSum, fSum2;			// sum of pixel values, sum of squares of pixel values.
	float	fMean;					// average value of pixels in region
	float	fStdDev;				// standard deviation of pixel values in region
	float	fMin;					// minimum pixel value in region
	int		iMinRow, iMinCol;		// location of minimum pixel value in region
	float	fMax;					// maximum pixel value in region
	int		iMaxRow, iMaxCol;		// location of maximum pixel value in region
	float	fMedian;				// median pixel value in region
	float	fMode;					// mode (most common) pixel value in region
	int		nHistogram[256];		// histogram of pixel values in region
}
A3ImageStats;

// Parameters used by A3ImageFindObjects to detect objects in an image
	
typedef struct A3FindObjectParams
{
	int		left, top;			// left, top corner of image rectangle to search, in pixels; objects outside this rectangle are ignored.
	int		width, height;		// width, height of image rectangle to search, in pixels; objects outside this rectangle are ignored.
	float	minAlt;				// minimum altitude above horizon to search, in radians; ignored if minAlt == -INFINITY or image center RA/Dec/PA are INFINITY.
	int		objectRadius;		// radius to use for object flux measurement, in pixels.
	int		backgroundRadius;	// radius to use for background measurement, in pixels.
	float	peakSigma;			// peak pixel value inside object radius must be this many standard deviations above background.
	float	meanSigma;			// mean pixel value inside object radius must be this many standard deviations above background.
	float	maxMajorAxis;		// maximum allowed major axis in pixels; INFINITY finds objects of any size.
	float	maxAxisRatio;		// maximum allowed ratio between major and minor axes; 0 or INFINITY finds objects of any shape.
}
A3FindObjectParams;
	
// Inputs and results of sky image object matching function - see A3ImageMatchObjects() and A3ImageIdentifyObjects()

typedef struct A3MatchParams
{
	int		left, top;		// INPUT: left, top corner of image rectangle to match, in pixels; references outside this rectangle are ignored.
	int		width, height;	// INPUT: width, height of image rectangle to match, in pixels; references outside this rectangle are ignored.
	int		tolerance;		// INPUT: maximum allowed distance from reference object's predicted position (x,y) to image object centroid, in pixels.
	float	timeout;		// INPUT: number of seconds to run computation before giving up if no match is found.  Ignored if zero.
	int		nObjsMatched;	// OUTPUT: number of reference sky objects that were actually matched with an image object.
	float	meanError;		// OUTPUT: average distance from image object to matching reference object, in radians.
}
A3MatchParams;

/*** A2RADecToXiEta ***********************************************************
 
	Converts Right Ascension and Declination to standard coordinates.
	
	void A2RADecToXiEta ( double a, double d, double a0, double d0,
	double *xi, double *eta )
 
	(a,d): R.A. and Dec. coordinates, in radians.
	(a0,d0): R.A. and Dec. of adopted reference point, in radians.
	(xi,eta): receives standard coordinates, in radians.
	
	Reference: Marsden, B.G., "How to Reduce Plate Measurements",
	Sky & Telescope, Sept 1982, p. 284.
	
*******************************************************************************/

void ASTROLIBDLL_API A2RADecToXiEta ( double, double, double, double, double *, double * );

/*** A2XiEtaToRADec ************************************************************
 
	Converts dimensionless "standard coordinates" to R.A./Dec.
	
	void A2XiEtaToRADec ( double xi, double eta, double ra0, double dec0,
	double *ra, double *dec )
 
	(xi,eta): standard coordinates
	(ra0,dec0): adopted R.A./Dec. of image center
	(ra,dec): receives R.A./Dec. corresponding to (xi,eta)
 
	Reference: Marsden, B.G., "How to Reduce Plate Measurements",
	Sky & Telescope, Sept 1982, p. 284.
	
**********************************************************************************/

void ASTROLIBDLL_API A2XiEtaToRADec ( double, double, double, double, double *, double * );

/*** A2XiEtaToXY *****************************************************************
 
	Converts dimensionless sky coordinates (xi,eta) to physical image coordinates.
	
	void A2XiEtaToXY ( double xi, double eta, double **coeffs, double *x, double *y )
 
	(xi,eta): sky coordinates, dimensionless units.
	(coeffs): matrix containing (xi,eta) -> (x,y) transformation coefficients.
	(x,y): receives physical image coordinates corresponding to (xi,eta).
	
***********************************************************************************/

void ASTROLIBDLL_API A2XiEtaToXY ( double, double, double **, double *, double * );

/**** A2XYToXiEta ****************************************************************
 
	Converts physical image coordinates to dimensionless sky coordinates (xi,eta).
	
	void A2XYToXiEta ( double x, double y, double **coeffs, double *xi, double *eta )
 
	(x,y): physical image coordinates.
	(coeffs): matrix containing (x,y) -> (xi,eta) transformation coefficients.
	(xi,eta): receives sky coordinates in dimensionless units.
 
**********************************************************************************/

void ASTROLIBDLL_API A2XYToXiEta ( double, double, double **, double *, double * );

/*** A2XYToRADec *****************************************************************
 
	Converts physical image coordinates to R.A./Dec. coordinates.
	
	void A2XYToRADec ( double x, double y, double **coefs, double ra0, double dec0,
	double *ra, double *dec )
 
	(x,y): physical image coordinates
	(coeffs): matrix containing (x,y) -> (xi,eta) transformation coefficients.
	(ra0,dec0): adopted coordinates of image center.
	(ra,dec): receives R.A./Dec. coordinates corresponding to (x,y).
 
***********************************************************************************/

void ASTROLIBDLL_API A2XYToRADec ( double, double, double **, double, double, double *, double * );

/*** RADecToXY ******************************************************************
 
	Converts R.A./Dec. coordinates to physical image coordinates.
	
	void RADecToXY ( double ra, double dec, double ra0, double dec0,
	double **coeffs, double *x, double *y )
 
	(ra,dec): receives R.A./Dec. coordinates corresponding to (x,y).
	(ra0,dec0): adopted coordinates of image center.
	(coeffs): matrix containing (xi,eta) -> (x,y) transformation coefficients.
	(x,y): receives physical image coordinates.
 
***********************************************************************************/

void ASTROLIBDLL_API A2RADecToXY ( double, double, double, double, double **, double *, double * );

/*** A2NewAstrometricSolution ****************************************************
 
	Allocates and initializes matrices for a new astrometric coefficient solution.
 
	int A2NewAstrometricSolution ( double ***a, double ***b, double ***c, double ***d )
 
	(a): recieves pointer to (xi,eta) -> (x,y) transformation covariance matrix.
	(b): receives pointer to (xi,eta) -> (x,y) transformation coefficient matrix.
	(c): recieves pointer to (x,y) -> (xi,eta) transformation covariance matrix.
	(d): receives pointer to (x,y) -> (xi,eta) transformation coefficient matrix.
	
	The function returns TRUE if successful and FALSE on failure.
	
	If don't want to allocate all of the matrices, pass NULL instead of a
	matrix pointer for any which you do not wish to allocate.
 
**********************************************************************************/

int ASTROLIBDLL_API A2NewAstrometricSolution ( double ***, double ***, double ***, double *** );

/*** A2DeleteAstrometricSolution **************************************************
 
	Frees memory for astrometric coefficient solution matrices.
	
	int A2DeleteAstrometricSolution ( double **a, double **b, double **c, double **d )
	
	(a): (xi,eta) -> (x,y) transformation covariance matrix.
	(b): (xi,eta) -> (x,y) transformation coefficient matrix.
	(c): (x,y) -> (xi,eta) transformation covariance matrix.
	(d): (x,y) -> (xi,eta) transformation coefficient matrix.
 
	If you don't want to delete all of the matrices, pass NULL instead of a
	matrix pointer for any matrix which you don't wish to de-allocate.
	
***********************************************************************************/

void ASTROLIBDLL_API A2DeleteAstrometricSolution ( double **, double **, double **, double ** );

/*** A2CopyAstrometricSolution **************************************************
 
	Copies astrometric coefficients and/or covariances from one matrix to another.
	
	void A2CopyAstrometricSolution ( double **a, double **b, double **c, double **d )
	
	(a): destination covariance matrix.
	(b): destination coefficient matrix.
	(c): source covariance matrix.
	(d): source coefficient matrix.
 
	All the matrices (a, b, c, and d) must have been created and initialized by
	the function NewAstrometricSolution().
	
	If you don't want to copy all of the matrices, pass NULL for any matrix which
	you don't wish to copy.  Note that to copy both sets of coefficient and
	covariance matrices allocated by NewAstrometricSolution, you must call this
	function twice.
	
***********************************************************************************/

void ASTROLIBDLL_API A2CopyAstrometricSolution ( double **, double **, double **, double ** );

/*** A2AugmentAstrometricSolution *************************************************
 
	Augments the normal equation matrices for an astrometric coefficient solution
	by adding position information for a reference star.
	
	int AugmentAstrometricSolution ( double **a, double **b, double **c, double **d,
	short order, double ra, double dec, double ra0, double dec0, double x, double y )
	
	(a): (xi,eta) -> (x,y) transformation covariance matrix.
	(b): (xi,eta) -> (x,y) transformation coefficient matrix.
	(c): (x,y) -> (xi,eta) transformation covariance matrix.
	(d): (x,y) -> (xi,eta) transformation coefficient matrix.
	(order): desired astrometric solution order (1-3).
	(ra0): right ascention of reference point, in radians.
	(dec0): declination of reference point, in radians.
	(ra): known right ascention of reference object, in radians.
	(dec): known declination of reference object, in radians.
	(x): measured X coordinate of reference object on image.
	(y): measured Y coordinate of reference object on image.
 
	The matrices (a, b, c, and d) must have been created and initialized
	by the function NewAstrometricSolution() prior to calling this function.
	
	Once you have created and initialized the astrometric solution matrices,
	you should call this function once for each known reference object whose
	position information you wish to add to the solution.  When you are done
	adding reference objects, call FitAstrometricSolution() to solve for the
	(x,y) -> (xi,eta) and (xi,eta) -> (x,y) transformation coefficients and
	associated covariance matrices.
	
*********************************************************************************/

int ASTROLIBDLL_API A2InitializeAstrometricSolution ( double **, double **, double **, double ** );

void ASTROLIBDLL_API A2AugmentAstrometricSolution ( double **, double **, double **, double **,
												 short, double, double, double, double, double, double );

/*** A2FitAstrometricSolution ****************************************************
 
	Computes astrometric solution coefficients by least-squares best fit.
	
	int A2FitAstrometricSolution ( double **a, double **b, double **c, double **d,
	short order )
	
	(a): (xi,eta) -> (x,y) transformation covariance matrix.
	(b): (xi,eta) -> (x,y) transformation coefficient matrix.
	(c): (x,y) -> (xi,eta) transformation covariance matrix.
	(d): (x,y) -> (xi,eta) transformation coefficient matrix.
	(order): desired astrometric solution order (1-3).
 
	The function returns TRUE if it is able to successfully determine a
	least-squares best-fit, and FALSE otherwise.
	
	The input normal equation matrices (a, b, c, d) must have been created and
	initialized by NewAstrometricSolution() prior to calling this function, and
	information for all desired all reference stars must have been added by calls
	to the function AugmentAstrometricSolution().
	
**********************************************************************************/

int ASTROLIBDLL_API A2FitAstrometricSolution ( double **, double **, double **, double **, short );

/*** A2SetAstrometricSolution *******************************************************
 
	Computes "ideal" coefficients for the (x,y) -> (xi,eta) and (xi,eta) -> (x,y)
	transformations, given an image reference pixel, field rotation, and scale.
	
	void A2SetAstrometricSolution ( double x, double y, double a, double s,
	double **b, double **d )
	
	(x): X coordinate of image reference pixel.
	(y): Y coordinate of image reference pixel.
	(a): field rotation angle, in radians.
	(s): image scale, in radians/pixel.
	(b): matrix to receive (xi,eta) -> (x,y) transformation coefficients.
	(d): matrix to receive (x,y) -> (xi,eta) transformation coefficients.
 
	The matrices (b) and (d) must have been previously created and initialized
	by the function NewAstrometricSolution().
	
	If you do not wish to compute one or the other transformation coefficient
	matrices, pass NULL instead of a pointer to it.
	
**************************************************************************************/

void ASTROLIBDLL_API A2SetAstrometricSolution ( double, double, double, double, double **, double ** );

/*** A2GetAstrometricSolution *********************************************************
 
	Computes field reference pixel, rotation, and scale from astrometric coefficients.
	
	void A2GetAstrometricSolution ( double **b, double **d, double *x, double *y,
	double *a, double *s )
 
	(b): matrix containing (xi,eta) -> (x,y) transformation coefficients.
	(d): matrix containing (x,y) -> (xi,eta) transformation coefficients.
	(x): receieves X coordinate of image reference pixel.
	(y): receieves Y coordinate of image reference pixel.
	(a): receives field rotation angle, in radians.
	(s): recieves image scale, in radians/pixel.
 
	If you do not wish to compute field parameters from one or the other
	transformation coefficient matrices, simply pass NULL instead of a pointer to it.
	
**************************************************************************************/

void ASTROLIBDLL_API A2GetAstrometricSolution ( double **, double **, double *, double *, double *, double * );

/*** A3ImageInit *********************************************************************
	 
	Initializes A3Image structure with time, location, field of view (if known)
 	and raster image data.
 
	void A3ImageInit ( A3Image *pImage, double lon, double lat, double jd,
		 double ra, double dec, double angle,
 		 int width, int height, int depth, void *data )

 	(pImage): pointer to A3Image structure to initialize
 	(lon):    longitude from which image was taken, radians; east positive
 	(lat):    latitude from which image was taken, radians; north positive
 	(jd):     julian date on which image was taken
 	(ra):     right ascension of approximate image center, radians
 	(dec):    declination of approximate image center, radians
 	(angle):  field-of-view angle corresponding to larger image dimension, radians
 	(width):  raster image width, pixels
 	(height): raster image height, pixels
 	(format): raster image data format
 	(data):   pointer to buffer containing raster image data.
 
	If the location and/or time the image was taken are unknown, pass INFINITY
	for lon, lat, and/or jd.
 
 	If the approximate image center and/or field-of-view are unknown, pass INFINITY
 	for ra, dec, and/or angle.  If you provide finite values for the approximate image
 	center, they should be within one field-of-view angle of the true image center.
 	If you provide a field-of-view angle, it should be within 10% of the true value.
 
 	The approximate values you provide will be used to eliminate references
 	that cannot possibly appear in the image, and to reject triangles of objects
	that are the wrong size to correspond to any references.  After identifying
 	objects, your approximate values will be replaced by values calculated from
 	the plate solution. Reasonable approximate values greatly speeds up object
 	identification and helps prevent false matches, but are not strictly needed.

	The image's object and reference arrays are deleted, and array counts set to zero.

	The raster data you provide is copied into the image; you may change or delete
	the original. The raster data pixel format must be one of the constants below.
 
**************************************************************************************/

#define A3IMAGE_MONO8	1	// monochrome, 8 bits per pixel
#define A3IMAGE_MONO16	2	// monochrome, 16 bits per pixel
#define A3IMAGE_RGB8	3	// red-green-blue, 8 bits per color channel
#define A3IMAGE_RGBA8	4	// red-green-blue-alpha, 8 bits per color channel
#define A3IMAGE_RGB16	6	// red-green-blue, 16 bits per color channel
#define A3IMAGE_RGBA16	8	// red-green-blue-alpha, 16 bits per color channel

void ASTROLIBDLL_API A3ImageInit ( A3Image *pImage, double lon, double lat, double jd,
					 double ra, double dec, double angle,
					 int width, int height, int depth, void *data );

/*** A3ImageSetLonLatJD *************************************************************
 
	Sets the location and time at which an image was taken.

	A3ImageSetLonLatJD ( A3Image *pImage, double lon, double lat, double jd )

	(pImage): pointer to A3Image structure to initialize.
	(lon):    longitude from which image was taken, radians; east positive
	(lat):    latitude from which image was taken, radians; north positive
	(jd):     julian date on which image was taken

	This function also computes the image's precession and horizon coordinate
	transformation matrices.
 
*************************************************************************************/
	
void ASTROLIBDLL_API A3ImageSetLonLatJD ( A3Image *pImage, double lon, double lat, double jd );
	
double ASTROLIBDLL_API A3AngleToScale ( double angle, double width );
double ASTROLIBDLL_API A3ScaleToAngle ( double scale, double width );

/*** A3ImageStatsInRectangle ********************************************************

	Analyzes statisics within a rectangular region of an image.
 
	void A3ImageStatsInRectangle ( A3Image *pImage, A3ImageStats *pStats,
         int iLeft, int iTop, int iWidth, int iHeight )

	(pImage): pointer to sky image structure with initialized data array.
 	(pStats): pointer to structure which receives statistical results.
 	(iLeft,iTop): top left corner of rectangle, in pixels, from (0,0).
 	(iWidth,iHeight): size of rectangle in pixels.
 
	Returns statistical analysis in the structure pointed to by pStats.

**************************************************************************************/

void ASTROLIBDLL_API A3ImageStatsInRectangle ( A3Image *pImage, A3ImageStats *pStats,
					 int iLeft, int iTop, int iWidth, int iHeight );

/*** A3ImageStatsInEllipse ***********************************************************
	 
	Analyzes statisics within a elliptical region of an image.
	 
	void A3ImageStatsInCircle ( A3Image *pImage, A3ImageStats *pStats,
	     iCol0, int iRow0, int iColRad, int iRowRad )
	 
	(pImage): pointer to sky image structure with initialized data array.
 	(pStats): pointer to structure which receives statistical results.
 	(iCol0,iRow0): coordinates of central pixel, from (0,0) = (left,top).
 	(iColRad): horizontal radius of ellipse in pixels.
	(iRowRad): vertical radius of circle in pixels.
 
	Returns statistical analysis in the structure pointed to by pStats.
	 
**************************************************************************************/
	
void ASTROLIBDLL_API A3ImageStatsInEllipse ( A3Image *pSkyImage, A3ImageStats *pStats,
                     int iCol0, int iRow0, int iColRad, int iRowRad );

/*** A3ImageStatsInPolygon ***********************************************************
 
	Analyzes statisics within an arbitrary polygonal region of an image.
 
	void A3ImageStatsInPolygon ( A3Image *pImage, A3ImageStats *pStats,
	     int nPoints, int pPointX[], int pPointY[] )
 
	(pImage): pointer to sky image structure with initialized data array.
	(pStats): pointer to structure which receives statistical results.
	(nPoints): number of pointa making up polygon.
	(pPointX): array of horizontal coordinates of polygon points.
	(pPointY): array of vertical coordinates of polygon points.

	Returns statistical analysis in the structure pointed to by pStats.
 
**************************************************************************************/

void ASTROLIBDLL_API A3ImageStatsInPolygon ( A3Image *pSkyImage, A3ImageStats *pStats,
					 int nPoints, int pPointX[], int pPointY[] );

/*** A3ImageXYToUVW ******************************************************************

	Converts 2-dimensional image coordinates to 3-dimensional unit vector in image
 	coordinate system.

	void A3ImageXYToUVW ( A3Image *image, float x, float y, double uvw[3] )

 	(image): contains sky image center coordinates, rotation, and scale.
 	(x,y):   2D input coordinates on sky image, in pixels.
 	(uvw):   receives output 3D unit vector
 
	The u axis runs through the image center perpendicular to the image plane;
	the v axis is parallel to X on the image; the w axis is parallel to Z on the image.
	The image center, where (uvw) = (1,0,0), is located at (image->x0,image->y0);
	The image scale at that point (image->scale) is in radians per pixel.

***************************************************************************************/
	
void ASTROLIBDLL_API A3ImageXYToUVW ( A3Image *image, float x, float y, double uvw[3] );

/*** A3ImageUVWToXY ******************************************************************
 
	Converts 3-dimensional unit vector in image coordinate system to 2-dimensional
 	image coordinates.
 
	void A3ImageUVWToXY ( A3Image *image, double uvw[3], float *x, float *y )
 
	(image): contains sky image center coordinates, rotation, and scale.
	(uvw):   input 3D unit vector in image coordinate system.
	(x,y):   receives 2D coordinates on sky image, in pixels.
 
	The u axis runs through the image center perpendicular to the image plane;
	the v axis is parallel to X on the image; the w axis is parallel to Z on the image.
	The image center, where (uvw) = (1,0,0), is located at (image->x0,image->y0);
	The image scale at that point (image->scale) is in radians per pixel.
 
***************************************************************************************/

void ASTROLIBDLL_API A3ImageUVWToXY ( A3Image *image, double uvw[3], float *x, float *y );

/*** A3ImageUVWToRADec ****************************************************************

	void A3ImageUVWToRADec ( A3Image *pImage, double uvw[3], double *ra, double *dec )

	(image):  contains sky image center coordinates, rotation, and scale.
	(uvw):    input 3D unit vector in image coordinate system.
	(ra,dec): receives J2000 equatorial coordinates, in radians.

	Converts 3-dimensional image vector uvw[] in the reference frame described above
	to spherical celestial coordinates (ra,dec) in radians.
 
	Image matrix must have been previously initialized or solved for.

***************************************************************************************/

void ASTROLIBDLL_API A3ImageUVWToRADec ( A3Image *pImage, double uvw[3], double *ra, double *dec );

/*** A3ImageRADecToUVW ****************************************************************
 
	void A3ImageRADecToUVW ( A3Image *pImage, double ra, double dec, double uvw[3] )
 
	(image):  contains sky image center coordinates, rotation, and scale.
	(ra,dec): input J2000 equatorial coordinates, in radians.
	(uvw):    receives output 3D unit vector
 
	Converts spherical celestial coordinates (ra,dec) in radians
	to 3-dimensional image vector uvw[] in the reference frame described above.
 
	Image matrix must have been previously initialized or solved for.
 
***************************************************************************************/

void ASTROLIBDLL_API A3ImageRADecToUVW ( A3Image *pImage, double ra, double dec, double uvw[3] );

/*** A3ImageXYToRADec ******************************************************************
 
	void A3ImageXYToRADec ( A3Image *pImage, float x, float y, double *ra, double *dec )
 
	(image):  contains sky image center coordinates, rotation, and scale.
	(x,y):    2D input coordinates on sky image, in pixels.
	(ra,dec): receives J2000 equatorial coordinates, in radians.
 
	Converts 2-dimensional image coordinates (x,y) in pixels
	to spherical celestial coordinates (ra,dec) in radians.

	Image matrix must have been previously initialized or solved for.
 
***************************************************************************************/

void ASTROLIBDLL_API A3ImageXYToRADec ( A3Image *pImage, float x, float y, double *ra, double *dec );

/*** A3ImageRADecToXY ******************************************************************
	 
	void A3ImageRADecToXY ( A3Image *pImage, double ra, double dec, float *x, float *y )
 
	(image):  contains sky image center coordinates, rotation, and scale.
	(ra,dec): input J2000 equatorial coordinates, in radians.
	(x,y):    receives 2D coordinates on sky image, in pixels.

	Converts spherical celestial coordinates (ra,dec) in radians
	to 2-dimensional image coordinates (x,y) in pixels.
 
	Image matrix and scale must have been previously initialized or solved for.

***************************************************************************************/

void ASTROLIBDLL_API A3ImageRADecToXY ( A3Image *pImage, double ra, double dec, float *x, float *y );

/*** A3ImageRADecToAzAlt **************************************************************
 
	void A3ImageRADecToAzAlt ( A3Image *pImage, double ra, double dec, float *az, float *alt )
 
	(image):  pointer to sky image with initialized coordinate transforms.
	(ra,dec): input J2000 equatorial coordinates, in radians.
	(az,alt): receives current local horizon coordinates, in radians.
 
	Converts J2000 equatorial coordinates (Right Ascension, Declination)
	to current local horizon coordinates (Azimuth, Altitude) in radians.
 
	Image longitude, latitude, and Julian Date must have been previously initialized.
 
***************************************************************************************/

void ASTROLIBDLL_API A3ImageRADecToAzAlt ( A3Image *pImage, double ra, double dec, double *az, double *alt );

/*** A3ImageAzAltToRADec ***************************************************************
 
	void A3ImageAzAltToRADec ( A3Image *pImage, double az, double alt, double *ra, double *dec )
 
	(image):  pointer to sky image with initialized coordinate transforms.
	(az,alt): input current local horizon coordinates, in radians.
	(ra,dec): receives J2000 equatorial coordinates, in radians.
 
	Converts current local horizon coordinates (Azimuth, Altitude)
	to J2000 equatorial coordinates (Right Ascension, Declination).
 
	Image longitude, latitude, and Julian Date must have been previously initialized.
 
***************************************************************************************/

void ASTROLIBDLL_API A3ImageAzAltToRADec ( A3Image *pImage, double az, double alt, double *ra, double *dec );

/*** A3ImageAzAltToXY ******************************************************************
 
	void A3ImageAzAltToXY ( A3Image *pImage, double az, double alt, float *x, float *y )
 
	(image):  contains sky image center coordinates, rotation, and scale.
	(az,alt): input current local horizon coordinates, in radians.
	(x,y):    receives 2D coordinates on sky image, in pixels.
 
	Converts current local horizon coordinates (Azimuth, Altitude)
	to 2-dimensional image coordinates (x,y) in pixels.
 
	Image longitude, latitude, julian date, matrix, and scale must have been
	previously initialized or solved for.
 
***************************************************************************************/

void ASTROLIBDLL_API A3ImageXYToAzAlt ( A3Image *pImage, float x, float y, double *az, double *alt );

/*** A3ImageRADecToXY ******************************************************************
 
	void A3ImageAzAltToXY ( A3Image *pImage, double az, double alt, float *x, float *y )
 
	(image):  contains sky image center coordinates, rotation, and scale.
	(az,alt): input current local horizon coordinates, in radians.
	(x,y):    receives 2D coordinates on sky image, in pixels.
 
	Converts current local horizon coordinates (Azimuth, Altitude) in radians
	to 2-dimensional image coordinates (x,y) in pixels.
 
	Image matrix and scale must have been previously initialized or solved for.
 
***************************************************************************************/

void ASTROLIBDLL_API A3ImageAzAltToXY ( A3Image *pImage, double az, double alt, float *x, float *y );

/*** A3ImageSetCenterRADec *************************************************************

	Sets an image's center J2000 equatorial coordinates, position angle, and rotation
	matrix.
 
	void A3ImageSetCenterRADec ( A3Image *pImage, double ra, double dec, double pa )

	(ra,dec): J2000 equatorial coordinates of image center, in radians.
	(pa):     position angle of vertical axis relative to north at image center, in radians.

	The function computes the image's J2000-equatorial to image coordinate (UVW)
	transformation matrix (pImage->imgmat).
 
****************************************************************************************/

void ASTROLIBDLL_API A3ImageSetCenterRADec ( A3Image *pImage, double ra, double dec, double pa );

/*** A3ImageSetCenterAzAlt *************************************************************
	 
	Sets an image's center current local horizon coordinates, position angle, and rotation
	matrix.
	 
	void A3ImageSetCenterAzAlt ( A3Image *pImage, double az, double alt, double pa )
	 
	(az,alt): current local horizon coordinates of image center, in radians.
	(pa):     position angle of vertical axis relative to zenith at image center, in radians.
	 
	The function computes the J2000 Right Ascension and Declination corresponding to the
	provided azimuth and altitude, sets the image center coordinates to the result, and
 	recomputes the image's J2000 equatorial to image coordinate (UVW) transformation matrix.
 
****************************************************************************************/

void ASTROLIBDLL_API A3ImageSetCenterAzAlt ( A3Image *pImage, double az, double alt, double pa );

/*** A3ImageSetMatrixFromCenter ********************************************************

	Computes the matrix which converts from celestial XYZ coordinates
	to image UVW coordinates.

	void A3ImageSetMatrixFromCenter ( A3Image *image, double matrix[3][3] )

	(image):  contains sky image center coordinates, rotation.
 	(matrix): receives matrix which converts from XYZ to UVW coordinates.

	The image center right ascension and declination (image->ra0,image->dec0)
 	and position angle of image vertical axis at image center (image->pa0)
 	are all in radians.
 
*****************************************************************************************/
	
void ASTROLIBDLL_API A3ImageSetMatrixFromCenter ( A3Image *image, double matrix[3][3] );

/*** A3ImageGetCenterFromMatrix ********************************************************
 
	Computes image center right and position angle, given the matrix which converts
	from celestial XYZ coordinates to image UVW coordinates
 
	void A3ImageGetCenterFromMatrix ( double matrix[3][3], A3Image *image )
 
	(matrix): matrix which converts from XYZ to UVW coordinates.
	(image):  receives image center coordinates and rotation.
 
	The image center right ascension and declination (image->ra0,image->dec0)
	and position angle of image vertical axis at image center (image->pa0)
	are all returned in radians.
 
*****************************************************************************************/

void ASTROLIBDLL_API A3ImageGetCenterFromMatrix ( double matrix[3][3], A3Image *image );

/*** A3ImageFindObjects *****************************************************************

	Find objects an image, with a tweakable set of object-detection parameters.
 
	int A3ImageFindObjects ( A3Image *pImage, A3FindObjectParams *pParams )
 
	(pImage): pointer to image with initialized raster data array
 	(pParams): pointer to image detection parameter set
 
 	The function returns the number of objects detected in the image.
 
	As currently implemented, this function will not find any objects within
 	pParams->backgroundRadius pixels of the outer margins of the image. No object
 	will be found within pParams->objectRadius pixels of another object's center.
 
 	The objects are added to the image's pObjects array, and the image's
 	nObjects count is incremented by the number of objects added.
 
	Use A3ImageFreeObjects() to release memory for the object array.
 
*****************************************************************************************/

int ASTROLIBDLL_API A3ImageFindObjects ( A3Image *pImage, A3FindObjectParams *pParams );

/*** A3ImageAddObject *******************************************************************

	Adds an object to an A3Image's array of objects at a specific index.
 
	int A3ImageAddObject ( A3Image *pImage, A3Object *pObject, int index )
 
 	(pImage): pointer to sky image to receive object.
	(pObject): pointer to object to add to image.
	(index):   index number at which to add object (from zero to pImage->nObjects).
 
	If the index is >= pImage->nObjects, inserts the new object at the end of array.
 	Otherwise, overwrites existing object at the specified index.
 
	Returns total number of entries in image object array.

***************************************************************************************/
	
int ASTROLIBDLL_API A3ImageAddObject ( A3Image *pImage, A3Object *pObject, int index );

/*** A3ImageAddReference ***************************************************************
	 
	Adds a reference to an A3Image's array of references.
	 
	int A3ImageAddReference ( A3Image *pImage, A3Reference *pReference, double minAlt )
	 
	(pImage):     pointer to sky image to receive reference.
	(pReference): pointer to reference to add to image.
	(minAlt):     minimum altitude above horizon (radians); pass -HUGE_VAL to ignore.
	
	If the index is >= pImage->nReferences, inserts the new reference at the end of array.
 	Otherwise, overwrites existing reference at the specified index.
	
	If the image location and time are known, i.e. pImage->lon, lat, and jd all have finite
 	values, then the reference is not added if it appeared below the minimum altitude
 	where/when the image was taken. In this case the function returns -1.

 	If the image center and scale are known, i.e. pImage->ra0, dec0, and scale all have
 	finite values, then the reference will is not added if it appears farther from the
 	image center than the diagonal angle across the image; the function returns -2.
 
	If the reference was added, the function returns the total number of entries
 	in the image reference array.
	 
*******************************************************************************************/
	
int ASTROLIBDLL_API A3ImageAddReference ( A3Image *pImage, A3Reference *pReference, double minAlt );

/*** A3ImageFindObjectIndex *****************************************************************
	 
	Finds objects near a 2D (x,y) position in an image.
	 
	int A3ImageFindObjectIndex ( A3Image *pImage, float x, float y, float radius )
	 
 	(pImage): pointer to sky image contianing objects.
	(x,y):    2D coordinates in sky image at which to search for object.
	(radius): radius in pixels within which to search for object.
	 
	Returns the index (from 0 to pImage->nbjects - 1 )of the first object within
	radius pixels of position (x,y), or -1 if no objects are found within that radius.

**************************************************************************************/

int ASTROLIBDLL_API A3ImageFindObjectIndex ( A3Image *pImage, float x, float y, float radius );

/*** A3ImageFreeObjects ****************************************************************
 
	Frees memory for a sky image's object array.
 
	void A3ImageFreeObjects ( A3Image *pImage )
 
	(pImage): pointer to sky image containing objects.
 
**************************************************************************************/

void ASTROLIBDLL_API A3ImageFreeObjects ( A3Image *pImage );

/*** A3ImageFreeReferences ***********************************************************
 
	Frees memory for a sky image's reference array.
 
	void A3ImageFreeReferences ( A3Image *pImage )
 
	(pImage): pointer to sky image containing references.
 
**************************************************************************************/

void ASTROLIBDLL_API A3ImageFreeReferences ( A3Image *pImage );

/*** A3ImageObjectXYToUVW ************************************************************
	 
	 Converts all image object coordinates from 2D (x,y) to 3D (u,v,w)
	 using the image center and scale values currently stored in the image.
 
	 void A3ImageObjectXYToUVW ( A3Image *pImage )
	 
	 (pImage): pointer to sky image containing objects.
	 
**************************************************************************************/

void ASTROLIBDLL_API A3ImageObjectXYToUVW ( A3Image *pImage );

/*** A3ImageObjectXYToUVW ************************************************************

	Returns square of distance between two XYZ vectors by pythagorean theorem.

	double A3VectorDistance2 ( double xyz1[3], double xyz2[3] )

	(xyz1,xyz2): two XYZ vectors
 
	Returns square of distance between vector endpoints.

**************************************************************************************/

double ASTROLIBDLL_API A3VectorDistance2 ( double xyz1[3], double xyz2[3] );

/*** A3SolveA3SolveXYZToUVWMatrix *****************************************************

	Given a source triangle defined by the vectors xyz1[3], xyz2[3], and xyz3[3],
	and a destination triangle defined by the vectors uvw1[3], uvw2[3] and uvw3[3],
	this function computes the linear transformation matrix needed to convert
	an arbitrary point (x,y,z) in source space to the corresponding point (u,v,w)
	in destination space.
 
	int A3SolveXYZToUVWMatrix  ( double xyz1[3], double xyz2[3], double xyz3[3],
 	    double uvw1[3], double uvw2[3], double uvw3[3], double matrix[3][3] );

	If the transformation matrix can be successfully computed, the function returns
	TRUE; otherwise, it returns FALSE. On success, the result is stored in matrix[3][3].
	
	To perform the conversion from (x,y,z) to (u,v,w) with the matrix:

	u = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
	v = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
	w = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;

***************************************************************************************/

int ASTROLIBDLL_API A3SolveXYZToUVWMatrix  ( double xyz1[3], double xyz2[3], double xyz3[3],
					double uvw1[3], double uvw2[3], double uvw3[3], double matrix[3][3] );

/*** A3ImageTripletMatch ****************************************************************
 
	Tests whether a triplet of objects in a sky image matches a corresponding triplet
 	of reference stars, and if so, computes the image scale.
 
	double A3ImageTripletsMatch ( A3Image *pA3Image, int i1, int i2, int i3,
	       int j1, int j2, int j3, double tolerance );

	(pImage):   sky image, conaining arrays of objects and references.
	(i1,i2,i3): indices in pImage's object array specifying the object triplet.
	(j1,j2,j3); indices in pImage's reference array specifying the correponding reference triplet.
 	(tolerance): factor by which the legs of the triangles are allowed to differ.

	If the angular distances between the objects matches the angular distances between
 	the corresponding references to within the specified tolerance, the function returns
 	TRUE. It returns FALSE otherwise.
 
	If the function returns TRUE (indicating a successful match), then the image scale will
 	be recomputed to the value (radians per pixel) implied by the triplet.  If the function
 	returns FALSE, the image scale will not be affected.
 
 	Note that a triplet match may be a "false positive" - the triangle of objects may have
 	the same shape as the triangle of reference stars (to within tolerance), but still may be
 	wrong.  Passing a known image scale in pImage beforehand greatly reduces this possibility.
 	Use A3ImageMatchObjects() to attempt a more complete match of references to objects based
 	on a successful triplet match from this function.
 
	// NOTE - TOLERANCE IS CURRENTLY IGNORED; HARD CODED TOLERANCE OF 0.01 (1%) IS USED INSTEAD.

************************************************************************************************/
	
int ASTROLIBDLL_API A3ImageTripletsMatch ( A3Image *pA3Image, int i1, int i2, int i3,
											 int j1, int j2, int j3, double tolerance );

/*** A3ImageSolveMatrixFromTriplet ***************************************************************

	Tests whether a triplet of objects in a sky image matches a corresponding triplet of reference
 	stars, and if so, computes the image scale and solves for the matrix which converts from celestial
 	coordinates (x,y,z) to image coordinates (u,v,w).

	int A3ImageSolveMatrixFromTriplet ( A3Image *pImage, int i1, int i2, int i3,
 	    int j1, int j2, int j3 );

	(pImage):   sky image, conaining arrays of objects and references.
	(i1,i2,i3): indices in pImage's object array specifying the object triplet.
	(j1,j2,j3); indices in pImage's reference array specifying the correponding reference triplet.

	This function wraps A3ImageTripletsMatch() with a hard-coded 1% tolerance. In addition, if the
 	triplets match, this function also solves the 3x3 transformation matrix which converts (x,y,z)
 	from RA/Dec celestial coordinate space to (u,v,w) in image coordinate space.
 
***************************************************************************************************/
	
int ASTROLIBDLL_API A3ImageSolveMatrixFromTriplet ( A3Image *pImage, int i1, int i2, int i3,
												   int j1, int j2, int j3 );

/*** A3ImageMatchObjects ************************************************************************

	Matches objects in an image with known reference sky objects using the image's matrix and scale.

	void A3ImageMatchObjects ( A3Image *pImage, A3MatchParams *pParams, int maxObjects, int maxRefs );

	(pImage): pointer to sky image containing arrays of objects and references.
 	(pParams): pointer to object-matching parameter set.
 
	For each reference that falls within the image (excluding a margin around its boundaries),
	a match is made with the closest object (if any) within tolerance (in radians) of the reference's
	celestial coordinates.
 
	The tolerance and margin are input in the A3MatchParams structure.  On return, objects matched
 	to a reference object will have their reference pointer set to that reference. Objects not
 	matched with any reference will have their reference pointer set to NULL.
 
	Statistical results are returned in the A3MatchParams structure.

***************************************************************************************************/

void ASTROLIBDLL_API A3ImageMatchObjects ( A3Image *pImage, A3MatchParams *pParams, int maxObjects, int maxRefs );

/*** A3ImageUnmatchObjects **********************************************************************
 
	Nullifies any existing matches between image objects and reference objects.
 
	void A3ImageUnmatchObjects ( A3Image *pImage );
 
*************************************************************************************************/

void ASTROLIBDLL_API A3ImageUnmatchObjects ( A3Image *pImage );
	
/*** A3ImageIdentifyObjects *************************************************************************

	Identifies objects in a sky image using its array of reference objects. If successful,
 	computes the image coordinate transformation matrix, scale, and center coordinates.
 
	void A3ImageIdentifyObjects ( A3Image *pImage, A3MatchParams *pParams );
 
	(pImage): pointer to sky image containing arrays of objects and references.
 
	The function returns TRUE if it can successfully match enough references with image objects,
	where "enough" is the value specifed in pParams->nRefMatchMin.
 
 	This is a high-level wrapper around A3ImageMatchObjects(). It works by comparing triplets of objects
 	with triplets of references.  When it finds a matching triplet, it uses the triplet to identiffy
 	additional references match other objects.
 
	The image's object and refernce arrays should both be sorted brightest to faintest.
	Additional inputs and outputs in the A3MatchParams structure are described above.
 
	Objects which have been matched to a known reference object will have their reference pointer
	set to that reference.  Objects not identified will have their reference pointer set to NULL.
 
	If more than half the references match an object, the identification is considered successful.
	If successful, this function also determines the image scale, center RA/Dec, position angle,
	and it solves the 3x3 transformation matrix which converts (x,y,z) from RA/Dec celestial
	coordinate space to (u,v,w) in image coordinate space.

*****************************************************************************************************/

int ASTROLIBDLL_API A3ImageIdentifyObjects ( A3Image *pImage, A3MatchParams *pParams );

void ASTROLIBDLL_API A3ImagePrint ( A3Image *pImage, FILE *pFile );
void ASTROLIBDLL_API A3ImagePrintObjects ( A3Image *pImage, FILE *pFile );
void ASTROLIBDLL_API A3ImagePrintReferences ( A3Image *pImage, FILE *pFile );
void ASTROLIBDLL_API A3ImagePrintIdentifiedObjects ( A3Image *pImage, int onlyIdentified, FILE *pFile );
int ASTROLIBDLL_API A3ImageSolveMatrix ( A3Image *pImage, double matrix[3][3] );

int ASTROLIBDLL_API A3ImageMaskBelowAltitude ( A3Image *pImage, double dMaskAlt );
int ASTROLIBDLL_API A3ImageMeasureObject ( A3Image *pImage, int iCol0, int iRow0, A3FindObjectParams *pParams, A3Object *pObject );
void ASTROLIBDLL_API A3ImageMeasureBackground ( A3Image *pImage, int iCol0, int iRow0, A3FindObjectParams *pParams, A3Object *pObject );

#ifdef __cplusplus
}
#endif

#endif	// ASTROMETRY_H

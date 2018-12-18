/*** COPYRIGHT NOTICE *********************************************************
 
	Copyright (c) 2016 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

#ifndef BRIGHT_STARS_H
#define BRIGHT_STARS_H

// star type codes - make sure these don't conflict with solar system
// and deep sky object type codes in SolarSystem.h and DeepSkyObjects.h!

#define STAR					20	// Single star
#define VARIABLE_STAR			21	// Variable star (has entry in General Catalog of Variable Stars)
#define DOUBLE_STAR				22	// Double star (has entry in Washington Double Star Catalog)
#define VARIABLE_DOUBLE_STAR	23	// Variable and double star; has entries in both GCVS and WDS.

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BrightStar
{
	char		type;		// type code: STAR ... VARIABLE_DOUBLE_STAR
	float		ra;			// J2000 Right Ascension in decimal hours
	float		dec;		// J2000 Declination in decimal degrees
	float		mag;		// Visual magnitude
	float		sep;		// Separation in arcseconds for double stars.
	float		dist;		// Distance in parsecs; zero if unknown.
	short		hr;			// Harvard Revised catalog number
	const char *bayer;		// Pointer to Bayer designation as ASCII NUL-terminated name string
	const char *name;		// Pointer to common name as ASCII NUL-terminated name string
	const char *spectype;	// Pointer to spectral type as ASCII NUL-terminated string.
}
BrightStar;

extern BrightStar	gBrightStars[];		// pointer to global array of BrightStar records
extern int			gNumBrightStars;	// number of BrightStar records in global array

#ifdef __cplusplus
}
#endif

#endif // BRIGHT_STARS_H

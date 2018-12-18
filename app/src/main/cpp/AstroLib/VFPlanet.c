/*** COPYRIGHT NOTICE ********************************************************
 
	Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

#include "AstroLib.h"

#define SIN 1
#define COS 2

#define MAX_NUM_ARGS 8

struct VFPTerm
{
  float amp;
  signed char tcoeff;
  signed char argfunc;
  signed char argcoeffs[MAX_NUM_ARGS];
};

static double SumVFPTerms ( double[], double, VFPTerm[], short, short );

/*******************************  SumVFPTerms  ***************************************/

double SumVFPTerms ( double args[], double t, VFPTerm terms[], short first, short last )
{
  double term, sum = 0.0;
  short k;
  char i;
  
  for ( k = first; k <= last; k++ )
  {
    term = 0.0;

    for ( i = 0; i < MAX_NUM_ARGS; i++ )
      if ( terms[k].argcoeffs[i] )
        term += terms[k].argcoeffs[i] * args[i];

    if ( terms[k].argfunc == SIN )
      term = sin ( term );
    else
      term = cos ( term );
 
    for ( i = 0; i < terms[k].tcoeff; i++ )
      term *= t;
 
    sum += term * terms[k].amp;
  }
  
  return ( sum );
}

/******************************  VFPSun  **********************************/

void VFPSun ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  {
	{ 0.606434, 0.03660110129 },
	{ 0.347343,-0.00014709391 },
	{ 0.779072, 0.00273790931 },
	{ 0.993126, 0.00273777850 },
	{ 0.140023, 0.00445036173 },
	{ 0.053856, 0.00145561327 },
	{ 0.056531, 0.00023080893 },
	{ 0.000000, 0.00000000000 }
  };
  
  static VFPTerm terms[16] = {
  {  6910, 0, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {    72, 0, SIN, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {   -17, 1, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {    -7, 0, COS, { 0, 0, 0, 1, 0, 0,-1, 0 } },
  {     6, 0, SIN, { 1, 0,-1, 0, 0, 0, 0, 0 } },
  {     5, 0, SIN, { 0, 0, 0, 4, 0,-8, 3, 0 } },
  {    -5, 0, COS, { 0, 0, 0, 2,-2, 0, 0, 0 } },
  {    -4, 0, SIN, { 0, 0, 0, 1,-1, 0, 0, 0 } },
  {     4, 0, COS, { 0, 0, 0, 4, 0,-8, 3, 0 } },
  {     3, 0, SIN, { 0, 0, 0, 2,-2, 0, 0, 0 } },
  {    -3, 0, SIN, { 0, 0, 0, 0, 0, 0, 1, 0 } },
  {    -3, 0, SIN, { 0, 0, 0, 2, 0, 0,-2, 0 } },
  
  {     0, 0, SIN, { 0, 0, 0, 0, 0, 0, 0, 0 } },

  { 100014, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {  -1675, 0, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {    -14, 0, COS, { 0, 0, 0, 2, 0, 0, 0, 0 } } };

  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[2] + SumVFPTerms ( args, t, terms, 0, 11 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 12, 12 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 13, 15 ) / 1.0e5;
}

/******************************  VFPMoon  **********************************/
	 
void VFPMoon ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] =
  {
    { 0.606434, 0.03660110129 },
    { 0.374897, 0.03629164709 },
    { 0.259091, 0.03674819520 },
    { 0.827362, 0.03386319198 },
    { 0.347343,-0.00014709391 },
    { 0.779072, 0.00273790931 },
    { 0.993126, 0.00273777850 },
    { 0.505498, 0.00445046867 }
  };

  static VFPTerm terms[123] = {
  { 22640, 0, SIN, { 0, 1, 0, 0, 0, 0, 0, 0 } },
  { -4586, 0, SIN, { 0, 1, 0,-2, 0, 0, 0, 0 } },
  {  2370, 0, SIN, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {   769, 0, SIN, { 0, 2, 0, 0, 0, 0, 0, 0 } },
  {  -668, 0, SIN, { 0, 0, 0, 0, 0, 0, 1, 0 } },
  {  -412, 0, SIN, { 0, 0, 2, 0, 0, 0, 0, 0 } },
  {  -212, 0, SIN, { 0, 2, 0,-2, 0, 0, 0, 0 } },
  {  -206, 0, SIN, { 0, 1, 0,-2, 0, 0, 1, 0 } },
  {   192, 0, SIN, { 0, 1, 0, 2, 0, 0, 0, 0 } },
  {   165, 0, SIN, { 0, 0, 0, 2, 0, 0,-1, 0 } },
  {   148, 0, SIN, { 0, 1, 0, 0, 0, 0,-1, 0 } },
  {  -125, 0, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {  -110, 0, SIN, { 0, 1, 0, 0, 0, 0, 1, 0 } },
  {   -55, 0, SIN, { 0, 0, 2,-2, 0, 0, 0, 0 } },
  {   -45, 0, SIN, { 0, 1, 2, 0, 0, 0, 0, 0 } },
  {    40, 0, SIN, { 0, 1,-2, 0, 0, 0, 0, 0 } },
  {   -38, 0, SIN, { 0, 1, 0,-4, 0, 0, 0, 0 } },
  {    36, 0, SIN, { 0, 3, 0, 0, 0, 0, 0, 0 } },
  {   -31, 0, SIN, { 0, 2, 0,-4, 0, 0, 0, 0 } },
  {    28, 0, SIN, { 0, 1, 0,-2, 0, 0,-1, 0 } },
  {   -24, 0, SIN, { 0, 0, 0, 2, 0, 0, 1, 0 } },
  {    19, 0, SIN, { 0, 1, 0,-1, 0, 0, 0, 0 } },
  {    18, 0, SIN, { 0, 0, 0, 1, 0, 0, 1, 0 } },
  {    15, 0, SIN, { 0, 1, 0, 2, 0, 0,-1, 0 } },
  {    14, 0, SIN, { 0, 2, 0, 2, 0, 0, 0, 0 } },
  {    14, 0, SIN, { 0, 0, 0, 4, 0, 0, 0, 0 } },
  {   -13, 0, SIN, { 0, 3, 0,-2, 0, 0, 0, 0 } },
  {   -11, 0, SIN, { 0, 1, 0, 0, 0,16, 0,-18} },
  {    10, 0, SIN, { 0, 2, 0, 0, 0, 0,-1, 0 } },
  {     9, 0, SIN, { 0, 1,-2,-2, 0, 0, 0, 0 } },
  {     9, 0, COS, { 0, 1, 0, 0, 0,16, 0,-18} },
  {    -9, 0, SIN, { 0, 2, 0,-2, 0, 0, 1, 0 } },
  {    -8, 0, SIN, { 0, 1, 0, 1, 0, 0, 0, 0 } },
  {     8, 0, SIN, { 0, 0, 0, 2, 0, 0,-2, 0 } },
  {    -8, 0, SIN, { 0, 2, 0, 0, 0, 0, 1, 0 } },
  {    -7, 0, SIN, { 0, 0, 0, 0, 0, 0, 2, 0 } },
  {    -7, 0, SIN, { 0, 1, 0,-2, 0, 0, 2, 0 } },
  {     7, 0, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    -6, 0, SIN, { 0, 1,-2, 2, 0, 0, 0, 0 } },
  {    -6, 0, SIN, { 0, 0, 2, 2, 0, 0, 0, 0 } },
  {    -4, 0, SIN, { 0, 1, 0,-4, 0, 0, 1, 0 } },
  {     4, 1, COS, { 0, 1, 0, 0, 0,16, 0,-18} },
  {    -4, 0, SIN, { 0, 2, 2, 0, 0, 0, 0, 0 } },
  {     4, 1, SIN, { 0, 1, 0, 0, 0,16, 0,-18} },
  {     3, 0, SIN, { 0, 1, 0,-3, 0, 0, 0, 0 } },
  {    -3, 0, SIN, { 0, 1, 0, 2, 0, 0, 1, 0 } },
  {    -3, 0, SIN, { 0, 2, 0,-4, 0, 0, 1, 0 } },
  {     3, 0, SIN, { 0, 1, 0, 0, 0, 0,-2, 0 } },
  {     3, 0, SIN, { 0, 1, 0,-2, 0, 0,-2, 0 } },
  {    -2, 0, SIN, { 0, 2, 0,-2, 0, 0,-1, 0 } },
  {    -2, 0, SIN, { 0, 0, 2,-2, 0, 0, 1, 0 } },
  {     2, 0, SIN, { 0, 1, 0, 4, 0, 0, 0, 0 } },
  {     2, 0, SIN, { 0, 4, 0, 0, 0, 0, 0, 0 } },
  {     2, 0, SIN, { 0, 0, 0, 4, 0, 0,-1, 0 } },
  {     2, 0, SIN, { 0, 2, 0,-1, 0, 0, 0, 0 } },
	     
  { 18461, 0, SIN, { 0, 0, 1, 0, 0, 0, 0, 0 } },
  {  1010, 0, SIN, { 0, 1, 1, 0, 0, 0, 0, 0 } },
  {  1000, 0, SIN, { 0, 1,-1, 0, 0, 0, 0, 0 } },
  {  -624, 0, SIN, { 0, 0, 1,-2, 0, 0, 0, 0 } },
  {  -199, 0, SIN, { 0, 1,-1,-2, 0, 0, 0, 0 } },
  {  -167, 0, SIN, { 0, 1, 1,-2, 0, 0, 0, 0 } },
  {   117, 0, SIN, { 0, 0, 1, 2, 0, 0, 0, 0 } },
  {    62, 0, SIN, { 0, 2, 1, 0, 0, 0, 0, 0 } },
  {    33, 0, SIN, { 0, 1,-1, 2, 0, 0, 0, 0 } },
  {    32, 0, SIN, { 0, 2,-1, 0, 0, 0, 0, 0 } },
  {   -30, 0, SIN, { 0, 0, 1,-2, 0, 0, 1, 0 } },
  {   -16, 0, SIN, { 0, 2, 1,-2, 0, 0, 0, 0 } },
  {    15, 0, SIN, { 0, 1, 1, 2, 0, 0, 0, 0 } },
  {    12, 0, SIN, { 0, 0, 1,-2, 0, 0,-1, 0 } },
  {    -9, 0, SIN, { 0, 1,-1,-2, 0, 0, 1, 0 } },
  {    -8, 0, SIN, { 0, 0, 1, 0, 1, 0, 0, 0 } },
  {     8, 0, SIN, { 0, 0, 1, 2, 0, 0,-1, 0 } },
  {    -7, 0, SIN, { 0, 1, 1,-2, 0, 0, 1, 0 } },
  {     7, 0, SIN, { 0, 1, 1, 0, 0, 0,-1, 0 } },
  {    -7, 0, SIN, { 0, 1, 1,-4, 0, 0, 0, 0 } },
  {    -6, 0, SIN, { 0, 0, 1, 0, 0, 0, 1, 0 } },
  {    -6, 0, SIN, { 0, 0, 3, 0, 0, 0, 0, 0 } },
  {     6, 0, SIN, { 0, 1,-1, 0, 0, 0,-1, 0 } },
  {    -5, 0, SIN, { 0, 0, 1, 1, 0, 0, 0, 0 } },
  {    -5, 0, SIN, { 0, 1, 1, 0, 0, 0, 1, 0 } },
  {    -5, 0, SIN, { 0, 1,-1, 0, 0, 0, 1, 0 } },
  {     5, 0, SIN, { 0, 0, 1, 0, 0, 0,-1, 0 } },
  {     5, 0, SIN, { 0, 0, 1,-1, 0, 0, 0, 0 } },
  {     4, 0, SIN, { 0, 3, 1, 0, 0, 0, 0, 0 } },
  {    -4, 0, SIN, { 0, 0, 1,-4, 0, 0, 0, 0 } },
  {    -3, 0, SIN, { 0, 1,-1,-4, 0, 0, 0, 0 } },
  {     3, 0, SIN, { 0, 1,-3, 0, 0, 0, 0, 0 } },
  {    -2, 0, SIN, { 0, 2,-1,-4, 0, 0, 0, 0 } },
  {    -2, 0, SIN, { 0, 0, 3,-2, 0, 0, 0, 0 } },
  {     2, 0, SIN, { 0, 2,-1, 2, 0, 0, 0, 0 } },
  {     2, 0, SIN, { 0, 1,-1, 2, 0, 0,-1, 0 } },
  {     2, 0, SIN, { 0, 2,-1,-2, 0, 0, 0, 0 } },
  {     2, 0, SIN, { 0, 3,-1, 0, 0, 0, 0, 0 } },

  {6036298, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {-327746, 0, COS, { 0, 1, 0, 0, 0, 0, 0, 0 } },
  { -57994, 0, COS, { 0, 1, 0,-2, 0, 0, 0, 0 } },
  { -46357, 0, COS, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {  -8904, 0, COS, { 0, 2, 0, 0, 0, 0, 0, 0 } },
  {   3865, 0, COS, { 0, 2, 0,-2, 0, 0, 0, 0 } },
  {  -3237, 0, COS, { 0, 0, 0, 2, 0, 0,-1, 0 } },
  {  -2688, 0, COS, { 0, 1, 0, 2, 0, 0, 0, 0 } },
  {  -2358, 0, COS, { 0, 1, 0,-2, 0, 0, 1, 0 } },
  {  -2030, 0, COS, { 0, 1, 0, 0, 0, 0,-1, 0 } },
  {   1719, 0, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {   1671, 0, COS, { 0, 1, 0, 0, 0, 0, 1, 0 } },
  {   1247, 0, COS, { 0, 1,-2, 0, 0, 0, 0, 0 } },
  {    704, 0, COS, { 0, 0, 0, 0, 0, 0, 1, 0 } },
  {    529, 0, COS, { 0, 0, 0, 2, 0, 0, 1, 0 } },
  {   -524, 0, COS, { 0, 1, 0,-4, 0, 0, 0, 0 } },
  {    398, 0, COS, { 0, 1, 0,-2, 0, 0,-1, 0 } },
  {   -366, 0, COS, { 0, 3, 0, 0, 0, 0, 0, 0 } },
  {   -295, 0, COS, { 0, 2, 0,-4, 0, 0, 0, 0 } },
  {   -263, 0, COS, { 0, 0, 0, 1, 0, 0, 1, 0 } },
  {    249, 0, COS, { 0, 3, 0,-2, 0, 0, 0, 0 } },
  {   -221, 0, COS, { 0, 1, 0, 2, 0, 0,-1, 0 } },
  {    185, 0, COS, { 0, 0, 2,-2, 0, 0, 0, 0 } },
  {   -161, 0, COS, { 0, 0, 0, 2, 0, 0,-2, 0 } },
  {    147, 0, COS, { 0, 1, 2,-2, 0, 0, 0, 0 } },
  {   -142, 0, COS, { 0, 0, 0, 4, 0, 0, 0, 0 } },
  {    139, 0, COS, { 0, 2, 0,-2, 0, 0, 1, 0 } },
  {   -118, 0, COS, { 0, 1, 0,-4, 0, 0, 1, 0 } },
  {   -116, 0, COS, { 0, 2, 0, 2, 0, 0, 0, 0 } },
  {   -110, 0, COS, { 0, 2, 0, 0, 0, 0,-1, 0 } } };

  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[0] + SumVFPTerms ( args, t, terms, 0, 54 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 55, 92 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 93, 122 ) / 1.0e5;  
}

/******************************  VFPMercury  **********************************/

void VFPMercury ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;
  
  static double arguments[MAX_NUM_ARGS][2] = 
  {
    { 0.606434, 0.03660110129 },
    { 0.347343,-0.00014709391 },
    { 0.779072, 0.00273790931 },
    { 0.993126, 0.00273777850 },
    { 0.700695, 0.01136771400 },
    { 0.485541, 0.01136759566 },
    { 0.566441, 0.01136762384 },
    { 0.140023, 0.00445036173 },
  };
    
  static VFPTerm terms[31] = {
  { 84378, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  { 10733, 0, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {  1892, 0, SIN, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {  -646, 0, SIN, { 0, 0, 0, 0, 0, 0, 2, 0 } },
  {   381, 0, SIN, { 0, 0, 0, 0, 0, 4, 0, 0 } },
  {  -306, 0, SIN, { 0, 0, 0, 0, 0, 1,-2, 0 } },
  {  -274, 0, SIN, { 0, 0, 0, 0, 0, 1, 2, 0 } },
  {   -92, 0, SIN, { 0, 0, 0, 0, 0, 2, 2, 0 } },
  {    83, 0, SIN, { 0, 0, 0, 0, 0, 5, 0, 0 } },
  {   -28, 0, SIN, { 0, 0, 0, 0, 0, 3, 2, 0 } },
  {    25, 0, SIN, { 0, 0, 0, 0, 0, 2,-2, 0 } },
  {    19, 0, SIN, { 0, 0, 0, 0, 0, 6, 0, 0 } },
  {    -9, 0, SIN, { 0, 0, 0, 0, 0, 4, 2, 0 } },
  {     8, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {     7, 0, COS, { 0, 0, 0, 0, 0, 2, 0,-5 } },

  { 24134, 0, SIN, { 0, 0, 0, 0, 0, 0, 1, 0 } },
  {  5180, 0, SIN, { 0, 0, 0, 0, 0, 1,-1, 0 } },
  {  4910, 0, SIN, { 0, 0, 0, 0, 0, 1, 1, 0 } },
  {  1124, 0, SIN, { 0, 0, 0, 0, 0, 2, 1, 0 } },
  {   271, 0, SIN, { 0, 0, 0, 0, 0, 3, 1, 0 } },
  {   132, 0, SIN, { 0, 0, 0, 0, 0, 2,-1, 0 } },
  {    67, 0, SIN, { 0, 0, 0, 0, 0, 4, 1, 0 } },
  {    18, 0, SIN, { 0, 0, 0, 0, 0, 3,-1, 0 } },
  {    17, 0, SIN, { 0, 0, 0, 0, 0, 5, 1, 0 } },
  {   -10, 0, SIN, { 0, 0, 0, 0, 0, 0, 3, 0 } },
  {    -9, 0, SIN, { 0, 0, 0, 0, 0, 1,-3, 0 } },
  
  {   39528, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {   -7834, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    -795, 0, COS, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {    -121, 0, COS, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {     -22, 0, COS, { 0, 0, 0, 0, 0, 4, 0, 0 } } };
  
  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[4] + SumVFPTerms ( args, t, terms, 0, 14 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 15, 25 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 26, 30 ) / 1.0e5;
}

/******************************  VFPVenus  **********************************/

void VFPVenus ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  { 0.606434, 0.03660110129,
    0.347343,-0.00014709391,
    0.779072, 0.00273790931,
    0.993126, 0.00273777850,
    0.505498, 0.00445046867,
    0.140023, 0.00445036173,
    0.292498, 0.00445040017,
    0.053856, 0.00145561327,
  };
      	
  static VFPTerm terms[11] = {
  {  2814, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {  -181, 0, SIN, { 0, 0, 0, 0, 0, 0, 2, 0 } },
  {   -20, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    12, 0, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {   -10, 0, COS, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {     7, 0, COS, { 0, 0, 0, 3, 0,-3, 0, 0 } },
  
  { 12215, 0, SIN, { 0, 0, 0, 0, 0, 0, 1, 0 } },
  {    83, 0, SIN, { 0, 0, 0, 0, 0, 1, 1, 0 } },
  {    83, 0, SIN, { 0, 0, 0, 0, 0, 1,-1, 0 } },
  
  {  72335, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {   -493, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } } };
  
  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[4] + SumVFPTerms ( args, t, terms, 0, 5 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 6, 8 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 9, 10 ) / 1.0e5;
}


/*****************************  VFPEarth  **********************************/

void VFPEarth ( double jd, double *lon, double *lat, double *rp )
{
	VFPSun ( jd, lon, lat, rp );

	*lon += PI;
	if ( *lon > TWO_PI )
		*lon -= TWO_PI;

	*lat = - *lat;
}

/******************************  VFPMars  **********************************/

void VFPMars ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  { 0.347343,-0.00014709391,
    0.779072, 0.00273790931,
    0.993126, 0.00273777850,
    0.140023, 0.00445036173,
    0.987353, 0.00145575328,
    0.053856, 0.00145561327,
    0.849694, 0.00145569465,
    0.056531, 0.00023080893
  };
  
  static VFPTerm terms[31] = {
  { 38451, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {  2238, 0, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {   181, 0, SIN, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {   -52, 0, SIN, { 0, 0, 0, 0, 0, 0, 2, 0 } },
  {    37, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   -22, 0, COS, { 0, 0, 0, 0, 0, 1, 0,-2 } },
  {   -19, 0, SIN, { 0, 0, 0, 0, 0, 1, 0,-1 } },
  {    17, 0, COS, { 0, 0, 0, 0, 0, 1, 0,-1 } },
  {    17, 0, SIN, { 0, 0, 0, 0, 0, 4, 0, 0 } },
  {   -16, 0, COS, { 0, 0, 0, 0, 0, 2, 0,-2 } },
  {    13, 0, COS, { 0, 0, 1, 0, 0,-2, 0, 0 } },
  {   -10, 0, SIN, { 0, 0, 0, 0, 0, 1,-2, 0 } },
  {   -10, 0, SIN, { 0, 0, 0, 0, 0, 1, 2, 0 } },
  {     7, 0, COS, { 0, 0, 1, 0, 0,-1, 0, 0 } },
  {    -7, 0, COS, { 0, 0, 2, 0, 0,-3, 0, 0 } },
  {    -5, 0, SIN, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {    -5, 0, SIN, { 0, 0, 1, 0, 0,-1, 0, 0 } },
  {    -5, 0, SIN, { 0, 0, 1, 0, 0,-2, 0, 0 } },
  {    -4, 0, COS, { 0, 0, 2, 0, 0,-4, 0, 0 } },
  {     4, 1, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {     4, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 1 } },
  {     3, 0, COS, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {     3, 0, SIN, { 0, 0, 0, 0, 0, 2, 0,-2 } },
  
  {  6603, 0, SIN, { 0, 0, 0, 0, 0, 0, 1, 0 } },
  {   622, 0, SIN, { 0, 0, 0, 0, 0, 1,-1, 0 } },
  {   615, 0, SIN, { 0, 0, 0, 0, 0, 1, 1, 0 } },
  {    64, 0, SIN, { 0, 0, 0, 0, 0, 2, 1, 0 } },
  
  { 153031, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  { -14170, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   -660, 0, COS, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {    -47, 0, COS, { 0, 0, 0, 0, 0, 3, 0, 0 } } };

  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[4] + SumVFPTerms ( args, t, terms, 0, 22 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 23, 26 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 27, 30 ) / 1.0e5;
}

/******************************  VFPJupiter  **********************************/

void VFPJupiter ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  { 0.347343,-0.00014709391,
    0.779072, 0.00273790931,
    0.993126, 0.00273777850,
    0.089608, 0.00023080893,
    0.056531, 0.00023080893,
    0.882987, 0.00009294371,
    0.400589, 0.00003269438,
    0.000000, 0.00000000000 };
  
  static VFPTerm terms[84] = {
  { 19934, 0, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {  5023, 1, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {  2511, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {  1093, 0, COS, { 0, 0, 0, 0, 2,-5, 0, 0 } },
  {   601, 0, SIN, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {  -479, 0, SIN, { 0, 0, 0, 0, 2,-5, 0, 0 } },
  {  -185, 0, SIN, { 0, 0, 0, 0, 2,-2, 0, 0 } },
  {   137, 0, SIN, { 0, 0, 0, 0, 3,-5, 0, 0 } },
  {  -131, 0, SIN, { 0, 0, 0, 0, 1,-2, 0, 0 } },
  {    79, 0, COS, { 0, 0, 0, 0, 1,-1, 0, 0 } },
  {   -76, 0, COS, { 0, 0, 0, 0, 2,-2, 0, 0 } },
  {   -74, 1, COS, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    68, 1, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    66, 0, COS, { 0, 0, 0, 0, 2,-3, 0, 0 } },
  {    63, 0, COS, { 0, 0, 0, 0, 3,-5, 0, 0 } },
  {    53, 0, COS, { 0, 0, 0, 0, 1,-5, 0, 0 } },
  {    49, 0, SIN, { 0, 0, 0, 0, 2,-3, 0, 0 } },
  {   -43, 1, SIN, { 0, 0, 0, 0, 2,-5, 0, 0 } },
  {   -37, 0, COS, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    25, 0, SIN, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {    25, 0, SIN, { 0, 0, 0, 0, 3, 0, 0, 0 } },
  {   -23, 0, SIN, { 0, 0, 0, 0, 1,-5, 0, 0 } },
  {   -19, 1, COS, { 0, 0, 0, 0, 2,-5, 0, 0 } },
  {    17, 0, COS, { 0, 0, 0, 0, 2,-4, 0, 0 } },
  {    17, 0, COS, { 0, 0, 0, 0, 3,-3, 0, 0 } },
  {   -14, 0, SIN, { 0, 0, 0, 0, 1,-1, 0, 0 } },
  {   -13, 0, SIN, { 0, 0, 0, 0, 3,-4, 0, 0 } },
  {    -9, 0, COS, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {     9, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    -9, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    -9, 0, SIN, { 0, 0, 0, 0, 3,-2, 0, 0 } },
  {     9, 0, SIN, { 0, 0, 0, 0, 4,-5, 0, 0 } },
  {     9, 0, SIN, { 0, 0, 0, 0, 2,-6, 3, 0 } },
  {    -8, 0, COS, { 0, 0, 0, 0, 4,-10, 0, 0 } },
  {     7, 0, COS, { 0, 0, 0, 0, 3,-4, 0, 0 } },
  {    -7, 0, COS, { 0, 0, 0, 0, 1,-3, 0, 0 } },
  {    -7, 0, SIN, { 0, 0, 0, 0, 4,-10, 0, 0 } },
  {    -7, 0, SIN, { 0, 0, 0, 0, 1,-3, 0, 0 } },
  {     6, 0, COS, { 0, 0, 0, 0, 4,-5, 0, 0 } },
  {    -6, 0, SIN, { 0, 0, 0, 0, 3,-3, 0, 0 } },
  {     5, 0, COS, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {    -4, 0, SIN, { 0, 0, 0, 0, 4,-4, 0, 0 } },
  {    -4, 0, COS, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {     4, 0, COS, { 0, 0, 0, 0, 2,-1, 0, 0 } },
  {    -4, 0, COS, { 0, 0, 0, 0, 3,-2, 0, 0 } },
  {    -4, 1, COS, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {     3, 1, SIN, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {     3, 0, COS, { 0, 0, 0, 0, 0, 5, 0, 0 } },
  {     3, 0, COS, { 0, 0, 0, 0, 5,-10, 0, 0 } },
  {     3, 0, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {    -2, 0, SIN, { 0, 0, 0, 2,-1, 0, 0, 0 } },
  {     2, 0, SIN, { 0, 0, 0, 2, 1, 0, 0, 0 } },
  {    -2, 1, SIN, { 0, 0, 0, 0, 3,-5, 0, 0 } },
  {    -2, 1, SIN, { 0, 0, 0, 0, 1,-5, 0, 0 } },
  
  { -4692, 0, COS, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {   259, 0, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {   227, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {  -227, 0, COS, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {    30, 1, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    21, 1, COS, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    16, 0, SIN, { 0, 0, 0, 0, 3,-5, 0, 0 } },
  {   -13, 0, SIN, { 0, 0, 0, 0, 1,-5, 0, 0 } },
  {   -12, 0, COS, { 0, 0, 0, 0, 3, 0, 0, 0 } },
  {    12, 0, SIN, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {     7, 0, COS, { 0, 0, 0, 0, 3,-5, 0, 0 } },
  {    -5, 0, COS, { 0, 0, 0, 0, 1,-5, 0, 0 } },
  	     
  { 520883, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  { -25122, 0, COS, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {   -604, 0, COS, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {    260, 0, COS, { 0, 0, 0, 0, 2,-2, 0, 0 } },
  {   -170, 0, COS, { 0, 0, 0, 0, 3,-5, 0, 0 } },
  {   -106, 0, SIN, { 0, 0, 0, 0, 2,-2, 0, 0 } },
  {    -91, 1, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    -84, 1, COS, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {     69, 0, SIN, { 0, 0, 0, 0, 2,-3, 0, 0 } },
  {    -67, 0, SIN, { 0, 0, 0, 0, 1,-5, 0, 0 } },
  {     66, 0, SIN, { 0, 0, 0, 0, 3,-5, 0, 0 } },
  {     63, 0, SIN, { 0, 0, 0, 0, 1,-1, 0, 0 } },
  {    -51, 0, COS, { 0, 0, 0, 0, 2,-3, 0, 0 } },
  {    -46, 0, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {    -29, 0, COS, { 0, 0, 0, 0, 1,-5, 0, 0 } },
  {     27, 0, COS, { 0, 0, 0, 0, 1,-2, 0, 0 } },
  {    -22, 0, COS, { 0, 0, 0, 0, 3, 0, 0, 0 } },
  {    -21, 0, SIN, { 0, 0, 0, 0, 2,-5, 0, 0 } } };

  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[3] + SumVFPTerms ( args, t, terms, 0, 53 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 54, 65 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 66, 83 ) / 1.0e5;
}

/******************************  VFPSaturn  **********************************/

void VFPSaturn ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  {
    0.347343,-0.00014709391,
    0.779072, 0.00273790931,
    0.993126, 0.00273777850,
    0.056531, 0.00023080893,
    0.133295, 0.00009294371,
    0.882987, 0.00009294371,
    0.821218, 0.00009294371,
    0.400589, 0.00003269438
  };

  static VFPTerm terms[131] = {
  { 23045, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {  5014, 1, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  { -2689, 0, COS, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {  2507, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {  1177, 0, SIN, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {  -826, 0, COS, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {   802, 0, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {   425, 0, SIN, { 0, 0, 0, 1, 0,-2, 0, 0 } },
  {  -229, 1, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {  -153, 0, COS, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {  -142, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {  -114, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   101, 1, SIN, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {   -70, 0, COS, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {    67, 0, SIN, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {    66, 0, SIN, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {    60, 1, COS, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {    41, 0, SIN, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {    39, 0, SIN, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {    31, 0, SIN, { 0, 0, 0, 1, 0,-1, 0, 0 } },
  {    31, 0, SIN, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {   -29, 0, COS, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {   -28, 0, SIN, { 0, 0, 0, 2, 0,-6, 0, 3 } },
  {    28, 0, COS, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {    22, 1, SIN, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {   -22, 0, SIN, { 0, 0, 0, 0, 0, 1, 0,-3 } },
  {    20, 0, SIN, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {    20, 0, COS, { 0, 0, 0, 4, 0,-10,0, 0 } },
  {    19, 0, COS, { 0, 0, 0, 0, 0, 2, 0,-3 } },
  {    19, 0, SIN, { 0, 0, 0, 4, 0,-10,0, 0 } },
  {   -17, 1, COS, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {   -16, 0, COS, { 0, 0, 0, 0, 0, 1, 0,-3 } },
  {   -12, 0, SIN, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {    12, 0, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {   -12, 0, SIN, { 0, 0, 0, 0, 0, 2, 0,-2 } },
  {   -11, 1, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {   -11, 0, COS, { 0, 0, 0, 2, 0,-7, 0, 0 } },
  {    10, 0, SIN, { 0, 0, 0, 0, 0, 2, 0,-3 } },
  {    10, 0, COS, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {     9, 0, SIN, { 0, 0, 0, 4, 0,-9, 0, 0 } },
  {    -8, 0, SIN, { 0, 0, 0, 0, 0, 1, 0,-2 } },
  {    -8, 0, COS, { 0, 0, 0, 0, 2, 1, 0, 0 } },
  {     8, 0, COS, { 0, 0, 0, 0, 2,-1, 0, 0 } },
  {     8, 0, COS, { 0, 0, 0, 0, 0, 1, 0,-1 } },
  {    -8, 0, SIN, { 0, 0, 0, 0, 2,-1, 0, 0 } },
  {     7, 0, SIN, { 0, 0, 0, 0, 2, 1, 0, 0 } },
  {    -7, 0, COS, { 0, 0, 0, 1, 0,-2, 0, 0 } },
  {    -7, 0, COS, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {    -6, 1, SIN, { 0, 0, 0, 4, 0,-10,0, 0 } },
  {     6, 1, COS, { 0, 0, 0, 4, 0,-10,0, 0 } },
  {     6, 1, SIN, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {    -5, 0, SIN, { 0, 0, 0, 3, 0,-7, 0, 0 } },
  {    -5, 0, COS, { 0, 0, 0, 3, 0,-3, 0, 0 } },
  {    -5, 0, COS, { 0, 0, 0, 0, 0, 2, 0,-2 } },
  {     5, 0, SIN, { 0, 0, 0, 3, 0,-4, 0, 0 } },
  {     5, 0, SIN, { 0, 0, 0, 2, 0,-7, 0, 0 } },
  {     4, 0, SIN, { 0, 0, 0, 3, 0,-3, 0, 0 } },
  {     4, 0, SIN, { 0, 0, 0, 3, 0,-5, 0, 0 } },
  {     4, 1, COS, { 0, 0, 0, 1, 0,-2, 0, 0 } },
  {     3, 1, COS, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {     3, 0, COS, { 0, 0, 0, 2, 0,-6, 0, 3 } },
  {    -3, 1, SIN, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {     3, 1, COS, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {    -3, 1, COS, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {     3, 0, COS, { 0, 0, 0, 3, 0,-7, 0, 0 } },
  {     3, 0, COS, { 0, 0, 0, 4, 0,-9, 0, 0 } },
  {     3, 0, SIN, { 0, 0, 0, 3, 0,-6, 0, 0 } },
  {     3, 0, SIN, { 0, 0, 0, 2, 0,-1, 0, 0 } },
  {     3, 0, SIN, { 0, 0, 0, 1, 0,-4, 0, 0 } },
  {     2, 0, COS, { 0, 0, 0, 0, 0, 3, 0,-3 } },
  {     2, 1, SIN, { 0, 0, 0, 1, 0,-2, 0, 0 } },
  {     2, 0, SIN, { 0, 0, 0, 0, 0, 4, 0, 0 } },
  {    -2, 0, COS, { 0, 0, 0, 3, 0,-4, 0, 0 } },
  {    -2, 0, COS, { 0, 0, 0, 2, 0,-1, 0, 0 } },
  {    -2, 0, SIN, { 0, 0, 0, 2, 0,-7, 0, 3 } },
  {     2, 0, COS, { 0, 0, 0, 1, 0,-4, 0, 0 } },
  {     2, 0, COS, { 0, 0, 0, 4, 0,-11,0, 0 } },
  {    -2, 0, SIN, { 0, 0, 0, 0, 0, 1, 0,-1 } },

  {  8297, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  { -3346, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   462, 0, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {  -189, 0, COS, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {   185, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {    79, 1, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   -71, 0, COS, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {    46, 0, SIN, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {   -45, 0, COS, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {    29, 0, SIN, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {   -20, 0, COS, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {    18, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   -14, 0, COS, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {   -11, 0, COS, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {   -10, 1, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {     9, 0, SIN, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {     8, 0, SIN, { 0, 0, 0, 1, 0,-1, 0, 0 } },
  {    -6, 0, SIN, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {     5, 0, SIN, { 0, 0, 0, 2, 0,-7, 0, 0 } },
  {    -5, 0, COS, { 0, 0, 0, 2, 0,-7, 0, 0 } },
  {     4, 0, SIN, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {    -4, 1, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {    -3, 0, COS, { 0, 0, 0, 1, 0,-1, 0, 0 } },
  {     3, 0, COS, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {     3, 1, SIN, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {     3, 0, SIN, { 0, 0, 0, 1, 0,-2, 0, 0 } },
  {     2, 0, SIN, { 0, 0, 0, 0, 0, 4, 0, 0 } },
  {    -2, 0, COS, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  
  { 955774, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  { -53252, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {  -1878, 0, SIN, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {  -1482, 0, COS, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {    817, 0, SIN, { 0, 0, 0, 1, 0,-1, 0, 0 } },
  {   -539, 0, COS, { 0, 0, 0, 1, 0,-2, 0, 0 } },
  {   -524, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    349, 0, SIN, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {    347, 0, SIN, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {    328, 1, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   -225, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    149, 0, COS, { 0, 0, 0, 2, 0,-6, 0, 0 } },
  {   -126, 0, COS, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {    104, 0, COS, { 0, 0, 0, 1, 0,-1, 0, 0 } },
  {    101, 0, COS, { 0, 0, 0, 2, 0,-5, 0, 0 } },
  {     98, 0, COS, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {    -73, 0, COS, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {    -62, 0, COS, { 0, 0, 0, 0, 0, 3, 0, 0 } },
  {     42, 0, SIN, { 0, 0, 0, 0, 0, 2, 0,-3 } },
  {     41, 0, SIN, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {    -40, 0, SIN, { 0, 0, 0, 1, 0,-3, 0, 0 } },
  {     40, 0, COS, { 0, 0, 0, 2, 0,-4, 0, 0 } },
  {    -28, 1, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {    -23, 0, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {     20, 0, SIN, { 0, 0, 0, 2, 0,-7, 0, 0 } } };

  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[4] + SumVFPTerms ( args, t, terms, 0, 77 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 78, 105 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 106, 130 ) / 1.0e5;
}

/******************************  VFPUranus  **********************************/

void VFPUranus ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  { 0.056531, 0.00023080893,
    0.882987, 0.00009294371,
    0.870169, 0.00003269438,
    0.400589, 0.00003269438,
    0.664614, 0.00003265562,
    0.725368, 0.00001672092,
    0.000000, 0.00000000000,
    0.000000, 0.00000000000,
  };
  	
  static VFPTerm terms[54] = {
  { 19397, 0, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {   570, 0, SIN, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {  -536, 1, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {   143, 0, SIN, { 0, 1, 0,-2, 0, 0, 0, 0 } },
  {   110, 1, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {   102, 0, SIN, { 0, 1, 0,-3, 0, 0, 0, 0 } },
  {    76, 0, COS, { 0, 1, 0,-3, 0, 0, 0, 0 } },
  {   -49, 0, SIN, { 1, 0, 0,-1, 0, 0, 0, 0 } },
  {    32, 2, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {   -30, 1, COS, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {    29, 0, SIN, { 2,-6, 0, 3, 0, 0, 0, 0 } },
  {    29, 0, COS, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {   -28, 0, COS, { 0, 0, 0, 1, 0,-1, 0, 0 } },
  {    23, 0, SIN, { 0, 0, 0, 3, 0, 0, 0, 0 } },
  {   -21, 0, COS, { 1, 0, 0,-1, 0, 0, 0, 0 } },
  {    20, 0, SIN, { 0, 0, 0, 1, 0,-1, 0, 0 } },
  {    20, 0, COS, { 0, 1, 0,-2, 0, 0, 0, 0 } },
  {   -19, 0, COS, { 0, 1, 0,-1, 0, 0, 0, 0 } },
  {    17, 0, SIN, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {    14, 0, SIN, { 0, 0, 0, 3, 0,-3, 0, 0 } },
  {    13, 0, SIN, { 0, 1, 0,-1, 0, 0, 0, 0 } },
  {   -12, 2, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {   -12, 0, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {    10, 0, SIN, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {    -9, 0, SIN, { 0, 0, 0, 0, 2, 0, 0, 0 } },
  {    -9, 2, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {     9, 0, COS, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {     8, 1, COS, { 0, 1, 0,-2, 0, 0, 0, 0 } },
  {     7, 1, COS, { 0, 1, 0,-3, 0, 0, 0, 0 } },
  {    -7, 1, SIN, { 0, 1, 0,-3, 0, 0, 0, 0 } },
  {     7, 1, SIN, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {     6, 0, SIN, { 2,-6, 0, 2, 0, 0, 0, 0 } },
  {     6, 0, COS, { 2,-6, 0, 2, 0, 0, 0, 0 } },
  {     5, 0, SIN, { 0, 1, 0,-4, 0, 0, 0, 0 } },
  {    -4, 0, SIN, { 0, 0, 0, 3, 0,-4, 0, 0 } },
  {     4, 0, COS, { 0, 0, 0, 3, 0,-3, 0, 0 } },
  {    -3, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    -2, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },

  {  2775, 0, SIN, { 0, 0, 0, 0, 1, 0, 0, 0 } },
  {   131, 0, SIN, { 0, 0, 0, 1,-1, 0, 0, 0 } },
  {   130, 0, SIN, { 0, 0, 0, 1, 1, 0, 0, 0 } },

  {1921216, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  { -90154, 0, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {  -2488, 1, SIN, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {  -2121, 0, COS, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {   -585, 0, COS, { 0, 1, 0,-2, 0, 0, 0, 0 } },
  {   -508, 1, COS, { 0, 0, 0, 1, 0, 0, 0, 0 } },
  {   -451, 0, COS, { 1, 0, 0,-1, 0, 0, 0, 0 } },
  {    336, 0, SIN, { 0, 1, 0,-1, 0, 0, 0, 0 } },
  {    198, 0, SIN, { 1, 0, 0,-1, 0, 0, 0, 0 } },
  {    118, 0, COS, { 0, 1, 0,-3, 0, 0, 0, 0 } },
  {    107, 0, SIN, { 0, 1, 0,-2, 0, 0, 0, 0 } },
  {   -103, 1, SIN, { 0, 0, 0, 2, 0, 0, 0, 0 } },
  {    -81, 0, COS, { 0, 0, 0, 3, 0,-3, 0, 0 } } };
  
  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[2] + SumVFPTerms ( args, t, terms, 0, 37 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 38, 40 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 41, 53 ) / 1.0e5;
}

/******************************  VFPNeptune  **********************************/

void VFPNeptune ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  { 0.056531, 0.00023080893,
    0.882987, 0.00009294371,
    0.870169, 0.00003269438,
    0.400589, 0.00003269438,
    0.846912, 0.00001672092,
    0.725368, 0.00001672092,
    0.480856, 0.00001663715,
    0.000000, 0.00000000000
  };

  static VFPTerm terms[26] = {
  {  3523, 0, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   -50, 0, SIN, { 0, 0, 0, 0, 0, 0, 2, 0 } },
  {   -43, 1, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    29, 0, SIN, { 1, 0, 0, 0, 0,-1, 0, 0 } },
  {    19, 0, SIN, { 0, 0, 0, 0, 0, 2, 0, 0 } },
  {   -18, 0, COS, { 1, 0, 0, 0, 0,-1, 0, 0 } },
  {    13, 0, COS, { 0, 1, 0, 0, 0,-1, 0, 0 } },
  {    13, 0, SIN, { 0, 1, 0, 0, 0,-1, 0, 0 } },
  {    -9, 0, SIN, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {     9, 0, COS, { 0, 0, 0, 2, 0,-2, 0, 0 } },
  {    -5, 0, COS, { 0, 0, 0, 2, 0,-3, 0, 0 } },
  {    -4, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {     4, 0, COS, { 0, 0, 0, 1, 0,-2, 0, 0 } },
  {     4, 2, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
	 
  {  6404, 0, SIN, { 0, 0, 0, 0, 0, 0, 1, 0 } },
  {    55, 0, SIN, { 0, 0, 0, 0, 0, 1, 1, 0 } },
  {    55, 0, SIN, { 0, 0, 0, 0, 0, 1,-1, 0 } },
  {   -33, 1, SIN, { 0, 0, 0, 0, 0, 0, 1, 0 } },
	 
  {3007175, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  { -25701, 0, COS, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {   -787, 0, COS, { 0, 0, 2,-1,-2, 0, 0, 0 } },
  {    409, 0, COS, { 1, 0, 0, 0, 0,-1, 0, 0 } },
  {   -314, 1, SIN, { 0, 0, 0, 0, 0, 1, 0, 0 } },
  {    250, 0, SIN, { 1, 0, 0, 0, 0,-1, 0, 0 } },
  {   -194, 0, SIN, { 0, 1, 0, 0, 0,-1, 0, 0 } },
  {    185, 0, COS, { 0, 1, 0, 0, 0,-1, 0, 0 } } };

  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[4] + SumVFPTerms ( args, t, terms, 0, 13 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 14, 17 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 18, 25 ) / 1.0e5;
}

/******************************  VFPPluto  **********************************/

void VFPPluto ( double jd, double *lon, double *lat, double *rp )
{
  double t, args[MAX_NUM_ARGS];
  char i;

  static double arguments[MAX_NUM_ARGS][2] = 
  { 0.663854, 0.00001115482,
    0.041020, 0.00001104864,
    0.357355, 0.00001104864,
    0.779072, 0.00273790931,
    0.000000, 0.00000000000,
    0.000000, 0.00000000000,
    0.000000, 0.00000000000,
    0.000000, 0.00000000000,
  };
  
  static VFPTerm terms[28] = {
  {101577, 0, SIN, { 0, 1, 0, 0, 0, 0, 0, 0 } },
  { 15517, 0, SIN, { 0, 2, 0, 0, 0, 0, 0, 0 } },
  { -3593, 0, SIN, { 0, 0, 2, 0, 0, 0, 0, 0 } },
  {  3414, 0, SIN, { 0, 3, 0, 0, 0, 0, 0, 0 } },
  { -2201, 0, SIN, { 0, 1,-2, 0, 0, 0, 0, 0 } },
  { -1871, 0, SIN, { 0, 1, 2, 0, 0, 0, 0, 0 } },
  {   839, 0, SIN, { 0, 4, 0, 0, 0, 0, 0, 0 } },
  {  -757, 0, SIN, { 0, 2, 2, 0, 0, 0, 0, 0 } },
  {  -285, 0, SIN, { 0, 3, 2, 0, 0, 0, 0, 0 } },
  {   227, 2, SIN, { 0, 1, 0, 0, 0, 0, 0, 0 } },
  {   218, 0, SIN, { 0, 2,-2, 0, 0, 0, 0, 0 } },
  {   200, 1, SIN, { 0, 1, 0, 0, 0, 0, 0, 0 } },
  
  { 57726, 0, SIN, { 0, 0, 1, 0, 0, 0, 0, 0 } },
  { 15257, 0, SIN, { 0, 1,-1, 0, 0, 0, 0, 0 } },
  { 14102, 0, SIN, { 0, 1, 1, 0, 0, 0, 0, 0 } },
  {  3870, 0, SIN, { 0, 2, 1, 0, 0, 0, 0, 0 } },
  {  1138, 0, SIN, { 0, 3, 1, 0, 0, 0, 0, 0 } },
  {   472, 0, SIN, { 0, 2,-1, 0, 0, 0, 0, 0 } },
  {   353, 0, SIN, { 0, 4, 1, 0, 0, 0, 0, 0 } },
  {  -144, 0, SIN, { 0, 1,-3, 0, 0, 0, 0, 0 } },
  {  -119, 0, SIN, { 0, 0, 3, 0, 0, 0, 0, 0 } },
  {  -111, 0, SIN, { 0, 1, 3, 0, 0, 0, 0, 0 } },
  
  {4074638, 0, COS, { 0, 0, 0, 0, 0, 0, 0, 0 } },
  {-958235, 0, COS, { 0, 1, 0, 0, 0, 0, 0, 0 } },
  {-116703, 0, COS, { 0, 2, 0, 0, 0, 0, 0, 0 } },
  { -22649, 0, COS, { 0, 3, 0, 0, 0, 0, 0, 0 } },
  {  -4996, 0, COS, { 0, 4, 0, 0, 0, 0, 0, 0 } } };

  t = jd - 2451545.0;
  
  for ( i = 0; i < MAX_NUM_ARGS; i++ )
    args[i] = TWO_PI * fmod ( arguments[i][0] + arguments[i][1] * t, 1.0 );
  
  t = t / 36525.0 + 1.0;
  
  *lon = fmod ( args[0] + SumVFPTerms ( args, t, terms, 0, 11 ) / ARCSEC_PER_RAD, TWO_PI );
  *lat = SumVFPTerms ( args, t, terms, 12, 21 ) / ARCSEC_PER_RAD;
  *rp = SumVFPTerms ( args, t, terms, 22, 27 ) / 1.0e5;
}

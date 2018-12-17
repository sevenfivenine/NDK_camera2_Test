/*** COPYRIGHT NOTICE *********************************************************
 
	Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

// These headers needed for AACurrentUTC() function below.

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "AstroLib.h"

#define floor(x)	floor((double)(x))	// to keep Visual C++ 2012 happy

/*** AALocalJD *************************************************************/

double AALocalJD ( void )
{
	time_t		timer;
	struct tm	t;
	double		jd;
	
	time ( &timer );
	t = *localtime ( &timer );
	jd = AADateTimeToJD ( t.tm_year + 1900, t.tm_mon + 1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec, 0.0, TRUE );

	return ( jd );
}

/*** AACurrentJD *************************************************************/

double AACurrentJD ( void )
{
	time_t		timer;
	struct tm	t;
	double		jd;
	
	time ( &timer );
	t = *gmtime ( &timer );
	jd = AADateTimeToJD ( t.tm_year + 1900, t.tm_mon + 1, t.tm_mday, t.tm_hour, t.tm_min, t.tm_sec, 0.0, TRUE );
	
	return ( jd );
}

/*** AACurrentUTC *************************************************************/

double AACurrentUTC ( void )
{
	double		jd = 0.0;
#ifdef WIN32
	SYSTEMTIME	systime = { 0 };

	GetSystemTime ( &systime );
	jd = AADateTimeToJD ( systime.wYear, systime.wMonth, systime.wDay, systime.wHour, systime.wMinute, systime.wSecond + systime.wMilliseconds / 1000.0, 0.0, TRUE );
#else
	struct timeval	tv = { 0 };
	
	gettimeofday ( &tv, NULL );
	jd = J1970 + ( tv.tv_sec + tv.tv_usec / 1000000.0 ) / 86400.0;
#endif
	return ( jd );
}

/*** AACurrentTimeZone **********************************************************/

double AACurrentTimeZone ( int *isDST )
{
	time_t		t = time ( NULL );
	struct tm	*tp = localtime ( &t );
	double		zone = tp->tm_gmtoff / 86400.0 - ( tp->tm_isdst ? 1.0 / 24.0 : 0.0 );
	
	if ( isDST )
		*isDST = tp->tm_isdst;
	
	return ( zone );
}

/*** AADateTimeToJD *********************************************************/
   
double AADateTimeToJD ( int year, short month, double day,
short hour, short min, double sec, double zone, short calendar )
{
	double	jd;
	int 	a, b;

	day += hour / 24.0 + min / 1440.0 + sec / 86400.0 - zone;

	if ( calendar == AA_CALENDAR_MAYAN )
	{
		jd = AAMayanLongCountToJD ( year );
		jd += day - floor ( day );
	}
	else if ( calendar == AA_CALENDAR_INDIAN )
	{
		jd = AAIndianToJD ( year, month, day );
	}
	else if ( calendar == AA_CALENDAR_PERSIAN )
	{
		jd = AAPersianToJD ( year, month, day );
	}
	else if ( calendar == AA_CALENDAR_ISLAMIC )
	{
		jd = AAIslamicToJD ( year, month, day );
	}
	else if ( calendar == AA_CALENDAR_HEBREW )
	{
		jd = AAHebrewToJD ( year, month, day );
	}
	else if ( calendar == AA_CALENDAR_JULIAN )
	{
		jd = AAJulianToJD ( year, month, day );
	}
	else if ( calendar == AA_CALENDAR_GREGORIAN )
	{
		jd = AAGregorianToJD ( year, month, day );
	}
	else
	{
		if ( month < 3 )
		{
			year -= 1;
			month += 12;
		}

		// Use Gregorian calendar after 4 October 1582
		
		if ( ( year > 1582 ) || ( year == 1582 && month > 10 ) || ( year == 1582 && month == 10 && day >= 5.0 ) )
		{
			a = year / 100;
			b = 2 - a + a / 4;
		}
		else
			b = 0;

		jd = floor ( 365.25 * ( year + 4716 ) ) + floor ( 30.6001 * ( month + 1 ) ) +
			day + b - 1524.5;
	}
	
	return ( jd );
}

/*** AAJDToDateTime **************************************************************/

void AAJDToDateTime ( double jd, double zone, int *year, short *month, double *day,
short *hour, short *min, double *sec, short calendar )
{
	int		z, a, b, c, d, e;
	double	f;

	jd += zone;
	
	if ( calendar == AA_CALENDAR_MAYAN )
	{
		*year = AAJDToMayanLongCount ( jd );
		AAJDToMayanHaab ( jd, month, day );
	}
	else if ( calendar == AA_CALENDAR_INDIAN )
	{
		AAJDToIndian ( jd, year, month, day );
	}
	else if ( calendar == AA_CALENDAR_PERSIAN )
	{
		AAJDToPersian ( jd, year, month, day );
	}
	else if ( calendar == AA_CALENDAR_ISLAMIC )
	{
		AAJDToIslamic ( jd, year, month, day );
	}
	else if ( calendar == AA_CALENDAR_HEBREW )
	{
		AAJDToHebrew ( jd, year, month, day );
	}
	else if ( calendar == AA_CALENDAR_JULIAN )
	{
		AAJDToJulian ( jd, year, month, day );
	}
	else if ( calendar == AA_CALENDAR_GREGORIAN )
	{
		AAJDToGregorian ( jd, year, month, day );
	}
	else
	{
		jd += 0.5;
		z = floor ( jd );
		f = jd - z;

		// Use Gregorian calendar after 4 October 1582

		if ( jd >= 2299161.0 )
		{
			a = ( z - 1867216.25 ) / 36524.25;
			a = z + 1 + a - a / 4;
		}
		else
			a = z;

		b = a + 1524;
		c = floor ( ( b - 122.1 ) / 365.25 );
		d = floor ( 365.25 * c );
		e = ( b - d ) / 30.6001;

		if ( e < 14 )
			*month = e - 1;
		else
			*month = e - 13;

		if ( *month > 2 )
			*year = c - 4716;
		else
			*year = c - 4715;

		*day = b - d - floor ( 30.6001 * e ) + f;
	}

	f = *day - floor ( *day );

	*hour = 24.0 * f;
	*min = 1440.0 * f - 60.0 * *hour;
	*sec = 86400.0 * f - 3600.0 * *hour - 60.0 * *min;
  
	/*** Correct floating-point roundoff error ***/
  
	if ( *sec <= 0.0 )
		*sec = 0.0;
}

/*** AAJulianYearToJD *****************************************************/

double AAJulianYearToJD ( double year )
{
	return ( J2000 + DAY_PER_JYEAR * ( year - 2000.0 ) );
}

/*** AABesselianYearToJD  ***************************************************/

double AABesselianYearToJD ( double year )
{
	return ( B1900 + DAY_PER_BYEAR * ( year - 1900.0 ) );
}

/*** AALocalWeekDay *******************************************************/

int AALocalWeekDay ( double jd, double zone )
{
	int	d = floor ( jd + zone + 1.5 );
	
	d = d % 7;
	if ( d < 0 )
		d += 7;
		
	return ( d );
}

/*** AADaylightSavingsTime *************************************************/

int AADaylightSavingsTime ( double jd, double zone, int rule, double *start, double *end )
{
	int		year;
	short	month, hour, min, weekday;
	double	day, sec;
	
	// Degenerate cases
	
	if ( rule == AA_DLST_ALWAYS )
	{
		*start = -HUGE_VAL;
		*end = HUGE_VAL;
		return ( TRUE );
	}
	else if ( rule == AA_DLST_NONE )
	{
		*start = 0.0;
		*end = 0.0;
		return ( FALSE );
	}
	
	AAJDToDateTime ( jd, zone, &year, &month, &day, &hour, &min, &sec, AA_CALENDAR_GREGORIAN );
	
	// USA, Canada, and Mexico: before 2007, Daylight Saving Time starts at 2 AM
	// on the first Sunday in April and ends at 2 AM on the last Sunday in October.
	// Beginning in 2007, Daylight Saving Time starts at 2 AM on the second Sunday
	// in March and ends at 2 AM on the first Sunday in November in the USA and
	// Canada; Mexico still observes the old algorithm.
		
	if ( rule == AA_DLST_USA_CANADA && year >= 2007 )
	{
		*start = AADateTimeToJD ( year, 3, 8.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		if ( weekday > 0 )
			*start = *start + ( 7 - weekday );
		
		*end = AADateTimeToJD ( year, 11, 1.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		if ( weekday > 0 )
			*end = *end + ( 7 - weekday );
	}
	else if ( rule == AA_DLST_MEXICO || ( rule == AA_DLST_USA_CANADA && year < 2007 ) )
	{
		*start = AADateTimeToJD ( year, 4, 1.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		if ( weekday > 0 )
			*start = *start + ( 7 - weekday );
		
		*end = AADateTimeToJD ( year, 10, 31.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		*end = *end - weekday;
	}
	
	// Brazil (southern states): Daylight Saving Time starts at local midnight on the first
	// Sunday in October and ends at local midnight on the third Sunday in February.
	
	else if ( rule == AA_DLST_BRAZIL )
	{
		*start = AADateTimeToJD ( year, 10, 1.0, 0, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		if ( weekday > 0 )
			*start = *start + ( 7 - weekday );
		
		*end = AADateTimeToJD ( year, 2, 15.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		if ( weekday > 0 )
			*end = *end + ( 7 - weekday );
	}
	
	// Chile: Daylight Saving Time starts at local midnight on the second
	// Saturday in October and ends at local on the second Saturday of March.
	
	else if ( rule == AA_DLST_CHILE )
	{
		*start = AADateTimeToJD ( year, 10, 8.0, 0, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		if ( weekday != 6 )
			*start = *start + ( 6 - weekday );
		
		*end = AADateTimeToJD ( year, 3, 8.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		if ( weekday != 6 )
			*end = *end + ( 6 - weekday );
	}

	// European Union: Daylight Saving Time starts at 1 AM UTC on the last
	// Sunday in March and ends at 1 AM UTC on the last Sunday in October.

	else if ( rule == AA_DLST_EUROPE )
	{
		*start = AADateTimeToJD ( year, 3, 31.0, 1, 0, 0.0, 0.0, 1 );
		weekday = AALocalWeekDay ( *start, 0.0 );
		*start = *start - weekday;
	
		*end = AADateTimeToJD ( year, 10, 31.0, 1, 0, 0.0, 0.0, 1 );
		weekday = AALocalWeekDay ( *end, 0.0 );
		*end = *end - weekday;
	}
	
	// Russia (and most former USSR republics): Daylight Saving Time starts at
	// 2 AM local time on the last Sunday in March and ends at 2 AM local time
	// on the last Sunday in October.  Russia abolished DLST in 2011.

	else if ( rule == AA_DLST_RUSSIA && year < 2011 )
	{
		*start = AADateTimeToJD ( year, 3, 31.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		*start = *start - weekday;
	
		*end = AADateTimeToJD ( year, 10, 31.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		*end = *end - weekday;
	}
	
	// Iran: Daylight Saving Time starts at midnight local time on March 22
	// and ends at midnight on September 22.
	
	else if ( rule == AA_DLST_IRAN )
	{
		*start = AADateTimeToJD ( year, 3, 22.0, 0, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		*start = *start - weekday;
		
		*end = AADateTimeToJD ( year, 9, 22.0, 0, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		*end = *end - weekday;
	}

	// Australia: Daylight Saving Time starts at 2 AM local time on the first
	// Sunday in October and ends at 2 AM local time on the last Sunday in March
	// (changed to first Sunday in April in 2008).

	else if ( rule == AA_DLST_AUSTRALIA )
	{
		*start = AADateTimeToJD ( year, 10, 1.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		*start = *start + ( 7 - weekday );

		if ( year >= 2008 )
		{
			*end = AADateTimeToJD ( year, 4, 1.0, 2, 0, 0.0, zone, 1 );
			weekday = AALocalWeekDay ( *end, zone );
			if ( weekday > 0 )
				*end = *end + ( 7 - weekday );
		}
		else
		{
			*end = AADateTimeToJD ( year, 3, 31.0, 2, 0, 0.0, zone, 1 );
			weekday = AALocalWeekDay ( *end, zone );
			*end = *end - weekday;
		}
	}

	// New Zealand: Daylight Saving Time starts at 2 AM local time on the last
	// Sunday in September and ends at 2 AM local time on the first Sunday in April.

	else if ( rule == AA_DLST_NEW_ZEALAND )
	{
		*start = AADateTimeToJD ( year, 9, 30.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		*start = *start - weekday;

		*end = AADateTimeToJD ( year, 4, 1.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		if ( weekday > 0 )
			*end = *end + ( 7 - weekday );
	}
	
	// Fiji: Daylight Saving Time starts at 2 AM local time on the last
	// Sunday in October and ends at 2 AM local time on the third Sunday in January.
	
	else if ( rule == AA_DLST_FIJI )
	{
		*start = AADateTimeToJD ( year, 10, 31.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *start, zone );
		*start = *start - weekday;
		
		*end = AADateTimeToJD ( year, 1, 15.0, 2, 0, 0.0, zone, 1 );
		weekday = AALocalWeekDay ( *end, zone );
		if ( weekday > 0 )
			*end = *end + ( 7 - weekday );
	}

	// unrecognized DLST rule - return FALSE and don't attempt
	// to compute DLST start/end dates.
	
	else
	{
		return ( FALSE );
	}
	
	// If within a particular year DLST starts before it ends
	// (i.e. the Northern hemisphere), then it is DLST if the
	// current date is after the start and before the end of
	// DLST.  If DLST ends before it starts in a particular
	// year (i.e. the Southern hemisphere) then it is DLST if
	// the current date is before the end or after the start
	// of DLST.
	
	if ( *start < *end )
	{
		if ( jd > *start && jd < *end )
			return ( TRUE );
		else
			return ( FALSE );
	}
	else
	{
		if ( jd < *end || jd > *start )
			return ( TRUE );
		else
			return ( FALSE );
	}
}

/*** AAJDToJulianYear ****************************************************/

double AAJDToJulianYear ( double jd )
{
	return ( ( jd - J2000 ) / DAY_PER_JYEAR + 2000.0 );
}

/*** AAJDToBesselianYear **************************************************/

double AAJDToBesselianYear ( double jd )
{
	return ( ( jd - B1900 ) / DAY_PER_BYEAR + 1900.0 );
}

/*** AAGregorianToJD *********************************************************/

double AAGregorianToJD ( int y, int m, double d )
{
	double	jd;
	
	jd = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4
	   + ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12
	   - ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4
	   + d - 32075.5;
	   
	return ( jd );
}

/*** AAJDToGregorian *********************************************************/

void AAJDToGregorian ( double jd, int *y, short *m, double *d )
{
	int		l, n, i, j;
	double	f;
	
	jd = jd + 0.5;
	j = floor ( jd );
	f = jd - j;
	
	l = j + 68569;
	n = floor ( ( 4 * l ) / 146097.0 );
	l = l - ( 146097 * n + 3 ) / 4;
	i = ( 4000 * ( l + 1 ) ) / 1461001;
	l = l - ( 1461 * i ) / 4 + 31;
	j = ( 80 * l ) / 2447;
	
	*d = l - ( 2447 * j ) / 80 + f;
	l = j / 11;
	*m = j + 2 - 12 * l;
	*y = 100 * ( n - 49 ) + i + l;
}

/*** AAGregorianToJD **********************************************************

// Alternate implementation adapted from:
// http://www.fourmilab.ch/documents/calendar/

#define GREGORIAN_EPOCH		1721425.5

int IsGregorianLeapYear ( int year )
{
    int leap;
	
	leap = ( year % 4 == 0 ) && ! ( ( year % 100 == 0 ) && ( year % 400 != 0 ) );
	
	return ( leap );
}

double AAGregorianToJD ( int year, int month, double day )
{
	double jd;
	
    jd = GREGORIAN_EPOCH - 1 +
         365 * ( year - 1 ) +
         floor ( ( year - 1 ) / 4.0 ) -
         floor ( ( year - 1 ) / 100.0 ) +
         floor ( ( year - 1 ) / 400.0 ) +
         floor ( ( 367 * month - 362 ) / 12.0 );
		 
	if ( month > 2 )
		jd -= IsGregorianLeapYear ( year ) ? 1 : 2;

	jd += day;
	
	return ( jd );
}

******************************************************************************/

/*** AAJDToGregorian **********************************************************

// Alternate implementation adapted from:
// http://www.fourmilab.ch/documents/calendar/

void AAJDToGregorian ( double jd, int *year, short *month, double *day )
{
    double wjd, depoch, quadricent, dqc, cent, dcent, quad, dquad,
        yindex, yearday, leapadj;

    wjd = floor ( jd - 0.5 ) + 0.5;
    depoch = wjd - GREGORIAN_EPOCH;
    quadricent = floor ( depoch / 146097.0 );
    dqc = fmod ( depoch, 146097.0 );
    cent = floor ( dqc / 36524.0 );
    dcent = fmod ( dqc, 36524.0 );
    quad = floor ( dcent / 1461.0 );
    dquad = fmod ( dcent, 1461.0 );
    yindex = floor ( dquad / 365.0 );
    *year = ( quadricent * 400 ) + ( cent * 100 ) + ( quad * 4 ) + yindex;
    if ( ! ( cent == 4 || yindex == 4 ) )
        (*year)++;
		
    yearday = wjd - AAGregorianToJD ( *year, 1, 1 );
	
	if ( wjd >= AAGregorianToJD ( year, 3, 1 ) )
		leapadj = IsGregorianLeapYear ( *year ) ? 1 : 2;
	else
		leapadj = 0.0;
		
    *month = floor ( ( ( ( yearday + leapadj ) * 12 ) + 373 ) / 367 );
    *day = jd - AAGregorianToJD ( *year, *month, 1.0 ) + 1.0;
}

***************************************************************************/

/*** AAJulianToJD *********************************************************/

double AAJulianToJD ( int y, int m, double d )
{
	double jd;
	
	jd = 367 * y
	   - ( 7 * ( y + 5001 + ( m - 9 ) / 7 ) ) / 4
	   + ( 275 * m ) / 9 + d + 1729776.5;
	   
	return ( jd );
}

/*** AAJDToJulian *********************************************************/

void AAJDToJulian ( double jd, int *y, short *m, double *d )
{
	int		j, k, l, n, i;
	double	f;
	
	jd = jd + 0.5;
	j = floor ( jd );
	f = jd - j;

	j = j + 1402;
	k = floor ( ( j - 1 ) / 1461.0 );
	l = j - 1461 * k;
	n = ( l - 1 ) / 365 - l / 1461;
	i = l - 365 * n + 30;
	j = ( 80 * i ) / 2447;
	
	*d = i - ( 2447 * j ) / 80 + f;
	i = j / 11;
	*m = j + 2 - 12 * i;
	*y = 4 * k + n + i - 4716;
}

/*** AAHebrewToJD **********************************************************/

int mod ( int l, int n )
{
	int m = l % n;
	
	if ( m < 0 )
		m += n;
		
	return ( m );
}

int hebrew_leap ( int year )
{
	return ( mod ( ( year * 7 ) + 1, 19 ) < 7 ? TRUE : FALSE );
}

//  How many months are there in a Hebrew year (12 = normal, 13 = leap)

int hebrew_year_months ( int year )
{
	return ( hebrew_leap ( year ) ? 13 : 12 );
}

//  Test for delay of start of new year and to avoid
//  Sunday, Wednesday, and Friday as start of the new year.

int hebrew_delay_1 ( int year )
{
    int months, parts, day;

    months = floor ( ( 235 * year - 234) / 19 );
    parts = 12084 + ( 13753 * months );
    day = months * 29 + floor ( parts / 25920.0 );

    if ( mod ( ( 3 * ( day + 1 ) ), 7 ) < 3 )
        day++;

    return ( day );
}

//  Check for delay in start of new year due to length of adjacent years

int hebrew_delay_2 ( int year )
{
    int last, present, next;

    last = hebrew_delay_1 ( year - 1 );
    present = hebrew_delay_1 ( year );
    next = hebrew_delay_1 ( year + 1 );

    if ( ( next - present) == 356 )
		return ( 2 );
	else if ( ( present - last ) == 382 )
		return ( 1 );
	else
		return ( 0 );
}

//  How many days are in a Hebrew year ?

int hebrew_year_days ( int year )
{
    return ( AAHebrewToJD ( year + 1, 7, 1 ) - AAHebrewToJD ( year, 7, 1 ) );
}

//  How many days are in a given month of a given year

int hebrew_month_days ( int year, short month )
{
    //  First of all, dispose of fixed-length 29 day months

    if ( month == 2 || month == 4 || month == 6 || month == 10 || month == 13 )
		return ( 29 );

    //  If it's not a leap year, Adar has 29 days

    if ( month == 12 && ! hebrew_leap ( year ) )
        return ( 29 );

    //  If it's Heshvan, days depend on length of year

    if ( month == 8 && ! ( mod ( hebrew_year_days ( year ), 10 ) == 5 ) )
		return ( 29 );

    //  Similarly, Kislev varies with the length of year

    if ( month == 9 && ( mod ( hebrew_year_days ( year ), 10 ) == 3 ) )
        return ( 29 );

    //  Nope, it's a 30 day month

    return ( 30 );
}

//  Finally, wrap it all up into...

#define HEBREW_EPOCH 347995.5

double AAHebrewToJD ( int year, short month, double day )
{
	double	jd = 0;
	short	mon, months;

    months = hebrew_year_months ( year );
    jd = HEBREW_EPOCH + hebrew_delay_1 ( year ) +
         hebrew_delay_2 ( year ) + day + 1;

    if ( month < 7 )
	{
        for ( mon = 7; mon <= months; mon++ )
            jd += hebrew_month_days ( year, mon );

        for ( mon = 1; mon < month; mon++ )
            jd += hebrew_month_days ( year, mon );

    }
	else
	{
        for ( mon = 7; mon < month; mon++ )
            jd += hebrew_month_days ( year, mon );
    }

	return ( jd );
}

/*** AAHebrewToJD **********************************************************/

void AAJDToHebrew ( double jd, int *year, short *month, double *day )
{
	double	jd0;
	int		i, count, first;

	jd0 = floor ( jd - 0.5 ) + 0.5;
    count = floor ( ( ( jd0 - HEBREW_EPOCH ) * 98496.0 ) / 35975351.0 );
    *year = count - 1;
    for ( i = count; jd0 >= AAHebrewToJD ( i, 7, 1 ); i++ )
		(*year)++;

    first = jd0 < AAHebrewToJD ( *year, 1, 1 ) ? 7 : 1;
	*month = first;
    for ( i = first; jd0 > AAHebrewToJD ( *year, i, hebrew_month_days ( *year, i ) ); i++ )
        (*month)++;

	*day = jd - AAHebrewToJD ( *year, *month, 1 ) + 1.0;
}

/*** AAIslamicToJD *********************************************************/

double AAIslamicToJD ( int y, int m, double d )
{
	int		jd0;
	double	jd;
	
	jd0 = 1948440; // : 1948939;
	
	jd = floor ( ( 11 * y + 3 ) / 30.0 )
	   + 354 * y
	   + 30 * m
	   - ( m - 1 ) / 2
	   + d + jd0 - 385.5;
	   
	return ( jd );
}

/*** AAJDToIslamic **********************************************************/

void AAJDToIslamic ( double jd, int *y, short *m, double *d )
{
	int		jd0, l, n, j;
	double	f;
	
	jd = jd + 0.5;
	j = floor ( jd );
	f = jd - j;

	jd0 = 1948440; // alternately 1948939;
	
	l = j - jd0 + 10632;
	n = floor ( ( l - 1 ) / 10631.0 );
	l = l - 10631 * n + 354;
	j = ( ( 10985 - l ) / 5316 ) * ( ( 50 * l ) / 17719 )
	  + ( l / 5670 ) * ( ( 43 * l ) / 15238 );
	l = l - ( ( 30 - j ) / 15 ) * ( ( 17719 * j ) / 50 )
	  - ( j / 16 ) * ( ( 15238 * j ) / 43 ) + 29;
	  
	*m = ( 24 * l ) / 709;
	*d = l - ( *m * 709 ) / 24 + f;
	*y = 30 * n + j - 30;
}

/*** AAIslamicToJD **********************************************************

// Alternate implementation adapted from:
// http://www.fourmilab.ch/documents/calendar/

#define ISLAMIC_EPOCH	1948439.5

double AAIslamicToJD ( int year, int month, double day )
{
    double jd;
	
	jd = day + 
	     ceil ( 29.5 * ( month - 1 ) ) +
         ( year - 1 ) * 354 +
         floor ( ( 3 + ( 11 * year ) ) / 30.0 ) +
         ISLAMIC_EPOCH - 1.0;
		 
	return ( jd );
}

*****************************************************************************/

/*** AAJDToIslamic **********************************************************

// Alternate implementation adapted from:
// http://www.fourmilab.ch/documents/calendar/

void AAJDToIslamic ( double jd, int *year, short *month, double *day )
{
	double	jd0;

    jd0 = floor ( jd - 0.5 ) + 0.5;
    *year = floor ( ( ( 30.0 * ( jd0 - ISLAMIC_EPOCH ) ) + 10646.0 ) / 10631.0 );
    *month = ceil ( ( jd0 - ( 29.0 + AAIslamicToJD ( *year, 1, 1 ) ) ) / 29.5 ) + 1.0;
	if ( *month > 12 )
		*month = 12;
		
    *day = ( jd - AAIslamicToJD ( *year, *month, 1 ) ) + 1.0;
}

****************************************************************************/

/*** AAIndianToJD  *********************************************************/

double AAIndianToJD ( int y, int m, double d )
{
	double jd;
	
	jd = 365 * y
	   + floor ( ( y + 78 - 1 / m ) / 4.0 )
	   + 31 * m - ( m + 9 ) / 11
	   - ( m / 7 ) * ( m - 7 )
	   - floor ( ( 3 * ( floor ( ( y + 78 - 1 / m ) / 100.0 ) + 1 ) ) / 4.0 )
	   + d + 1749578.5;
	   
	return ( jd );
}

/*** AAJDToIndian ***********************************************************/

void AAJDToIndian ( double jd, int *y, short *m, double *d )
{
	int		l, n, i, j;
	double	f;
	
	jd = jd + 0.5;
	j = floor ( jd );
	f = jd - j;

	l = j + 68518;
	n = floor ( ( 4 * l ) / 146097.0 );
	l = l - ( 146097 * n + 3 ) / 4;
	i = ( 4000 * ( l + 1 ) ) / 1461001;
	l = l - ( 1461 * i ) / 4 + 1;
	j = ( ( l - 1 ) / 31 ) * ( 1 - l / 185 )
	  + ( l / 185 ) * ( ( l - 156 ) / 30 + 5 ) - l / 366;
	  
	*d = l - 31 * j + ( ( j + 2 ) / 8 ) * ( j - 5 ) + f;
	l = j / 11;
	*m = j + 2 - 12 * l;
	*y = 100 * ( n - 49 ) + l + i - 78;
}

/*** AAPersianToJD ***********************************************************/

#define PERSIAN_EPOCH	1948320.5

double AAPersianToJD ( int year, short month, double day )
{
	double	jd;
	int	epbase, epyear;

	epbase = year - 474;
    epyear = 474 + mod ( epbase, 2820 );

    jd = day +
       ( month <= 7 ? ( month - 1 ) * 31 : ( month - 1 ) * 30 + 6 ) +
       floor ( ( epyear * 682 - 110 ) / 2816.0 ) + ( epyear - 1 ) * 365 +
       floor ( epbase / 2820.0 ) * 1029983 +
	   PERSIAN_EPOCH - 1;
			
	return ( jd );
}

/*** AAJDToPersian ***********************************************************/

void AAJDToPersian ( double jd, int *year, short *month, double *day )
{
	double	jd0;
	int	depoch, cycle, cyear, ycycle, aux1, aux2, yday;
	
	jd0 = floor ( jd - 0.5 ) + 0.5;

    depoch = jd0 - AAPersianToJD ( 475, 1, 1 );
    cycle = floor ( depoch / 1029983.0 );
    cyear = mod ( depoch, 1029983 );
    if ( cyear == 1029982 )
	{
        ycycle = 2820;
    }
	else
	{
        aux1 = floor ( cyear / 366.0 );
        aux2 = mod ( cyear, 366 );
        ycycle = floor ( ( ( 2134 * aux1 ) + ( 2816 * aux2 ) + 2815 ) / 1028522.0 ) +
                 aux1 + 1;
    }
	
    *year = ycycle + ( 2820 * cycle ) + 474;
    yday = ( jd0 - AAPersianToJD ( *year, 1, 1 ) ) + 1;
    *month = ( yday <= 186 ) ? ceil ( yday / 31.0 ) : ceil ( ( yday - 6 ) / 30.0 );
    *day = jd - AAPersianToJD ( *year, *month, 1 ) + 1;
}

/*** AAMayanLongCountToJD *****************************************************/

#define MAYAN_EPOCH	584282.5

double AAMayanLongCountToJD ( int lcount )
{
	return ( lcount + MAYAN_EPOCH );
}

int AAJDToMayanLongCount ( double jd )
{
	return ( floor ( jd - MAYAN_EPOCH ) );
}

double AAMayanToJD ( short baktun, short katun, short tun, short uinal, short kin )
{
	double jd;
	
	jd = MAYAN_EPOCH
       + baktun * 144000L
	   + katun * 7200L
	   + tun * 360L
	   + uinal * 20L
	   + kin;
	
	return ( jd );
}

void AAJDToMayan ( double jd, short *baktun, short *katun, short *tun, short *uinal, short *kin )
{
	double	jd0;
    int		d;

    jd0 = floor ( jd - 0.5 ) + 0.5;
    d = jd0 - MAYAN_EPOCH;
    *baktun = floor ( d / 144000.0 );
    d = mod ( d, 144000 );
    *katun = floor ( d / 7200.0 );
    d = mod ( d, 7200 );
    *tun = floor ( d / 360.0 );
    d = mod ( d, 360 );
    *uinal = floor ( d / 20.0 );
    *kin = mod ( d, 20 );
}

void AAJDToMayanHaab ( double jd, short *haab, double *day )
{
	double	lcount;

    lcount = jd - MAYAN_EPOCH;
    *day = fmod ( lcount + 8 + ( ( 18 - 1 ) * 20 ), 365.0 );
	if ( *day < 0.0 )
		*day += 365.0;
		
	*haab = floor ( *day / 20.0 ) + 1;
	*day = fmod ( *day, 20.0 );
}

/*** AADeltaT **************************************************************/

double AADeltaT ( double jd )
{
	double y = AAJDToJulianYear ( jd ) - 0.5 / 12.0;
	double u, u2, u3, u4, u5, u6;
	double t, t2, t3, t4, t5, t6, t7;
	double dt = 0;
	
	if ( y < -500.0 )
	{
		u = ( y - 1820.0 ) / 100.0;
		dt = -20.0 + 32.0 * u * u;
	}
	else if ( y < 500.0 )
	{
		u = y / 100.0;
		u2 = u * u;
		u3 = u2 * u;
		u4 = u3 * u;
		u5 = u4 * u;
		u6 = u5 * u;
		dt = 10538.6 - 1014.41 * u + 33.78311 * u2 - 5.952053 * u3
		- 0.1798452 * u4 + 0.022174192 * u5 + 0.0090316521 * u6;
	}
	else if ( y < 1600.0 )
	{
		u = ( y - 1000.0 ) / 100.0;
		u2 = u * u;
		u3 = u2 * u;
		u4 = u3 * u;
		u5 = u4 * u;
		u6 = u5 * u;
		dt = 1574.2 - 556.01 * u + 71.23472 * u2 + 0.319781 * u3
		- 0.8503463 * u4 - 0.005050998 * u5 + 0.0083572073 * u6;
	}
	else if ( y < 1700.0 )
	{
		t = y - 1600.0;
		t2 = t * t;
		t3 = t2 * t;
		dt = 120.0 - 0.9808 * t - 0.01532 * t2 + t3 / 7129.0;
	}
	else if ( y < 1800.0 )
	{
		t = y - 1700.0;
		t2 = t * t;
		t3 = t2 * t;
		t4 = t3 * t;
		dt = 8.83 + 0.1603 * t - 0.0059285 * t2 + 0.00013336 * t3 - t4 / 1174000.0;
	}
	else if ( y < 1860.0 )
	{
		t = y - 1800.0;
		t2 = t * t;
		t3 = t2 * t;
		t4 = t3 * t;
		t5 = t4 * t;
		t6 = t5 * t;
		t7 = t6 * t;
		dt = 13.72 - 0.332447 * t + 0.0068612 * t2 + 0.0041116 * t3 - 0.00037436 * t4
		+ 0.0000121272 * t5 - 0.0000001699 * t6 + 0.000000000875 * t7;
	}
	else if ( y < 1900.0 )
	{
		t = y - 1860.0;
		t2 = t * t;
		t3 = t2 * t;
		t4 = t3 * t;
		t5 = t4 * t;
		dt = 7.62 + 0.5737 * t - 0.251754 * t2 + 0.01680668 * t3
		- 0.0004473624 * t4 + t5 / 233174.0;
	}
	else if ( y < 1920.0 )
	{
		t = y - 1900.0;
		t2 = t * t;
		t3 = t2 * t;
		t4 = t3 * t;
		dt = -2.79 + 1.494119 * t - 0.0598939 * t2 + 0.0061966 * t3 - 0.000197 * t4;
	}
	else if ( y < 1940.0 )
	{
		t = y - 1920.0;
		t2 = t * t;
		t3 = t2 * t;
		dt = 21.20 + 0.84493 * t - 0.076100 * t2 + 0.0020936 * t3;
	}
	else if ( y < 1960.0 )
	{
		t = y - 1950.0;
		t2 = t * t;
		t3 = t2 * t;
		dt = 29.07 + 0.407 * t - t2 / 233.0 + t3 / 2547.0;
	}
	else if ( y < 1985.0 )
	{
		t = y - 1975.0;
		t2 = t * t;
		t3 = t2 * t;
		dt = 45.45 + 1.067 * t - t2 / 260.0 - t3 / 718.0;
	}
	else if ( y < 2005.0 )
	{
		t = y - 2000.0;
		t2 = t * t;
		t3 = t2 * t;
		t4 = t3 * t;
		t5 = t4 * t;
		dt = 63.86 + 0.3345 * t - 0.060374 * t2 + 0.0017275 * t3 + 0.000651814 * t4
		+ 0.00002373599 * t5;
	}
	else if ( y < 2050.0 )
	{
		t = y - 2000.0;
		t2 = t * t;
		// 		The new formula below the commented-out original better fits actual
		//		published Delta T data from 2000 to 2015, and still maintains the
		//		projected value of 93 seconds at year 2050.
		//		dt = 62.92 + 0.32217 * t + 0.005589 * t2;
		dt = 63.83 + 0.1102 * t + 0.009464 * t2;
	}
	else if ( y < 2150.0 )
	{
		u = ( y - 1820.0 ) / 100.0;
		dt = -20 + 32 * u * u - 0.5628 * ( 2150.0 - y );
	}
	else
	{
		u = ( y - 1820.0 ) / 100.0;
		dt = -20.0 + 32.0 * u * u;
	}
	
	return ( dt );
}

/*** AAGreenwichMeanSiderealTime *******************************************************/

double AAGreenwichMeanSiderealTime ( double jd )
{
	double jd0, gmst, t, t1, t2, t3;
	
	jd0 = floor ( jd - 0.5 ) + 0.5;
	
	t = ( jd0 - 2451545.0 ) / 36525.0;
	if ( t > 100.0 )
		t1 = 100.0;
	else if ( t < -100.0 )
		t1 = -100.0;
	else
		t1 = t;
	
	t2 = t1 * t1;
	t3 = t2 * t1;
	
	gmst = 24110.54841 + 8640184.812866 * t + 0.093104 * t2 - 0.0000062 * t3;
	gmst /= 86400.0;
	gmst = gmst - floor ( gmst ) + ( jd - jd0 ) * SIDEREAL_PER_SOLAR;
	
	return ( Mod2Pi ( gmst * TWO_PI ) );
}

/*** AALocalMeanSiderealTime *******************************************************/

double AALocalMeanSiderealTime ( double jd, double lon )
{
	return ( Mod2Pi ( AAGreenwichMeanSiderealTime ( jd ) + lon ) );
}

/*** AASemiDiurnalArc *****************************************************/
   
double AASemiDiurnalArc ( double a, double d, double f )
{
	double ch;
	
	ch = ( sin ( a ) - sin ( d ) * sin ( f ) ) / ( cos ( d ) * cos ( f ) );
	
	if ( ch >= 1.0 )
		return ( 0.0 );
		
	if ( ch <= -1.0 )
		return ( PI );

 	return ( acos ( ch ) );
}

/*** AARiseSetTime *****************************************************/

double AARiseSetTime ( double ra, double dec, double jd, int sign,
double lon, double lat, double alt )
{
	double ha, lst, theta;
	
	/*** Compute the object's hour angle when it is on the horizon.
	     This is the difference between the object's right ascension
	     and the line of right ascension which is on the meridian
	     at the time the object appears on the horizon. ***/

	ha = AASemiDiurnalArc ( alt, dec, lat );

	/*** If the object never sets, return infinity; if it never rises,
	     return negative infinity. ***/
	     
	if ( ha == PI && sign != 0 )
		return ( INFINITY );
	
	if ( ha == 0.0 )
		return ( -INFINITY );
	
	/*** Compute the local sidereal time.  This is the same as the circle
	     of right ascension which is currently on the local meridian. ***/
	
	lst = AALocalMeanSiderealTime ( jd, lon );
	
	/*** Now compute the angular distance that the earth needs to turn
	     through to make the star appear on the horizon.  If the angle
	     is outside the range -PI to +PI, reduce it to that range. ***/
	
	theta = ( ra - lst + sign * ha );
	
	while ( theta > PI )
		theta = theta - TWO_PI;
		
	while ( theta < -PI )
		theta = theta + TWO_PI;
	
	/*** Obtain the time of rising or setting by adding the amount of time
	     the earth takes to rotate through the angle calculated above to the
	     current time. ***/ 
		
	jd = jd + theta / TWO_PI / SIDEREAL_PER_SOLAR;

	return ( jd );
}

/*** AARiseSetTimeSearch *****************************************************/

double AARiseSetTimeSearch ( AARiseSetSearchProcPtr proc, void *data, double jd,
int sign, double lon, double lat, double alt, double precision, int imax )
{
	double		ra, dec, oldjd;
	int			i = 0;
		     
	do
	{
		oldjd = jd;
		(*proc) ( jd, &ra, &dec, data );
		jd = AARiseSetTime ( ra, dec, jd, sign, lon, lat, alt );
		i = i + 1;
	}
	while ( fabs ( jd - oldjd ) > precision && fabs ( jd ) != INFINITY && i < imax );
			
	return ( jd );
}

/*** AADailyRiseSetTimeSearch *************************************************/

double AADailyRiseSetTimeSearch ( AARiseSetSearchProcPtr proc, void *data,
double jd, int sign, double zone, double lon, double lat, double alt,
double precis, int imax )
{
	double	start, end;
	
	/*** Find the julian dates that correspond to the start and end
	     of the local day. ***/
	     
	start = floor ( jd - 0.5 + zone ) + 0.5 - zone;
	end = start + 1.0;

	/*** Search for the object's exact rise/set time, starting from
	     the middle of the local day. ***/
	     
	jd = AARiseSetTimeSearch ( proc, data, start + 0.5, sign,
	     lon, lat, alt, precis, imax );
	
	/*** If the resulting rise/set time is after the end of the day,
	     start searching again from the middle of the previous day;
	     similarly, if the resulting rise/set time is before the start
	     of the current day, start searching again from the middle
	     of the following day. ***/

	if ( jd > end )
		jd = AARiseSetTimeSearch ( proc, data, start - 0.5, sign,
		     lon, lat, alt, precis, imax );
	else if ( jd < start )
		jd = AARiseSetTimeSearch ( proc, data, end + 0.5, sign,
		     lon, lat, alt, precis, imax );

	/*** If the resulting rise/set time is still before the beginning or
	     after the end of the local day, the object does not rise or set
	     at all on that day.  Return positive infinity to indicate the
	     object does not set or negative infinity to indicate that it
	     does not rise. ***/
	
	if ( jd > end || jd < start )
	{
		if ( sign == -1 )
			jd = -INFINITY;
		else
			jd = INFINITY;
	}
				
	return ( jd );
}

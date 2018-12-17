/*** COPYRIGHT NOTICE *********************************************************
 
	Copyright (c) 1992-2016 Southern Stars Group, LLC.  All Rights Reserved.
 
******************************************************************************/

#include "NMatrix.h"

/****************************  NCreateVector  *********************************/

void *NCreateVector ( size_t size, int n )
{
	void *vector;
	
	vector = (void *) malloc ( size * n );
	if ( vector == NULL )
		return ( NULL );
		
	memset ( vector, 0, size * n );
	return ( vector );
}

/****************************  NCreateMatrix  *********************************/

void **NCreateMatrix ( size_t size, int m, int n )
{
	int	i;
	void	**matrix;
	
	if ( m == 0 || n == 0 )
		return ( NULL );
	
	matrix = (void **) NCreateVector ( sizeof ( void * ), m + 1 );
	if ( matrix == NULL )
		return ( NULL );
			
	for ( i = 0; i < m; i++ )
	{
		matrix[i] = NCreateVector ( size, n );
		if ( matrix[i] == NULL )
		{
			NDestroyMatrix ( matrix );
			return ( NULL );
		}
	}
	
	return ( matrix );
}

/****************************  NDestroyVector  *******************************/

void NDestroyVector ( void *vector )
{
	free ( vector );
}

/***************************  NDestroyMatrix  *********************************/

void NDestroyMatrix ( void **matrix )
{
	int i;
	
	if ( matrix )
		for ( i = 0; matrix[i]; i++ )
			free ( matrix[i] );
		
	free ( matrix );
}

/*****************************  NCopyVector  ********************************/

void NCopyVector ( double *b, double *a, int n )
{
	int i;
	
	for ( i = 0; i < n; i++ )
		b[i] = a[i];
}

/*****************************  NCopyMatrix  ********************************/

void NCopyMatrix ( double **b, double **a, int m, int n )
{
	int i, j;
	
	for ( i = 0; i < m; i++ )
		for ( j = 0; j < n; j++ )
			b[i][j] = a[i][j];
}

/**************************  NInitializeMatrix  ******************************/

void NInitializeMatrix ( double **a, int m, int n, double x )
{
	int i, j;
	
	for ( i = 0; i < m; i++ )
		for ( j = 0; j < n; j++ )
			a[i][j] = x;
}

/**************************  NTransposeMatrix  *******************************/

void NTransposeMatrix ( double **b, double **a, int m, int n )
{
	int i, j;
	
	for ( i = 0; i < m; i++ )
		for ( j = 0; j < n; j++ )
			b[j][i] = a[i][j];
}

/****************************  NMultiplyMatrix  *****************************/

void NMultiplyMatrix ( double **c, double **a, double **b, int m, int n, int p )
{
	int i, j, k;
	
	for ( i = 0; i < m; i++ )
		for ( k = 0; k < p; k++ )
			for ( c[i][k] = 0.0, j = 0; j < n; j++ )
				c[i][k] += a[i][j] * b[j][k];
}

/****************************  NIdentityMatrix  ******************************/

void NIdentityMatrix ( double **a, int m, int n )
{
	int i, j;
	
	for ( i = 0; i < m; i++ )
		for ( j = 0; j < n; j++ )
			if ( i == j )
				a[i][j] = 1.0;
			else
				a[i][j] = 0.0;
}

/****************************  NCreateNormalEqns  *****************************/

int NCreateNormalEqns ( int m, int n, double ***a, double ***b )
{
	int i, j;
	
	*a = NMatrix ( double, m, m );
	if ( *a == NULL )
		return ( FALSE );
		
	*b = NMatrix ( double, m, n );
	if ( *b == NULL )
	{
		NDestroyMatrix ( (void **) *a );
		return ( FALSE );
	}
		
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < m; j++ )
			(*a)[i][j] = 0.0;
		
		for ( j = 0; j < n; j++ )
			(*b)[i][j] = 0.0;
	}
	
	return ( TRUE );	
}

/************************  NAugmentNormalEqns  ****************************/

void NAugmentNormalEqns ( int m, int n, double *x, double *y, double **a, double **b )
{
	int i, j;
	
	for ( i = 0; i < m; i++ )
	{
		for ( j = 0; j < m; j++ )
			a[i][j] += x[i] * x[j];
		
		for ( j = 0; j < n; j++ )
			b[i][j] += x[i] * y[j];
	}
}

/***********************  NGaussJordanSolveMatrixEqn  **********************/
	
int NGaussJordanSolveMatrixEqn ( double **a, int m, double **b, int n )
{
	int		*indxc, *indxr, *ipiv;
	int		i, j, k, l, ll, irow = 0, icol = 0;
	double	big, temp, pivinv;
	
	/*** these arrays are used for bookkeeping on the pivoting ***/
		
	if ( ( indxr = NVector ( int, m ) ) == NULL )
		return ( FALSE );
		
	if ( ( indxc = NVector ( int, m ) ) == NULL )
	{
		free ( indxr );
		return ( FALSE );
	}
	
	if ( ( ipiv = NVector ( int, m ) ) == NULL )
	{
		free ( indxr );
		free ( indxc );
		return ( FALSE );
	}
	
	for ( j = 0; j < m; j++ )
		ipiv[j] = 0;

	/*** this is the main loop over the columns to be reduced. ***/
	
	for ( i = 0; i < m; i++ )
	{	
		big = 0.0;
		
		/*** this is the outer loop of the search for a pivot element ***/
		
		for ( j = 0; j < m; j++ )
		{
			if ( ipiv[j] != 1 )
			{
				for ( k = 0; k < m; k++ )
				{
					if ( ipiv[k] == 0 )
					{
						if ( ( temp = fabs ( a[j][k] ) >= big ) != 0.0 )
						{
							big = temp;
							irow = j;
							icol = k;
						}
					}
					else
					{
						if ( ipiv[k] > 1 )
						{
							free ( indxc );
							free ( indxr );
							free ( ipiv );
							return ( FALSE );
						}
					}
				}
			}
		}
		
		++( ipiv[icol] );
		
		/*** We now have the pivot element, so we interchange rows, if needed, to put
			the pivot element on the diagonal.  The columns are not physically inter-
			changed, only relabeled: indxc[i], the column of the ith pivot element, is
			the ith column that is reduced, while indxr[i] is the row in which that
			pivot element was originally located.  If indxr[i] ­ indxc[i] there is an
			implied column interchange.  With this form of bookkeeping, the solution
			matrix will end up the correct order, and the inverse matrix will be
			scrambled by columns. ***/
			
		if ( irow != icol )
		{
			for ( l = 0; l < m; l++ )
			{
				temp = a[irow][l];
				a[irow][l] = a[icol][l];
				a[icol][l] = temp;
			}
			
			for ( l = 0; l < n; l++ )
			{
				temp = b[irow][l];
				b[irow][l] = b[icol][l];
				b[icol][l] = temp;
			}
		}
			
		/*** We are now ready to divide the pivot row by the pivot element, located
			at irow and icol ***/
		
		indxr[i] = irow;
		indxc[i] = icol;
		
		if ( ( temp = a[icol][icol] ) == 0.0 )
		{
			free ( indxc );
			free ( indxr );
			free ( ipiv );
			return ( FALSE );
		}
		
		pivinv = 1.0 / temp;
		a[icol][icol] = 1.0;
		
		for ( l = 0; l < m; l++ )
			a[icol][l] *= pivinv;
			
		for ( l = 0; l < n; l++ )
			b[icol][l] *= pivinv;
			
		/*** next, we reduce the rows, except for the pivot one... ***/
		
		for ( ll = 0; ll < m; ll++ )
		{
			if ( ll != icol )
			{
				temp = a[ll][icol];
				a[ll][icol] = 0.0;
				
				for ( l = 0; l < m; l++ )
					a[ll][l] -= a[icol][l] * temp;
					
				for ( l = 0; l < n; l++ )
					b[ll][l] -= b[icol][l] * temp;
			}
		}
	}
	
	/*** This is the end of the main loop over columns of the reduction.  It only remains
		to unscrable the solution in view of the column interchanges.  We do this by
		interchanging pairs of columns in the reverse order that the permutation was
		built up. ***/
		
	for ( l = m - 1; l >= 0; l-- )
	{
		if ( indxr[l] != indxc[l] )
			for ( k = 1; k < m; k++ )
			{
				temp = a[k][indxr[l]];
				a[k][indxr[l]] = a[k][indxc[l]];
				a[k][indxc[l]] = temp;
			}
	}
	
	/*** And we are done. ***/
	
	free ( indxr );
	free ( indxc );
	free ( ipiv );
	
	return ( TRUE );
}


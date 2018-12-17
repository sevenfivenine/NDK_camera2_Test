#ifndef NMATRIX_H
#define NMATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef TRUE
#define TRUE	1
#endif
	
#ifndef FALSE
#define	FALSE	0
#endif

/*** NCreateMatrix **********************************************************

	Allocates memory for a new matrix of arbitrary dimensions.
	
	void **NCreateMatrix ( size_t size, long m, long n )
	
	(size): size of an individual matrix element in bytes.
	   (m): desired number of rows for matrix.
	   (n): desired number of columns for matrix.
	
	The function returns an array of pointer to the matrix rows, or NULL
	on failure.  Memory for each matrix row is allocated separately, and
	all elements of each row are initialized to zero.

	Because this function can create a matrix with elements of any size,
	you will need to cast the matrix row pointer array it returns to the
	proper type.
	
	For example, the following code creates a 5 x 10 matrix with short
	integer elements, and another with double-precision floating-point
	elements:
	
	short **matrix1 = (short **) NewMatrix ( sizeof ( short ), 5, 10 );
	double **matrix2 = (double **) NewMatrix ( sizeof ( double ), 5, 10 );
	
	To refer to an individual element in the ith row and jth column, use
	the notation matrix[i][j].  Both the row and column number should be
	counted from zero, following standard C array notation.

	Reference: Numerical Recipes in C, pp. 16-20.
	
***************************************************************************/

void *NCreateVector ( size_t, int );
void **NCreateMatrix ( size_t, int, int );

#define NVector(type,n)		(type *)  NCreateVector ( sizeof ( type ), n )
#define NMatrix(type,m,n)	(type **) NCreateMatrix ( sizeof ( type ), m, n )

/*** NDestroyMatrix ********************************************************

	Frees memory for a matrix allocated with NewMatrix().
	
	void NDestroyMatrix ( void **matrix )

	(matrix): array of pointers to rows of matrix to free.
	
	This function returns nothing.  After this funtion is called, the
	matrix row pointer array given (matrix) will be invalid.  You should
	only use this function to free matrices created with NewMatrix()!!!

	Reference: Numerical Recipes in C, pp. 16-20.

*****************************************************************************/

void NDestroyVector ( void * );
void NDestroyMatrix ( void ** );

void NCopyVector ( double *, double *, int );
void NCopyMatrix ( double **, double **, int, int );

void NInitializeMatrix ( double **, int, int, double );
void NIdentityMatrix ( double **, int, int );

void NTransposeMatrix ( double **, double **, int, int );
void NMultiplyMatrix ( double **, double **, double **, int, int, int );

/*** NCreateNormalEqns ***************************************************
 
	Allocates and initializes normal equation matrices for a least-squares
	solution to an overdetermined linear system of equations.
	
	int NCreateNormalEqns ( int m, int n, double ***a, double ***b )
	
	(m): number of independent or left-hand-side (LHS) variables.
	(n): number of dependent or right-hand-side (RHS) variables.
	(a): recieves pointer to LHS normal equation matrix.
	(b): recieves pointer to RHS normal equation matrix.
	
	The function returns TRUE if successful and FALSE on failure.
	
	On a successful return, the variable (a) will point to a matrix of
	(m) rows by (m) columns, and (b) will point to a matrix of (m) rows
	by (n) columns.  Both matrices will be initialized to zero.
	
	The values required for the parameters (m) and (n) depend on the
	the system of equations whose least-squares solution you wish to
	determine.  The variable (m) represents the number of independent
	(LHS) variables in the equation, and (n) represents the number of
	dependent (RHS) variables.
 
	In general, a linear system of equations with M independent variables
	and N dependent variables takes the following form:
 
	y1 = a11 * x1 + a12 * x2 + ... + a1M * xM
	
	y2 = a21 * x2 + a22 * x2 + ... + a2M * xM
	
	.        .          .      .         .
	.        .          .       .        .
	.        .          .        .       .
	
	yN = aN1 * xN + aN2 * x2 + ... + aNM * xM
 
	Here, x1...xM are the independent variables, y1...yN are the dependent
	variables, and a11...aNM are the unknown coefficients whose values will
	be determined by least-squares solution.
 
	We can write this in vector-matrix form as y = x * A, where:
	
	y = | y1  y2  ... yN  | is the vector of dependent variables;
	
	x = | x1  x2  ... xM  | is the vector of independent variables;
	
	    | a11 a21 ... aN1 |
	    |                 |
	    | a12 a22 ... aN2 |
	    |                 |
	A = |  .   .  .    .  | is the matrix of unknown coefficients.
	    |  .   .   .   .  |
	    |  .   .    .  .  |
	    |                 |
	    | a1M a2M ... aNM |
 
	In order to solve for the matrix of unknown coefficients, you must have
	a set of "data points" for which the values of both the dependent and
	independent variables are known.  The system of equations with the
	values of the data points "plugged in" can then be written as:
	
	Y = X * A
 
	Here, each row of matrix X contains the values of the dependent variables
	x1...xM at a particular data point, and the corresponding row of matrix Y
	contains the values of the dependent variables y1...yN at that same point.
	If the number of known data points is P, then matrix X will contain P rows
	and M columns, while matrix Y will contain P rows and N columns.
	
	In order to solve for the matrix A, the number of data points (P) must be
	equal to or greater than the number of rows of the coefficient matrix (M),
	so that the number of known values equals or exceeds the number of unknown
	coefficients.
	
	 T       T
	X * Y = X * X * A
 
	  T    -1    T
	(X * X)  * (X * Y) = A
	
****************************************************************************/

int NCreateNormalEqns ( int, int, double ***, double *** );

/*** NAugmentNormalEqns ****************************************************
	 
	 Adds the values of a known data point to the normal equations.
	 
	 void NAugmentNormalEnqs ( int m, int n, double *x, double *y, double **a, double **b )
	 
	 (m): number of independent or right-hand-side (RHS) variables.
	 (n): number of dependent or left-hand-side (LHS) variables.
	 (x): values of the independent variables at a known data point.
	 (y): values of the dependent variables at a known data point.
	 (a): pointer to RHS normal equation matrix.
	 (b): pointer to LHS normal equation matrix.
	 
	 This function returns nothing.  It assumes that (a) points to a matrix
	 of (m) rows by (m) columns, and that (b) points to a matrix of (m) rows
	 by (n) columns.  Both matrices should be zero-initialized before calling
	 this function.  You can ensure that this has been done correctly by
	 using the function NewNormalEqns() to create both matrices.  After you
	 have done so, you should call AugmentNormalEqns() once for each known
	 data point you wish to include in the least-squares solution.
	 
	 After you have finished calling AugmentNormalEqns() for all of the
	 known data points, the normal equation matrices will contain the sums
	 of the products of the values of the dependent and independent variables
	 at each of the known data points.
	 
	 The RHS normal equation matrix (a) contains the sums:
	 
	 sum(x1*x1) sum(x1*x2) ... sum(x1*xM)
	 
	 sum(x2*x1) sum(x2*x2) ... sum(x2*xM)
	 
	 .          .      .       .
	 .          .       .      .
	 .          .        .     .
	 
	 sum(xM*x1) sum(xM*x2)     sum(xM*xM)
	 
	 The LHS normal equation matrix (b) contains the sums:
	 
	 sum(x1*y1) sum(x1*y2) ... sum(x1*yN)
	 
	 sum(x2*y1) sum(x2*y2) ... sum(x2*yN)
	 
	 .          .      .       .
	 .          .       .      .
	 .          .        .     .
	 
	 sum(xM*y1) sum(xM*y2) ... sum(xM*yN)
	 
	 Here, x1 ... xM are the independent variables, y1 ... yN are the
	 dependent variables, and sum(...) means "sum over all of the known
	 data points".
	 
	 At this point the original least-squares problem can be rewritten as
	 
	 Y = X * A
	 
	 where X is the RHS normal equation matrix (a), Y is the LHS normal
	 equation matrix (b).
	 
	 To solve the normal equations and obtain a least-squares best-fit
	 solution, you may call the routine NGaussJordanSolveMatrixEqn().
	 
***************************************************************************/

void NAugmentNormalEqns ( int, int, double *, double *, double **, double ** );

/*** NGaussJordanSolveMatrixEqn ********************************************

	Linear equation solution by Gauss-Jordan elimination.
	
	void NGaussJordanSolveMatrixEqn ( double **a, int m, double **b, int n )
	
	(a): input matrix of (m) by (m) elemnts.
	(b): input matrix of (m) rows and (n) columns.
	
	The function returns TRUE if it is successfully able to solve the
	matrix equation, and FALSE otherwise.
	
	On output, (a) is replaced by its inverse, and (b) is replaced by the
	corresponding set of solution vectors.
	
	Reference: Numerical Recipes in C, pp. 32-36
	
****************************************************************************/

int NGaussJordanSolveMatrixEqn ( double **, int, double **, int );

#ifdef __cplusplus
}
#endif
#endif // NMATRIX_H

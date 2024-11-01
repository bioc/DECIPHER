/****************************************************************************
 *                 Distance Between Encoded Integer Vectors                 *
 *                           Author: Erik Wright                            *
 ****************************************************************************/

/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

// for math functions
#include <math.h>

// for div_t div
#include <stdlib.h>

// DECIPHER header file
#include "DECIPHER.h"

// decode integer vector
// NOTE:  Zeros the elements of SEXP x
SEXP intDist(SEXP x, SEXP levels, SEXP bins, SEXP maxBins, SEXP numRows, SEXP totRows, SEXP power)
{	
	int i, j, k, pos;
	int *X = INTEGER(x);
	int l = asInteger(levels);
	int b = asInteger(bins);
	int mB = asInteger(maxBins);
	int n = asInteger(numRows);
	int cols = length(x)/n;
	int t = asInteger(totRows);
	double p = asReal(power);
	div_t Y;
	
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, b, n));
	int *rans;
	rans = INTEGER(ans);
	
	// zero memory
	for (i = 0; i < b*n; i++)
		rans[i] = 0;
	
	// convert encoded vector to integers
	for (i = 0; i < n; i++) {
		for (j = 0; j < cols; j++) {
			pos = mB*j + b*i;
			while (X[j*n + i] > 0) {
				Y = div(X[j*n + i], l);
				rans[pos] = Y.rem;
				X[j*n + i] = Y.quot;
				pos++;
			}
		}
	}
	
	// calculate the difference in area between integer vectors
	double avg = 0, avg_j, diff;
	int *skip = R_Calloc(n, int); // initialized to zero
	int weight;
	for (i = 0; i < n; i++) {
		if (skip[i] != 0)
			continue;
		weight = 1;
		for (j = i + 1; j < n; j++) {
			avg_j = 0;
			for (k = 0; k < b; k++) {
				diff = rans[i*b + k] - rans[j*b + k];
				if (diff > 0) {
					avg_j += diff;
				} else {
					avg_j -= diff;
				}
			}
			if (avg_j == 0) {
				skip[j] = 1;
				weight++;
			} else {
				if (p == 1) { // L1-norm
					avg += avg_j*weight;
				} else { // LP-norm
					avg += pow(avg_j/(double)(b*(l - 1)), p)*weight;
				}
			}
		}
	}
	R_Free(skip);
	
	SEXP dist;
	PROTECT(dist = allocVector(REALSXP, 1));
	double *d;
	d = REAL(dist);
	
	// calculate the weighted mean including missing dists
	// sum of distances divided by half triangle of matrix
	if (p == 1) { // average L1-norm
		d[0] = avg/(double)(b*(l - 1))/(double)((t*t - t)/2);
	} else { // average LP-norm
		d[0] = pow(avg, 1/p)/(double)((t*t - t)/2);
	}
	
	UNPROTECT(2);
	
	return dist;
}

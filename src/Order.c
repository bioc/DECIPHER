/****************************************************************************
 *                       Obtain Ordering of a Vector                        *
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

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

// for math functions
#include <math.h>

// for calloc/free
#include <stdlib.h>

// DECIPHER header file
#include "DECIPHER.h"

// order x (positive integers only)
SEXP radixOrder(SEXP x, SEXP ascending, SEXP keepNAs)
{
	int i, k, o;
	int l = length(x);
	int *v = INTEGER(x);
	int keep = asInteger(keepNAs); // whether to keep NAs when ordering
	int s = asInteger(ascending); // start of index
	
	int *order = (int *) malloc(l*sizeof(int)); // thread-safe on Windows
	int m = 1;
	int NAs = 0;
	if (keep) {
		for (i = 0; i < l; i++) {
			order[i] = i;
			if (v[i] > m)
				m = v[i];
		}
	} else {
		for (i = 0; i < l; i++) {
			order[i] = i;
			if (v[i] == NA_INTEGER) {
				NAs++;
			} else if (v[i] > m) {
				m = v[i];
			}
		}
	}
	m = (int)ceil(log2((double)(m + 1)));
	
	int R; // size of radix key
	i = 0;
	do {
		i++;
		R = (int)ceil((double)m/(double)i);
	} while(R > 8);
	m = i;
	int count = 1 << R; // 2^R
	
	unsigned int mask = 1;
	for (i = 1; i < R; i++)
		mask |= 1 << i; // R ones
	
	// least significant digit first
	for (k = 0; k < m; k++) {
		// subset bits in radix k
		int *counts = (int *) calloc(count, sizeof(int)); // initialized to zero (thread-safe on Windows)
		o = k*R;
		for (i = 0; i < l; i++)
			counts[(v[i] >> o) & mask]++;
		
		// cumulative sum from zero
		int one = 0;
		for (i = 1; i < count; i++) {
			counts[i] = counts[i - 1] + counts[i];
			int two = counts[i - 1];
			counts[i - 1] = one;
			one = two;
		}
		counts[count - 1] = one;
		
		// move orders
		int *temp = (int *) malloc(l*sizeof(int)); // thread-safe on Windows
		for (i = 0; i < l; i++)
			temp[counts[(v[order[i]] >> o) & mask]++] = order[i]; // slow when cache misses
		free(counts);
		
		// replace orders
		free(order);
		order = temp;
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l - NAs));
	int *rans = INTEGER(ans);
	if (keep) {
		if (s) {
			for (i = 0; i < l; i++)
				rans[i] = order[i] + 1;
		} else {
			for (i = 0; i < l; i++)
				rans[i] = order[l - i - 1] + 1;
		}
	} else {
		k = 0;
		if (s) {
			for (i = 0; i < l; i++)
				if (v[order[i]] != NA_INTEGER)
					rans[k++] = order[i] + 1;
		} else {
			for (i = 0; i < l; i++)
				if (v[order[l - i - 1]] != NA_INTEGER)
					rans[k++] = order[l - i - 1] + 1;
		}
	}
	free(order);
	
	UNPROTECT(1);
	
	return ans;
}

// dereplicate x using its ordering
SEXP dereplicate(SEXP x, SEXP o)
{
	int *X = INTEGER(x);
	int *O = INTEGER(o);
	int l = length(x);
	
	int *numbers = (int *) malloc(l*sizeof(int)); // thread-safe on Windows
	int *counts = (int *) calloc(l, sizeof(int)); // initialized to zero (thread-safe on Windows)
	
	int count = 1;
	int i = 0;
	int j = 0;
	int k = 1;
	int prev, curr;
	if (l > 0)
		prev = X[O[j] - 1];
	while (k < l) {
		curr = X[O[k] - 1];
		if (curr == prev) {
			count++;
		} else {
			numbers[i] = O[j];
			counts[i] = count;
			i++;
			count = 1;
			j = k;
			prev = curr;
		}
		k++;
	}
	if (l > 0) {
		numbers[i] = O[j];
		counts[i] = count;
		i++;
	}
	
	SEXP ans1, ans2;
	PROTECT(ans1 = allocVector(INTSXP, i));
	PROTECT(ans2 = allocVector(INTSXP, i));
	int *rans1 = INTEGER(ans1);
	int *rans2 = INTEGER(ans2);
	
	for (j = 0; j < i; j++) {
		rans1[j] = numbers[j];
		rans2[j] = counts[j];
	}
	free(numbers);
	free(counts);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	
	UNPROTECT(3);
	
	return ret_list;
}

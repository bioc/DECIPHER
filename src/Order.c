/****************************************************************************
 *                       Obtain Ordering of a Vector                        *
 *                           Author: Erik Wright                            *
 ****************************************************************************/

// for OpenMP parallel processing
#ifdef _OPENMP
#include <omp.h>
#undef match
#endif

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
SEXP radixOrder(SEXP x, SEXP ascending, SEXP keepNAs, SEXP nThreads)
{
	int p, b, j, o, key, m = 1;
	R_xlen_t i, k, *swap, l = xlength(x);
	int *v = INTEGER(x);
	int s = asInteger(ascending); // start of index
	int keep = asInteger(keepNAs); // whether to keep NAs when ordering
	int nthreads = asInteger(nThreads);
	
	R_xlen_t *order = (R_xlen_t *) malloc(l*sizeof(R_xlen_t)); // thread-safe on Windows
	if (keep) {
		k = l;
		for (i = 0; i < l; i++) {
			order[i] = i;
			if (v[i] > m)
				m = v[i];
		}
	} else {
		k = 0;
		for (i = 0; i < l; i++) {
			if (v[i] != NA_INTEGER) {
				order[k++] = i;
				if (v[i] > m)
					m = v[i];
			}
		}
	}
	
	// select key size
	if (l >= 8192) {
		key = 8;
	} else {
		key = 4;
	}
	
	m = (int)ceil(log2((double)(m + 1)));
	
	int R; // size of radix key
	j = 0;
	do {
		j++;
		R = (int)ceil((double)m/(double)j);
	} while(R > key);
	m = j;
	int count = 1 << R; // 2^R
	
	unsigned int mask = 1;
	for (j = 1; j < R; j++)
		mask |= 1 << j; // R ones
	
	order = (R_xlen_t *) realloc(order, k*sizeof(R_xlen_t)); // thread-safe on Windows
	R_xlen_t *temp = (R_xlen_t *) malloc(k*sizeof(R_xlen_t)); // thread-safe on Windows
	R_xlen_t *counts = (R_xlen_t *) malloc(count*sizeof(R_xlen_t)); // thread-safe on Windows
	
	// sort most significant key
	for (p = 0; p < count; p++)
		counts[p] = 0;
	o = (--m)*R;
	for (i = 0; i < k; i++)
		counts[(v[order[i]] >> o) & mask]++;
	
	R_xlen_t *bounds;
	if (m != 0) {
		// record bin boundaries
		bounds = (R_xlen_t *) malloc((count + 1)*sizeof(R_xlen_t)); // thread-safe on Windows
		bounds[0] = 0;
		for (p = 0; p < count; p++)
			bounds[p + 1] = bounds[p] + counts[p];
	}
	
	// cumulative sum from zero
	R_xlen_t one = 0;
	for (p = 1; p < count; p++) {
		counts[p] = counts[p - 1] + counts[p];
		R_xlen_t two = counts[p - 1];
		counts[p - 1] = one;
		one = two;
	}
	counts[count - 1] = one;
	
	if (m == 0) {
		// move orders
		for (i = 0; i < k; i++)
			temp[counts[(v[order[i]] >> o) & mask]++] = order[i];
		free(counts);
		swap = order;
		order = temp;
		temp = swap;
	} else {
		// move orders
		int *unsorted = (int *) calloc(count, sizeof(int)); // initialized to zero (thread-safe on Windows)
		int *last = (int *) calloc(count, sizeof(int)); // initialized to zero (thread-safe on Windows)
		for (i = 0; i < k; i++) {
			int val = v[order[i]];
			int index = (val >> o) & mask;
			temp[counts[index]++] = order[i];
			unsorted[index] = unsorted[index] || (val < last[index]);
			last[index] = val;
		}
		free(counts);
		free(last);
		
		// record bins that can be sorted
		int *bin = (int *) malloc(count*sizeof(int)); // thread-safe on Windows
		int bins = 0;
		for (b = 0; b < count; b++)
			if (unsorted[b])
				bin[bins++] = b;
		free(unsorted);
		
		if (bins) {
			// copy current order
			for (i = 0; i < k; i++)
				order[i] = temp[i];
			
			int NTHREADS;
			if (nthreads > 1) {
				NTHREADS = k/100000;
				if (NTHREADS < 1) {
					NTHREADS = 1;
				} else if (NTHREADS > nthreads) {
					NTHREADS = nthreads;
				}
			} else {
				NTHREADS = 1;
			}
			
			// sort least significant keys
			for (j = 0; j < m; j++) {
				o = j*R;
				#ifdef _OPENMP
				#pragma omp parallel num_threads(NTHREADS)
				{
				#endif
					R_xlen_t *counts = (R_xlen_t *) malloc(count*sizeof(R_xlen_t)); // thread-safe on Windows
					
					#ifdef _OPENMP
					#pragma omp for private(b,i) schedule(dynamic)
					#endif
					for (b = 0; b < bins; b++) { // each bin
						// count binned values
						for (p = 0; p < count; p++)
							counts[p] = 0;
						
						for (i = bounds[bin[b]]; i < bounds[bin[b] + 1]; i++)
							counts[(v[order[i]] >> o) & mask]++; // slow when cache misses
						
						// cumulative sum from zero
						R_xlen_t one = 0;
						for (p = 1; p < count; p++) {
							counts[p] = counts[p - 1] + counts[p];
							R_xlen_t two = counts[p - 1];
							counts[p - 1] = one;
							one = two;
						}
						counts[count - 1] = one;
						
						// shift start of bins
						for (p = 0; p < count; p++)
							counts[p] += bounds[bin[b]];
						
						// move orders
						for (i = bounds[bin[b]]; i < bounds[bin[b] + 1]; i++)
							temp[counts[(v[order[i]] >> o) & mask]++] = order[i]; // slow when cache misses
					}
					free(counts);
				#ifdef _OPENMP
				}
				#endif
				
				// swap orders
				swap = order;
				order = temp;
				temp = swap;
			}
		} else {
			swap = order;
			order = temp;
			temp = swap;
		}
		free(bounds);
		free(bin);
	}
	free(temp);
	
	SEXP ans;
	if (l > 2147483647) {
		PROTECT(ans = allocVector(REALSXP, k));
		double *rans = REAL(ans);
		if (s) {
			for (i = 0; i < k; i++)
				rans[i] = order[i] + 1;
		} else {
			for (i = 0; i < k; i++)
				rans[i] = order[k - i - 1] + 1;
		}
	} else {
		PROTECT(ans = allocVector(INTSXP, k));
		int *rans = INTEGER(ans);
		if (s) {
			for (i = 0; i < k; i++)
				rans[i] = order[i] + 1;
		} else {
			for (i = 0; i < k; i++)
				rans[i] = order[k - i - 1] + 1;
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

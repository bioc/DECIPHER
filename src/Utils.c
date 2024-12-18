/****************************************************************************
 *                       Low-level Utility Functions                        *
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

// for calloc/free and qsort
#include <stdlib.h>

// DECIPHER header file
#include "DECIPHER.h"

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// first matches of x[z...] == y[1]
SEXP multiMatch(SEXP x, SEXP y, SEXP z)
{
	int i, size = length(x);
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *u = INTEGER(z);
	int start = -1, stop = -1;
	
	for (i = *u - 1; i < size; i++) {
		if (v[i] == w[0]) {
			start = i;
			stop = i;
			for (i = start + 1; i < size; i++) {
				if (v[i] == w[0]) {
					stop = i;
				} else {
					break;
				}
			}
			break;
		}
	}
	
	SEXP ans;
	if (start != -1) {
		PROTECT(ans = allocVector(INTSXP, stop - start + 1));
		int *rans = INTEGER(ans);
		for (i = start; i <= stop; i++) {
			rans[i - start] = i + 1;
		}
	} else {
		PROTECT(ans = allocVector(INTSXP, 0));
	}
	UNPROTECT(1);
	
	return ans;
}

// first matches of x[z...] >= y[1]
SEXP multiMatchUpper(SEXP x, SEXP y, SEXP z)
{
	int i, size = length(x);
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *u = INTEGER(z);
	int start = -1, stop = -1;
	
	for (i = *u - 1; i < size; i++) {
		if (v[i] >= w[0]) {
			start = i;
			stop = i;
			for (i = start + 1; i < size; i++) {
				if (v[i] == v[start]) {
					stop = i;
				} else {
					break;
				}
			}
			break;
		}
	}
	
	SEXP ans;
	if (start != -1) {
		PROTECT(ans = allocVector(INTSXP, stop - start + 1));
		int *rans = INTEGER(ans);
		for (i = start; i <= stop; i++) {
			rans[i - start] = i + 1;
		}
	} else {
		PROTECT(ans = allocVector(INTSXP, 0));
	}
	UNPROTECT(1);
	
	return ans;
}

// first matches of x[...z] <= y[1]
SEXP multiMatchLower(SEXP x, SEXP y, SEXP z)
{
	int i, size = length(x);
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *u = INTEGER(z);
	int start = -1, stop = -1;
	
	for (i = *u - 1; i >= 0; i--) {
		if (v[i] <= w[0]) {
			start = i;
			stop = i;
			for (i = start + 1; i < size; i++) {
				if (v[i] == v[start]) {
					stop = i;
				} else {
					break;
				}
			}
			break;
		}
	}
	
	SEXP ans;
	if (start != -1) {
		PROTECT(ans = allocVector(INTSXP, stop - start + 1));
		int *rans = INTEGER(ans);
		for (i = start; i <= stop; i++) {
			rans[i - start] = i + 1;
		}
	} else {
		PROTECT(ans = allocVector(INTSXP, 0));
	}
	UNPROTECT(1);
	
	return ans;
}

// index of first non-NA elements
SEXP multiMatchCharNotNA(SEXP x)
{
	int i, size = length(x);
	int stop = 0;
	
	for (i = 0; i < size; i++) {
		if (STRING_ELT(x, i) != NA_STRING) {
			stop = i + 1;
		} else {
			break;
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, stop));
	int *rans = INTEGER(ans);
	for (i = 0; i < stop; i++) {
		rans[i] = i + 1;
	}
	
	UNPROTECT(1);
	
	return ans;
}

// same as x %in% y for ordered integer vectors
SEXP intMatch(SEXP x, SEXP y)
{
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int i, j;
	int size_x = length(x);
	int size_y = length(y);
	
	SEXP ans;
	PROTECT(ans = allocVector(LGLSXP, size_x));
	int *rans = INTEGER(ans);
	
	int s = 0;
	for (i = 0; i < size_x; i++) {
		rans[i] = 0;
		for (j = s; j < size_y; j++) {
			if (v[i] == w[j]) {
				rans[i] = 1;
				break;
			} else if (v[i] < w[j]) {
				break;
			}
		}
		s = j;
	}
	
	UNPROTECT(1);
	
	return ans;
}

// matrix of d[i, j] = length x[i] %in% y[j] / min(length)
// requires a list of ordered integers
SEXP matchListsDual(SEXP x, SEXP y, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	int i, j, size_x = length(x), size_y = length(y), before, v, *rPercentComplete;
	int o, p, start, lx, ly, *X, *Y, count;
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size_x, size_y));
	double *rans = REAL(ans);
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	int nthreads = asInteger(nThreads);
	
	if (v) { // initialize progress variables
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	int **py = (int **) malloc(size_y*sizeof(int *)); // thread-safe on Windows
	int *l = (int *) malloc(size_y*sizeof(int)); // thread-safe on Windows
	for (i = 0; i < size_y; i++) {
		py[i] = INTEGER(VECTOR_ELT(y, i));
		l[i] = length(VECTOR_ELT(y, i));
	}
	
	for (i = 0; i < size_x; i++) {
		X = INTEGER(VECTOR_ELT(x, i));
		lx = length(VECTOR_ELT(x, i));
		#ifdef _OPENMP
		#pragma omp parallel for private(j,o,p,start,count,Y,ly) schedule(guided) num_threads(nthreads)
		#endif
		for (j = 0; j < size_y; j++) {
			Y = py[j];
			ly = l[j];
			
			if (lx > 0 && ly > 0) {
				int first = -1;
				int last = -1;
				for (o = 0; o < lx; o++) {
					if (X[o] >= Y[0]) {
						first = o;
						break;
					}
				}
				if (first == -1) { // no overlap
					*(rans + j*size_x + i) = NA_REAL;
					continue;
				}
				
				for (o = lx - 1; o >= 0; o--) {
					if (X[o] <= Y[ly - 1]) {
						last = o;
						break;
					}
				}
				if (last == -1) { // no overlap
					*(rans + j*size_x + i) = NA_REAL;
					continue;
				}
				
				count = 0;
				start = 0;
				for (o = first; o <= last; o++) {
					for (p = start; p < ly; p++) {
						if (X[o] == Y[p]) {
							count++;
							start = p + 1;
							break;
						} else if (Y[p] > X[o]) {
							break;
						}
					}
				}
				
				if (lx > ly) {
					*(rans + j*size_x + i) = (double)count/(double)ly;
				} else {
					*(rans + j*size_x + i) = (double)count/(double)lx;
				}
			} else {
				*(rans + j*size_x + i) = NA_REAL;
			}
		}
		
		if (v) {
			// print the percent completed so far
			*rPercentComplete = floor(100*(double)(i + 1)/size_x);
			
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	free(py);
	free(l);
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}

// dist matrix = ordered matches (x[i], x[j]) / min(length)
// requires a list of unsorted integers retaining the order of x
SEXP matchOrder(SEXP x, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	R_xlen_t i, j, size_x = xlength(x);
	int before, v, *rPercentComplete;
	int lx, ly, *X, *Y;
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, size_x*(size_x - 1)/2));
	double *rans = REAL(ans);
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	int nthreads = asInteger(nThreads);
	
	if (v) { // initialize progress variables
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	int **px = (int **) malloc(size_x*sizeof(int *)); // thread-safe on Windows
	int *l = (int *) malloc(size_x*sizeof(int)); // thread-safe on Windows
	for (i = 0; i < size_x; i++) {
		px[i] = INTEGER(VECTOR_ELT(x, i));
		l[i] = length(VECTOR_ELT(x, i));
	}
	
	for (i = 0; i < size_x; i++) {
		X = px[i];
		lx = l[i];
		#ifdef _OPENMP
		#pragma omp parallel for private(j,Y,ly) schedule(guided) num_threads(nthreads)
		#endif
		for (j = i + 1; j < size_x; j++) {
			Y = px[j];
			ly = l[j];
			
			int link_x = -1, link_y = -1; // last established link between X and Y
			int matches = 0; // running number of matches
			int pos_x, pos_y, off; // current position
			int offset = 1; // distance offset from link
			int delta, forward;
			while (offset + link_x < lx && offset + link_y < ly) { // within sequence
				pos_y = link_y + 1;
				pos_x = link_x + offset;
				
				if (matches) {
					if (forward) {
						for (off = 1; off <= offset; off++, pos_x--, pos_y++) {
							if (X[pos_x] == Y[pos_y]) {
								link_x = pos_x;
								link_y = pos_y;
								offset = 0;
								matches++;
							}
						}
					} else { // backward
						for (off = 1; off <= offset; off++, pos_x--, pos_y++) {
							if (X[lx - pos_x - 1] == Y[ly - pos_y - 1]) {
								link_x = pos_x;
								link_y = pos_y;
								offset = 0;
								matches++;
							}
						}
					}
				} else {
					for (off = 1; off <= offset; off += delta, pos_x -= delta, pos_y += delta) {
						if (X[pos_x] == Y[pos_y]) {
							link_x = pos_x;
							link_y = pos_y;
							offset = 0;
							matches++;
							forward = 1;
							break;
						} else if (X[lx - pos_x - 1] == Y[ly - pos_y - 1]) {
							link_x = pos_x;
							link_y = pos_y;
							offset = 0;
							matches++;
							forward = 0; // backward
							break;
						}
						delta = (off < 10) ? 1 : (off/5);
					}
				}
				
				offset++;
			}
			
			if (lx < ly) {
				rans[size_x*i - i*(i + 1)/2 + j - i - 1] = 1 - (double)matches/(double)lx;
			} else {
				rans[size_x*i - i*(i + 1)/2 + j - i - 1] = 1 - (double)matches/(double)ly;
			}
		}
		
		if (v) {
			// print the percent completed so far
			//*rPercentComplete = floor(100*((double)i/((double)size_x - 1)));
			//*rPercentComplete = floor(100*(double)((i + 1)*size_x+(i + 1))/((size_x - 1)*size_x+(size_x - 1)));
			*rPercentComplete = floor(100*(double)(2*size_x - 2 - i)*(i + 1)/((size_x - 1)*size_x));
			
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	free(px);
	free(l);
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}

// returns shared ranges between pairs in two unordered (gapped) lists
SEXP matchRanges(SEXP x, SEXP y, SEXP wordSize, SEXP maxLength, SEXP threshold)
{
	int i, j, size_x = length(x), size_y = length(y), size;
	int lx, ly, *X, *Y, *X_pos, *Y_pos, wS, *rans;
	int l = asInteger(maxLength);
	double thresh = asReal(threshold);
	wS = asInteger(wordSize);
	int *bits = R_Calloc(l*2, int); // initialized to zero
	
	if (size_x > size_y) {
		size = size_x;
	} else {
		size = size_y;
	}
	
	for (i = 0; i < size; i += 2) {
		j = i;
		X = INTEGER(VECTOR_ELT(x, i));
		Y = INTEGER(VECTOR_ELT(y, j));
		X_pos = INTEGER(VECTOR_ELT(x, i + 1));
		Y_pos = INTEGER(VECTOR_ELT(y, j + 1));
		lx = length(VECTOR_ELT(x, i));
		ly = length(VECTOR_ELT(y, j));
		
		int link_x = -1, link_y = -1; // last established link between X and Y
		int pos_x, pos_y; // current position
		
		int offset = 1; // distance offset from link
		while (offset <= (lx - link_x + ly - link_y - 2)) { // within sequence
			for (int off = 1; off <= offset; off++) {
				pos_y = link_y + off;
				pos_x = link_x + offset - off + 1;
				if (pos_y < ly &&
					pos_x < lx &&
					X[pos_x] == Y[pos_y] &&
					X[pos_x] != NA_INTEGER) {
					if (pos_x == (link_x + 1) &&
						pos_y == (link_y + 1)) {
						if (*(bits + l + X_pos[pos_x] - 1) == 0) { // new anchor
							*(bits + l + X_pos[pos_x] - 1) = Y_pos[pos_y]; // anchor position
							*(bits + X_pos[pos_x] - 1) += 1; // increment anchoring
						} else if (*(bits + l + X_pos[pos_x] - 1) == Y_pos[pos_y]) { // previous anchor
							*(bits + X_pos[pos_x] - 1) += 1; // increment anchoring
						} else { // does not match previous anchoring
							*(bits + X_pos[pos_x] - 1) = 0; // reset inconsistent anchoring
							*(bits + l + X_pos[pos_x] - 1) = 0;
						}
					}
					link_x = pos_x;
					link_y = pos_y;
					offset = 0;
				}
			}
			offset++;
		}
		
		R_CheckUserInterrupt();
	}
	
	size /= 2;
	// calculate ranges of anchors
	int *temp = R_Calloc(l, int); // initialized to zero
	int match = 0, count = -1, last_end_x = -10000, last_end_y = -10000;
	for (i = 0; i < l; i++) {
		//Rprintf("\n%d start_x=%d start_y=%d end_x=%d end_y=%d last_end_x=%d last_end_y=%d", i, i - wS + 2, *(bits + i + l) - wS + 1, i+1, *(bits + i + l), last_end_x, last_end_y);
		if ((double)*(bits + i)/(double)size >= thresh) {
			if (match == 0) { // start of anchor range
				if (i - wS + 2 > last_end_x + 100 &&
					*(bits + i + l) - wS + 1 > last_end_y + 100) {
					match = 1;
					count++;
					last_end_x = i - 100 - wS;
					last_end_y = *(bits + i + l - 1) - 100 - wS;
					*(temp + count*4) = i + 1;//i - wS + 2;
					*(temp + count*4 + 2) = *(bits + i + l);//*(bits + i + l) - wS + 1;
					*(temp + count*4 + 1) = i + 1;
					*(temp + count*4 + 3) = *(bits + i + l);
				}
			} else if (i - wS + 2 > last_end_x + 1000 &&
				*(bits + i + l) - wS + 1 > last_end_y + 1000) {
				count++; // start new range
				last_end_x = i - 100 - wS;
				last_end_y = *(bits + i + l - 1) - 100 - wS;
				*(temp + count*4) = i + 1;//i - wS + 2;
				*(temp + count*4 + 2) = *(bits + i + l);//*(bits + i + l) - wS + 1;
				*(temp + count*4 + 1) = i + 1;
				*(temp + count*4 + 3) = *(bits + i + l);
			} else if (i + 1 > last_end_x + 100 &&
				*(bits + i + l) > last_end_y + 100 &&
				*(bits + i + l) > *(temp + count*4 + 3)) { // extend end
				*(temp + count*4 + 1) = i + 1;
				*(temp + count*4 + 3) = *(bits + i + l);
			}
		} else if (match == 1) { // set to last match
			match = 0;
			last_end_x = *(temp + count*4 + 1);
			last_end_y = *(temp + count*4 + 3);
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, 4, count + 1));
	rans = INTEGER(ans);
	for (i = 0; i <= count; i++) {
		*(rans + i*4) = *(temp + i*4);
		*(rans + i*4 + 1) = *(temp + i*4 + 1);
		*(rans + i*4 + 2) = *(temp + i*4 + 2);
		*(rans + i*4 + 3) = *(temp + i*4 + 3);
	}
	
	UNPROTECT(1);
	R_Free(bits);
	R_Free(temp);
	
	return ans;
}

// which (bl <= x <= bu)
// x must be in ascending order
SEXP boundedMatches(SEXP x, SEXP bl, SEXP bu)
{
	int i, mid, start = 0, end = length(x), count = 0, size_x = length(x);
	int lowBound = asInteger(bl);
	int upBound = asInteger(bu);
	int *v = INTEGER(x);
	int *buffer = (int *) R_alloc(size_x, sizeof(int));
	
	while (start < end) {
		mid = floor(start + (end - start)/2);
		if (v[mid] >= lowBound) {
			end = mid;
		} else if (start == mid) {
			break;
		} else {
			start = mid;
		}
	}
	
	for (i = end; i < size_x; i++) {
		if (v[i] >= lowBound && v[i] <= upBound) {
			buffer[count] = i + 1;
			count++;
		} else {
			break;
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, count));
	memcpy(INTEGER(ans), buffer, sizeof(int) * count);
	
	UNPROTECT(1);
	
	return ans;
}

// first unmatched occurrence of x in y for ascending order integer vectors
// requires NAs to be first in the ordering (treated as the most negative)
// similar to match(x, y, incomparables=NA) without repetition in the result
// performs the following reorderings of the inputs and output:
// ans[o1] <- o2[intMatchOnce(x[o1], y[o2])] + 1
// where o1 and o2 are indexed starting at zero, and ans starts at 1
SEXP intMatchOnce(SEXP x, SEXP y, SEXP o1, SEXP o2)
{
	int *v = INTEGER(x);
	int *w = INTEGER(y);
	int *p = INTEGER(o1);
	int *q = INTEGER(o2);
	
	int i, j, k, temp, start = 0;
	int size_x = length(x);
	int size_y = length(y);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, size_x));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < size_x; i++) {
		rans[p[i]] = NA_INTEGER;
		if (v[p[i]] != NA_INTEGER)
			break;
	}
	
	for (; i < size_x; i++) { // i = i
		temp = NA_INTEGER;
		for (j = start; j < size_y; j++) {
			if (v[p[i]] < w[q[j]]) {
				start = j;
				break;
			} else if (v[p[i]] == w[q[j]]) {
				k = j + 1;
				if (k < size_y &&
					w[q[j]] == w[q[k]]) {
					start = k; // prevent repeated matching
					temp = j;
				} else {
					start = j; // allow repeated matching
					temp = j;
				}
				break;
			}
		}
		
		if (temp == NA_INTEGER) {
			rans[p[i]] = NA_INTEGER;
		} else {
			rans[p[i]] = q[temp] + 1;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// first unmatched occurrence of x in x for ascending order integer vectors
// prevents matches between the same positions in x (otherwise ans = 1, 2, ...)
// requires NAs to be first in the ordering (treated as the most negative)
// performs the following reorderings of the inputs and output:
// ans[o1] <- o1[intMatchOnce(x[o1], x[o1])] + 1
// where o1 is indexed starting at zero, and ans starts at 1
SEXP intMatchSelfOnce(SEXP x, SEXP o1)
{
	int *v = INTEGER(x);
	int *p = INTEGER(o1);
	
	int i, j, k, temp;
	int size_x = length(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, size_x));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < size_x; i++) {
		rans[p[i]] = NA_INTEGER;
		if (v[p[i]] != NA_INTEGER)
			break;
	}
	
	for (; i < size_x; i++) { // i = i
		temp = NA_INTEGER;
		for (j = i + 1; j < size_x; j++) {
			if (v[p[i]] < v[p[j]]) {
				break;
			} else if (v[p[i]] == v[p[j]]) {
				k = j + 1;
				if (k < size_x &&
					v[p[j]] == v[p[k]]) {
					temp = j;
				} else {
					temp = j;
				}
				break;
			}
		}
		
		if (temp == NA_INTEGER) {
			rans[p[i]] = NA_INTEGER;
		} else {
			rans[p[i]] = p[temp] + 1;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// index of first maxes in x by z groups in y
SEXP groupMax(SEXP x, SEXP y, SEXP z)
{
	double *v = REAL(x); // values
	int *w = INTEGER(y); // groups
	int *u = INTEGER(z); // unique groups
	int l = length(x); // number of values
	int n = length(z); // number of groups
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, n));
	int *rans = INTEGER(ans);
	
	int curr = 0;
	for (int i = 0; i < n; i++) {
		double max = -1e53;
		while (curr < l && w[curr] == u[i]) {
			if (v[curr] > max) {
				rans[i] = curr + 1; // index starting at 1
				max = v[curr];
			}
			curr++;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// sum integers randomly projected within bins
SEXP sumBins(SEXP v, SEXP bins)
{
	int i, j;
	unsigned int val;
	int n = length(v);
	int b = asInteger(bins);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, b));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < b; i++)
		rans[i] = 0;
	
	for (i = 0; i < n; i++) {
		int *V = INTEGER(VECTOR_ELT(VECTOR_ELT(v, i), 0));
		int l = length(VECTOR_ELT(VECTOR_ELT(v, i), 0));
		
		int prev = -1;
		for (j = 0; j < l; j++) {
			val = V[j];
			if (val != prev) {
				prev = val;
				val ^= val << 13;
				val ^= val >> 17;
				val ^= val << 5;
				val = val % b;
				rans[val]++;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// performs Xorshift random number generation [1, base]
SEXP xorShift(SEXP seed, SEXP base)
{
	int i;
	unsigned int val;
	int b = asInteger(base);
	int *s = INTEGER(seed);
	int l = length(seed);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < l; i++) {
		val = s[i];
		val ^= val << 13;
		val ^= val >> 17;
		val ^= val << 5;
		rans[i] = (val % b) + 1;
	}
	
	UNPROTECT(1);
	
	return ans;
}

// selects groups from a grouping vector (stopping when reaching next)
SEXP selectGroups(SEXP ordering, SEXP initial, SEXP final, SEXP num, SEXP next)
{
	int i, g, k = 0;
	R_xlen_t j, end;
	
	int *o = INTEGER(ordering);
	int *fin = INTEGER(final);
	int N = asInteger(num);
	int l = length(initial);
	int n = asInteger(next);
	
	int M = 0;
	for (i = 0; i < l; i++)
		if (fin[i] != NA_INTEGER)
			M += fin[i];
	if (M > N)
		M = N;
	
	int *groups = (int *) malloc(M*sizeof(int)); // thread-safe on Windows
	
	i = 0;
	if (isInteger(initial)) {
		int *ini = INTEGER(initial);
		while (i < l && k < M) {
			if (fin[i] != NA_INTEGER) {
				j = (R_xlen_t)ini[i] - 1;
				end = j + fin[i];
				while (j < end && k < M) {
					g = o[j++];
					if (g == n)
						break;
					groups[k++] = g;
				}
			}
			i++;
		}
	} else {
		double *ini = REAL(initial);
		while (i < l && k < M) {
			if (fin[i] != NA_INTEGER) {
				j = (R_xlen_t)ini[i] - 1;
				end = j + fin[i];
				while (j < end && k < M) {
					g = o[j++];
					if (g == n)
						break;
					groups[k++] = g;
				}
			}
			i++;
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, k));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < k; i++)
		rans[i] = groups[i];
	free(groups);
	
	UNPROTECT(1);
	
	return ans;
}

// first hit where x[1] == y[...]
// all y values must be sorted ascending and distinct
SEXP firstHit(SEXP x, SEXP y)
{
	int size = length(y);
	int v = INTEGER(x)[0];
	int *w = INTEGER(y);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, 1));
	int *rans = INTEGER(ans);
	rans[0] = NA_INTEGER;
	
	int start = 0;
	int end = size - 1;
	int loc = end/2;
	
	if (size > 0) {
		if (v >= w[start] && v <= w[end]) {
			if (v == w[start]) {
				rans[0] = start + 1;
			} else if (v == w[end]) {
				rans[0] = end + 1;
			} else {
				while (loc > start) {
					if (w[loc] == v) {
						rans[0] = loc + 1;
						break;
					} else if (w[loc] < v) {
						start = loc;
					} else {
						end = loc;
					}
					loc = (end - start)/2 + start;
				}
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// uniques a sorted integer vector
// drops any NAs at start of vector
SEXP sortedUnique(SEXP v)
{
	int i, n = 0;
	int *V = INTEGER(v);
	int l = length(v);
	
	int *w = (int *) malloc(l*sizeof(int)); // thread-safe on Windows
	
	int prev = NA_INTEGER;
	for (i = 0; i < l; i++) {
		if (V[i] != prev) {
			w[n++] = i;
			prev = V[i];
		}
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, n));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < n; i++)
		rans[i] = V[w[i]];
	free(w);
	
	UNPROTECT(1);
	
	return ans;
}

// (maybe in-place) splitting of partitions
SEXP splitPartitions(SEXP order, SEXP partition, SEXP var, SEXP minSize, SEXP minSplit)
{
	int l = length(partition); // also length of order and var
	int *o = INTEGER(order);
	double *v = REAL(var);
	int m = asInteger(minSize);
	double s = asReal(minSplit);
	
	int shared;
	SEXP ans;
	if (MAYBE_SHARED(partition)) {
		shared = 1;
		PROTECT(ans = duplicate(partition));
	} else {
		shared = 0;
		ans = partition;
	}
	int *p = INTEGER(ans);
	
	int count = 0; // new partition number
	int prev = 0; // previous partition number
	int change; // index before start of partition
	for (int j = l - 1; j >= 0; j--) { // assume order is decreasing
		if (prev != p[o[j] - 1]) { // change in partition number
			count++;
			prev = p[o[j] - 1];
			change = j - 1;
		} else if (change - j >= m && // large enough partition
			v[o[j] - 1] <= s) { // split partition
			count++;
			change = j - 1;
		}
		p[o[j] - 1] = count;
	}
	
	if (shared)
		UNPROTECT(1);
	
	return ans;
}

// create a list of partition indices
// order of partition must be decreasing
SEXP indexPartitions(SEXP order, SEXP partition, SEXP keep)
{
	int i, j, m, n, *rans;
	int l = length(partition); // also length of order and keep
	int *o = INTEGER(order);
	int *p = INTEGER(partition);
	int *k = INTEGER(keep);
	
	int prev = 0; // previous partition number
	int count = 0; // number of partitions
	for (i = 0; i < l; i++) {
		if (prev != p[o[i] - 1]) {
			count++;
			prev = p[o[i] - 1];
		}
	}
	
	SEXP ans, ans1, ans2;
	PROTECT(ans1 = allocVector(VECSXP, count));
	PROTECT(ans2 = allocVector(LGLSXP, count));
	int *rans2 = INTEGER(ans2);
	
	j = l - 1; // start of partition
	prev = 0;
	count = 0;
	for (i = j; i >= -1; i--) { // start from first partition
		if (i == -1 || prev != p[o[i] - 1]) {
			if (prev != 0) {
				PROTECT(ans = allocVector(INTSXP, j - i));
				rans = INTEGER(ans);
				
				m = 0;
				n = 0;
				while (j > i) {
					n += k[o[j] - 1];
					rans[m++] = o[j--];
				}
				// j == i for next iteration
				
				SET_VECTOR_ELT(ans1, count, ans);
				UNPROTECT(1);
				if (n < 2) {
					rans2[count] = 0;
				} else {
					rans2[count] = 1;
				}
				count++;
			}
			if (i >= 0)
				prev = p[o[i] - 1];
		}
	}
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 2));
	
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	
	UNPROTECT(3);
	
	return ret_list;
}

// detect the number of available cores
SEXP detectCores()
{
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, 1));
	int *rans = INTEGER(ans);
	
	#ifdef _OPENMP
	rans[0] = omp_get_num_procs();
	#else
	rans[0] = 1;
	#endif
	
	UNPROTECT(1);
	
	return ans;
}

int comparePointers(const void *one, const void *two)
{
	const unsigned int **left  = (const unsigned int **)one;
	const unsigned int **right = (const unsigned int **)two;
	
	return (**left < **right) - (**right < **left);
}

int *orderInt(unsigned int *a, int n)
{
	int i;
	
	unsigned int **pointers = malloc(n*sizeof(unsigned int *));
	
	for (i = 0; i < n; i++)
		pointers[i] = a + i;
	
	qsort(pointers, n, sizeof(unsigned int *), comparePointers); // unstable sort
	
	int *indices = malloc(n*sizeof(int));
	for (i = 0; i < n; i++)
		indices[i] = pointers[i] - a;
	
	free(pointers);
	
	return indices;
}

// matching occurrence of columns in a set of sequences
// optionally flag uninformative (constant or missing) columns
SEXP matchColumns(SEXP x, SEXP letters)
{
	// initialize variables
	int i, j, k, l, p;
	
	XStringSet_holder x_set, l_set;
	Chars_holder x_i, l_i;
	x_set = hold_XStringSet(x);
	int n = get_length_from_XStringSet_holder(&x_set); // number of sequences
	l_set = hold_XStringSet(letters);
	l_i = get_elt_from_XStringSet_holder(&l_set, 0);
	int c = l_i.length; // number of states
	int *lkup;
	if (c > 0) {
		lkup = (int *) malloc(256*sizeof(int)); // thread-safe on Windows
		for (i = 0; i < 256; i++)
			lkup[i] = NA_INTEGER;
		for (i = 0; i < c; i++)
			lkup[(unsigned char)l_i.ptr[i]] = i;
	}
	
	// determine longest sequence
	l = 0;
	for (i = 0; i < n; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (x_i.length > l)
			l = x_i.length;
	}
	
	// hash the columns using xorshift
	unsigned int *hash = (unsigned int *) calloc(l, sizeof(unsigned int)); // initialized to zero (thread-safe on Windows)
	for (i = 0; i < n; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		for (j = 0; j < x_i.length; j++) {
			hash[j] ^= (unsigned int)x_i.ptr[j];
			hash[j] ^= hash[j] << 13;
			hash[j] ^= hash[j] >> 17;
			hash[j] ^= hash[j] << 5;
		}
	}
	
	int *o = orderInt(hash, l);
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l));
	int *rans = INTEGER(ans);
	
	// compare columns
	k = 0; // position in hash
	for (j = 0; j < l; j++) { // column in x
		rans[o[j]] = o[j] + 1; // no match yet
		while (k < j && hash[o[k]] > hash[o[j]])
			k++;
		p = k;
		while (p < j && hash[o[j]] == hash[o[p]]) {
			int match = 1;
			for (i = 0; i < n; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				if (x_i.length > o[j]) {
					if (x_i.length <= o[p] || x_i.ptr[o[j]] != x_i.ptr[o[p]]) {
						match = 0;
						break;
					}
				} else if (x_i.length > o[p]) {
					match = 0;
					break;
				}
			}
			if (match) {
				rans[o[j]] = rans[o[p]];
				break;
			}
			p++;
		}
		if (c > 0 && rans[o[j]] == o[j] + 1) { // filter uninformative sites
			int first = c;
			int second = c;
			for (i = 0; i < n; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				if (o[j] < x_i.length) {
					int val = lkup[(unsigned char)x_i.ptr[o[j]]];
					if (val != NA_INTEGER) {
						if (first == c) {
							first = val;
							second = val;
						} else if (val != first) {
							second = val;
							break;
						}
					}
				}
			}
			if (first == second) // missing or constant column
				rans[o[j]] = 0; // flag as uninformative
		}
	}
	
	free(hash);
	free(o);
	if (c > 0)
		free(lkup);
	
	UNPROTECT(1);
	
	return ans;
}

// hash list of integer vectors into an integer vector
SEXP hashList(SEXP x)
{
	int i, j, *X, lx;
	int l = length(x);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, l));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < l; i++) {
		X = INTEGER(VECTOR_ELT(x, i));
		lx = length(VECTOR_ELT(x, i));
		
		rans[i] = 0;
		for (j = 0; j < lx; j++) {
			rans[i] ^= (unsigned int)X[j];
			rans[i] ^= rans[i] << 13;
			rans[i] ^= rans[i] >> 17;
			rans[i] ^= rans[i] << 5;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// first row and value of two column integer matrix with value greater than row index
SEXP firstRow(SEXP x)
{
	int i = 0;
	int *X = INTEGER(x);
	int l = length(x);
	int j = l/2;
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, 2));
	int *rans = INTEGER(ans);
	rans[0] = 0;
	
	while (j < l) {
		if (X[i] > X[j]) {
			if (X[i] > i) {
				rans[0] = i + 1;
				rans[1] = X[i];
				break;
			}
		} else {
			if (X[j] > i) {
				rans[0] = i + 1;
				rans[1] = X[j];
				break;
			}
		}
		i++;
		j++;
	}
	
	UNPROTECT(1);
	
	return ans;
}

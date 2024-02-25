/****************************************************************************
 *                         Creates Distance Marix                           *
 *                           Author: Erik Wright                            *
 ****************************************************************************/
 
 // for OpenMP parallel processing
 #ifdef _OPENMP
 #include <omp.h>
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

static double distance(const Chars_holder *P, const Chars_holder *S, int start, int end, int pGapLetters, int width, double coverage)
{
	double distance;
	int i, j, mismatches, gapGapMatches, gapLetterMatches, count, letters, state;
	const char *p, *s;
	
	distance = 0;
	gapGapMatches = 0;
	gapLetterMatches = 0;
	mismatches = 0;
	count = 0;
	letters = 0;
	state = 0;
	
	// walk along the sequence from (position start + 1) to (length - end - 1)
	for (i = start, j = start, p = (P->ptr + start), s = (S->ptr + start);
		(i < (P->length - end)) && (i < (S->length - end));
		i++, j++, p++, s++)
	{
		if (!((*p) & 0x20 || (*s) & 0x20)) { // not masked
			count++; // increment the length covered
			if (!((*p) & (*s))) { // sequences are not equal
				if (((*p) & 0x40 && (*s) & 0x10) || ((*p) & 0x10 && (*s) & 0x40)) { // gap-gap match
					gapGapMatches++; // don't include gap-gap matches in length
				} else if ((*p) & 0x10 || (*s) & 0x10 || (*p) & 0x40 || (*s) & 0x40) { // gap-letter match
					if (pGapLetters == 0) { // don't penalize gap-letter matches
						gapLetterMatches++; // don't include gap-letter matches in length
					} else if (pGapLetters == 1) { // penalize gap-letter matches
						mismatches++; // count gap-letter matches as mis-matches
					} else { // penalize state changes
						if ((*p) & 0x10 || (*p) & 0x40) {
							if (state != 1) {
								mismatches++;
							} else {
								gapLetterMatches++;
							}
							state = 1;
						} else {
							if (state != 2) {
								mismatches++;
							} else {
								gapLetterMatches++;
							}
							state = 2;
						}
					}
				} else {
					mismatches++; // mis-match
					letters++;
					state = 0;
				}
			} else { // sequences are equal
				if (((*p) & 0x10 && (*s) & 0x10) || ((*p) & 0x40 && (*s) & 0x40)) { // gap-gap match
					gapGapMatches++; // don't include gap-gap matches in length
				} else {
					letters++;
					state = 0;
				}
			}
		}
	}
	
	//Rprintf("start%d end%d",start,end);
	//Rprintf("\nmismatches:%d gapGapMatches:%d gapLetterMatches:%d count:%d",mismatches,gapGapMatches,gapLetterMatches,count);
	
	// calculate distance as the percent mis-matches
	if (coverage > 0 && (double)letters/(double)width < coverage) {
		distance = NA_REAL;
	} else {
		distance = (double)mismatches/((double)count - (double)gapGapMatches - (double)gapLetterMatches);
	}
	return distance;
}

static double distanceAA(const Chars_holder *P, const Chars_holder *S, int start, int end, int pGapLetters, int width, double coverage)
{
	double distance;
	int i, j, mismatches, gapGapMatches, gapLetterMatches, count, letters, state;
	const char *p, *s;
	
	distance = 0;
	gapGapMatches = 0;
	gapLetterMatches = 0;
	mismatches = 0;
	count = 0;
	letters = 0;
	state = 0;
	
	// walk along the sequence from (position start + 1) to (length - end - 1)
	for (i = start, j = start, p = (P->ptr + start), s = (S->ptr + start);
		(i < (P->length - end)) && (i < (S->length - end));
		i++, j++, p++, s++)
	{
		if ((*p) ^ 0x2B && (*s) ^ 0x2B) { // not masked
			count++; // increment the length covered
			if ((*p) ^ (*s) && // sequences are not equal
				!(!((*p) ^ 0x58) && !(!((*s) ^ 0x2D) || !((*s) ^ 0x2B) || !((*s) ^ 0x2A))) && !(!((*s) ^ 0x58) && !(!((*p) ^ 0x2D) || !((*p) ^ 0x2B) || !((*p) ^ 0x2A))) && // not (X && !(non-letter))
				!(!((*p) ^ 0x42) && (!((*s) ^ 0x4E) || !((*s) ^ 0x44))) && !(!((*s) ^ 0x42) && (!((*p) ^ 0x4E) || !((*p) ^ 0x44))) && // not (B && (N or D))
				!(!((*p) ^ 0x4A) && (!((*s) ^ 0x49) || !((*s) ^ 0x4C))) && !(!((*s) ^ 0x4A) && (!((*p) ^ 0x49) || !((*p) ^ 0x4C))) && // not (J && (I or L))
				!(!((*p) ^ 0x5A) && (!((*s) ^ 0x51) || !((*s) ^ 0x45))) && !(!((*s) ^ 0x5A) && (!((*p) ^ 0x51) || !((*p) ^ 0x45)))) { // not (Z && (Q or E))
				if ((!((*p) ^ 0x2D) && !((*s) ^ 0x2E)) || (!((*p) ^ 0x2E) && !((*s) ^ 0x2D))) { // gap-gap match
					gapGapMatches++; // don't include gap-gap matches in length
				} else if (!((*p) ^ 0x2D) || !((*s) ^ 0x2D) || !((*p) ^ 0x2E) || !((*s) ^ 0x2E)) { // gap-letter match
					if (pGapLetters == 0) { // don't penalize gap-letter matches
						gapLetterMatches++; // don't include gap-letter matches in length
					} else if (pGapLetters == 1) { // penalize gap-letter matches
						mismatches++; // count gap-letter matches as mis-matches
					} else { // penalize state changes
						if ((*p) ^ 0x2D && (*p) ^ 0x2E) {
							if (state != 1) {
								mismatches++;
							} else {
								gapLetterMatches++;
							}
							state = 1;
						} else {
							if (state != 2) {
								mismatches++;
							} else {
								gapLetterMatches++;
							}
							state = 2;
						}
					}
				} else {
					mismatches++; // mis-match
					letters++;
					state = 0;
				}
			} else { // sequences are equal
				if ((!((*p) ^ 0x2D) && !((*s) ^ 0x2D)) || (!((*p) ^ 0x2E) && !((*s) ^ 0x2E))) { // gap-gap match
					gapGapMatches++; // don't include gap-gap matches in length
				} else {
					letters++;
					state = 0;
				}
			}
		}
	}
	
	//Rprintf("start%d end%d",start,end);
	//Rprintf("\nmismatches:%d gapGapMatches:%d gapLetterMatches:%d count:%d",mismatches,gapGapMatches,gapLetterMatches,count);
	
	// calculate distance as the percent mis-matches
	if (coverage > 0 && (double)letters/(double)width < coverage) {
		distance = NA_REAL;
	} else {
		distance = (double)mismatches/((double)count - (double)gapGapMatches - (double)gapLetterMatches);
	}
	return distance;
}

static int frontTerminalGaps(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
		i < P->length;
		i++, p++)
	{
		if ((*p) & 0x10 || (*p) & 0x40) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGaps(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
		i >= 0;
		i--, p--)
	{
		if ((*p) & 0x10 || (*p) & 0x40) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int totalGaps(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
		i < P->length;
		i++, p++)
	{
		if ((*p) & 0x10 || (*p) & 0x40) // gap character
			gaps++; // count gaps
	}
	return gaps;
}

static int frontTerminalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
		i < P->length;
		i++, p++)
	{
		if (!((*p) ^ 0x2D) || !((*p) ^ 0x2E)) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int endTerminalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
		i >= 0;
		i--, p--)
	{
		if (!((*p) ^ 0x2D) || !((*p) ^ 0x2E)) { // gap character
			gaps++; // count gaps
		} else { // not a gap
			return gaps;
		}
	}
	return gaps;
}

static int totalGapsAA(const Chars_holder *P)
{
	int i, gaps;
	const char *p;
	gaps = 0;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
		i < P->length;
		i++, p++)
	{
		if (!((*p) ^ 0x2D) || !((*p) ^ 0x2E)) // gap character
			gaps++; // count gaps
	}
	return gaps;
}

SEXP distMatrix(SEXP x, SEXP t, SEXP terminalGaps, SEXP penalizeGapLetters, SEXP fullMatrix, SEXP output, SEXP e, SEXP minCoverage, SEXP method, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	XStringSet_holder x_set;
	Chars_holder x_i, x_j;
	int start, end, seqLength_i, seqLength_j, index, width, useMax;
	double E = asReal(e);
	R_xlen_t x_length, i, j, last;
	int pGapLetters, tGaps, fM = asLogical(fullMatrix);
	int before, v, *rPercentComplete;
	double *rans, soFar;
	int nthreads = asInteger(nThreads);
	int o = asInteger(output);
	double coverage = asReal(minCoverage);
	if (coverage >= 0) {
		useMax = 0;
	} else {
		useMax = 1;
		coverage *= -1;
	}
	int mode = asInteger(method);
	SEXP ans, percentComplete, utilsPackage;
	v = asLogical(verbose);
	if (v) { // percent complete variables
		soFar = 0;
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	if (x_length < 2) { // there is only one sequence
		PROTECT(ans = NEW_INTEGER(0));
	} else {
		if (o == 1) { // type is "matrix"
			if (fM) {
				last = x_length - 1;
				PROTECT(ans = allocMatrix(REALSXP, x_length, x_length));
			} else {
				last = 1;
				PROTECT(ans = allocMatrix(REALSXP, 1, x_length));
			}
		} else { // type is "dist"
			last = x_length - 1;
			PROTECT(ans = allocVector(REALSXP, x_length*(x_length - 1)/2));
		}
		rans = REAL(ans);
		
		// find the terminal gap lengths
		// always needed to identify no-overlap
		int *gapLengths[2];
		gapLengths[0] = (int *) calloc(x_length, sizeof(int)); // initialized to zero (thread-safe on Windows)
		gapLengths[1] = (int *) calloc(x_length, sizeof(int)); // initialized to zero (thread-safe on Windows)
		if (asInteger(t) == 3) { // AAStringSet
			for (i = 0; i < x_length; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				gapLengths[0][i] = frontTerminalGapsAA(&x_i);
				gapLengths[1][i] = endTerminalGapsAA(&x_i);
				//Rprintf("\nstart:%dend:%d",gapLengths[0][i],gapLengths[1][i]);
			}
		} else { // DNAStringSet or RNAStringSet
			for (i = 0; i < x_length; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				gapLengths[0][i] = frontTerminalGaps(&x_i);
				gapLengths[1][i] = endTerminalGaps(&x_i);
				//Rprintf("\nstart:%dend:%d",gapLengths[0][i],gapLengths[1][i]);
			}
		}
		
		// find the sequence lengths
		int *seqLengths = (int *) calloc(x_length, sizeof(int)); // initialized to zero (thread-safe on Windows)
		if (asInteger(t) == 3) { // AAStringSet
			for (i = 0; i < x_length; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				seqLengths[i] = x_i.length - totalGapsAA(&x_i);
			}
		} else { // DNAStringSet or RNAStringSet
			for (i = 0; i < x_length; i++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, i);
				seqLengths[i] = x_i.length - totalGaps(&x_i);
			}
		}
		
		tGaps = asLogical(terminalGaps);
		pGapLetters = asLogical(penalizeGapLetters);
		for (i = 0; i < last; i++) {
			// extract each ith DNAString from the DNAStringSet
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			seqLength_i = x_i.length;
			
			#ifdef _OPENMP
			#pragma omp parallel for private(j,x_j,seqLength_j,start,end,index,width) schedule(guided) num_threads(nthreads)
			#endif
			for (j = (i+1); j < x_length; j++) {
				// extract each jth DNAString from the DNAStringSet
				x_j = get_elt_from_XStringSet_holder(&x_set, j);
				seqLength_j = x_j.length;
				
				if (o == 1) { // matrix
					index = j + x_length*i;
				} else { // dist
					index = x_length*i - i*(i + 1)/2 + j - i - 1;
				}
				
				// find the distance for each row of the matrix
				if ((seqLength_i - gapLengths[1][i]) <= gapLengths[0][j] ||
					gapLengths[0][i] >= (seqLength_j - gapLengths[1][j])) {
					// no overlap between sequences
					if (tGaps && // include terminal gaps
						coverage == 0) { // do not require coverage
						rans[index] = 1;
					} else {
						rans[index] = NA_REAL;
					}
				} else {
					if (!tGaps) { // don't include terminal gaps
						// find the intersection of both string's ranges
						// to shorten the sequence comparison for speed
						if (gapLengths[0][i] >= gapLengths[0][j]) {
							start = gapLengths[0][i];
						} else {
							start = gapLengths[0][j];
						}
						if ((seqLength_i - gapLengths[1][i]) <= (seqLength_j - gapLengths[1][j])) {
							end = gapLengths[1][i];
						} else {
							end = gapLengths[1][j];
						}
						if (seqLengths[i] <= seqLengths[j]) {
							if (useMax == 0) {
								width = seqLengths[i];
							} else {
								width = seqLengths[j];
							}
						} else {
							if (useMax == 0) {
								width = seqLengths[j];
							} else {
								width = seqLengths[i];
							}
						}
					} else { // use whole sequence including terminal gaps
						if (mode == 1) { // overlap
							start = 0;
							end = 0;
							if (seqLengths[i] <= seqLengths[j]) {
								if (useMax == 0) {
									width = seqLengths[i];
								} else {
									width = seqLengths[j];
								}
							} else {
								if (useMax == 0) {
									width = seqLengths[j];
								} else {
									width = seqLengths[i];
								}
							}
						} else if (mode == 2) { // shortest
							if (seqLengths[i] <= seqLengths[j]) {
								if (useMax == 0) {
									width = seqLengths[i];
								} else {
									width = seqLengths[j];
								}
								start = gapLengths[0][i];
								end = gapLengths[1][i];
							} else {
								if (useMax == 0) {
									width = seqLengths[j];
								} else {
									width = seqLengths[i];
								}
								start = gapLengths[0][j];
								end = gapLengths[1][j];
							}
						} else { // longest
							if (seqLengths[i] >= seqLengths[j]) {
								if (useMax == 0) {
									width = seqLengths[j];
								} else {
									width = seqLengths[i];
								}
								start = gapLengths[0][i];
								end = gapLengths[1][i];
							} else {
								if (useMax == 0) {
									width = seqLengths[i];
								} else {
									width = seqLengths[j];
								}
								start = gapLengths[0][j];
								end = gapLengths[1][j];
							}
						}
					}
					if (asInteger(t) == 3) { // AAStringSet
						rans[index] = distanceAA(&x_i, &x_j, start, end, pGapLetters, width, coverage);
					} else {
						rans[index] = distance(&x_i, &x_j, start, end, pGapLetters, width, coverage);
					}
					if (E > 0) {
						if (rans[index] >= E) {
							rans[index] = R_PosInf;
						} else {
							rans[index] = 1 - rans[index]/E;
							rans[index] = -E*log(rans[index]);
						}
					}
				}
			}
			if (fM && o == 1) // make the matrix symetrical
				for (j = (i+1); j < x_length; j++)
					rans[i + x_length*j] = rans[j + x_length*i];
			
			if (o == 1) // set the matrix diagonal to zero distance
				rans[i*x_length + i] = 0;
			
			if (v) { // print the percent completed so far
				soFar = (2*last - i)*(i + 1);
				*rPercentComplete = floor(100*soFar/(last*(last + 1)));
				if (*rPercentComplete > before) { // when the percent has changed
					// tell the progress bar to update in the R console
					eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
					before = *rPercentComplete;
				}
			} else {
				R_CheckUserInterrupt();
			}
		}
		if (fM && o == 1) // set the last element of the diagonal to zero
			rans[(x_length - 1)*x_length + (x_length - 1)] = 0;
		
		free(gapLengths[0]);
		free(gapLengths[1]);
		free(seqLengths);
	}
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;	
}

//.Call("gaps", myDNAStringSet, 1L, PACKAGE="DECIPHER")
SEXP gaps(SEXP x, SEXP t)
{
	XStringSet_holder x_set;
	Chars_holder x_i;
	int x_length, i;
	double *rans;
	SEXP ans;
	
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);

	PROTECT(ans = allocMatrix(REALSXP, x_length, 3));
	rans = REAL(ans);
	
	// find the three lengths for each sequence
	if (asInteger(t) == 3) { // AAStringSet
		for (i = 0; i < x_length; i++) {
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			rans[i + x_length*0] = frontTerminalGapsAA(&x_i);
			rans[i + x_length*1] = endTerminalGapsAA(&x_i);
			rans[i + x_length*2] = x_i.length - rans[i + x_length*1] - rans[i + x_length*0];
		}
	} else { // DNAStringSet or RNAStringSet
		for (i = 0; i < x_length; i++) {
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			rans[i + x_length*0] = frontTerminalGaps(&x_i);
			rans[i + x_length*1] = endTerminalGaps(&x_i);
			rans[i + x_length*2] = x_i.length - rans[i + x_length*1] - rans[i + x_length*0];
		}
	}
	
	UNPROTECT(1);
	return ans;
}

//.Call("firstSeqsEqual", dna1, dna2, start1, end1, start2, end2, first1, first2, PACKAGE="DECIPHER")
SEXP firstSeqsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP first_x, SEXP first_y)
{	
	int i, j;
	XStringSet_holder x_set;
	XStringSet_holder y_set;
	Chars_holder x_i, y_i;
	int sx = asInteger(start_x);
	int ex = asInteger(end_x);
	int fx = asInteger(first_x);
	int sy = asInteger(start_y);
	int ey = asInteger(end_y);
	int fy = asInteger(first_y);
	
	
	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1));
	int *rans;
	rans = INTEGER(ans);
	*(rans) = 1; // equal
	if ((sx - ex) != (sy - ey)) { // different sequence lengths
		*(rans) = 0; // not equal
	} else {
		x_set = hold_XStringSet(x);
		y_set = hold_XStringSet(y);
		x_i = get_elt_from_XStringSet_holder(&x_set, fx - 1);
		y_i = get_elt_from_XStringSet_holder(&y_set, fy - 1);
		for (i = sx - 1, j = sy - 1;
			 i < ex; // i <= ex - 1 covers j <= ey - 1 because equal length
			 i++, j++) {
			if (x_i.ptr[i] != y_i.ptr[j]) {
				*(rans) = 0; // not equal
				break;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP firstSeqsGapsEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP t, SEXP first_x, SEXP first_y)
{	
	int i, j;
	XStringSet_holder x_set;
	XStringSet_holder y_set;
	Chars_holder x_i, y_i;
	int sx = asInteger(start_x);
	int ex = asInteger(end_x);
	int fx = asInteger(first_x);
	int sy = asInteger(start_y);
	int ey = asInteger(end_y);
	int fy = asInteger(first_y);
	
	
	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1));
	int *rans;
	rans = INTEGER(ans);
	*(rans) = 1; // equal
	if ((sx - ex) != (sy - ey)) { // different sequence lengths
		*(rans) = 0; // not equal
	} else {
		x_set = hold_XStringSet(x);
		y_set = hold_XStringSet(y);
		x_i = get_elt_from_XStringSet_holder(&x_set, fx - 1);
		y_i = get_elt_from_XStringSet_holder(&y_set, fy - 1);
		if (asInteger(t) == 3) { // AAStringSet
			for (i = sx - 1, j = sy - 1;
				 i < ex; // i <= ex - 1 covers j <= ey - 1 because equal length
				 i++, j++) {
				if ((!((*((char *)x_i.ptr + i)) ^ 0x2D) || !((*((char *)x_i.ptr + i)) ^ 0x2E)) ^
					(!((*((char *)y_i.ptr + j)) ^ 0x2D) || !((*((char *)y_i.ptr + j)) ^ 0x2E))) {
						*(rans) = 0; // not equal
						break;
					}
			}
		} else { // DNAStringSet or RNAStringSet
			for (i = sx - 1, j = sy - 1;
				 i < ex; // i <= ex - 1 covers j <= ey - 1 because equal length
				 i++, j++) {
				if (((*((char *)x_i.ptr + i)) & 0x10 || (*((char *)x_i.ptr + i)) & 0x40) ^
					((*((char *)y_i.ptr + j)) & 0x10 || (*((char *)y_i.ptr + j)) & 0x40)) {
					*(rans) = 0; // not equal
					break;
				}
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP firstSeqsPosEqual(SEXP x, SEXP y, SEXP start_x, SEXP end_x, SEXP start_y, SEXP end_y, SEXP t, SEXP first_x, SEXP first_y)
{
	XStringSet_holder x_set;
	XStringSet_holder y_set;
	Chars_holder x_i, y_i;
	int sx = asInteger(start_x);
	int ex = asInteger(end_x);
	int fx = asInteger(first_x);
	int sy = asInteger(start_y);
	int ey = asInteger(end_y);
	int fy = asInteger(first_y);
	int type = asInteger(t);
	
	int sizex = 100, sizey = 100;
	int *nx = Calloc(sizex, int); // number of gaps in x
	int *px = Calloc(sizex, int); // position of gaps in x
	int *ny = Calloc(sizey, int); // number of gaps in y
	int *py = Calloc(sizey, int); // position of gaps in y
	
	int i = sx - 1, j = sy - 1; // position in x or y
	int cx = 0, cy = 0; // number of sites
	int ci = 1, cj = 1; // include site in count
	int gx = -1, gy = -1; // index of current gap
	int itx = 0, ity = 0; // iterate position
	
	x_set = hold_XStringSet(x);
	y_set = hold_XStringSet(y);
	x_i = get_elt_from_XStringSet_holder(&x_set, fx - 1);
	y_i = get_elt_from_XStringSet_holder(&y_set, fy - 1);
	
	while (i < ex && j < ey) {
		if (ci) {
			if (type == 3) { // AAStringSet
				if (!(!((*((char *)x_i.ptr + i)) ^ 0x2D) || !((*((char *)x_i.ptr + i)) ^ 0x2E))) {
					cx++; // site
				}
			} else { // DNAStringSet or RNAStringSet
				if (!((*((char *)x_i.ptr + i)) & 0x10 || (*((char *)x_i.ptr + i)) & 0x40)) {
					cx++; // site
				}
			}
		}
		if (cj) {
			if (type == 3) { // AAStringSet
				if (!(!((*((char *)y_i.ptr + j)) ^ 0x2D) || !((*((char *)y_i.ptr + j)) ^ 0x2E))) {
					cy++; // site
				}
			} else { // DNAStringSet or RNAStringSet
				if (!((*((char *)y_i.ptr + j)) & 0x10 || (*((char *)y_i.ptr + j)) & 0x40)) {
					cy++; // site
				}
			}
		}
		
		if (cx > cy) {
			j++;
			if (itx) {
				nx[gx]++;
			} else {
				ci = 0;
				itx = 1;
				gx++;
				
				if (gx >= (sizex - 1)) {
					sizex += 100;
					nx = Realloc(nx, sizex, int);
					px = Realloc(px, sizex, int);
				}
				
				nx[gx] = 1;
				px[gx] = i + 1;
			}
		} else if (cx < cy) {
			i++;
			if (ity) {
				ny[gy]++;
			} else {
				cj = 0;
				ity = 1;
				gy++;
				
				if (gy >= (sizey - 1)) {
					sizey += 100;
					ny = Realloc(ny, sizey, int);
					py = Realloc(py, sizey, int);
				}
				
				ny[gy] = 1;
				py[gy] = j + 1;
			}
		} else {
			ci = 1;
			cj = 1;
			i++;
			j++;
			itx = 0;
			ity = 0;
		}
	}
	
	if (i < ex) {
		gy++;
		ny[gy] = ex - i;
		py[gy] = ey + 1;
		
		if (!ci)
			i++;
			
			while (i < ex) {
				if (type == 3) { // AAStringSet
					if (!(!((*((char *)x_i.ptr + i)) ^ 0x2D) || !((*((char *)x_i.ptr + i)) ^ 0x2E))) {
						cx++; // site
					}
				} else { // DNAStringSet or RNAStringSet
					if (!((*((char *)x_i.ptr + i)) & 0x10 || (*((char *)x_i.ptr + i)) & 0x40)) {
						cx++; // site
					}
				}
				i++;
			}
	}
	
	if (j < ey) {
		gx++;
		nx[gx] = ey - j;
		px[gx] = ex + 1;
		
		if (!cj)
			j++;
			
			while (j < ey) {
				if (type == 3) { // AAStringSet
					if (!(!((*((char *)y_i.ptr + j)) ^ 0x2D) || !((*((char *)y_i.ptr + j)) ^ 0x2E))) {
						cy++; // site
					}
				} else { // DNAStringSet or RNAStringSet
					if (!((*((char *)y_i.ptr + j)) & 0x10 || (*((char *)y_i.ptr + j)) & 0x40)) {
						cy++; // site
					}
				}
				j++;
			}
	}
	
	SEXP ret_list, ans;
	if (cx == cy) { // same number of sites
		PROTECT(ret_list = allocVector(VECSXP, 4));
		int *rans;
		
		PROTECT(ans = allocVector(INTSXP, gx + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gx; i++) {
			*(rans + i) = px[i];
		}
		SET_VECTOR_ELT(ret_list, 0, ans);
		PROTECT(ans = allocVector(INTSXP, gx + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gx; i++) {
			*(rans + i) = nx[i];
		}
		SET_VECTOR_ELT(ret_list, 1, ans);
		PROTECT(ans = allocVector(INTSXP, gy + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gy; i++) {
			*(rans + i) = py[i];
		}
		SET_VECTOR_ELT(ret_list, 2, ans);
		PROTECT(ans = allocVector(INTSXP, gy + 1));
		rans = INTEGER(ans);
		for (i = 0; i <= gy; i++) {
			*(rans + i) = ny[i];
		}
		SET_VECTOR_ELT(ret_list, 3, ans);
	}
	
	Free(nx);
	Free(px);
	Free(ny);
	Free(py);
	
	if (cx == cy) { // same number of sites
		UNPROTECT(5);
		
		return ret_list;
	} else {
		return R_NilValue;
	}
}

// compute overlap between indicies x (length == 1) and y (length >= 1)
SEXP computeOverlap(SEXP x, SEXP y, SEXP v, SEXP wordSize, SEXP mult, SEXP entropy, SEXP maxAlpha, SEXP u1, SEXP u2, SEXP seqs, SEXP dropScore, SEXP subMatrix, SEXP letters, SEXP terminalGaps, SEXP penalizeGapLetters, SEXP minCoverage, SEXP method, SEXP nThreads)
{
	int i, j, k, l, n, p, d, lx, new, *t, *keep, *I, *X, *OX, pos, p1, p2, t1, t2, ov, OV, off, g1, g2, g, o, count, useMax;
	
	int i1 = asInteger(x) - 1;
	I = INTEGER(y);
	X = INTEGER(VECTOR_ELT(VECTOR_ELT(v, i1), 0));
	OX = INTEGER(VECTOR_ELT(VECTOR_ELT(v, i1), 1));
	lx = length(VECTOR_ELT(VECTOR_ELT(v, i1), 0));
	int wS = asInteger(wordSize);
	double logN = log2(asReal(mult));
	double logE = log2(asReal(entropy));
	int maxA = asInteger(maxAlpha);
	int tGaps = asLogical(terminalGaps);
	int pGapLetters = asLogical(penalizeGapLetters);
	double coverage = asReal(minCoverage);
	if (coverage >= 0) {
		useMax = 0;
	} else {
		useMax = 1;
		coverage *= -1;
	}
	int mode = asInteger(method);
	int global = tGaps != 0 && pGapLetters != 0;
	int nthreads = asInteger(nThreads);
	l = length(y);
	
	XStringSet_holder s_set = hold_XStringSet(seqs);
	Chars_holder s_x = get_elt_from_XStringSet_holder(&s_set, asInteger(u1) - 1);
	int *u = INTEGER(u2); // sequences indices
	int w1 = s_x.length;
	
	// fill the vector of divisors
	int divisors[wS];
	divisors[0] = 1;
	for (i = 1; i < wS; i++)
		divisors[i] = divisors[i - 1]*maxA;
	
	double *sM = REAL(subMatrix);
	double dS = asReal(dropScore);
	XStringSet_holder l_set = hold_XStringSet(letters);
	Chars_holder l_i = get_elt_from_XStringSet_holder(&l_set, 0);
	int *lkup_row = (int *) malloc(256*sizeof(int)); // thread-safe on Windows
	int *lkup_col = (int *) malloc(256*sizeof(int)); // thread-safe on Windows
	for (i = 0; i < 256; i++) {
		lkup_row[i] = NA_INTEGER;
		lkup_col[i] = NA_INTEGER;
	}
	for (i = 0; i < l_i.length; i++) {
		lkup_row[(unsigned char)l_i.ptr[i]] = i;
		lkup_col[(unsigned char)l_i.ptr[i]] = i*l_i.length;
	}
	
	SEXP sims;
	PROTECT(sims = allocVector(REALSXP, l));
	double *sim = REAL(sims);
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, l + 1));
	
	int **pY = (int **) malloc(l*sizeof(int *)); // thread-safe on Windows
	int **pOY = (int **) malloc(l*sizeof(int *)); // thread-safe on Windows
	int *ly = (int *) malloc(l*sizeof(int)); // thread-safe on Windows
	
	for (i = 0; i < l; i++) {
		pY[i] = INTEGER(VECTOR_ELT(VECTOR_ELT(v, I[i] - 1), 0));
		pOY[i] = INTEGER(VECTOR_ELT(VECTOR_ELT(v, I[i] - 1), 1));
		ly[i] = length(VECTOR_ELT(VECTOR_ELT(v, I[i] - 1), 0));
	}
	
	int maxX = 0;
	int *ox = (int *) malloc(lx*sizeof(int)); // thread-safe on Windows
	for (i = 0; i < lx; i++) {
		if (OX[i] > maxX)
			maxX = OX[i];
		ox[i] = OX[i] - 1;
	}
	
	int **pt = (int **) malloc(l*sizeof(int *)); // thread-safe on Windows
	int *N = (int *) malloc(l*sizeof(int)); // thread-safe on Windows
	int **keeps = (int **) malloc(l*sizeof(int *)); // thread-safe on Windows
	
	#ifdef _OPENMP
	#pragma omp parallel num_threads(nthreads)
	{
	#endif
		int *m = (int *) calloc(maxX, sizeof(int)); // initialized to zero (thread-safe on Windows)
		
		#ifdef _OPENMP
		#pragma omp for private(i,j,k,n,p,new,t,keep,d,pos,p1,p2,t1,t2,ov,OV,off,g1,g2,g,o,count)
		#endif
		for (i = 0; i < l; i++) {
			int *Y = pY[i];
			int *OY = pOY[i];
			Chars_holder s_y = get_elt_from_XStringSet_holder(&s_set, u[i] - 1);
			int w2 = s_y.length;
			
			int goalSize = (int)ceil((logN + log2((w1 > w2) ? w1 : w2))/logE);
			int divisor;
			if (goalSize > wS) {
				goalSize = wS;
				divisor = 1;
			} else if (goalSize < 1) {
				goalSize = 1;
				divisor = divisors[wS - 1];
			} else { // 1 <= goalSize <= wordSize
				divisor = divisors[wS - goalSize];
			}
			
			j = 0;
			k = 0;
			if (divisor == 1) { // goalSize == wordSize
				while (j < lx && k < ly[i]) {
					if (X[j] == Y[k]) {
						m[ox[j]] = OY[k];
						j++;
						k++;
					} else if (X[j] < Y[k]) {
						do {
							j++;
						} while (j < lx && X[j] < Y[k]);
					} else {
						do {
							k++;
						} while (k < ly[i] && Y[k] < X[j]);
					}
				}
			} else {
				while (j < lx && k < ly[i]) {
					if (X[j]/divisor == Y[k]/divisor) {
						m[ox[j]] = OY[k];
						j++;
						k++;
					} else if (X[j]/divisor < Y[k]/divisor) {
						do {
							j++;
						} while (j < lx && X[j]/divisor < Y[k]/divisor);
					} else {
						do {
							k++;
						} while (k < ly[i] && Y[k]/divisor < X[j]/divisor);
					}
				}
			}
			
			// count the number of anchors
			n = 0;
			new = 1;
			int one, two;
			for (j = 0; j < maxX; j++) {
				two = m[j];
				if (two > 0) {
					if (new) {
						new = 0;
						n++;
					} else if (two != one) {
						n++;
					}
					one = two + 1;
				} else {
					new = 1;
				}
			}
			
			// record anchor ranges
			t = (int *) malloc(3*n*sizeof(int)); // thread-safe on Windows
			k = -1;
			new = 1;
			p = -1;
			for (j = 0; j < maxX; j++) {
				two = m[j];
				if (two > 0) {
					if (new ||
						two != one) {
						new = 0;
						k++;
						t[++p] = j + 1;
						t[++p] = m[j];
						t[++p] = goalSize;
					} else {
						t[p]++;
					}
					one = two + 1;
					m[j] = 0; // clear memory for next i
				} else {
					if (k == n - 1)
						break;
					new = 1;
				}
			}
			
			// separate overlapping regions
			j = 1;
			while (j < n) {
				k = j - 1;
				while (k >= 0 && k > j - goalSize) {
					int d1 = t[0 + 3*j] - t[0 + 3*k] - t[2 + 3*k];
					int d2 = t[1 + 3*j] - t[1 + 3*k] - t[2 + 3*k];
					if (d1 < 0 || d2 < 0) {
						if (d1 < d2) {
							if (t[2 + 3*k] > -d1)
								t[2 + 3*k] += d1;
						} else if (t[2 + 3*k] > -d2) {
							t[2 + 3*k] += d2;
						}
					}
					k--;
				}
				j++;
			}
			
			N[i] = n;
			pt[i] = t;
			
			// chain anchors
			int *b = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
			int *s = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
			for (j = 0; j < n; j++) {
				b[j] = -1;
				s[j] = t[2 + 3*j];
			}
			j = 1;
			p = 0;
			while (j < n) {
				for (k = 0; k < j; k++) {
					if (t[1 + 3*k] < t[1 + 3*j] && // starts within bounds
						((t[1 + 3*k] + t[2 + 3*k] <= t[1 + 3*j] &&
						t[0 + 3*k] + t[2 + 3*k] <= t[0 + 3*j]) ||
						t[0 + 3*j] - t[0 + 3*k] == t[1 + 3*j] - t[1 + 3*k])) { // ends are compatible
						if (s[k] + t[2 + 3*j] > s[j]) { // extend chain
							s[j] = s[k] + t[2 + 3*j];
							b[j] = k;
						}
					}
				}
				if (s[j] > s[p]) // higher score
					p = j;
				j++;
			}
			free(s);
			
			// rectify anchors
			keep = (int *) calloc(n, sizeof(int)); // initialized to zero (thread-safe on Windows)
			keeps[i] = keep;
			if (n > 0) {
				while (p >= 0) { // traceback
					keep[p] = 1;
					p = b[p];
				}
			}
			free(b);
			k = -1;
			j = 0;
			pos = 0; // number of matching positions in anchors
			g1 = 0;
			g2 = 0;
			count = 0;
			double tempScore;
			int matches, c1, c2, lkup1, lkup2, bound1 = 0, bound2 = 0;
			while (j < n) {
				if (keep[j]) {
					if (k >= 0 &&
						t[0 + 3*k] + t[2 + 3*k] >= t[0 + 3*j] &&
						t[0 + 3*j] - t[0 + 3*k] == t[1 + 3*j] - t[1 + 3*k]) {
						// merge anchors
						keep[j] = 0;
						pos -= t[2 + 3*k];
						t[2 + 3*k] = t[2 + 3*j] + t[0 + 3*j] - t[0 + 3*k];
						pos += t[2 + 3*k];
					} else {
						if (k >= 0) { // previous range
							// remove overlap
							d = t[0 + 3*k] + t[2 + 3*k] - t[0 + 3*j];
							if (d > 0) {
								t[0 + 3*j] += d;
								t[1 + 3*j] += d;
								t[2 + 3*j] -= d;
							}
							d = t[1 + 3*k] + t[2 + 3*k] - t[1 + 3*j];
							if (d > 0) {
								t[0 + 3*j] += d;
								t[1 + 3*j] += d;
								t[2 + 3*j] -= d;
							}
							
							// record gaps
							g = t[1 + 3*j]; // seq2 current start
							g -= t[1 + 3*k] + t[2 + 3*k] - 1; // seq2 previous end
							g -= t[0 + 3*j]; // seq1 current start
							g += t[0 + 3*k] + t[2 + 3*k] - 1; // seq1 previous end
							if (g > 0) {
								g2 -= g;
								count++;
							} else if (g < 0) {
								g1 += g;
								count++;
							} // else g = 0
						} else {
							p1 = t[0 + 3*j]; // starting position in seq1
							p2 = t[1 + 3*j]; // starting position in seq2
						}
						pos += t[2 + 3*j];
						
						// try extending left
						c1 = t[0 + 3*j] - 1;
						c2 = t[1 + 3*j] - 1;
						tempScore = 0;
						matches = 0;
						while (c1 > bound1 && c2 > bound2 && tempScore > dS) {
							c1--;
							c2--;
							lkup1 = lkup_row[(unsigned char)s_x.ptr[c1]];
							lkup2 = lkup_col[(unsigned char)s_y.ptr[c2]];
							if (lkup1 != NA_INTEGER && lkup2 != NA_INTEGER) {
								tempScore += sM[lkup1 + lkup2];
								if (s_x.ptr[c1] == s_y.ptr[c2])
									matches++;
							} else { // encountered missing character
								break;
							}
							if (tempScore > 0) {
								tempScore = 0;
								pos += matches;
								matches = 0;
								t[2 + 3*j] += t[0 + 3*j] - 1 - c1;
								t[0 + 3*j] = c1 + 1;
								t[1 + 3*j] = c2 + 1;
							}
						}
						
						if (k >= 0) { // try extending right
							bound1 = t[0 + 3*j] - 2;
							bound2 = t[1 + 3*j] - 2;
							c1 = t[0 + 3*k] + t[2 + 3*k] - 2;
							c2 = t[1 + 3*k] + t[2 + 3*k] - 2;
							tempScore = 0;
							matches = 0;
							while (c1 < bound1 && c2 < bound2 && tempScore > dS) {
								c1++;
								c2++;
								lkup1 = lkup_row[(unsigned char)s_x.ptr[c1]];
								lkup2 = lkup_col[(unsigned char)s_y.ptr[c2]];
								if (lkup1 != NA_INTEGER && lkup2 != NA_INTEGER) {
									tempScore += sM[lkup1 + lkup2];
									if (s_x.ptr[c1] == s_y.ptr[c2])
										matches++;
								} else { // encountered missing character
									break;
								}
								if (tempScore > 0) {
									tempScore = 0;
									pos += matches;
									matches = 0;
									t[2 + 3*k] = c1 - t[0 + 3*k] + 2;
								}
							}
						}
						
						k = j;
					}
					bound1 = t[0 + 3*k] + t[2 + 3*k] - 1;
					bound2 = t[1 + 3*k] + t[2 + 3*k] - 1;
				}
				j++;
			}
			
			// compute similarity
			if (k >= 0) { // at least one range
				bound1 = w1 - 1;
				bound2 = w2 - 1;
				c1 = t[0 + 3*k] + t[2 + 3*k] - 2;
				c2 = t[1 + 3*k] + t[2 + 3*k] - 2;
				tempScore = 0;
				matches = 0;
				while (c1 < bound1 && c2 < bound2 && tempScore > dS) {
					c1++;
					c2++;
					lkup1 = lkup_row[(unsigned char)s_x.ptr[c1]];
					lkup2 = lkup_col[(unsigned char)s_y.ptr[c2]];
					if (lkup1 != NA_INTEGER && lkup2 != NA_INTEGER) {
						tempScore += sM[lkup1 + lkup2];
						if (s_x.ptr[c1] == s_y.ptr[c2])
							matches++;
					} else { // encountered missing character
						break;
					}
					if (tempScore > 0) {
						tempScore = 0;
						pos += matches;
						matches = 0;
						t[2 + 3*k] = c1 - t[0 + 3*k] + 2;
					}
				}
				
				t1 = w1 - (t[0 + 3*k] + t[2 + 3*k] - 1); // remaining positions in seq1
				t2 = w2 - (t[1 + 3*k] + t[2 + 3*k] - 1); // remaining positions in seq2
				
				if (global &&
					pGapLetters == 1) {
					if (mode == 1) { // overlap
						if (p1 >= p2) {
							o = 1;
							if (t1 >= t2) { // 2 within 1
								ov = w1;
							} else { // end of 1 overlaps start of 2
								ov = w1 + t2 - t1;
							}
						} else {
							o = 2;
							if (t2 >= t1) { // 1 within 2
								ov = w2;
							} else { // end of 2 overlaps start of 1
								ov = w2 + t1 - t2;
							}
						}
					} else if (mode == 2) { // shortest
						if (w1 < w2) {
							ov = w1;
							o = 1;
						} else {
							ov = w2;
							o = 2;
						}
					} else { // longest
						if (w1 < w2) {
							ov = w2;
							o = 2;
						} else {
							ov = w1;
							o = 1;
						}
					}
				}
				
				if (coverage != 0 ||
					global == 0 ||
					pGapLetters != 1) {
					if (t1 > t2) {
						t1 = t1 - t2;
						t2 = 0;
					} else if (t2 > t1) {
						t2 = t2 - t1;
						t1 = 0;
					} else {
						t1 = 0;
						t2 = 0;
					}
					if (global && pGapLetters != 0) // pGapLetters = NA
						count += (t1 != t2) + (p1 != p2);
					if (p1 <= p2 && t1 <= t2) { // 1 within 2
						OV = w1;
						off = g1;
						if (global == 0 || pGapLetters != 1)
							o = 1;
					} else if (p2 <= p1 && t2 <= t1) { // 2 within 1
						OV = w2;
						off = g2;
						if (global == 0 || pGapLetters != 1)
							o = 2;
					} else if (p1 > p2) { // end of 1 overlaps start of 2
						if (w1 - p1 + p2 > w2 - t2) {
							OV = w1 - p1 + p2;
							off = g1;
							if (global == 0 || pGapLetters != 1)
								o = 1;
						} else {
							OV = w2 - t2;
							off = g2;
							if (global == 0 || pGapLetters != 1)
								o = 2;
						}
					} else { // end of 2 overlaps start of 1
						if (w2 - p2 + p1 > w1 - t1) {
							OV = w2 - p2 + p1;
							off = g2;
							if (global == 0 || pGapLetters != 1)
								o = 2;
						} else {
							OV = w1 - t1;
							off = g1;
							if (global == 0 || pGapLetters != 1)
								o = 1;
						}
					}
					
					if (global == 0 || pGapLetters != 1)
						ov = OV;
				}
				
				if ((useMax == 0 &&
					((w2 <= w1 && (double)(OV + off)/(double)w2 < coverage) ||
					(w1 <= w2 && (double)(OV + off)/(double)w1 < coverage))) ||
					(useMax == 1 &&
					((w2 <= w1 && (double)(OV + off)/(double)w1 < coverage) ||
					(w1 <= w2 && (double)(OV + off)/(double)w2 < coverage)))) {
					sim[i] = 0;
				} else {
					if (pGapLetters == 1) {
						if (o == 1) {
							sim[i] = (double)pos/((double)(ov - g2));
						} else {
							sim[i] = (double)pos/((double)(ov - g1));
						}
					} else if (pGapLetters == 0) {
						if (o == 1) {
							sim[i] = (double)pos/((double)(ov + g1));
						} else {
							sim[i] = (double)pos/((double)(ov + g2));
						}
					} else {
						if (o == 1) {
							sim[i] = (double)pos/((double)(ov + count + g1));
						} else {
							sim[i] = (double)pos/((double)(ov + count + g2));
						}
					}
				}
			} else {
				sim[i] = 0;
			}
		}
		
		free(m);
	#ifdef _OPENMP
	}
	#endif
	free(ox);
	free(pY);
	free(pOY);
	free(ly);
	
	// second loop not tread-safe
	for (i = 0; i < l; i++) {
		t = pt[i];
		n = N[i];
		keep = keeps[i];
		
		// record final anchors
		k = 0;
		for (j = 0; j < n; j++)
			if (keep[j])
				k++;
		SEXP ans;
		PROTECT(ans = allocMatrix(INTSXP, 4, k));
		int *rans = INTEGER(ans);
		k = -1;
		for (j = 0; j < n; j++) {
			if (keep[j]) {
				k++;
				rans[0 + 4*k] = t[0 + 3*j];
				rans[1 + 4*k] = t[0 + 3*j] + t[2 + 3*j] - 1;
				rans[2 + 4*k] = t[1 + 3*j];
				rans[3 + 4*k] = t[1 + 3*j] + t[2 + 3*j] - 1;
			}
		}
		free(keep);
		free(t);
		
		SET_VECTOR_ELT(ret_list, i, ans);
		UNPROTECT(1);
	}
	free(pt);
	free(N);
	free(keeps);
	
	SET_VECTOR_ELT(ret_list, l, sims);
	UNPROTECT(2);
	
	return ret_list;
}

// difference in overlap between pairs
SEXP overlap(SEXP res, SEXP widths1, SEXP widths2)
{
	int i, l, *x;
	int w1 = asInteger(widths1);
	int *w2 = INTEGER(widths2);
	int n = length(res);
	
	SEXP ans;
	PROTECT(ans = allocVector(INTSXP, n));
	int *rans = INTEGER(ans);
	
	for (i = 0; i < n; i++) {
		x = INTEGER(VECTOR_ELT(res, i)); // matrix of anchor ranges
		l = length(VECTOR_ELT(res, i)); // a multiple of 4
		
		rans[i] = 1;
		if (l == 0) {
			rans[i] += w1 + w2[i];
		} else {
			if (x[0] > x[2]) {
				rans[i] += x[0] - x[2];
			} else {
				rans[i] += x[2] - x[0];
			}
			int one = w1 - x[l - 3];
			int two = w2[i] - x[l - 1];
			if (two > one) {
				rans[i] += two - one;
			} else {
				rans[i] += one - two;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// in-place addition of cophenetic distances
SEXP cophenetic(SEXP Index1, SEXP N, SEXP D, SEXP H)
{
	int i, j, val;
	int *I = INTEGER(Index1);
	int l1 = length(Index1);
	int n = asInteger(N);
	double *d = REAL(D);
	double h = asReal(H);
	
	char *t = Calloc(n, char);
	for (i = 0; i < l1; i++)
		t[I[i] - 1] = 1;
	int l2 = n;
	for (i = 0; i < n; i++)
		if (t[i])
			l2--;
	int *J = Calloc(l2, int);
	j = 0;
	for (i = 0; i < n; i++)
		if (!t[i])
			J[j++] = i + 1;
	Free(t);
	
	for (i = 0; i < l1; i++) {
		for (j = 0; j < l2; j++) {
			if (I[i] < J[j]) {
				val = n*(I[i] - 1) - I[i]*(I[i] - 1)/2 + J[j] - I[i] - 1;
			} else {
				val = n*(J[j] - 1) - J[j]*(J[j] - 1)/2 + I[i] - J[j] - 1;
			}
			d[val] += h;
		}
	}
	Free(J);
	
	return D;
}

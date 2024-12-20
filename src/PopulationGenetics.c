/****************************************************************************
 *                      Population Genetics Functions                       *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

SEXP correlationProfile(SEXP x, SEXP readingFrame, SEXP maxN, SEXP verbose, SEXP pBar)
{
	int i, j, k, p, I, J, v, before, *rPercentComplete;
	double weight, soFar, tot;
	int *rF = INTEGER(readingFrame);
	int N = asInteger(maxN);
	v = asLogical(verbose);
	
	XStringSet_holder x_set;
	R_xlen_t x_length;
	Chars_holder x_i, x_j;
	
	SEXP percentComplete, utilsPackage;
	if (v) { // percent complete variables
		soFar = 0;
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
		
		tot = 0;
		for (i = 0; i < length(x); i++) {
			x_set = hold_XStringSet(VECTOR_ELT(x, i));
			x_length = get_length_from_XStringSet_holder(&x_set);
			if (x_length > 0)
				tot += ((double)x_length*(double)(x_length - 1))/2;
		}
	}
	
	// initialize results
	SEXP ans1, ans2, ans3;
	PROTECT(ans1 = allocVector(REALSXP, N));
	double *P = REAL(ans1);
	PROTECT(ans2 = allocVector(REALSXP, N));
	double *C = REAL(ans2);
	PROTECT(ans3 = allocVector(REALSXP, 2));
	double *subs = REAL(ans3);
	
	// clear memory
	subs[0] = 0;
	subs[1] = 0;
	for (i = 0; i < N; i++) {
		P[i] = 0;
		C[i] = 0;
	}
	
	for (i = 0; i < length(x); i++) {
		x_set = hold_XStringSet(VECTOR_ELT(x, i));
		x_length = get_length_from_XStringSet_holder(&x_set);
		
		weight = 1/((double)x_length - 1);
		I = 0;
		while (I < x_length - 1) {
			x_i = get_elt_from_XStringSet_holder(&x_set, I);
			J = I + 1;
			while (J < x_length) {
				x_j = get_elt_from_XStringSet_holder(&x_set, J);
				
				for (j = 0; j < x_i.length && j < x_j.length; j++) {
					if (x_i.ptr[j] < 16) { // base
						subs[1] += weight;
						
						if (!((x_i.ptr[j]) & (x_j.ptr[j]))) { // unequal
							subs[0] += weight;
							
							k = j + 1;
							if (rF[i] == NA_INTEGER) {
								p = 0;
							} else {
								p = ((k - rF[i] + 1) % 3);
							}
							while (p < N && k < x_i.length && k < x_j.length) {
								if ((x_i.ptr[k] == 16 || x_i.ptr[k] == 64) &&
									(x_j.ptr[k] == 16 || x_j.ptr[k] == 64)) { // both gap
									k += 3; // maintain reading frame
									continue;
								}
								if (x_i.ptr[k] >= 16 || x_j.ptr[k] >= 16) // non-base
									break;
								if (!((x_i.ptr[k]) & (x_j.ptr[k]))) // unequal
									P[p] += weight;
								C[p] += weight;
								k++;
								p++;
							}
						}
					}
				}
				J++;
			}
			I++;
			
			if (v) { // print the percent completed so far
				soFar += (double)(J - I);
				*rPercentComplete = floor(100*soFar/tot);
				if (*rPercentComplete > before) { // when the percent has changed
					// tell the progress bar to update in the R console
					eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
					before = *rPercentComplete;
				}
			} else {
				R_CheckUserInterrupt();
			}
		}
	}
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	SET_VECTOR_ELT(ret_list, 2, ans3);
	
	if (v) {
		UNPROTECT(6);
	} else {
		UNPROTECT(4);
	}
	
	return ret_list;
}

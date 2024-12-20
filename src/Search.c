/****************************************************************************
 *               Utilities for Searching an Inverted Index                  *
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

// for calloc/free
#include <stdlib.h>

// for math functions
#include <math.h>

// for time and difftime
#include <time.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

void heapSelect(double *score, int *A, int l, int k)
{
	int p, c, temp;
	int end = l;
	int start = end/2;
	while (end > l - k) {
		if (start > 0) {
			start--;
		} else {
			end--;
			// swap
			temp = A[0];
			A[0] = A[end];
			A[end] = temp;
		}
		p = start;
		c = 2*p + 1;
		while (c < end) {
			if (c + 1 < end && score[A[c]] < score[A[c + 1]])
				c++;
			if (score[A[p]] < score[A[c]]) {
				// swap
				temp = A[p];
				A[p] = A[c];
				A[c] = temp;
				p = c;
			} else {
				break;
			}
			c = 2*p + 1;
		}
	}
}

// returns hits between queries and targets in an inverted index
SEXP searchIndex(SEXP query, SEXP wordSize, SEXP stepSize, SEXP logFreqs, SEXP count, SEXP location, SEXP index, SEXP positions, SEXP sepC, SEXP gapC, SEXP total, SEXP minScore, SEXP scoreOnly, SEXP pattern, SEXP subject, SEXP subMatrix, SEXP letters, SEXP dropScore, SEXP limitTarget, SEXP limitQuery, SEXP alphabet, SEXP correction, SEXP background, SEXP iterations, SEXP threshold, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	int i, j, k, p, c;
	int n = length(query); // number of sequences in the query
	int K = asInteger(wordSize); // k-mer length
	int step = asInteger(stepSize); // separation between k-mers
	double *freqs = REAL(logFreqs); // -log of normalized letter frequencies
	int levels = length(logFreqs); // number of letters in the reduced alphabet
	int *num = INTEGER(count); // count of each k-mer in target
	int *loc = INTEGER(location); // location of k-mer in target
	int *ind = INTEGER(index); // index of target sequence
	int *pos = INTEGER(positions); // number of matchable positions in target
	int L = length(count); // number of possible of k-mers
	double tot = asReal(total); // total size of target database
	double minS = asReal(minScore); // minimum score or NA to calculate
	int sO = asInteger(scoreOnly); // FALSE to output anchor positions
	int limitT = asInteger(limitTarget); // maximum number of matches per subject
	int limitQ = asInteger(limitQuery); // maximum number of matches per pattern
	int *alpha = INTEGER(alphabet); // (possibly reduced) letter mapping to freqs
	int size = length(alphabet); // number of letters in standard alphabet
	double correct = asReal(correction); // weight of subMatrix background correction
	double *bg = REAL(background); // background frequencies of standard alphabet
	int maxIt = asInteger(iterations); // number of profile search iterations
	double thresh = asReal(threshold); // expect-value for profile inclusion
	int nthreads = asInteger(nThreads);
	
	// if subMatrix provided then query/target sequences present
	XStringSet_holder p_set, s_set, l_set;
	Chars_holder l_i;
	double *sM, dS;
	int *lkup_row, *lkup_col, *pwv;
	double *scores, *addScores;
	int sizeM = length(subMatrix);
	if (sizeM > 0) {
		l_set = hold_XStringSet(letters);
		l_i = get_elt_from_XStringSet_holder(&l_set, 0);
		if (l_i.length < size)
			error("Incorrect size of substitutionMatrix.");
		if (l_i.length*l_i.length != sizeM)
			error("Incorrect size of substitutionMatrix.");
		sM = REAL(subMatrix);
		dS = asReal(dropScore);
		p_set = hold_XStringSet(pattern);
		s_set = hold_XStringSet(subject);
		lkup_row = (int *) malloc(256*sizeof(int)); // thread-safe on Windows
		lkup_col = (int *) malloc(256*sizeof(int)); // thread-safe on Windows
		for (i = 0; i < 256; i++) {
			lkup_row[i] = NA_INTEGER;
			lkup_col[i] = NA_INTEGER;
		}
		for (i = 0; i < l_i.length; i++) {
			lkup_row[(unsigned char)l_i.ptr[i]] = i;
			lkup_col[(unsigned char)l_i.ptr[i]] = i*l_i.length;
		}
		if (maxIt > 0) {
			pwv = (int *) malloc(K*sizeof(int)); // thread-safe on Windows
			pwv[0] = 1;
			for (i = 1; i < K; i++)
				pwv[i] = pwv[i - 1]*levels;
		}
	} else {
		// calculate -log(expected k-mer frequency)
		scores = (double *) calloc(L, sizeof(double)); // initialized to zero (thread-safe on Windows)
		addScores = (double *) calloc(L, sizeof(double)); // initialized to zero (thread-safe on Windows)
		c = 1;
		for (i = 1; i <= K; i++) {
			k = -1;
			j = 0;
			p = 0;
			while (j < L) {
				if (j == p) {
					p += c;
					if (k == levels - 1) {
						k = 0;
					} else {
						k++;
					}
				}
				scores[j] += freqs[k];
				if (i > K - step)
					addScores[j] += freqs[k];
				j++;
			}
			c *= levels;
		}
	}
	
	// set up a timer to minimize interrupt checks
	time_t start, end;
	double elapsed;
	time(&start);
	
	int before, v, *rPercentComplete;
	double N, soFar;
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	if (v) { // percent complete variables
		N = 0; // completed iterations
		soFar = 0;
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	// determine cost for the distance between k-mer matches
	double sC = asReal(sepC); // cost for separation between k-mers
	double gC = asReal(gapC); // cost for minimum gaps between k-mers
	int maxSep = (int)sqrt((double)L); // one match expected by chance every maxSep unmasked k-mers
	double *sepCost = (double *) malloc((maxSep + 1)*sizeof(double)); // thread-safe on Windows
	double *gapCost = (double *) malloc((maxSep + 1)*sizeof(double)); // thread-safe on Windows
	for (i = 0; i <= maxSep; i++) {
		sepCost[i] = sqrt(i);
		gapCost[i] = gC*sepCost[i];
		sepCost[i] *= sC;
	}
	
	// determine the offset for each k-mer
	R_xlen_t *offset = (R_xlen_t *) malloc(L*sizeof(R_xlen_t)); // thread-safe on Windows
	offset[0] = 0;
	k = 0;
	for (i = 1; i < L; i++) {
		offset[i] = offset[k] + num[k];
		k = i;
	}
	
	// build threadsafe vectors for inputs and outputs
	double **vecs = (double **) malloc(n*sizeof(double *)); // thread-safe on Windows
	int **ptrs = (int **) malloc(n*sizeof(int *)); // thread-safe on Windows
	int *l = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	for (i = 0; i < n; i++) {
		ptrs[i] = INTEGER(VECTOR_ELT(query, i)); // query k-mers
		l[i] = length(VECTOR_ELT(query, i)); // query length
	}
	int ***matrices;
	if (sO == 0)
		matrices = (int ***) malloc(n*sizeof(int **)); // thread-safe on Windows
	
	int negK = -1*K;
	int abort = 0;
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,k,p,c) schedule(dynamic) num_threads(nthreads)
	#endif
	for (i = 0; i < n; i++) { // each sequence
		if (abort == 0) {
			int *w = ptrs[i]; // k-mers
			
			// record the number of unmasked positions
			int width = 0; // positions in the query
			k = -2; // last unmasked position
			for (j = 0; j < l[i]; j++) { // each k-mer
				if (w[j] != NA_INTEGER) { // unmasked
					if (k == j - 1) {
						width++; // query is always staggered by one position
					} else {
						// all masked (NA) positions are at least K long
						width += K; // new k-mer
					}
					k = j;
				}
			}
			
			// count target occurrences of each query k-mer
			int s = 0; // total number of k-mers shared with targets
			int *counts = (int *) malloc(l[i]*sizeof(int)); // thread-safe on Windows
			for (j = 0; j < l[i]; j++) { // each query k-mer
				if (w[j] == NA_INTEGER) {
					counts[j] = 0;
				} else {
					counts[j] = num[w[j]];
					s += counts[j];
					if (s < 0) { // signed integer overflow
						abort = i + 1;
						break;
					}
				}
			}
			if (abort != 0 // exit
				|| width == 0 // no query
				|| s == 0) { // no targets
				free(counts);
				
				int *set = (int *) malloc(0*sizeof(int)); // thread-safe on Windows
				double *score = (double *) malloc(0*sizeof(double)); // thread-safe on Windows
				ptrs[i] = set;
				vecs[i] = score;
				l[i] = 0; // number of results
				if (sO == 0) {
					int **anchors = (int **) malloc(0*sizeof(int *)); // thread-safe on Windows
					matrices[i] = anchors;
				}
				continue;
			}
			
			// record target occurrences of each query k-mer
			int *posQuery = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
			int *posTarget = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
			int *set = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
			int *set2;
			double *score, *addScore, *score2, *addScore2;
			if (sizeM > 0) {
				s = 0;
				for (j = 0; j < l[i]; j++) { // each query k-mer
					if (counts[j] > 0) { // w[j] != NA_INTEGER
						R_xlen_t P = offset[w[j]]; // position of k-mer in loc and ind
						for (k = 0; k < counts[j]; k++) { // each target instance of k-mer
							posQuery[s] = j + 1;
							posTarget[s] = loc[P];
							set[s] = ind[P];
							P++;
							s++;
						}
					}
				}
			} else {
				score = (double *) malloc(s*sizeof(double)); // thread-safe on Windows
				addScore = (double *) malloc(s*sizeof(double)); // thread-safe on Windows
				s = 0;
				for (j = 0; j < l[i]; j++) { // each query k-mer
					if (counts[j] > 0) { // w[j] != NA_INTEGER
						R_xlen_t P = offset[w[j]]; // position of k-mer in loc and ind
						for (k = 0; k < counts[j]; k++) { // each target instance of k-mer
							posQuery[s] = j + 1;
							posTarget[s] = loc[P];
							set[s] = ind[P];
							score[s] = scores[w[j]];
							addScore[s] = addScores[w[j]];
							P++;
							s++;
						}
					}
				}
			}
			
			double *subScores; // substitution scores (if maxIt > 0)
			int lkup; // index of letter
			if (sizeM > 0) {
				Chars_holder p_i;
				p_i = get_elt_from_XStringSet_holder(&p_set, i);
				
				double *sMcorr; // background corrected substitution matrix
				if (correct > 0) { // apply composition correction to substitution matrix
					double *fg = (double *) malloc(size*sizeof(double)); // thread-safe on Windows
					for (j = 0; j < size; j++)
						fg[j] = correct*bg[j]; // initialize to pseudocount
					for (j = 0; j < p_i.length; j++) {
						lkup = lkup_row[(unsigned char)p_i.ptr[j]];
						if (lkup >= 0 && lkup < size)
							fg[lkup]++;
					}
					double weight = 0;
					for (j = 0; j < size; j++)
						weight += fg[j];
					for (j = 0; j < size; j++)
						fg[j] /= weight;
					sMcorr = (double *) malloc(sizeM*sizeof(double)); // thread-safe on Windows
					for (j = 0; j < sizeM; j++)
						sMcorr[j] = sM[j];
					for (j = 0; j < size; j++)
						for (k = 0; k < size; k++)
							sMcorr[j*l_i.length + k] += log(bg[j]*bg[k]/(fg[j]*fg[k]));
					free(fg);
				} else {
					sMcorr = sM;
				}
				
				// calculate substitution scores
				subScores = (double *) malloc(l_i.length*p_i.length*sizeof(double)); // thread-safe on Windows
				for (j = 0; j < p_i.length; j++) {
					if (j < l[i] && w[j] == NA_INTEGER) {
						lkup = NA_INTEGER; // masked position
					} else {
						lkup = lkup_col[(unsigned char)p_i.ptr[j]];
					}
					if (lkup == NA_INTEGER) {
						for (k = 0; k < l_i.length; k++)
							subScores[j*l_i.length + k] = NA_REAL;
					} else {
						for (k = 0; k < l_i.length; k++)
							subScores[j*l_i.length + k] = sMcorr[lkup + k];
					}
				}
				if (correct > 0)
					free(sMcorr);
			}
			
			int *chain; // set of joined hits per match
			int *len; // length of each hit
			int *res; // indices of results
			c = l[i]; // current number of batches
			for (int it = 0; it <= maxIt; it++) { // each search iteration
				// merge sort on set then posTarget
				int *o1 = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				for (j = 0; j < s; j++)
					o1[j] = j; // initial order
				int *o2; // pointer to previous order
				for (j = 1; j < c; j++)
					counts[j] += counts[j - 1]; // batch size
				int group1, group2;
				while (c > 1) {
					o2 = o1; // store previous order
					o1 = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
					p = 0; // position in o1
					j = 0; // position in group 1
					group1 = 0; // group 1 in counts
					group2 = 1; // group 2 in counts
					// sort pairs of groups
					while (group2 < c) {
						k = counts[group1]; // position in group 2
						while (j < counts[group1] && k < counts[group2]) {
							if (set[o2[j]] == set[o2[k]]) { // apply tiebreaker
								if (posTarget[o2[j]] <= posTarget[o2[k]]) {
									o1[p++] = o2[j++];
								} else {
									o1[p++] = o2[k++];
								}
							} else if (set[o2[j]] < set[o2[k]]) {
								o1[p++] = o2[j++];
							} else {
								o1[p++] = o2[k++];
							}
						}
						while (j < counts[group1])
							o1[p++] = o2[j++];
						while (k < counts[group2])
							o1[p++] = o2[k++];
						group1 = group2 + 1;
						group2 = group1 + 1;
						j = k;
					}
					while (p < s) {
						o1[p] = o2[p];
						p++;
					}
					
					// merge pairs of groups
					p = 0;
					k = 0;
					while (k < c - 1) {
						k += 2;
						counts[p++] = counts[k - 1];
					}
					if (k < c)
						counts[p++] = counts[k++];
					c = p;
					free(o2);
				}
				free(counts);
				
				// reorder all vectors
				int *posQuery2 = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				int *posTarget2 = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				set2 = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				if (sizeM > 0) {
					for (j = 0; j < s; j++) {
						posQuery2[j] = posQuery[o1[j]];
						posTarget2[j] = posTarget[o1[j]];
						set2[j] = set[o1[j]];
					}
				} else {
					score2 = (double *) malloc(s*sizeof(double)); // thread-safe on Windows
					addScore2 = (double *) malloc(s*sizeof(double)); // thread-safe on Windows
					for (j = 0; j < s; j++) {
						posQuery2[j] = posQuery[o1[j]];
						posTarget2[j] = posTarget[o1[j]];
						set2[j] = set[o1[j]];
						score2[j] = score[o1[j]];
						addScore2[j] = addScore[o1[j]];
					}
					free(score);
					free(addScore);
					score = score2;
					addScore = addScore2;
				}
				free(o1);
				free(posQuery);
				free(posTarget);
				free(set);
				posQuery = posQuery2;
				posTarget = posTarget2;
				set = set2;
				
				// collapse adjacent hits
				len = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				int *origin = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				char *keep = (char *) malloc(s*sizeof(char)); // thread-safe on Windows
				for (j = 0; j < s; j++) {
					len[j] = K;
					origin[j] = j;
					keep[j] = 1;
				}
				c = 0; // current position
				k = 0; // previous position
				p = s; // remaining vector length
				int deltaQuery, deltaTarget;
				while (c < s - 1) {
					c++; // advance current position
					if (set[k] == set[c]) { // same index
						j = k;
						while (j < c) {
							deltaTarget = posTarget[c] - posTarget[j];
							if (deltaTarget > step) {
								k = j + 1; // advance previous position
							} else if (deltaTarget == step) {
								deltaQuery = posQuery[c] - posQuery[j];
								if (deltaQuery == step) { // merge positions
									p--;
									keep[c] = 0;
									origin[c] = origin[j];
									len[origin[j]] = len[origin[j]] + step;
									if (sizeM == 0)
										score[origin[j]] = score[origin[j]] + addScore[c];
									break; // done merging
								}
							} else { // deltaTarget == 0L
								break; // reached same target position
							}
							j++;
						}
					} else {
						k = c;
					}
				}
				if (sizeM == 0)
					free(addScore);
				free(origin);
				
				// resize all vectors
				posQuery2 = (int *) malloc(p*sizeof(int)); // thread-safe on Windows
				posTarget2 = (int *) malloc(p*sizeof(int)); // thread-safe on Windows
				set2 = (int *) malloc(p*sizeof(int)); // thread-safe on Windows
				int *len2 = (int *) malloc(p*sizeof(int)); // thread-safe on Windows
				score2 = (double *) malloc(p*sizeof(double)); // thread-safe on Windows
				p = 0;
				if (sizeM > 0) {
					for (j = 0; j < s; j++) {
						if (keep[j]) {
							posQuery2[p] = posQuery[j];
							posTarget2[p] = posTarget[j];
							set2[p] = set[j];
							len2[p] = len[j];
							p++;
						}
					}
				} else {
					for (j = 0; j < s; j++) {
						if (keep[j]) {
							posQuery2[p] = posQuery[j];
							posTarget2[p] = posTarget[j];
							set2[p] = set[j];
							len2[p] = len[j];
							score2[p] = score[j];
							p++;
						}
					}
					free(score);
				}
				s = p;
				free(keep);
				free(posQuery);
				free(posTarget);
				free(set);
				free(len);
				posQuery = posQuery2;
				posTarget = posTarget2;
				set = set2;
				len = len2;
				score = score2;
				
				double prevScore, tempScore;
				int gap, sep, temp;
				if (sizeM > 0) { // rescore and extend hits
					Chars_holder p_i, s_j;
					p_i = get_elt_from_XStringSet_holder(&p_set, i);
					int p1, p2; // position in query or target
					int prev = 0; // index previous sequence in target
					int bound; // boundary in target sequence
					int maxLen; // maximum observed length
					
					// score and extend hits
					for (j = 0; j < s; j++) {
						if (prev != set[j]) {
							s_j = get_elt_from_XStringSet_holder(&s_set, set[j] - 1);
							prev = set[j];
							maxLen = 0;
						}
						
						// score hit
						score[j] = 0;
						p1 = posQuery[j] + len[j] - 1;
						p2 = posTarget[j] + len[j] - 1;
						while (p1 >= posQuery[j]) {
							p1--;
							p2--;
							lkup = lkup_row[(unsigned char)s_j.ptr[p2]];
							if (lkup != NA_INTEGER && subScores[p1*l_i.length + lkup] != NA_REAL)
								score[j] += subScores[p1*l_i.length + lkup];
						}
						
						// try extending left
						k = j - 1;
						bound = 0; // left bound
						while (k >= 0 && set[k] == set[j]) {
							deltaTarget = posTarget[j] - posTarget[k] - len[k];
							if (deltaTarget >= maxSep + maxLen) {
								break;
							} else if (deltaTarget >= 0) {
								deltaQuery = posQuery[j] - posQuery[k] - len[k];
								if (deltaQuery >= 0 && deltaQuery <= maxSep) {
									if (deltaQuery > deltaTarget) {
										gap = deltaQuery - deltaTarget;
									} else {
										gap = deltaTarget - deltaQuery;
									}
									if (gap == 0) { // intersecting
										bound = posTarget[k] + len[k];
										break;
									}
								}
							}
							k--;
						}
						tempScore = 0;
						while (p1 > 0 && p2 > bound && tempScore > dS) {
							p1--;
							p2--;
							lkup = lkup_row[(unsigned char)s_j.ptr[p2]];
							if (lkup != NA_INTEGER && subScores[p1*l_i.length + lkup] != NA_REAL) {
								tempScore += subScores[p1*l_i.length + lkup];
							} else { // encountered missing character
								break;
							}
							if (tempScore > 0) { // new max score
								score[j] += tempScore;
								tempScore = 0;
								len[j] += posQuery[j] - p1 - 1;
								posQuery[j] = p1 + 1;
								posTarget[j] = p2 + 1;
							}
						}
						
						// try extending right
						k = j + 1;
						bound = s_j.length - 1; // right bound
						while (k < s && set[k] == set[j]) {
							deltaTarget = posTarget[k] - posTarget[j] - len[j];
							if (deltaTarget >= maxSep) {
								break;
							} else if (deltaTarget >= 0) {
								deltaQuery = posQuery[k] - posQuery[j] - len[j];
								if (deltaQuery >= 0 && deltaQuery <= maxSep) {
									if (deltaQuery > deltaTarget) {
										gap = deltaQuery - deltaTarget;
									} else {
										gap = deltaTarget - deltaQuery;
									}
									if (gap == 0) { // intersecting
										bound = posTarget[k] - 2;
										break;
									}
								}
							}
							k++;
						}
						tempScore = 0;
						p1 = posQuery[j] + len[j] - 1;
						p2 = posTarget[j] + len[j] - 1;
						while (p1 < p_i.length && p2 <= bound && tempScore > dS) {
							lkup = lkup_row[(unsigned char)s_j.ptr[p2]];
							if (lkup != NA_INTEGER && subScores[p1*l_i.length + lkup] != NA_REAL) {
								tempScore += subScores[p1*l_i.length + lkup];
							} else { // encountered missing character
								break;
							}
							if (tempScore > 0) { // new max score
								score[j] += tempScore;
								tempScore = 0;
								len[j] = p1 - posQuery[j] + 2;
							}
							p1++;
							p2++;
						}
						
						if (len[j] > maxLen)
							maxLen = len[j];
						
						// (re)order by posTarget
						p2 = j;
						p1 = p2 - 1;
						while (p2 > 0 && // within bounds
							set[p1] == set[p2] && // same index
							posTarget[p2] < posTarget[p1]) { // out of order
							temp = len[p1];
							len[p1] = len[p2];
							len[p2] = temp;
							temp = posQuery[p1];
							posQuery[p1] = posQuery[p2];
							posQuery[p2] = temp;
							temp = posTarget[p1];
							posTarget[p1] = posTarget[p2];
							posTarget[p2] = temp;
							tempScore = score[p1];
							score[p1] = score[p2];
							score[p2] = tempScore;
							p2 = p1;
							p1--;
						}
					}
					
					// pull back overlap
					for (j = 0; j < s; j++) {
						if (prev != set[j]) {
							s_j = get_elt_from_XStringSet_holder(&s_set, set[j] - 1);
							prev = set[j];
						}
						
						prevScore = 0; // best score after chaining
						double bestScore = 0; // score difference from overlap
						int bestLen = 0; // length difference from overlap
						k = j + 1;
						while (k < s && set[k] == set[j]) {
							deltaTarget = posTarget[j] + len[j] - posTarget[k]; // overlap
							if (deltaTarget > maxSep) { // cannot chain
								break;
							} else if (posQuery[k] > posQuery[j] && // can chain
								deltaTarget < len[k] && deltaTarget < len[j]) { // not fully overlapping
								deltaQuery = posQuery[j] + len[j] - posQuery[k]; // overlap
								if (deltaQuery <= maxSep && // can chain
									deltaQuery < len[k] && deltaQuery < len[j]) { // not fully overlapping
									if (deltaQuery > deltaTarget) {
										sep = deltaQuery;
										gap = deltaQuery - deltaTarget;
									} else {
										sep = deltaTarget;
										gap = deltaTarget - deltaQuery;
									}
									if (gap <= maxSep && // can chain
										score[j] + gapCost[gap] > 0 && // can chain
										score[k] + gapCost[gap] > 0) { // can chain
										if (sep > 0) { // score chaining after pulling back overlap
											p1 = posQuery[j] + len[j] - 1;
											p2 = posTarget[j] + len[j] - 1;
											p = p1 - sep;
											tempScore = 0;
											while (p1 > p) {
												p1--;
												p2--;
												lkup = lkup_row[(unsigned char)s_j.ptr[p2]];
												tempScore += subScores[p1*l_i.length + lkup];
											}
											if (score[k] + gapCost[gap] - tempScore > prevScore) {
												prevScore = score[k] + gapCost[gap] - tempScore;
												bestLen = sep;
												bestScore = tempScore;
											}
										} else { // score chaining without overlap
											sep *= -1;
											if (sep <= maxSep && // can chain
												score[k] + gapCost[gap] + sepCost[sep] > prevScore) {
												prevScore = score[k] + gapCost[gap] + sepCost[sep];
												bestLen = 0;
												bestScore = 0;
											}
										}
									}
								}
							}
							k++;
						}
						len[j] -= bestLen;
						score[j] -= bestScore;
					}
				}
				
				// determine significance of occurrences
				chain = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				origin = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				int *cov = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				for (j = 0; j < s; j++) {
					chain[j] = j;
					origin[j] = j;
					cov[j] = len[j] - 1;
				}
				j = 0; // last hit
				k = 1; // current hit
				int delta;
				while (k < s) {
					if (set[k] != set[j]) { // switched index
						j = k;
					} else {
						prevScore = score[k];
						p = j;
						while (p < k) {
							deltaTarget = posTarget[k] - posTarget[p] - len[p];
							if (deltaTarget > maxSep) {
								if (p == j)
									j = p + 1; // limit search space
							} else if (deltaTarget > negK) {
								deltaQuery = posQuery[k] - posQuery[p] - len[p];
								if (deltaQuery > negK && deltaQuery <= maxSep) {
									if (deltaTarget < 0 || deltaQuery < 0) { // overlap due to indels
										if (deltaQuery < deltaTarget) {
											delta = deltaQuery;
										} else {
											delta = deltaTarget;
										}
										deltaQuery -= delta;
										deltaTarget -= delta;
										
										if (deltaQuery <= maxSep && deltaTarget <= maxSep) {
											// deduct the approximate score of removed positions
											tempScore = score[p] + prevScore*(double)(len[k] + delta)/(double)len[k];
											if (tempScore > score[k]) {
												if (deltaQuery > deltaTarget) {
													gap = deltaQuery - deltaTarget;
													sep = deltaTarget;
												} else {
													gap = deltaTarget - deltaQuery;
													sep = deltaQuery;
												}
												tempScore = tempScore + gapCost[gap];
												tempScore = tempScore + sepCost[sep];
												if (tempScore > score[k]) {
													score[k] = tempScore;
													chain[k] = p;
													origin[k] = origin[p];
													cov[k] = len[k] - 1 + cov[p] + delta;
												}
											}
										}
									} else {
										tempScore = score[p] + prevScore;
										if (tempScore > score[k]) {
											if (deltaQuery > deltaTarget) {
												gap = deltaQuery - deltaTarget;
												sep = deltaTarget;
											} else {
												gap = deltaTarget - deltaQuery;
												sep = deltaQuery;
											}
											tempScore = tempScore + gapCost[gap];
											tempScore = tempScore + sepCost[sep];
											if (tempScore > score[k]) {
												score[k] = tempScore;
												chain[k] = p;
												origin[k] = origin[p];
												cov[k] = len[k] - 1 + cov[p];
											}
										}
									}
								}
							}
							p++;
						}
					}
					k++;
				}
				
				// correct for size of k-mer search space
				for (j = 0; j < s; j++) {
					temp = pos[set[j] - 1] - cov[j];
					if (temp > step) // correct for multiple target k-mers
						score[j] -= log(((double)temp)/((double)step));
					if (width > cov[j]) // correct for multiple query k-mers
						score[j] -= log((double)(width - cov[j]));
				}
				free(cov);
				
				// eliminate lower scoring chains with same origin
				int *maxOrigin; // pointer to maximum origin
				maxOrigin = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
				for (j = 0; j < s; j++) {
					if (origin[j] == j) {
						maxOrigin[j] = j;
					} else if (score[maxOrigin[origin[j]]] < score[j]) {
						maxOrigin[origin[j]] = j;
						maxOrigin[j] = NA_INTEGER;
					} else {
						maxOrigin[j] = NA_INTEGER;
					}
				}
				free(origin);
				keep = calloc(s, sizeof(char)); // initialized to zero (thread-safe on Windows)
				c = 0; // count of results
				for (j = 0; j < s; j++) {
					if (maxOrigin[j] != NA_INTEGER) {
						keep[maxOrigin[j]] = 1;
						c++;
					}
				}
				free(maxOrigin);
				res = (int *) malloc(c*sizeof(int)); // thread-safe on Windows
				p = 0;
				for (j = 0; j < s; j++)
					if (keep[j])
						res[p++] = j;
				free(keep);
				
				if (it < maxIt) { // repeat search with more k-mers
					Chars_holder p_i, s_j;
					p_i = get_elt_from_XStringSet_holder(&p_set, i);
					int p1, p2; // position in query or target
					int prev = 0; // index previous sequence in target
					double cutoff, weight;
					
					// calculate sequence profile
					double *profile = (double *) calloc(l_i.length*p_i.length, sizeof(double)); // initialized to zero (thread-safe on Windows)
					for (j = 0; j < c; j++) {
						cutoff = log((tot - (double)pos[set[res[j]] - 1])/(double)step/thresh);
						if (score[res[j]] > cutoff) {
							if (prev != set[res[j]]) {
								s_j = get_elt_from_XStringSet_holder(&s_set, set[res[j]] - 1);
								prev = set[res[j]];
							}
							weight = 1 - exp(cutoff - score[res[j]]); // probability of homology
							p = res[j];
							while (p != chain[p]) {
								p = chain[p];
								p1 = posQuery[p] + len[p] - 1;
								p2 = posTarget[p] + len[p] - 1;
								while (p1 >= posQuery[p]) {
									p1--;
									p2--;
									lkup = lkup_row[(unsigned char)s_j.ptr[p2]];
									if (lkup != NA_INTEGER)
										profile[p1*l_i.length + lkup] += weight;
								}
							}
						}
					}
					
					// convert profile to letters
					double temp_profile[l_i.length];
					int *patterns; // most frequent query k-mers
					int *patterns2; // alternative query k-mers
					patterns = (int *) calloc(p_i.length, sizeof(int)); // initialized to zero (thread-safe on Windows)
					patterns2 = (int *) calloc(p_i.length, sizeof(int)); // initialized to zero (thread-safe on Windows)
					for (j = 0; j < p_i.length; j++) {
						for (k = 0; k < levels; k++)
							temp_profile[k] = 0;
						weight = 0;
						for (k = 0; k < size; k++) {
							temp_profile[alpha[k]] += profile[j*l_i.length + k];
							weight += profile[j*l_i.length + k];
						}
						if (weight == 0) {
							patterns[j] = NA_INTEGER;
						} else {
							for (k = 0; k < levels; k++) {
								if (temp_profile[k] > 0) {
									temp_profile[k] /= weight; // normalize profile
									temp_profile[k] = log(temp_profile[k]);
									temp_profile[k] += freqs[k]; // convert to log-odds
									if (temp_profile[k] > temp_profile[patterns[j]]) {
										patterns2[j] = patterns[j];
										patterns[j] = k;
									}
								}
							}
						}
					}
					free(profile);
					
					// convert profile to k-mers
					s = 0; // unmasked k-mers in original sequence
					p = 0; // number of alternative pattern k-mers
					for (j = 0; j < l[i]; j++) {
						if (w[j] != NA_INTEGER) { // query k-mer is not masked
							s++;
							if (patterns[j] != NA_INTEGER) { // enriched letter
								for (k = 1; k < K; k++) {
									if (patterns[j + k] == NA_INTEGER) {
										patterns[j] = NA_INTEGER;
										break;
									}
									patterns[j] += pwv[k]*patterns[j + k];
									patterns2[j] += pwv[k]*patterns[j + k];
								}
								if (patterns[j] != NA_INTEGER) { // enriched letters
									if (patterns[j] != w[j]) { // most frequent pattern k-mer
										p++;
									} else if (patterns[j] != patterns2[j]) { // alternative pattern k-mer
										patterns[j] = patterns2[j];
										p++;
									}
								}
							}
						}
					}
					free(patterns2);
					
					if (p == 0) { // no alternative pattern k-mers
						free(patterns);
						break;
					}
					
					// count target occurrences of each query k-mer
					counts = (int *) malloc((s + p)*sizeof(int)); // thread-safe on Windows
					s = 0; // total number of k-mers shared with targets
					p = 0; // position in counts
					for (j = 0; j < l[i]; j++) { // each query k-mer
						if (w[j] != NA_INTEGER) {
							counts[p] = num[w[j]];
							s += counts[p++];
							if (s < 0) // signed integer overflow
								break;
							if (patterns[j] != NA_INTEGER && // enriched letters
								patterns[j] != w[j]) { // alternative pattern k-mer
								counts[p] = num[patterns[j]];
								s += counts[p++];
								if (s < 0) // signed integer overflow
									break;
							}
						}
					}
					
					if (s < 0) { // revert to previous results
						free(patterns);
						free(counts);
						break;
					}
					
					free(set);
					free(score);
					free(posQuery);
					free(posTarget);
					free(len);
					free(chain);
					free(res);
					
					// record target occurrences of each query k-mer
					posQuery = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
					posTarget = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
					set = (int *) malloc(s*sizeof(int)); // thread-safe on Windows
					s = 0; // target k-mer
					c = 0; // position in counts
					for (j = 0; j < l[i]; j++) { // each query k-mer
						if (w[j] != NA_INTEGER) {
							R_xlen_t P = offset[w[j]]; // position of k-mer in loc and ind
							for (k = 0; k < counts[c]; k++) { // each target instance of k-mer
								posQuery[s] = j + 1;
								posTarget[s] = loc[P];
								set[s] = ind[P];
								P++;
								s++;
							}
							c++;
							if (patterns[j] != NA_INTEGER && // enriched letters
								patterns[j] != w[j]) { // alternative pattern k-mer
								R_xlen_t P = offset[patterns[j]]; // position of k-mer in loc and ind
								for (k = 0; k < counts[c]; k++) { // each target instance of k-mer
									posQuery[s] = j + 1;
									posTarget[s] = loc[P];
									set[s] = ind[P];
									P++;
									s++;
								}
								c++;
							}
						}
					}
					free(patterns);
				}
			}
			if (sizeM > 0)
				free(subScores);
			
			// eliminate chains below the minimum score
			k = 0;
			if (ISNA(minS)) { // determine the minimum score per target
				double mS; // minimum score
				for (j = 0; j < c; j++) {
					mS = log((tot - (double)pos[set[res[j]] - 1])/(double)step);
					if (score[res[j]] >= mS)
						res[k++] = res[j];
				}
			} else { // apply minimum score threshold
				for (j = 0; j < c; j++) {
					if (score[res[j]] >= minS)
						res[k++] = res[j];
				}
			}
			c = k; // count passing minScore
			
			// select top scores within each set
			int *res2;
			if (limitT > 0) {
				res2 = (int *) malloc(c*sizeof(int)); // thread-safe on Windows
				k = 0;
				p = 0;
				for (j = 1; j < c; j++) {
					if (set[res[j]] != set[res[p]]) {
						s = j - p; // number to select
						if (s > limitT) {
							heapSelect(score, res + p, s, limitT);
							s = limitT;
						}
						p = 0;
						while (p < s) {
							p++;
							res2[k++] = res[j - p];
						}
						p = j;
					}
				}
				s = c - p; // number to select
				if (s > limitT) {
					heapSelect(score, res + p, s, limitT);
					s = limitT;
				}
				p = 0;
				while (p < s) {
					p++;
					res2[k++] = res[c - p];
				}
				c = k;
				free(res);
				res = res2;
			}
			
			// select top scores overall
			if (limitQ > 0 && c > limitQ) {
				res2 = (int *) malloc(c*sizeof(int)); // thread-safe on Windows
				heapSelect(score, res, c, limitQ);
				k = 0;
				p = 0;
				while (p < limitQ) {
					p++;
					res2[k++] = res[c - p];
				}
				c = limitQ;
				free(res);
				res = res2;
			}
			
			set2 = (int *) malloc(c*sizeof(int)); // thread-safe on Windows
			score2 = (double *) malloc(c*sizeof(double)); // thread-safe on Windows
			for (j = 0; j < c; j++) {
				set2[j] = set[res[j]];
				score2[j] = score[res[j]];
			}
			free(set);
			free(score);
			set = set2;
			score = score2;
			
			ptrs[i] = set;
			vecs[i] = score;
			l[i] = c; // number of results
			
			if (sO == 0) { // include anchor positions with output
				int **anchors = (int **) malloc(c*sizeof(int *)); // thread-safe on Windows
				for (j = 0; j < c; j++) {
					// measure length of chain
					p = res[j];
					k = 1; // length of chain
					while (p != chain[p]) {
						p = chain[p]; // traceback
						k++; // lengthen chain
					}
					
					// record anchor positions
					int *anchor = (int *) malloc((4*k + 1)*sizeof(int)); // thread-safe on Windows
					anchor[0] = k; // store number of anchors in first position
					p = res[j];
					k = 1;
					anchor[k++] = posQuery[p];
					anchor[k++] = posQuery[p] + len[p] - 1;
					anchor[k++] = posTarget[p];
					anchor[k++] = posTarget[p] + len[p] - 1;
					while (p != chain[p]) {
						p = chain[p]; // traceback
						anchor[k++] = posQuery[p];
						anchor[k++] = posQuery[p] + len[p] - 1;
						anchor[k++] = posTarget[p];
						anchor[k++] = posTarget[p] + len[p] - 1;
					}
					anchors[j] = anchor;
				}
				matrices[i] = anchors;
			}
			
			free(posQuery);
			free(posTarget);
			free(len);
			free(chain);
			free(res);
			
			if (v) {
				#ifdef _OPENMP
				#pragma omp critical
				{
					N++;
				}
				#else
				N++;
				#endif
			}
			
			#ifdef _OPENMP
			int master = omp_get_thread_num();
			#else
			int master = 0;
			#endif
			
			if (master == 0) { // master thread
				time(&end);
				elapsed = difftime(end, start);
				
				if (elapsed >= 1) { // at least 1 second has elapsed
					start = end;
					abort = checkInterrupt(); // call back to R
					if (abort == 0 && v) {
						soFar = N/n;
						*rPercentComplete = floor(100*soFar);
						if (*rPercentComplete > before) { // when the percent has changed
							// tell the progress bar to update in the R console
							eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
							before = *rPercentComplete;
						}
					}
				}
			}
		} else {
			int *set = (int *) malloc(0*sizeof(int)); // thread-safe on Windows
			double *score = (double *) malloc(0*sizeof(double)); // thread-safe on Windows
			ptrs[i] = set;
			vecs[i] = score;
			l[i] = 0; // number of results
			if (sO == 0) {
				int **anchors = (int **) malloc(0*sizeof(int *)); // thread-safe on Windows
				matrices[i] = anchors;
			}
		}
	}
	
	free(offset);
	free(sepCost);
	free(gapCost);
	if (sizeM > 0) {
		free(lkup_row);
		free(lkup_col);
		if (maxIt > 0)
			free(pwv);
	} else {
		free(addScores);
		free(scores);
	}
	
	int **anchors; // pointers to anchor positions
	if (abort != 0) { // release memory
		for (i = 0; i < n; i++) {
			int *set = ptrs[i];
			double *score = vecs[i];
			if (sO == 0) {
				anchors = matrices[i];
				for (j = 0; j < l[i]; j++) {
					int *anchor = anchors[j];
					free(anchor);
				}
				free(anchors);
			}
			free(set);
			free(score);
		}
		
		free(vecs);
		free(ptrs);
		free(l);
		if (sO == 0)
			free(matrices);
		
		if (abort < 0) {
			error("Received user interrupt.");
		} else {
			error("Too many target k-mer hits for myXStringSet[%d].", abort);
		}
	}
	
	c = 0;
	for (i = 0; i < n; i++)
		c += l[i];
	
	SEXP ans, ans0, ans1, ans2, ans3, ret_list;
	PROTECT(ans0 = allocVector(INTSXP, c));
	int *rans0 = INTEGER(ans0);
	PROTECT(ans1 = allocVector(INTSXP, c));
	int *rans1 = INTEGER(ans1);
	PROTECT(ans2 = allocVector(REALSXP, c));
	double *rans2 = REAL(ans2);
	if (sO == 0)
		PROTECT(ans3 = allocVector(VECSXP, c));
	
	k = 0; // position in output vectors
	for (i = 0; i < n; i++) {
		int *set = ptrs[i];
		double *score = vecs[i];
		if (sO == 0)
			anchors = matrices[i];
		for (j = 0; j < l[i]; j++) {
			rans0[k] = i + 1;
			rans1[k] = set[j];
			rans2[k] = score[j];
			if (sO == 0) {
				int *anchor = anchors[j];
				c = anchor[0]; // number of anchors
				PROTECT(ans = allocMatrix(INTSXP, 4, c));
				int *rans = INTEGER(ans);
				p = 0;
				int e1 = 0, e2 = 0; // ends
				while (c > 0) { // reverse anchor order
					// pull-back anchor overlap due to indels
					int d1 = anchor[4*c - 3] - e1;
					int d2 = anchor[4*c - 1] - e2;
					if (d2 < d1)
						d1 = d2;
					if (d1 <= 0) {
						d1 = 1 - d1;
					} else {
						d1 = 0;
					}
					rans[p++] = anchor[4*c - 3] + d1;
					e1 = anchor[4*c - 2];
					rans[p++] = e1;
					rans[p++] = anchor[4*c - 1] + d1;
					e2 = anchor[4*c];
					rans[p++] = e2;
					c--;
				}
				free(anchor);
				SET_VECTOR_ELT(ans3, k, ans);
				UNPROTECT(1);
			}
			k++;
		}
		free(set);
		free(score);
		if (sO == 0)
			free(anchors);
	}
	
	free(vecs);
	free(ptrs);
	free(l);
	if (sO == 0) {
		PROTECT(ret_list = allocVector(VECSXP, 4));
		free(matrices);
	} else { // score only
		PROTECT(ret_list = allocVector(VECSXP, 3));
	}
	
	SET_VECTOR_ELT(ret_list, 0, ans0);
	SET_VECTOR_ELT(ret_list, 1, ans1);
	SET_VECTOR_ELT(ret_list, 2, ans2);
	
	if (sO == 0) {
		SET_VECTOR_ELT(ret_list, 3, ans3);
		UNPROTECT(5);
	} else { // score only
		UNPROTECT(4);
	}
	if (v)
		UNPROTECT(2);
	
	return ret_list;
}

// update count of k-mers in a list of k-mers
SEXP countIndex(SEXP num, SEXP query, SEXP step)
{
	if (MAYBE_SHARED(num))
		error(".Call function 'countIndex' called in incorrect context.");
	
	int i, j, n, *kmers;
	int *rans = INTEGER(num);
	int l = length(query); // number of sequences
	int s = asInteger(step); // step size between k-mers
	
	// set up a timer to minimize interrupt checks
	time_t start, end;
	double elapsed;
	time(&start);
	
	for (i = 0; i < l; i++) {
		kmers = INTEGER(VECTOR_ELT(query, i)); // query k-mers
		n = length(VECTOR_ELT(query, i)); // query length
		for (j = 0; j < n; j += s)
			if (kmers[j] != NA_INTEGER) // unmasked position
				rans[kmers[j]]++; // increment count
		
		time(&end);
		elapsed = difftime(end, start);
		if (elapsed >= 1) { // at least 1 second has elapsed
			start = end;
			int abort = checkInterrupt(); // call back to R
			
			if (abort != 0)
				error("Received user interrupt.");
		}
	}
	
	return R_NilValue;
}

// update record of k-mers in a list of k-mers
SEXP updateIndex(SEXP offset, SEXP query, SEXP wordSize, SEXP step, SEXP location, SEXP index, SEXP positions, SEXP count)
{
	if (MAYBE_SHARED(offset))
		error(".Call function 'updateIndex' called in incorrect context.");
	if (MAYBE_SHARED(location))
		error(".Call function 'updateIndex' called in incorrect context.");
	if (MAYBE_SHARED(index))
		error(".Call function 'updateIndex' called in incorrect context.");
	if (MAYBE_SHARED(positions))
		error(".Call function 'updateIndex' called in incorrect context.");
	
	int i, j, k, n, *kmers;
	double *off = REAL(offset);
	int l = length(query); // number of sequences
	int K = asInteger(wordSize); // k-mer length
	int s = asInteger(step); // step size between k-mers
	int *loc = INTEGER(location);
	int *ind = INTEGER(index);
	int *pos = INTEGER(positions);
	int c = asInteger(count); // processed count
	pos += c; // offset starting position by count
	
	// set up a timer to minimize interrupt checks
	time_t start, end;
	double elapsed;
	time(&start);
	
	for (i = 0; i < l; i++) {
		c++;
		kmers = INTEGER(VECTOR_ELT(query, i)); // query k-mers
		n = length(VECTOR_ELT(query, i)); // query length
		
		k = -1*s - 1; // last hit
		for (j = 0; j < n; j += s) {
			if (kmers[j] != NA_INTEGER) {
				// record the number of unmasked positions
				if (k == j - s) {
					pos[i] += s;
				} else {
					pos[i] += K;
				}
				k = j;
				
				// update the inverted index
				ind[(R_xlen_t)off[kmers[j]]] = c;
				loc[(R_xlen_t)off[kmers[j]]] = j + 1;
				off[kmers[j]]++;
			}
		}
		
		time(&end);
		elapsed = difftime(end, start);
		if (elapsed >= 1) { // at least 1 second has elapsed
			start = end;
			int abort = checkInterrupt(); // call back to R
			
			if (abort != 0)
				error("Received user interrupt.");
		}
	}
	
	return R_NilValue;
}

// determine offset and approximate letter frequencies
SEXP approxFreqs(SEXP offset, SEXP freqs, SEXP count)
{
	if (MAYBE_SHARED(offset))
		error(".Call function 'approxFreqs' called in incorrect context.");
	if (MAYBE_SHARED(freqs))
		error(".Call function 'approxFreqs' called in incorrect context.");
	
	int l = length(offset);
	double *off = REAL(offset);
	double *freq = REAL(freqs);
	int s = length(freqs);
	int *num = INTEGER(count);
	int binSize = l/s;
	int n = binSize;
	
	int j;
	int i = 0; // index in num
	int k = 0; // index in freqs
	freq[0] = num[0];
	while (i < l - 1) {
		j = i + 1;
		off[j] = off[i] + (double)num[i];
		i = j;
		if (i >= n) {
			k++;
			n += binSize;
		}
		freq[k] += num[i];
	}
	
	return R_NilValue;
}

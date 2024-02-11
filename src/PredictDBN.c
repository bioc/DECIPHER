/****************************************************************************
 *                Predicts 3 State Secondary RNA Structure                  *
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

static int firstpos(const Chars_holder *P)
{
	int i;
	const char *p;
	
	// start from the beginning of the sequence
	for (i = 0, p = P->ptr;
		 i < P->length;
		 i++, p++)
	{
		if (!((*p) & 0x10 || (*p) & 0x40)) {
			return i;
		}
	}
	return i;
}

static int lastpos(const Chars_holder *P)
{
	int i;
	const char *p;
	
	// start from the end of the sequence
	for (i = (P->length - 1), p = (P->ptr + P->length - 1);
		 i >= 0;
		 i--, p--)
	{
		if (!((*p) & 0x10 || (*p) & 0x40)) {
			return i;
		}
	}
	return i;
}

void Traceback(double *MI, int tot, int *unpaired, int *pos, char *states, char leftSymbol, char rightSymbol, int i, int j) {
	while (j > (i + 1)) { // prevent < 2 base hairpins
//		Rprintf("\ni = %d j = %d MI[%d, %d] = %1.0f", i + 1, j + 1, unpaired[i] + 1, unpaired[j] + 1, MI[unpaired[j]*tot + unpaired[i]]);
		if (MI[unpaired[j]*tot + unpaired[i]] > 1e9) { // bifurcation
			Traceback(MI, tot, unpaired, pos, states, leftSymbol, rightSymbol, MI[unpaired[j]*tot + unpaired[i]] - 1e9 + 1, j);
			Traceback(MI, tot, unpaired, pos, states, leftSymbol, rightSymbol, i, MI[unpaired[j]*tot + unpaired[i]] - 1e9);
			break;
		} else if (MI[unpaired[j]*tot + unpaired[i]] < 0 && MI[unpaired[j]*tot + unpaired[i]] > -1e9) {
			i -= MI[unpaired[j]*tot + unpaired[i]];
		} else if (MI[unpaired[j]*tot + unpaired[i]] < -1e9) {
			j += MI[unpaired[j]*tot + unpaired[i]] + 1e9;
		} else { // base pairing
//			if (states[pos[unpaired[i]]] != '.' || states[pos[unpaired[j]]] != '.')
//				error("crossed-over twice");
			states[pos[unpaired[i]]] = leftSymbol;
			states[pos[unpaired[j]]] = rightSymbol;
			i++;
			j--;
		}
	}
}

double Choose(double N, double K) {
	double result = 1;
	for (int i = 1; i <= K; i++) {
		result *= N - K + i;
		result /= i;
	}
	return result;
}

double pNoRun(double N, double K, double p) {
	double prob1 = 0, prob2 = 0;
	double i = 0;
	double top, c;
	while (1) {
		top = N - (i + 1)*K;
		if (top < i)
			break;
		c = Choose(top, i);
		prob1 += c*pow(-1*(1 - p)*pow(p, K), i);
		i++;
	}
	prob1 = prob1*pow(p, K);
	
	i = 1;
	while (1) {
		top = N - i*K;
		if (top < i)
			break;
		c = Choose(top, i);
		prob2 += c*pow(-1*(1 - p)*pow(p, K), i);
		i++;
	}
	
	prob1 = 1 - prob1 + prob2;
	if (isnan(prob1) || prob1 > 1 || prob1 < 0) {
		return 0;
	} else {
		return prob1;
	}
}

SEXP predictDBN(SEXP x, SEXP output, SEXP minOccupancy, SEXP impact, SEXP avgProdCorr, SEXP slope, SEXP shift, SEXP weights, SEXP pseudoknots, SEXP threshold, SEXP patterns, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	int i, j, k, p, s, d, l;
	
	// type of output
	// 1 = DBN
	// 2 = pairs
	// 3 = scores
	// 4 = structures
	// 5 = search
	// 6 = evidence
	int o = asInteger(output);
	double minO = asReal(minOccupancy);
	double *coef = REAL(impact);
	double APC = asReal(avgProdCorr);
	double sl = asReal(slope);
	double sh = asReal(shift);
	double *w = REAL(weights);
	int pseudo = asInteger(pseudoknots);
	double thresh = asReal(threshold);
	int before, v, last, *rPercentComplete;
	double *rans, soFar;
	int nthreads = asInteger(nThreads);
	SEXP ans, ans_s, percentComplete, utilsPackage;
	v = asLogical(verbose);
	if (v) { // percent complete variables
		soFar = 0;
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	XStringSet_holder x_set;
	x_set = hold_XStringSet(x);
	int x_length = get_length_from_XStringSet_holder(&x_set);
	Chars_holder x_s;
	x_s = get_elt_from_XStringSet_holder(&x_set, 0);
	int width = x_s.length;
	
	// initialize an array of terminal gap lengths
	int *endpoints = Calloc(2*x_length, int); // initialized to zero
	
	for (i = 0; i < x_length; i++) {
		x_s = get_elt_from_XStringSet_holder(&x_set, i);
		endpoints[i] = firstpos(&x_s);
		endpoints[x_length + i] = lastpos(&x_s);
	}
	
	// initialize an array of letter (A, C, G, T/U, other) frequencies
	double *counts = Calloc(5*width, double); // initialized to zero
	
	for (s = 0; s < x_length; s++) {
		x_s = get_elt_from_XStringSet_holder(&x_set, s);
		// assume x_s.length is equal to width
		for (p = endpoints[s]; p <= endpoints[x_length + s]; p++) {
			switch (x_s.ptr[p]) {
				case 1: // A
					counts[5*p] += w[s];
					break;
				case 2: // C
					counts[5*p + 1] += w[s];
					break;
				case 4: // G
					counts[5*p + 2] += w[s];
					break;
				case 8: // T/U
					counts[5*p + 3] += w[s];
					break;
				default: // other
					counts[5*p + 4] += w[s];
					break;
			}
		}
	}
	
	// initialize an array of positions
	int *pos = Calloc(width, int); // initialized to zero
	double sum; // sum of weights
	int tot = 0; // total number of positions >= minO
	for (p = 0; p < width; p++) {
		sum = counts[5*p] + counts[5*p + 1] + counts[5*p + 2] + counts[5*p + 3];
		if ((sum/x_length) >= minO) { // at least minOccupancy
			// normalize to sum of weights
			sum += counts[5*p + 4];
			counts[5*p] /= sum;
			counts[5*p + 1] /= sum;
			counts[5*p + 2] /= sum;
			counts[5*p + 3] /= sum;
			counts[5*p + 4] /= sum;
			
			pos[tot] = p;
			tot++;
			//Rprintf("\np = %d A=%1.2f C=%1.2f G=%1.2f U=%1.2f other=%1.2f", p + 1, counts[5*p], counts[5*p + 1], counts[5*p + 2], counts[5*p + 3], counts[5*p + 4]);
		}
	}
	
	// initialize an array of mutual information
	double *MI = Calloc(tot*tot, double); // initialized to zero
	double *rowMeans = Calloc(tot, double); // initialized to zero
	
	last = tot - 1;
	for (i = 0; i < (tot - 1); i++) {
		#ifdef _OPENMP
		#pragma omp parallel for private(j,l,s,x_s) schedule(guided) num_threads(nthreads)
		#endif
		for (j = i + 1; j < tot; j++) {
			double AU = 0, UA = 0, GC = 0, CG = 0, GU = 0, UG = 0, other = 0, temp = 0, bg;
			l = 0;
			
			for (s = 0; s < x_length; s++) {
				if (pos[i] < endpoints[s] || pos[j] > endpoints[x_length + s])
					continue;
				
				x_s = get_elt_from_XStringSet_holder(&x_set, s);
				
				int p1, p2;
				switch (x_s.ptr[pos[i]]) {
					case 1: // A
						p1 = 1;
						break;
					case 2: // C
						p1 = 2;
						break;
					case 4: // G
						p1 = 3;
						break;
					case 8: // T/U
						p1 = 4;
						break;
					default: // other
						p1 = 5;
						break;
				}
				switch (x_s.ptr[pos[j]]) {
					case 1: // A
						p2 = 1;
						break;
					case 2: // C
						p2 = 2;
						break;
					case 4: // G
						p2 = 3;
						break;
					case 8: // T/U
						p2 = 4;
						break;
					default: // other
						p2 = 5;
						break;
				}
				
				if (p1 == 5) { // -
					other += w[s];
				} else if (p2 == 5) { // -
					other += w[s];
				} else if (p1 == 1) { // A
					if (p2 == 4) { // T/U
						AU += w[s];
					} else {
						other += w[s];
					}
				} else if (p1 == 2) { // C
					if (p2 == 3) { // G
						CG += w[s];
					} else {
						other += w[s];
					}
				} else if (p1 == 3) { // G
					if (p2 == 2) { // C
						GC += w[s];
					} else if (p2 == 4) { // T/U
						GU += w[s];
					} else {
						other += w[s];
					}
				} else if (p1 == 4) { // T/U
					if (p2 == 1) { // A
						UA += w[s];
					} else if (p2 == 3) { // G
						UG += w[s];
					} else {
						other += w[s];
					}
				}
			}
			
			// normalize to x_length
			AU /= x_length;
			UA /= x_length;
			GC /= x_length;
			CG /= x_length;
			GU /= x_length;
			UG /= x_length;
			other /= x_length; // other does not include terminal gaps
			
			bg = counts[5*pos[i]]*counts[5*pos[j] + 3]; // AU
			if (bg > 0 && AU > 0)
				temp += AU*log2(AU/bg)*coef[0];
			bg = counts[5*pos[i] + 3]*counts[5*pos[j]]; // UA
			if (bg > 0 && UA > 0)
				temp += UA*log2(UA/bg)*coef[0];
			bg = counts[5*pos[i] + 2]*counts[5*pos[j] + 1]; // GC
			if (bg > 0 && GC > 0)
				temp += GC*log2(GC/bg)*coef[1];
			bg = counts[5*pos[i] + 1]*counts[5*pos[j] + 2]; // CG
			if (bg > 0 && CG > 0)
				temp += CG*log2(CG/bg)*coef[1];
			bg = counts[5*pos[i] + 2]*counts[5*pos[j] + 3]; // GU
			if (bg > 0 && GU > 0)
				temp += GU*log2(GU/bg)*coef[2];
			bg = counts[5*pos[i] + 3]*counts[5*pos[j] + 2]; // UG
			if (bg > 0 && UG > 0)
				temp += UG*log2(UG/bg)*coef[2];
			
			temp += other*coef[3];
			MI[i*tot + j] = temp;
			//Rprintf("\ni = %d j = %d MI = %1.2f", pos[i] + 1, pos[j] + 1, MI[i*tot + j]);
			
			rowMeans[j] += temp;
		}
		for (j = i + 1; j < tot; j++)
			rowMeans[i] += MI[i*tot + j];
		
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
	Free(endpoints);
	Free(counts);
	
	// apply Average Product Correction (APC)
	double avg;
	if (APC != 0) {
		avg = 0;
		for (i = 0; i < tot; i++) {
			rowMeans[i] /= tot - 1;
			avg += rowMeans[i];
		}
		avg /= tot;
		
		for (i = 0; i < (tot - 1); i++)
			for (j = i + 1; j < tot; j++)
				MI[i*tot + j] -= APC*rowMeans[i]*rowMeans[j]/avg;
	}
	Free(rowMeans);
	
//	for (i = 0; i < tot; i++) {
//		Rprintf("\n");
//		for (j = 0; j < tot; j++)
//			Rprintf("%1.2f ", MI[j*tot + i]);
//	}
	
	// determine the nth highest value
	int n = (int)(tot/2); // max number of MI values that could be paired
	double *vals = Calloc(n, double); // initialized to zero
	for (i = 0; i < n; i++)
		vals[i] = -1e12;
	double minVal = -1e12; // minimum value in vals
	int index = 0; // index of minimum in vals
	// determine the nth largest value in MI
	for (i = 0; i < (tot - 1); i++) {
		for (j = i + 1; j < tot; j++) {
			if (MI[i*tot + j] > minVal) {
				// replace minVal
				vals[index] = MI[i*tot + j];
				minVal = vals[index];
				
				// find the new minVal
				for (p = 0; p < n; p++) {
					if (vals[p] < minVal) {
						minVal = vals[p];
						index = p;
					}
				}
			}
		}
	}
	Free(vals);
	
	// apply sigmoidal transformation
	sh *= minVal; // shift
	sl /= -1*(sh - minVal); // -slope
	n = 0; // number of positions above threshold
	for (i = 0; i < (tot - 1); i++) {
		for (j = i + 1; j < tot; j++) {
			MI[i*tot + j] = 1/(1 + exp(sl*log(MI[i*tot + j]/sh)));
			if (ISNAN(MI[i*tot + j])) {
				MI[i*tot + j] = 0;
			} else if (o < 3 && MI[i*tot + j] < thresh) {
				MI[i*tot + j] = 0;
			}
			if (o == 6 &&  MI[i*tot + j] >= thresh)
				n++;
			//Rprintf("\ni = %d j = %d MI = %1.2f", pos[i] + 1, pos[j] + 1, MI[i*tot + j]);
		}
	}
	//Rprintf("\nshift = %1.2f slope = %1.2f minVal = %1.2f", sh, -1*sl, minVal);
	
	if (patterns != R_NilValue) {
		SEXP P1 = VECTOR_ELT(patterns, 0);
		SEXP P2 = VECTOR_ELT(patterns, 1);
		for (i = 0; i < (tot - 1); i++) {
			SEXP P11 = VECTOR_ELT(P1, pos[i]);
			SEXP P21 = VECTOR_ELT(P2, pos[i]);
			int L = length(P11);
			if (L == 0)
				continue;
			int *J = INTEGER(P11);
			double *W = REAL(P21);
			int count = 0;
			j = 0;
			while (count < L && j < tot) {
				if (W[count] < thresh) {
					count++;
				} else {
					while (j < tot) {
						if (J[count] - 1 > pos[j]) {
							j++;
						} else if (J[count] - 1 == pos[j]) {
							if (MI[i*tot + j] < W[count])
								MI[i*tot + j] = W[count];
							count++;
							j++;
							break;
						} else {
							count++;
							break;
						}
					}
				}
			}
		}
	}
	
	char *states;
	int *unpaired, q, *leftMax, *rightMax;
	double *MI2, *rowMax, *colMax;
	double match, left, right, prevL, prevR;
	char leftSymbol, rightSymbol;
	
	if (o < 3) { // perform traceback
		l = tot; // number remaining unpaired
		q = 10; // block size
		states = Calloc(width + 1, char);
		for (i = 0; i < width; i++)
			states[i] = '-';
		for (i = 0; i < tot; i++)
			states[pos[i]] = '.';
		states[width] = '\0'; // end (null terminate) the string
		unpaired = Calloc(tot, int); // initialized to zero
		for (i = 1; i < tot; i++)
			unpaired[i] = i;
		for (p = 0; p <= pseudo; p++) {
			if (p < pseudo) { // copy MI
				MI2 = Calloc(tot*tot, double); // initialized to zero
				for (i = 0; i < (tot - 1); i++)
					for (j = i + 1; j < tot; j++)
						MI2[i*tot + j] = MI[i*tot + j];
			}
			
			n = ceil((double)l/(double)q);
			rowMax = Calloc(l*n, double); // initialized to zero
			colMax = Calloc(l*n, double); // initialized to zero
			for (d = 2; d <= l; d++) {
				#ifdef _OPENMP
				#pragma omp parallel for private(i,j,k,match,left,right,prevL,prevR) schedule(guided) num_threads(nthreads)
				#endif
				for (i = 0; i < (l - d + 1); i++) {
					j = i + d - 1; // i <= j
					
					if (d > 2) { // not along first diagonals
						match = MI[unpaired[i]*tot + unpaired[j]];
						match += MI[(unpaired[i + 1])*tot + unpaired[j - 1]];
					} else {
						match = 0;
					}
					
					left = MI[(unpaired[i + 1])*tot + unpaired[j]];
					if (i < l) {
						prevL = MI[unpaired[j]*tot + unpaired[i + 1]];
						if (prevL < 0 && prevL > -1e9) {
							prevL -= 1;
						} else {
							prevL = -1;
						}
					} else {
						prevL = -1;
					}
					
					right = MI[unpaired[i]*tot + unpaired[j - 1]];
					if (j > 1) {
						prevR = MI[(unpaired[j - 1])*tot + unpaired[i]];
						if (prevR < -1e9) {
							prevR -= 1;
						} else {
							prevR = -1e9 - 1;
						}
					} else {
						prevR = -1e9 - 1;
					}
					
					if (match > left && match > right) {
						MI[unpaired[i]*tot + unpaired[j]] = match;
					} else if (left > right) {
						MI[unpaired[i]*tot + unpaired[j]] = left;
						MI[unpaired[j]*tot + unpaired[i]] = prevL;
					} else {
						MI[unpaired[i]*tot + unpaired[j]] = right;
						MI[unpaired[j]*tot + unpaired[i]] = prevR;
					}
					
					// bifurcation
					k = i + 3;
					while (k <= (j - 4)) {
						if ((k % q) == 0) {
							if ((rowMax[i*n + k/q] + colMax[j*n + k/q]) <= MI[unpaired[i]*tot + unpaired[j]]) {
								k += q;
								continue;
							}
						}
						
						if ((MI[unpaired[i]*tot + unpaired[k]] + MI[(unpaired[k + 1])*tot + unpaired[j]]) > MI[unpaired[i]*tot + unpaired[j]]) {
							MI[unpaired[i]*tot + unpaired[j]] = MI[unpaired[i]*tot + unpaired[k]] + MI[(unpaired[k + 1])*tot + unpaired[j]];
							MI[unpaired[j]*tot + unpaired[i]] = k + 1e9;
						}
						k++;
					}
					if (MI[unpaired[i]*tot + unpaired[j]] > rowMax[i*n + j/q])
						rowMax[i*n + j/q] = MI[unpaired[i]*tot + unpaired[j]];
					if (i > 0 && MI[unpaired[i]*tot + unpaired[j]] > colMax[j*n + (i - 1)/q])
						colMax[j*n + (i - 1)/q] = MI[unpaired[i]*tot + unpaired[j]];
				}
			}
			Free(rowMax);
			Free(colMax);
			
			if (p == 0) {
				leftSymbol = '(';
				rightSymbol = ')';
			} else if (p == 1) {
				leftSymbol = '[';
				rightSymbol = ']';
			} else if (p == 2) {
				leftSymbol = '{';
				rightSymbol = '}';
			} else if (p == 3) {
				leftSymbol = '<';
				rightSymbol = '>';
			}
			
			Traceback(MI, tot, unpaired, pos, states, leftSymbol, rightSymbol, 0, l - 1);
			
			if (p < pseudo) { // replace MI with the original
				Free(MI);
				MI = MI2;
				
				// excluded paired positions from unpaired
				l = 0;
				for (i = 0; i < tot; i++) {
					if (states[pos[i]] == '.')
						unpaired[l++] = i;
				}
				if (l == 0) // all positions paired
					break;
			}
		}
		Free(unpaired);
		
		PROTECT(ans = allocVector(STRSXP, 1));
		SET_STRING_ELT(ans, 0, mkChar(states));
		Free(states);
	} else if (o == 3) { // scores
		PROTECT(ans = allocMatrix(REALSXP, 3, width)); // [state][pos]
		rans = REAL(ans);
		for (i = 0; i < 3*width; i++)
			rans[i] = 0;
		
		// determine the largest value in each row/column
		for (i = 0; i < (tot - 1); i++) {
			for (j = i + 1; j < tot; j++) {
				if (MI[i*tot + j] > rans[3*pos[i] + 1])
					rans[3*pos[i] + 1] = MI[i*tot + j];
				if (MI[i*tot + j] > rans[3*pos[j] + 2])
					rans[3*pos[j] + 2] = MI[i*tot + j];
			}
		}
		
		// normalize the scores
		for (i = 0; i < tot; i++) {
			sum = rans[3*pos[i] + 1] + rans[3*pos[i] + 2];
			if (sum > 1) {
				rans[3*pos[i] + 1] /= sum;
				rans[3*pos[i] + 2] /= sum;
			} else {
				rans[3*pos[i]] = 1 - sum;
			}
		}
	} else if (o == 6) { // evidence
		PROTECT(ans = allocMatrix(REALSXP, n, 3));
		rans = REAL(ans);
		
		int count = 0;
		for (i = 0; i < (tot - 1); i++) {
			for (j = i + 1; j < tot; j++) {
				if (MI[i*tot + j] >= thresh) {
					rans[count] = pos[i] + 1; // left pair
					rans[count + n] = pos[j] + 1; // right pair
					rans[count + 2*n] = MI[i*tot + j]; // score
					count++;
					if (count == n)
						break;
				}
			}
			if (count == n)
				break;
		}
	} else { // structures or search
		// initialize constants for output type "structures"
		int minLoop = 3;
		double pMatch = 0.3535533905932737863687; // sqrt(32/256)
		double NN[256] = { // absolute value of nearest neighbor base pairing free energies
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.93, // AA
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.11, 0.00, // AC
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.35, 0.00, 1.27, // AG
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.33, 0.00, 1.00, 0.00, // AU
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.24, 0.00, 0.00, 0.00, 0.00, // CA
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 3.26, 0.00, 0.00, 0.00, 0.00, 0.00, // CC
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 3.42, 0.00, 2.51, 0.00, 0.00, 0.00, 0.00, // CG
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.35, 0.00, 1.53, 0.00, 0.00, 0.00, 0.00, 0.00, // CU
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.08, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.55, // GA
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.36, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.41, 0.00, // GC
			0.00, 0.00, 0.00, 0.00, 0.00, 3.26, 0.00, 2.11, 0.00, 0.00, 0.00, 0.00, 0.00, 1.53, 0.00, 0.00, // GG
			0.00, 0.00, 0.00, 0.00, 2.11, 0.00, 1.41, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, // GU
			0.00, 0.00, 0.00, 1.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.36, 0.00, 0.00, 0.00, 0.00, // UA
			0.00, 0.00, 2.08, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 2.11, 0.00, 0.00, 0.00, 0.00, 0.00, // UC
			0.00, 2.24, 0.00, 1.36, 0.00, 0.00, 0.00, 0.00, 0.00, 2.51, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, // UG
			0.93, 0.00, 0.55, 0.00, 0.00, 0.00, 0.00, 0.00, 1.27, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 // UU
		};
		
		leftMax = Calloc(tot, int); // initialized to zero
		rightMax = Calloc(tot, int); // initialized to zero
		
		// locate the largest value in each row/column
		for (i = 0; i < (tot - 1); i++) {
			for (j = i + minLoop + 1; j < tot; j++) {
				if (MI[i*tot + j] > MI[i*tot + leftMax[i]]) // same (; diff )
					leftMax[i] = j; // index = (; value = )
				if (MI[i*tot + j] > MI[rightMax[j]*tot + j]) // same ); diff (
					rightMax[j] = i; // index = ); value = (
			}
		}
		
		PROTECT(ans = allocVector(VECSXP, x_length));
		
		for (s = 0; s < x_length; s++) {
			//Rprintf("\n\nSEQUENCE = %d", s);
			x_s = get_elt_from_XStringSet_holder(&x_set, s);
			
			n = 0;
			int *nucs = Calloc(width, int); // initialized to zero
			int *nogaps = Calloc(width, int); // initialized to zero
			
			for (i = 0; i < width; i++) {
				if (x_s.ptr[i] != 16 && x_s.ptr[i] != 32 && x_s.ptr[i] != 64) {
					if (x_s.ptr[i] & 0x2) { // C
						nogaps[n] = 1;
					} else if (x_s.ptr[i] & 0x4) { // G
						nogaps[n] = 2;
					} else if (x_s.ptr[i] & 0x8) { // U
						nogaps[n] = 3;
					} // else leave as A
					nucs[i] = n++;
				} else {
					nucs[i] = -1; // no nucleotide at this alignment site
				}
			}
			
			PROTECT(ans_s = allocMatrix(REALSXP, 3, n)); // [state][pos]
			
			rans = REAL(ans_s);
			for (i = 0; i < 3*n; i++)
				rans[i] = 0;
			
			// record anchor positions that agree with the consensus structure
			int *anchor = Calloc(tot, int); // initialized to zero
			for (i = 0; i < tot; i++) {
				if (((x_s.ptr[pos[i]] & 0x1) && (x_s.ptr[pos[leftMax[i]]] & 0x8)) || // A/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[leftMax[i]]] & 0x1)) || // U/A
					((x_s.ptr[pos[i]] & 0x2) && (x_s.ptr[pos[leftMax[i]]] & 0x4)) || // C/G
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[leftMax[i]]] & 0x2)) || // G/C
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[leftMax[i]]] & 0x8)) || // G/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[leftMax[i]]] & 0x4))) { // U/G
					if (rans[3*nucs[pos[i]] + 1] < MI[i*tot + leftMax[i]])
						rans[3*nucs[pos[i]] + 1] = MI[i*tot + leftMax[i]];
					if (MI[i*tot + leftMax[i]] >= thresh)
						anchor[i] = leftMax[i]; // left anchor
				} else if (MI[i*tot + leftMax[i]] >= thresh) {
					anchor[i] = -1*leftMax[i]; // missing left anchor
				}
				if (((x_s.ptr[pos[i]] & 0x1) && (x_s.ptr[pos[rightMax[i]]] & 0x8)) || // A/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[rightMax[i]]] & 0x1)) || // U/A
					((x_s.ptr[pos[i]] & 0x2) && (x_s.ptr[pos[rightMax[i]]] & 0x4)) || // C/G
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[rightMax[i]]] & 0x2)) || // G/C
					((x_s.ptr[pos[i]] & 0x4) && (x_s.ptr[pos[rightMax[i]]] & 0x8)) || // G/U
					((x_s.ptr[pos[i]] & 0x8) && (x_s.ptr[pos[rightMax[i]]] & 0x4))) { // U/G
					if (rans[3*nucs[pos[i]] + 2] < MI[rightMax[i]*tot + i])
						rans[3*nucs[pos[i]] + 2] = MI[rightMax[i]*tot + i];
					if (MI[rightMax[i]*tot + i] >= thresh)
						anchor[i] = rightMax[i]; // right anchor
				} else if (MI[rightMax[i]*tot + i] >= thresh) {
					anchor[i] = -1*rightMax[i]; // missing right anchor
				}
			}
			
			if (o == 5) { // search
				// eliminate single anchors between two unanchored positions
				int last1[2] = {-1};
				int last2[2] = {-1};
				for (i = 0; i < tot; i++) {
					if (anchor[i] > 0) { // anchored
						// shift previous
						last2[0] = last1[0];
						last2[1] = last1[1];
						last1[0] = i;
						last1[1] = anchor[i];
					} else if (anchor[i] < 0) { // missing anchor
						if (last1[1] > 0 && last2[1] < 0) { // middle anchor
							// zero the score
							if (last1[0] < last1[1]) { // left
								rans[3*nucs[pos[last1[0]]] + 1] = 0;
							} else { // right
								rans[3*nucs[pos[last1[0]]] + 2] = 0;
							}
							last1[1] *= -1; // unanchor
							anchor[last1[0]] = last1[1];
						}
						// shift previous
						last2[0] = last1[0];
						last2[1] = last1[1];
						last1[0] = i;
						last1[1] = anchor[i];
					} // else unanchored
				}
				
				// flag regions between anchors to trigger search
				last1[0] = -1; // last left anchor
				last1[1] = -1; // last right anchor
				for (i = 0; i < tot; i++) {
					if (anchor[i] > 0) {
						if (i < anchor[i]) { // left anchored
							if (last1[0] < last1[1]) { // right anchored previously
								if (i - 2 > last1[1]) // space between anchors
									anchor[last1[1] + 1] = -1*(i - 1); // mark as missing anchor
								last1[1] = -1; // reset last right anchor
							}
							last1[0] = i;
						} else { // right anchored
							last1[1] = i;
						}
					} else if (anchor[i] == 0) { // unanchored
						if (last1[0] >= 0) {
							anchor[i] = -1*(anchor[last1[0]] - 1); // mark as missing anchor
							last1[0] = -1; // reset last left anchor
							last1[1] = -1; // reset last right anchor
						}
					} else { // already a missing anchor
						last1[0] = -1; // reset last left anchor
						last1[1] = -1; // reset last right anchor
					}
				}
				
				// normalize the scores
				for (i = 0; i < tot; i++) {
					if (nucs[pos[i]] >= 0) {
						sum = rans[3*nucs[pos[i]] + 1] + rans[3*nucs[pos[i]] + 2];
						if (sum > 1) {
							rans[3*nucs[pos[i]] + 1] /= sum;
							rans[3*nucs[pos[i]] + 2] /= sum;
						} else {
							rans[3*nucs[pos[i]]] = 1 - sum;
						}
					}
				}
				
				// fold regions containing missing anchors
				int range1[2] = {0};
				int range2[2] = {0};
				for (i = 0; i < tot; i++) {
					//Rprintf("\ni = %d anchor = %d nuc = %d nucanchor = %d pos1 = %d pos2 = %d outer loop", i, anchor[i], nucs[pos[i]], anchor[i] > 0 ? nucs[pos[anchor[i]]] : 0, pos[i], anchor[i] > 0 ? pos[anchor[i]] : (anchor[i] < 0 ? pos[-1*anchor[i]] : 0));
					if (anchor[i] > 0) { // anchored
						range1[0] = i;
						range2[1] = -1; // flag as new anchor
					} else if (anchor[i] < 0 && i < -1*anchor[i] && range2[1] == -1) { // new missing left anchor(s)
						// set right boundary in alignment
						range2[1] = pos[-1*anchor[i]];
						range2[0] = range2[1];
						
						// commence search for anchor boundaries
						for (i = i + 1; i < tot; i++) {
							//Rprintf("\ni = %d anchor = %d nuc = %d nucanchor = %d pos1 = %d pos2 = %d inner loop", i, anchor[i], nucs[pos[i]], anchor[i] > 0 ? nucs[pos[anchor[i]]] : 0, pos[i], anchor[i] > 0 ? pos[anchor[i]] : (anchor[i] < 0 ? pos[-1*anchor[i]] : 0));
							if (anchor[i] > 0) { // anchored
								range1[1] = i;
								
								// search for right boundaries in the sequence
								int prev, next = 0;
								for (j = range1[0]; j < tot; j++) {
									// find the anchor closest to the left-side of the right boundaries
									if (anchor[j] > 0 && pos[j] < range2[0]) { // anchor before range2[0]
										prev = j;
									} else if (pos[j] >= range2[0]) {
										break;
									}
								}
								for (; j < tot; j++) {
									// find the anchor closest to the right-side of the right boundaries
									if (anchor[j] > 0) {
										if (pos[j] > range2[1])
											next = j;
										break;
									}
								}
								if (next == 0) // pairing within the right boundaries
									break;
								
								range1[0] = nucs[pos[range1[0]]];// + 1;
								range1[1] = nucs[pos[range1[1]]];// - 1;
								range2[0] = nucs[pos[prev]];// + 1;
								range2[1] = nucs[pos[next]];// - 1;
								//Rprintf("\nBounds left: %d - %d right %d - %d", range1[0], range1[1], range2[0], range2[1]);
								
								// find stem loops in the region between range1 and range2
								int l1 = range1[1] - range1[0];
								int l2 = range2[1] - range2[0];
								if (l1 <= 0 || l2 <= 0)
									break;
								
								int *s1 = Calloc(l1, int); // initialized to zero
								for (j = 0; j < l1; j++)
									s1[j] = 4*nogaps[range1[0] + j] + nogaps[range1[0] + j + 1];
								int *s2 = Calloc(l2, int); // initialized to zero
								for (j = 0; j < l2; j++)
									s2[j] = nogaps[range2[0] + j] + 4*nogaps[range2[0] + j + 1];
								
								double *foldm = Calloc(l1*l2, double); // initialized to zero
								int *foldn = Calloc(l1*l2, int); // initialized to zero
								
								int maxd = 0; // longest possible diagonal traceback
								int first, pi, pj;
								double gapN, gapL, gapR; // no gap, gap in range1, gap in range2
								for (pi = l1 - 1; pi >= 0; pi--) {
									first = 1;
									for (pj = 0; pj < l2; pj++) {
										if (range1[0] + pi + minLoop + 1 < range2[0] + pj) { // require >= minLoop length loops
											if (first) {
												int di = l1 - pi;
												int dj = l2 - pj;
												if (di < dj) {
													if (dj > maxd)
														maxd = dj;
												} else {
													if (di > maxd)
														maxd = di;
												}
												first = 0;
											}
											
											if (pi < l1 - 1 && pj > 0) {
												gapN = foldm[l1*(pj - 1) + pi + 1] + NN[16*s2[pj] + s1[pi]];
											} else {
												gapN = NN[16*s2[pj] + s1[pi]];
											}
											if (pi < l1 - 1) { // i is unpaired
												gapL = foldm[l1*pj + pi + 1];
											} else {
												gapL = 0;
											}
											if (pj > 0) { // j is unpaired
												gapR = foldm[l1*(pj - 1) + pi];
											} else {
												gapR = 0;
											}
											
											if (gapL > gapN && gapL > gapR) {
												foldm[l1*pj + pi] = gapL;
												foldn[l1*pj + pi] = -2;
											} else if (gapR > gapN) {
												foldm[l1*pj + pi] = gapR;
												foldn[l1*pj + pi] = -1;
											} else {
												if (NN[16*s2[pj] + s1[pi]] == 0) {
													foldn[l1*pj + pi] = -3;
												} else if (pi < l1 - 1 && pj > 0 && foldn[l1*(pj - 1) + pi + 1] > 0) {
													foldn[l1*pj + pi] = foldn[l1*(pj - 1) + pi + 1] + 1;
												} else {
													foldn[l1*pj + pi] = 1;
												}
												foldm[l1*pj + pi] = gapN;
											}
										}
									}
								}
								
								int maxpos[2];
								maxpos[0] = 0;
								maxpos[1] = l2 - 1;
								for (pi = 1; pi < l1; pi++) {
									if (foldm[l1*maxpos[1] + pi] >= foldm[l1*maxpos[1] + maxpos[0]])
										maxpos[0] = pi;
								}
								for (pj = l2 - 2; pj >= 0; pj--) {
									if (foldm[l1*pj + maxpos[0]] >= foldm[l1*maxpos[1] + maxpos[0]])
										maxpos[1] = pj;
								}
								pi = maxpos[0];
								pj = maxpos[1];
								//Rprintf("\nmyMax = %d %d", pi, pj);
								
								double *prob = Calloc(maxd, double); // initialized to zero
								for (j = 0; j < maxd; j++)
									prob[j] = -1;
								
								while(range2[0] + pj > range1[0] + pi + minLoop + 1 && pi < l1 && pj >= 0) {
									if (foldn[l1*pj + pi] == -3) { // unpaired
										if (pi + 1 >= l1 || pj - 1 < 0) {
											break;
										} else if (foldn[l1*(pj - 1) + pi + 1] < 0) {
											pi += 2;
											pj -= 2;
										} else {
											pi++;
											pj--;
										}
									} else if (foldn[l1*pj + pi] > 0) { // paired
										int val = foldn[l1*pj + pi];
										if (prob[val - 1] < 0)
											prob[val - 1] = pNoRun((double)(maxd + 1), (double)(val + 1), pMatch);
										
										for (j = 0; j < val + 1; j++) { // add new pairs
											if (prob[val - 1] > rans[3*(range1[0] + pi + j) + 1]) { // new left pair
												rans[3*(range1[0] + pi + j) + 1] = prob[val - 1];
												sum = rans[3*(range1[0] + pi + j) + 1] + rans[3*(range1[0] + pi + j) + 2];
												if (sum > 1) {
													rans[3*(range1[0] + pi + j) + 1] /= sum;
													rans[3*(range1[0] + pi + j) + 2] /= sum;
												} else if (sum > 0) {
													rans[3*(range1[0] + pi + j)] = 1 - sum;
												}
											}
											
											if (prob[val - 1] > rans[3*(range2[0] + pj + 1 - j) + 2]) { // new right pair
												rans[3*(range2[0] + pj + 1 - j) + 2] = prob[val - 1];
												sum = rans[3*(range2[0] + pj + 1 - j) + 1] + rans[3*(range2[0] + pj + 1 - j) + 2];
												if (sum > 1) {
													rans[3*(range2[0] + pj + 1 - j) + 1] /= sum;
													rans[3*(range2[0] + pj + 1 - j) + 2] /= sum;
												} else if (sum > 0) {
													rans[3*(range2[0] + pj + 1 - j)] = 1 - sum;
												}
											}
										}
										
										//Rprintf("\npairing: L %d R %d pos %d", range1[0] + pi, range2[0] + pj + 1, val + 1);
										//Rprintf(" score = %1.2f", prob[val - 1]);
										
										if (pi + val >= l1 || pj - val < 0) {
											break;
										} else if (foldn[l1*(pj - val) + pi + val] < 0) {
											pi = pi + val + 1;
											pj = pj - val - 1;
										} else {
											pi += val;
											pj -= val;
										}
									} else if (foldn[l1*pj + pi] == -2) { // gap in L
										pi++;
									} else { // gap in R
										pj--;
									}
								}
								
								Free(s1);
								Free(s2);
								Free(foldm);
								Free(foldn);
								Free(prob);
								
								i--; // allow current anchor to be reused
								break;
							} else if (anchor[i] < 0) { // missing left anchor
								//if (pos[-1*anchor[i]] > range2[0])
								//	break; // inconsistent pairings
								if (pos[-1*anchor[i]] < range2[0]) // consistent pairing
									range2[0] = pos[-1*anchor[i]]; // right boundary in alignment
							}
						}
					} // else unanchored
				}
				
				Free(anchor);
				Free(nucs);
				Free(nogaps);
			} else {
				// normalize the scores
				for (i = 0; i < tot; i++) {
					if (nucs[pos[i]] >= 0) {
						sum = rans[3*nucs[pos[i]] + 1] + rans[3*nucs[pos[i]] + 2];
						if (sum > 1) {
							rans[3*nucs[pos[i]] + 1] /= sum;
							rans[3*nucs[pos[i]] + 2] /= sum;
						} else {
							rans[3*nucs[pos[i]]] = 1 - sum;
						}
					}
				}
			}
			
			SET_VECTOR_ELT(ans, s, ans_s);
			UNPROTECT(1); // ans_s
		}
		
		Free(leftMax);
		Free(rightMax);
	}
	
	Free(MI);
	Free(pos);
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}

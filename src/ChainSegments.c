/****************************************************************************
 *                         Chains Graph of Segments                         *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

// for math functions
#include <math.h>

/* for Calloc/Free */
#include <R_ext/RS.h>

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

// in-place fills missing sequence in the input
// replaces values in consecutive positions that are discontinuous
// e.g., 51703 51656 51705 would become 51703 51704 51705
SEXP fillOverlaps(SEXP m, SEXP n)
{
	if (MAYBE_SHARED(m))
		error(".Call function 'fillOverlaps' called in incorrect context.");
	
	int i;
	int *x = INTEGER(m);
	int l = length(m);
	int y = asInteger(n);
	
	int k = y - 1; // current position
	int last = k; // last non-incrementing position
	int p = 0; // last position
	
	while (k < l) {
		if (x[k] != NA_INTEGER &&
			x[p] != NA_INTEGER &&
			x[k] == (x[p] + y - 1)) {
			if (last > p) {
				for (i = p + 1; i < k; i++)
					x[i] = x[i - 1] + 1;
				last = p;
			}
		} else {
			last = k;
		}
		k++;
		p++;
	}
	
	return m;
}

// adjusts starts and ends to be within widths
SEXP indexByContig(SEXP starts, SEXP ends, SEXP order, SEXP index, SEXP widths)
{
	int j, k, p;
	int *o = INTEGER(order);
	int *w = INTEGER(widths);
	int *i = INTEGER(index);
	int l = length(starts);
	
	SEXP ans1, ans2, ans3;
	PROTECT(ans1 = allocVector(INTSXP, l));
	int *rans = INTEGER(ans1);
	PROTECT(ans2 = duplicate(starts));
	int *s = INTEGER(ans2);
	PROTECT(ans3 = duplicate(ends));
	int *e = INTEGER(ans3);
	
	// fill initial values
	for (j = 0; j < l; j++) {
		p = o[j] - 1;
		if (s[p] > w[0])
			break;
		rans[p] = i[0];
	}
	
	// index sequences by contig
	k = 1;
	for (; j < l; j++) { // j = j
		p = o[j] - 1;
		while (s[p] > w[k]) {
			k++;
		}
		s[p] = s[p] - w[k - 1];
		e[p] = e[p] - w[k - 1];
		rans[p] = i[k];
	}
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	SET_VECTOR_ELT(ret_list, 2, ans3);
	
	UNPROTECT(4);
	
	return ret_list;
}

SEXP chainSegments(SEXP x_s, SEXP x_e, SEXP x_i, SEXP x_f, SEXP y_s, SEXP y_e, SEXP y_i, SEXP y_f, SEXP weights, SEXP sepCost, SEXP sepPower, SEXP gapCost, SEXP gapPower, SEXP shiftCost, SEXP codingCost, SEXP maxSep, SEXP maxGap, SEXP ordering, SEXP minScore, SEXP maxW, SEXP allowOverlap)
{
	int i, k, max, dy, dx, sep, gap, xok, xoj;
	int *xs = INTEGER(x_s);
	int *xe = INTEGER(x_e);
	int *xi = INTEGER(x_i);
	int *xf = INTEGER(x_f);
	int *ys = INTEGER(y_s);
	int *ye = INTEGER(y_e);
	int *yi = INTEGER(y_i);
	int *yf = INTEGER(y_f);
	double *we = REAL(weights);
	double sepC = asReal(sepCost);
	double sepP = asReal(sepPower);
	double gapC = asReal(gapCost);
	double gapP = asReal(gapPower);
	double shiC = asReal(shiftCost);
	double codC = asReal(codingCost);
	double maxS = asReal(maxSep);
	double maxG = asReal(maxGap);
	double totW = asReal(maxW);
	int *xo = INTEGER(ordering);
	double minS = asReal(minScore);
	int aO = asInteger(allowOverlap);
	
	// initialize an array of sep multipliers
	double *SEPS = Calloc(maxS + 1, double); // initialized to zero
	// initialize an array of gap multipliers
	double *GAPS = Calloc(maxG + 1, double); // initialized to zero
	for (i = 1; i <= maxS; i++)
		SEPS[i] = pow(i, sepP)*sepC;
	for (i = 1; i <= maxG; i++)
		GAPS[i] = pow(i, gapP)*gapC;
	
	int l = length(x_s);
	i = 0;
	int j = -1;
	int prev = 0;
	double cost = log(totW); // cost of a single search
	minS += cost; // require additional cost removed at end
	double temp, score;
	
	// initialize an array of activities
	int *A = Calloc(l, int); // initialized to zero
	// initialize an array of scores
	double *S = Calloc(l, double); // initialized to zero
	// initialize an array of return indicies
	int *R = Calloc(l, int); // initialized to zero
	// initialize an array of origin indices
	int *O = Calloc(l, int); // initialized to zero
	
	while (i < l || j < (l - 1)) {
		if (i < l &&
			xs[i] <= xe[xo[j + 1]] &&
			xi[i] == xi[xo[j + 1]]) {
			// use i, left of rectangle
			if (i > 0 &&
				xi[i] != xi[i - 1])
				prev = 0; // deactivate left
			
			while (prev > 0 &&
				(xs[i] - xe[xo[j - prev + 1]]) > maxS) {
				prev--;
			}
			
			// find the highest scoring rectangle below the start
			max = -1;
			score = 0;
			for (k = j - prev + 1; k <= j; k++) {
				xok = xo[k];
				
				if (A[xok] != 1 ||
					yi[i] != yi[xok])
					continue;
				
				dy = ys[i] - ye[xok] - 1;
				if (dy <= 0)
					continue;
				
				dx = xs[i] - xe[xok] - 1;
				if (dx <= 0)
					continue;
				
				if (dx > dy) {
					sep = dy;
					gap = dx - dy;
				} else {
					sep = dx;
					gap = dy - dx;
				}
				if (sep > maxS ||
					gap > maxG)
					continue;
				
				// add cost for shifting reading frames
				if (xf[i] == xf[xok] &&
					(xf[i] == 0 || (dx % 3) == 0)) {
					temp = 0;
				} else {
					if (xf[i] == 0 || xf[xok] == 0) {
						temp = codC;
					} else {
						temp = shiC;
					}
				}
				if (!(yf[i] == yf[xok] &&
					(yf[i] == 0 || (dy % 3) == 0))) {
					if (yf[i] == 0 || yf[xok] == 0) {
						temp = codC;
					} else {
						temp = shiC;
					}
				}
				
				temp += SEPS[sep];
				temp += GAPS[gap];
				temp += S[xok];
				//temp += we[i];
				//if (sep > 0 && gap > 0) {
				//	temp -= log((double)(sep*gap));
				//} else if (sep > 0) {
				//	temp -= log((double)sep);
				//} else if (gap > 0) {
				//	temp -= log((double)gap);
				//}
				//if (xf[xok] == 0) { // nucleotide hit
				//	temp -= log((double)(sep + gap + 1)/(xe[xok] - xs[xok] + 1)*totW);
				//} else {
				//	temp -= log((double)(sep + gap + 1)/(xe[xok] - xs[xok] + 1)*totW/3);
				//}
				
				if (temp > score) {
					max = xok;
					score = temp;
				}
			}
			//Rprintf("\ni = %d score = %1.2f xo[i] = %d we[i] = %1.2f", i, score, xo[i], we[i]);
			if (score > 0) {
				S[i] = score + we[i];
				R[i] = max;// + 1;
				O[i] = O[max];// + 1;
			} else {
				S[i] = we[i];
				O[i] = i;// + 1;
				R[i] = -1;//0;
			}
			
			i++;
		} else {
			// use j, right of rectangle
			j++;
			xoj = xo[j];
			
			// find the nearest rectangle below the end
			max = -1;
			score = 1e53;
			for (k = j - prev; k < j; k++) {
				xok = xo[k];
				if (A[xok] != 1 ||
					yi[xoj] != yi[xok])
					continue;
				
				dy = ye[xoj] - ye[xok];
				if (dy >= 0 &&
					dy < score) {
					max = xok;
					score = dy;
				}
			}
			
			if (max != -1 &&
				S[xoj] > S[max]) {
				// find rectangles above the end with lower score and the same origin
				for (k = j - prev; k < j; k++) {
					xok = xo[k];
					if (ye[xok] >= ye[xoj] &&
						S[xok] < S[xoj] &&
						yi[xok] == yi[xoj] &&
						O[xok] == O[xoj])
						A[xok] = 0; // deactivate
				}
			}
			
			A[xoj] = 1;
			prev++;
		}
	}
	
	int count = 0;
	int *p;
	int n, min, overlap;
	int size = 1000;
	int **ptrs = Calloc(size, int *); // chains
	int *lens = Calloc(size, int); // length of each chain
	double *scores = Calloc(size, double); // score of each chain
	// start and end of the rectangle encompassing each chain
	int *rectXS = Calloc(size, int);
	int *rectXE = Calloc(size, int);
	int *rectYS = Calloc(size, int);
	int *rectYE = Calloc(size, int);
	int *rectXI = Calloc(size, int);
	int *rectYI = Calloc(size, int);
	
	for (i = 0; i < l; i++)
		A[i] = 1; // re-activate all
	
	while (1) {
		// find the highest active score
		max = 0;
		score = 0;
		for (i = 0; i < l; i++) {
			if (A[i] == 1 &&
				S[i] > score) {
				score = S[i];
				max = i;
			}
		}
		if (score < minS)
			break;
		
		j = 0; // length of chain
		// find the maximum length chain
		// with the same origin as max
		// and passing through max score
		for (i = l - 1; i >= O[max]; i--) {
			if (O[i] == O[max] &&
				A[i] == 1) {
				n = 1;
				k = i;
				if (S[k] == score) {
					overlap = 1;
				} else {
					overlap = 0;
				}
				A[i] = 0; // deactivate
				while (R[k] != -1) {
					k = R[k];
					if (A[k] != 0)
						A[k] = 0; // deactivate
					if (S[k] == score)
						overlap = 1;
					n++;
				}
				if (overlap == 1 &&
					n > j) {
					max = i;
					j = n;
				}
			}
		}
		
		// score must be increasing
		while (R[max] != -1 &&
			S[max] < S[R[max]]) {
			max = R[max];
		}
		
		// pull back overlapping end
		if (aO) { // allow overlaps in x or y
			while (max != -1 &&
				max >= O[max]) {
				overlap = 0;
				for (i = 0; i < count; i++) {
					if (xi[max] == rectXI[i] &&
						yi[max] == rectYI[i] &&
						xe[max] >= rectXS[i] &&
						xe[max] <= rectXE[i] &&
						ye[max] >= rectYS[i] &&
						ye[max] <= rectYE[i]) {
						overlap = 1;
						break;
					}
				}
				if (overlap == 0)
					break;
				max = R[max];
			}
		} else {
			while (max != -1 &&
				max >= O[max]) {
				overlap = 0;
				for (i = 0; i < count; i++) {
					if (xi[max] == rectXI[i] &&
						yi[max] == rectYI[i] &&
						((xe[max] >= rectXS[i] &&
						xe[max] <= rectXE[i]) ||
						(ye[max] >= rectYS[i] &&
						ye[max] <= rectYE[i]))) {
						overlap = 1;
						break;
					}
				}
				if (overlap == 0)
					break;
				max = R[max];
			}
		}
		
		if (max == -1) // no chain
			continue;
		
		// pull back overlapping start
		min = max;
		n = min;
		j = 0;
		if (aO) { // allow overlaps in x or y
			while (n >= O[max]) {
				overlap = 0;
				for (i = 0; i < count; i++) {
					if (xi[n] == rectXI[i] &&
						yi[n] == rectYI[i] &&
						xs[n] >= rectXS[i] &&
						xs[n] <= rectXE[i] &&
						ys[n] >= rectYS[i] &&
						ys[n] <= rectYE[i]) {
						overlap = 1;
						break;
					}
				}
				if (overlap == 1)
					break;
				j++;
				min = n; // prior value of min
				n = R[min];
			}
		} else {
			while (n >= O[max]) {
				overlap = 0;
				for (i = 0; i < count; i++) {
					if (xi[n] == rectXI[i] &&
						yi[n] == rectYI[i] &&
						((xs[n] >= rectXS[i] &&
						xs[n] <= rectXE[i]) ||
						(ys[n] >= rectYS[i] &&
						ys[n] <= rectYE[i]))) {
						overlap = 1;
						break;
					}
				}
				if (overlap == 1)
					break;
				j++;
				min = n; // prior value of min
				n = R[min];
			}
		}
		
		if (j == 0)
			continue; // fully overlapping
		
		// shorten to reach minScore
		score = S[max] - S[min] + we[min];
		while (score < minS && max > min) {
			max = R[max]; // shorten from end
			j--;
			score = S[max] - S[min] + we[min];
		}
		
		if (score >= minS) {
			// find the nearest rectangle
			int minDx = 2e9, minDy = 2e9, minX = -1, minY = -2, merge = 0, upX, upY;
			for (i = 0; i < count; i++) {
				if (xi[max] == rectXI[i] &&
					yi[max] == rectYI[i]) {
					dx = xs[min] - rectXE[i];
					if (dx < 0) {
						dx = rectXS[i] - xe[max];
					}
					if (dx < minDx) {
						minDx = dx;
						minX = i;
						upX = (xs[min] > rectXE[i]) ? 1 : 0;
					}
					dy = ys[min] - rectYE[i];
					if (dy < 0) {
						dy = rectYS[i] - ye[max];
					}
					if (dy < minDy) {
						minDy = dy;
						minY = i;
						upY = (ys[min] > rectYE[i]) ? 1 : 0;
					}
				}
			}
			if (aO) { // allow overlaps in x or y
				if (minDx < 0 && minDy < 0) // completely within rectangle
					continue;
			} else {
				if (minDx < 0 || minDy < 0) // completely overlapping in x or y
					continue;
			}
			
			sep = 1e9;
			gap = 1e9;
			if (minX == minY && upX == upY) {
				if (minDx < minDy) {
					sep = minDx;
					gap = minDy - minDx;
				} else {
					sep = minDy;
					gap = minDx - minDy;
				}
				if (gap <= maxG && sep <= maxS) {
					temp = scores[minX] + score + SEPS[sep] + GAPS[gap];
					if (temp >= minS &&
						minDx > 0 &&
						minDy > 0)
						merge = 1;
				}
			}
			
			if (merge) {
				ptrs[minX] = Realloc(ptrs[minX], j + lens[minX], int);
				p = ptrs[minX];
				
				if (upX) { // new chain is last
					j = j + lens[minX];
					lens[minX] = j;
					
					// add new chain at the end
					i = max;
					p[--j] = i;
					while (R[i] >= min) {
						i = R[i];
						p[--j] = i;
					}
					
					rectXE[minX] = xe[max];
					rectYE[minX] = ye[max];
				} else { // new chain is first
					// shift old chain to the end
					for (i = lens[minX] - 1; i >= 0; i--)
						p[i + j] = p[i];
					
					lens[minX] = j + lens[minX];
					
					// add new chain at the beginning
					i = max;
					p[--j] = i;
					while (R[i] >= min) {
						i = R[i];
						p[--j] = i;
					}
					
					rectXS[minX] = xs[min];
					rectYS[minX] = ys[min];
				}
				
				scores[count] = temp;
				
				continue;
			}
			
			if (count >= size) {
				size += 1000;
				ptrs = Realloc(ptrs, size, int *);
				lens = Realloc(lens, size, int);
				scores = Realloc(scores, size, double);
				rectXS = Realloc(rectXS, size, int);
				rectXE = Realloc(rectXE, size, int);
				rectYS = Realloc(rectYS, size, int);
				rectYE = Realloc(rectYE, size, int);
				rectXI = Realloc(rectXI, size, int);
				rectYI = Realloc(rectYI, size, int);
			}
			
			ptrs[count] = Calloc(j, int);
			lens[count] = j;
			p = ptrs[count];
			
			i = max;
			p[--j] = i;
			while (R[i] >= min) {
				i = R[i];
				p[--j] = i;
			}
			
			scores[count] = score;
			
			rectXI[count] = xi[max];
			rectYI[count] = yi[max];
			rectXS[count] = xs[min];
			rectXE[count] = xe[max];
			rectYS[count] = ys[min];
			rectYE[count] = ye[max];
			
			count++;
		}
	}
	
	Free(SEPS);
	Free(GAPS);
	Free(A);
	Free(S);
	Free(R);
	Free(O);
	Free(rectXS);
	Free(rectXE);
	Free(rectYS);
	Free(rectYE);
	Free(rectXI);
	Free(rectYI);
	
	SEXP ret, chains, chain, cs;
	int *pchain;
	double *pcs;
	PROTECT(cs = allocVector(REALSXP, count));
	pcs = REAL(cs);
	PROTECT(chains = allocVector(VECSXP, count));
	for (i = 0; i < count; i++) {
		PROTECT(chain = allocVector(INTSXP, lens[i]));
		pchain = INTEGER(chain);
		
		p = ptrs[i];
		for (j = 0; j < lens[i]; j++)
			pchain[j] = p[j] + 1;
		
		Free(p);
		
		SET_VECTOR_ELT(chains, i, chain);
		UNPROTECT(1);
		
		pcs[i] = scores[i] - cost; // apply cost for first search
	}
	
	Free(ptrs);
	Free(lens);
	Free(scores);
	
	PROTECT(ret = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ret, 0, chains);
	SET_VECTOR_ELT(ret, 1, cs);
	
	UNPROTECT(3);
	
	return ret;
}

SEXP extendMatches(SEXP X1, SEXP X2, SEXP starts1, SEXP ends1, SEXP index1, SEXP starts2, SEXP ends2, SEXP index2, SEXP width1, SEXP width2, SEXP subMatrix, SEXP letters, SEXP dropScore, SEXP nThreads)
{
	int i, p1, p2, lkup1, lkup2, boundL1, boundR1, boundL2, boundR2;
	
	XStringSet_holder x1_set = hold_XStringSet(X1);
	XStringSet_holder x2_set = hold_XStringSet(X2);
	Chars_holder x1 = get_elt_from_XStringSet_holder(&x1_set, 0);
	Chars_holder x2 = get_elt_from_XStringSet_holder(&x2_set, 0);
	
	int l = length(starts1);
	int *s1 = INTEGER(starts1);
	int *e1 = INTEGER(ends1);
	int *i1 = INTEGER(index1);
	int *s2 = INTEGER(starts2);
	int *e2 = INTEGER(ends2);
	int *i2 = INTEGER(index2);
	int *w1 = INTEGER(width1);
	int *w2 = INTEGER(width2);
	
	SEXP ans0, ans1, ans2, ans3, ans4;
	PROTECT(ans0 = allocVector(REALSXP, l));
	double *rans0 = REAL(ans0);
	PROTECT(ans1 = allocVector(INTSXP, l));
	int *rans1 = INTEGER(ans1);
	PROTECT(ans2 = allocVector(INTSXP, l));
	int *rans2 = INTEGER(ans2);
	PROTECT(ans3 = allocVector(INTSXP, l));
	int *rans3 = INTEGER(ans3);
	PROTECT(ans4 = allocVector(INTSXP, l));
	int *rans4 = INTEGER(ans4);
	
	double *sM = REAL(subMatrix);
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
	
	double dS = asReal(dropScore);
	int nthreads = asInteger(nThreads);
	
	#ifdef _OPENMP
	#pragma omp parallel for private(i,p1,p2,lkup1,lkup2,boundL1,boundL2,boundR1,boundR2) schedule(dynamic) num_threads(nthreads)
	#endif
	for (i = 0; i < l; i++) {
		double delta;
		
		rans0[i] = 0; // score
		rans1[i] = s1[i]; // starts1
		rans2[i] = e1[i]; // ends1
		rans3[i] = s2[i]; // starts2
		rans4[i] = e2[i]; // ends2
		if (i1[i] == 1) {
			boundL1 = 0;
		} else {
			boundL1 = w1[i1[i] - 2];
		}
		if (i2[i] == 1) {
			boundL2 = 0;
		} else {
			boundL2 = w2[i2[i] - 2];
		}
		boundR1 = w1[i1[i] - 1] - 1;
		boundR2 = w2[i2[i] - 1] - 1;
		
		// accumulate score for region
		p1 = s1[i] + boundL1 - 1;
		p2 = s2[i] + boundL2 - 1;
		while (p1 < e1[i] + boundL1) {
			lkup1 = lkup_row[(unsigned char)x1.ptr[p1]];
			lkup2 = lkup_col[(unsigned char)x2.ptr[p2]];
			if (lkup1 != NA_INTEGER && lkup2 != NA_INTEGER)
				rans0[i] += sM[lkup1 + lkup2];
			p1++;
			p2++;
		}
		
		// try extending left
		p1 = s1[i] + boundL1 - 1;
		p2 = s2[i] + boundL2 - 1;
		delta = 0;
		while (p1 > boundL1 && p2 > boundL2 && delta > dS) {
			p1--;
			p2--;
			lkup1 = lkup_row[(unsigned char)x1.ptr[p1]];
			lkup2 = lkup_col[(unsigned char)x2.ptr[p2]];
			if (lkup1 == NA_INTEGER || lkup2 == NA_INTEGER) {
				break;
			} else {
				delta += sM[lkup1 + lkup2];
			}
			if (delta > 0) {
				rans0[i] += delta;
				delta = 0;
				rans1[i] = p1 + 1 - boundL1;
				rans3[i] = p2 + 1 - boundL2;
			}
		}
		
		// try extending right
		p1 = e1[i] + boundL1 - 1;
		p2 = e2[i] + boundL2 - 1;
		delta = 0;
		while (p1 < boundR1 && p2 < boundR2 && delta > dS) {
			p1++;
			p2++;
			lkup1 = lkup_row[(unsigned char)x1.ptr[p1]];
			lkup2 = lkup_col[(unsigned char)x2.ptr[p2]];
			if (lkup1 == NA_INTEGER || lkup2 == NA_INTEGER) {
				break;
			} else {
				delta += sM[lkup1 + lkup2];
			}
			if (delta > 0) {
				rans0[i] += delta;
				delta = 0;
				rans2[i] = p1 + 1 - boundL1;
				rans4[i] = p2 + 1 - boundL2;
			}
		}
	}
	free(lkup_row);
	free(lkup_col);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 5));
	SET_VECTOR_ELT(ret_list, 0, ans0);
	SET_VECTOR_ELT(ret_list, 1, ans1);
	SET_VECTOR_ELT(ret_list, 2, ans2);
	SET_VECTOR_ELT(ret_list, 3, ans3);
	SET_VECTOR_ELT(ret_list, 4, ans4);
	
	UNPROTECT(6);
	
	return ret_list;
}

SEXP withdrawMatches(SEXP order, SEXP starts1, SEXP ends1, SEXP index1, SEXP starts2, SEXP ends2, SEXP index2, SEXP width1, SEXP width2, SEXP score, SEXP bufferSize)
{
	int i, j, p1, p2, boundL1, boundL2, temp1, temp2;
	
	int *o = INTEGER(order);
	int l = length(starts1);
	int *s1 = INTEGER(starts1);
	int *e1 = INTEGER(ends1);
	int *i1 = INTEGER(index1);
	int *s2 = INTEGER(starts2);
	int *e2 = INTEGER(ends2);
	int *i2 = INTEGER(index2);
	int *w1 = INTEGER(width1);
	int l1 = w1[length(width1) - 1];
	int *w2 = INTEGER(width2);
	int l2 = w2[length(width2) - 1];
	double *s = REAL(score);
	int buffer = asInteger(bufferSize);
	
	SEXP ans0, ans1, ans2, ans3, ans4;
	PROTECT(ans0 = allocVector(REALSXP, l));
	double *rans0 = REAL(ans0);
	PROTECT(ans1 = allocVector(INTSXP, l));
	int *rans1 = INTEGER(ans1);
	PROTECT(ans2 = allocVector(INTSXP, l));
	int *rans2 = INTEGER(ans2);
	PROTECT(ans3 = allocVector(INTSXP, l));
	int *rans3 = INTEGER(ans3);
	PROTECT(ans4 = allocVector(INTSXP, l));
	int *rans4 = INTEGER(ans4);
	
	int *pos1 = (int *) malloc(l1*sizeof(int)); // thread-safe on Windows
	int *pos2 = (int *) malloc(l2*sizeof(int)); // thread-safe on Windows
	
	for (i = 0; i < l1; i++)
		pos1[i] = buffer;
	for (i = 0; i < l2; i++)
		pos2[i] = buffer;
	
	for (j = l - 1; j >= 0; j--) {
		i = o[j] - 1; // index
		if (i1[i] == 1) {
			boundL1 = 0;
		} else {
			boundL1 = w1[i1[i] - 2];
		}
		if (i2[i] == 1) {
			boundL2 = 0;
		} else {
			boundL2 = w2[i2[i] - 2];
		}
		
		// pull back overlapping start
		p1 = s1[i] + boundL1 - 1;
		p2 = s2[i] + boundL2 - 1;
		while (p1 < e1[i] + boundL1 - 1 &&
			((p2 + 1 >= pos1[p1] + buffer && p2 + 1 <= pos1[p1] - buffer) ||
			(p1 + 1 >= pos2[p2] + buffer && p1 + 1 <= pos2[p2] - buffer))) {
			p1++;
			p2++;
		}
		rans1[i] = p1 + 1 - boundL1;
		rans3[i] = p2 + 1 - boundL2;
		
		// pull back overlapping end
		p1 = e1[i] + boundL1 - 1;
		p2 = e2[i] + boundL2 - 1;
		while (p1 >= rans1[i] + boundL1 &&
			((p2 + 1 >= pos1[p1] + buffer && p2 + 1 <= pos1[p1] - buffer) ||
			(p1 + 1 >= pos2[p2] + buffer && p1 + 1 <= pos2[p2] - buffer))) {
			p1--;
			p2--;
		}
		rans2[i] = p1 + 1 - boundL1;
		rans4[i] = p2 + 1 - boundL2;
		
		// mark positions as occupied
		p1 = rans2[i] + boundL1;
		p2 = rans4[i] + boundL2;
		while (p1 >= rans1[i] + boundL1) {
			temp1 = p1;
			temp2 = p2;
			p1--;
			p2--;
			if (pos1[p1] == buffer)
				pos1[p1] = temp2;
			if (pos2[p2] == buffer)
				pos2[p2] = temp1;
		}
		
		// lower score proportionally
		rans0[i] = s[i]*((double)(rans2[i] - rans1[i] + 1)/(double)(e1[i] - s1[i] + 1));
	}
	free(pos1);
	free(pos2);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 5));
	SET_VECTOR_ELT(ret_list, 0, ans0);
	SET_VECTOR_ELT(ret_list, 1, ans1);
	SET_VECTOR_ELT(ret_list, 2, ans2);
	SET_VECTOR_ELT(ret_list, 3, ans3);
	SET_VECTOR_ELT(ret_list, 4, ans4);
	
	UNPROTECT(6);
	
	return ret_list;
}

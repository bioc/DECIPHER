/****************************************************************************
 *                       Performs Pairwise Alignment                        *
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

// for calloc/free
#include <stdlib.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"

// strcpy
#include <string.h>

// for time and difftime
#include <time.h>

// DECIPHER header file
#include "DECIPHER.h"

static void integerEncode(const Chars_holder *P, int start, int end, int* v)
{
	int i, j;
	const char *p;
	
	for (i = 0, j = start, p = (P->ptr + start);
		j <= end;
		i++, j++, p++)
	{
		switch (*p) {
			case 1: // A
				v[i] = 0;
				break;
			case 2: // C
				v[i] = 1;
				break;
			case 4: // G
				v[i] = 2;
				break;
			case 8: // T
				v[i] = 3;
				break;
			case 6: // S
				v[i] = 1; // CG
				break;
			case 10: // Y
				v[i] = 1; // CT
				break;
			case 12: // K
				v[i] = 2; // GT
				break;
			case 14: // B
				v[i] = 1; // CGT
				break;
			default: // treat as A
				v[i] = 0;
				break;
		}
	}
}

static void integerEncodeAA(const Chars_holder *P, int start, int end, int* v)
{
	int i, j;
	const char *p;
	
	for (i = 0, j = start, p = (P->ptr + start);
		j <= end;
		i++, j++, p++)
	{
		switch (*p) {
			case 65: // A
				v[i] = 0;
				break;
			case 82: // R
				v[i] = 1;
				break;
			case 78: // N
				v[i] = 2;
				break;
			case 68: // D
				v[i] = 3;
				break;
			case 67: // C
				v[i] = 4;
				break;
			case 81: // Q
				v[i] = 5;
				break;
			case 69: // E
				v[i] = 6;
				break;
			case 71: // G
				v[i] = 7;
				break;
			case 72: // H
				v[i] = 8;
				break;
			case 73: // I
				v[i] = 9;
				break;
			case 76: // L
				v[i] = 10;
				break;
			case 75: // K
				v[i] = 11;
				break;
			case 77: // M
				v[i] = 12;
				break;
			case 70: // F
				v[i] = 13;
				break;
			case 80: // P
				v[i] = 14;
				break;
			case 83: // S
				v[i] = 15;
				break;
			case 84: // T
				v[i] = 16;
				break;
			case 87: // W
				v[i] = 17;
				break;
			case 89: // Y
				v[i] = 18;
				break;
			case 86: // V
				v[i] = 19;
				break;
			case 66: // B = N or D
				v[i] = 2;
				break;
			case 90: // Z = Q or E
				v[i] = 5;
				break;
			case 74: // J = I or L
				v[i] = 9;
				break;
			default: // treat as A
				v[i] = 0;
				break;
		}
	}
}

SEXP alignPair(SEXP x, SEXP y, SEXP s1, SEXP e1, SEXP s2, SEXP e2, SEXP go, SEXP ge, SEXP tg, SEXP maxLength, SEXP type, SEXP subMatrix, SEXP nThreads)
{
	int i, l1, l2, i1, j1, i2, j2, d, l, u;
	int *p1, *p2, *p3, *p4;
	int *S1 = INTEGER(s1);
	int *E1 = INTEGER(e1);
	int *S2 = INTEGER(s2);
	int *E2 = INTEGER(e2);
	int *SM = INTEGER(subMatrix);
	int GO = asInteger(go); // gap opening
	int GE = asInteger(ge); // gap extension
	int TG = asInteger(tg); // terminal gaps
	int ML = asInteger(maxLength); // maximum length to skip alignment of equal-length regions
	int t = asInteger(type); // type of XStringSet
	int nthreads = asInteger(nThreads);
	int n = length(s1); // number of regions
	
	XStringSet_holder x_set, ans_holder;
	Chars_holder X, Y;
	x_set = hold_XStringSet(x);
	X = get_elt_from_XStringSet_holder(&x_set, INTEGER(y)[0] - 1);
	Y = get_elt_from_XStringSet_holder(&x_set, INTEGER(y)[1] - 1);
	
	int **P1 = (int **) calloc(n, sizeof(int *)); // initialized to zero (thread-safe on Windows)
	int **P2 = (int **) calloc(n, sizeof(int *)); // initialized to zero (thread-safe on Windows)
	int **P3 = (int **) calloc(n, sizeof(int *)); // initialized to zero (thread-safe on Windows)
	int **P4 = (int **) calloc(n, sizeof(int *)); // initialized to zero (thread-safe on Windows)
	int *N1 = (int *) calloc(n, sizeof(int)); // initialized to zero (thread-safe on Windows)
	int *N2 = (int *) calloc(n, sizeof(int)); // initialized to zero (thread-safe on Windows)
	
	int *square;
	if (t == 3) { // AAStringSet
		square = (int *) malloc(20*sizeof(int)); // thread-safe on Windows
		for (i = 0; i < 20; i++)
			square[i] = i*20;
	} else { // DNAStringSet or RNAStringSet
		square = (int *) malloc(4*sizeof(int)); // thread-safe on Windows
		for (i = 0; i < 4; i++)
			square[i] = i*4;
	}
	
	#ifdef _OPENMP
	#pragma omp parallel for private(i,l1,l2,i1,j1,i2,j2,d,l,u,p1,p2,p3,p4) num_threads(nthreads)
	#endif
	for (i = 0; i < n; i++) { // each region
		l1 = E1[i] - S1[i];
		l2 = E2[i] - S2[i];
		if (l1 < 0 && l2 < 0) {
			continue; // nothing to align
		} else if (l1 > 0 && l2 < 0) {
			N1[i] = 1;
			p1 = (int *) malloc(sizeof(int)); // thread-safe on Windows
			p2 = (int *) malloc(sizeof(int)); // thread-safe on Windows
			p1[0] = 0;
			p2[0] = l1 + 1;
			P1[i] = p1;
			P2[i] = p2;
			continue;
		} else if (l2 > 0 && l1 < 0) {
			N2[i] = 1;
			p3 = (int *) malloc(sizeof(int)); // thread-safe on Windows
			p4 = (int *) malloc(sizeof(int)); // thread-safe on Windows
			p3[0] = 0;
			p4[0] = l2 + 1;
			P3[i] = p3;
			P4[i] = p4;
			continue;
		} else if (l1 == l2 && l1 < ML) {
			continue; // assume already aligned
		} // else need to align
		
		l1++;
		l2++;
		
		int *v1 = (int *) malloc(l1*sizeof(int)); // thread-safe on Windows
		int *v2 = (int *) malloc(l2*sizeof(int)); // thread-safe on Windows
		
		if (t == 3) { // AAStringSet
			integerEncodeAA(&X, S1[i] - 1, E1[i] - 1, v1);
			integerEncodeAA(&Y, S2[i] - 1, E2[i] - 1, v2);
		} else { // DNAStringSet or RNAStringSet
			integerEncode(&X, S1[i] - 1, E1[i] - 1, v1);
			integerEncode(&Y, S2[i] - 1, E2[i] - 1, v2);
		}
		
		int *index = (int *) malloc((l2 + 1)*sizeof(int));
		for (i1 = 0; i1 <= l2; i1++)
			index[i1] = i1*(l1 + 1);
		
		// initialize matrix
		int *m = (int *) malloc((l1 + 1)*(l2 + 1)*sizeof(int)); // thread-safe on Windows
		int *o = (int *) malloc((l1 + 1)*(l2 + 1)*sizeof(int)); // thread-safe on Windows
		m[0] = 0;
		o[0] = 0;
		
		// fill gap opening at beginning
		if (i == 0) {
			i1 = 0;
			j1 = 1;
			while (j1 <= l2) {
				m[i1 + index[j1]] = TG;
				o[i1 + index[j1]] = j1;
				j1++;
			}
			i1 = 1;
			j1 = 0;
			while (i1 <= l1) {
				m[i1 + index[j1]] = TG;
				o[i1 + index[j1]] = -1*i1;
				i1++;
			}
		} else {
			i1 = 0;
			j1 = 1;
			i2 = GO;
			while (j1 <= l2) {
				i2 += GE;
				m[i1 + index[j1]] = i2;
				o[i1 + index[j1]] = j1;
				j1++;
			}
			i1 = 1;
			j1 = 0;
			j2 = GO;
			while (i1 <= l1) {
				j2 += GE;
				m[i1 + index[j1]] = j2;
				o[i1 + index[j1]] = -1*i1;
				i1++;
			}
		}
		
		j2 = 1;
		j1 = 0;
		while (j2 <= l2) {
			i2 = 1;
			i1 = 0;
			while (i2 <= l1) {
				d = m[i1 + index[j1]] + SM[v1[i1] + square[v2[j1]]];
				if (o[i1 + index[j2]] < 0) {
					u = m[i1 + index[j2]] + GE;
				} else {
					u = m[i1 + index[j2]] + GO + GE;
				}
				if (o[i2 + index[j1]] > 0) {
					l = m[i2 + index[j1]] + GE;
				} else {
					l = m[i2 + index[j1]] + GO + GE;
				}
				if (d >= u && d >= l) { // diagonal
					o[i2 + index[j2]] = 0;
					m[i2 + index[j2]] = d;
				} else if (u >= l) { // up
					if (o[i1 + index[j2]] < 0) {
						o[i2 + index[j2]] = o[i1 + index[j2]] - 1;
					} else {
						o[i2 + index[j2]] = -1;
					}
					m[i2 + index[j2]] = u;
				} else { // left
					if (o[i2 + index[j1]] > 0) {
						o[i2 + index[j2]] = o[i2 + index[j1]] + 1;
					} else {
						o[i2 + index[j2]] = 1;
					}
					m[i2 + index[j2]] = l;
				}
				i1 = i2;
				i2++;
			}
			j1 = j2;
			j2++;
		}
		free(v1);
		free(v2);
		
		// fill gap closing at end
		if (i == n - 1) {
			i1 = l1 - 1;
			// j1 == l2
			while (i1 >= 0) {
				m[i1 + index[j1]] += TG;
				i1--;
			}
			i1 = l1;
			j1 = l2 - 1;
			while (j1 >= 0) {
				m[i1 + index[j1]] += TG;
				j1--;
			}
		} else {
			i1 = l1 - 1;
			j1 = l2;
			i2 = GO;
			while (i1 >= 0) {
				i2 += GE;
				m[i1 + index[j1]] += i2;
				i1--;
			}
			i1 = l1;
			j1 = l2 - 1;
			j2 = GO;
			while (j1 >= 0) {
				j2 += GE;
				m[i1 + index[j1]] += j2;
				j1--;
			}
		}
		
		// find the maximum score
		i2 = l1;
		j2 = l2;
		if (l2 > 0) {
			i1 = l1 - 1;
			j1 = l2;
			while (i1 >= 0) {
				if (m[i1 + index[j1]] > m[i2 + index[j2]]) {
					i2 = i1;
					j2 = j1;
				}
				i1--;
			}
		}
		if (l1 > 0) {
			i1 = l1;
			j1 = l2 - 1;
			while (j1 >= 0) {
				if (m[i1 + index[j1]] > m[i2 + index[j2]]) {
					i2 = i1;
					j2 = j1;
				}
				j1--;
			}
		}
		free(m);
		
		i1 = i2;
		j1 = j2;
		
		// perform traceback to count indels
		if (j1 < l2) {
			N2[i]++;
		} else if (i1 < l1) {
			N1[i]++;
		}
		while (i2 >= 0 && j2 >= 0) {
			if (o[i2 + index[j2]] == 0) {
				i2--;
				j2--;
			} else if (o[i2 + index[j2]] > 0) {
				N2[i]++;
				j2 -= o[i2 + index[j2]];
			} else {
				N1[i]++;
				i2 += o[i2 + index[j2]];
			}
		}
		
		if (N1[i] > 0 || N2[i] > 0) {
			i2 = i1;
			j2 = j1;
			
			if (N1[i] > 0) {
				p1 = (int *) malloc(N1[i]*sizeof(int)); // thread-safe on Windows
				p2 = (int *) malloc(N1[i]*sizeof(int)); // thread-safe on Windows
				N1[i] = 0; // reset count
			}
			if (N2[i] > 0) {
				p3 = (int *) malloc(N2[i]*sizeof(int)); // thread-safe on Windows
				p4 = (int *) malloc(N2[i]*sizeof(int)); // thread-safe on Windows
				N2[i] = 0; // reset count
			}
			
			// perform traceback
			if (j1 < l2) {
				p3[0] = l1;
				p4[0] = l2 - j1;
				N2[i]++;
			} else if (i1 < l1) {
				p1[0] = l2;
				p2[0] = l1 - i1;
				N1[i]++;
			}
			while (i2 >= 0 && j2 >= 0) {
				if (o[i2 + index[j2]] == 0) {
					i2--;
					j2--;
				} else if (o[i2 + index[j2]] > 0) {
					p3[N2[i]] = i2;
					p4[N2[i]] = o[i2 + index[j2]];
					N2[i]++;
					j2 -= o[i2 + index[j2]];
				} else {
					p1[N1[i]] = j2;
					p2[N1[i]] = -1*o[i2 + index[j2]];
					N1[i]++;
					i2 += o[i2 + index[j2]];
				}
			}
			
			if (N1[i] > 0) {
				P1[i] = p1;
				P2[i] = p2;
			}
			if (N2[i] > 0) {
				P3[i] = p3;
				P4[i] = p4;
			}
		}
		
		free(o);
		free(index);
	}
	free(square);
	
	int n1 = 0, n2 = 0;
	for (i = 0; i < n; i++) {
		n1 += N1[i];
		n2 += N2[i];
	}
	
	int *res1 = (int *) malloc(n1*sizeof(int)); // thread-safe on Windows
	int *res2 = (int *) malloc(n1*sizeof(int)); // thread-safe on Windows
	int *res3 = (int *) malloc(n2*sizeof(int)); // thread-safe on Windows
	int *res4 = (int *) malloc(n2*sizeof(int)); // thread-safe on Windows
	
	j1 = 0;
	j2 = 0;
	for (i = 0; i < n; i++) {
		if (N1[i] > 0) {
			p1 = P1[i];
			p2 = P2[i];
			
			for (i1 = N1[i] - 1; i1 >= 0; i1--) {
				res1[j1] = p1[i1] + S2[i];
				res2[j1] = p2[i1];
				j1++;
			}
			
			free(p1);
			free(p2);
		}
		
		if (N2[i] > 0) {
			p3 = P3[i];
			p4 = P4[i];
			
			for (i2 = N2[i] - 1; i2 >= 0; i2--) {
				res3[j2] = p3[i2] + S1[i];
				res4[j2] = p4[i2];
				j2++;
			}
			
			free(p3);
			free(p4);
		}
	}
	free(P1);
	free(P2);
	free(P3);
	free(P4);
	free(N1);
	free(N2);
	
	// insert gaps
	SEXP ans_width, ans;
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the widths of the aligned (equal width) XStringSet
	int sum = 0;
	for (i = 0; i < n2; i++)
		sum += res4[i];
	PROTECT(ans_width = NEW_INTEGER(2));
	int *width = INTEGER(ans_width);
	width[0] = X.length + sum;
	width[1] = width[0]; // same length after alignment
	
	// set the class of the XStringSet
	char ans_classname[40];
	if (t == 1) {
		strcpy(ans_classname, "DNAStringSet");
	} else if (t == 2) {
		strcpy(ans_classname, "RNAStringSet");
	} else { // t == 3
		strcpy(ans_classname, "AAStringSet");
	}
	
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	Chars_holder ans_elt_holder;
	
	// insert gaps in sequence X
	ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, 0);
	sum = 0; // position in ans_elt_holder.ptr
	int start = 0; // position in X
	for (i = 0; i < n2; i++) {
		if ((res3[i] - 1) > start) { // copy over sequence
			memcpy((char *) ans_elt_holder.ptr + sum, X.ptr + start, (res3[i] - 1 - start) * sizeof(char));
			sum += (res3[i] - 1 - start);
			start += (res3[i] - 1 - start);
		}
		if (res4[i] > 0) { // insert gaps
			if (t == 3) { // AAStringSet
				memset((char *) ans_elt_holder.ptr + sum, 45, res4[i] * sizeof(char));
			} else { // DNAStringSet or RNAStringSet
				memset((char *) ans_elt_holder.ptr + sum, 16, res4[i] * sizeof(char));
			}
			sum += res4[i];
		}
	}
	if (sum < ans_elt_holder.length) {
		memcpy((char *) ans_elt_holder.ptr + sum, X.ptr + start, (ans_elt_holder.length - sum) * sizeof(char));
	}
	
	// insert gaps in sequence Y
	ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, 1);
	sum = 0; // position in ans_elt_holder.ptr
	start = 0; // position in Y
	for (i = 0; i < n1; i++) {
		if ((res1[i] - 1) > start) { // copy over sequence
			memcpy((char *) ans_elt_holder.ptr + sum, Y.ptr + start, (res1[i] - 1 - start) * sizeof(char));
			sum += (res1[i] - 1 - start);
			start += (res1[i] - 1 - start);
		}
		if (res2[i] > 0) { // insert gaps
			if (t == 3) { // AAStringSet
				memset((char *) ans_elt_holder.ptr + sum, 45, res2[i] * sizeof(char));
			} else { // DNAStringSet or RNAStringSet
				memset((char *) ans_elt_holder.ptr + sum, 16, res2[i] * sizeof(char));
			}
			sum += res2[i];
		}
	}
	if (sum < ans_elt_holder.length) {
		memcpy((char *) ans_elt_holder.ptr + sum, Y.ptr + start, (ans_elt_holder.length - sum) * sizeof(char));
	}
	
	free(res1);
	free(res2);
	free(res3);
	free(res4);
	
	UNPROTECT(2);
	
	return ans;
}

static int alignRegion(const Chars_holder *s1, const Chars_holder *s2, int pos1, int pos2, int l1, int l2, int anchor, int bandWidth, double GO, double GE, double dropScore, double *subMatrix, int *lkup_row, int *lkup_col, int *results, int **indels1, int **lengths1, int **indels2, int **lengths2)
{
	// s1: pointer to pattern sequence
	// s2: pointer to subject sequence
	// pos1: start position in pattern sequence
	// pos2: start position in subject sequence
	// l1: number of positions to align in pattern
	// l2: number of positions to align in subject
	// anchor: 0 = global without terminal gap penalties; +/- 1 = local (+ = anchored at left side, - = anchored at right side); 2 = global with terminal gap penalties
	// bandWidth: width of adaptive band around max score in alignment matrix
	// GO: gap opening penalty
	// GE: (affine) gap extension penalty
	// dropScore: discontinue extension when score decreases by dropScore (applicable when anchor = +/- 1)
	// subMatrix: (square) substitution matrix
	// lkup_row: vector to convert characters to row indices in substitution matrix
	// lkup_col: vector to convert characters to column indices in substitution matrix
	
	// initialize variables
	int i, c1, c2, p1, p2, count, max_count, temp;
	double score, subScore;
	int I, J; // value at a sequence position
	const char *p, *s; // pattern and subject pointers
	p = s1->ptr - 1; // pointer to initial position in pattern
	s = s2->ptr - 1; // pointer to initial position in subject
	
	// determine allowable bandWidth
	if (l1 - pos1 + 1 < bandWidth)
		bandWidth = l1 - pos1 + 1;
	if (l2 - pos2 + 1 < bandWidth)
		bandWidth = l2 - pos2 + 1;
	if (bandWidth < 4)
		bandWidth = 4; // minimum bandWidth
	if (bandWidth % 2 == 1)
		bandWidth++; // bandWidth must be an even number
	
	// initialize size variables
	int half = bandWidth/2; // bandWidth is even so half stays even
	int tot = l1 - pos1 + l2 - pos2 + 1; // maximum number of diagonals
	int diags; // columns in initial alignment matrix
	int size = bandWidth + 2; // number of rows in o
	if (anchor == 0 || anchor == 2) {
		diags = tot; // initialize to maximum number of possible diagonals
	} else {
		diags = bandWidth; // anticipate at least bandWidth diagonals
	}
	
	// initialize score matrix [(bandWidth rows) x (3 columns)]
	double *m = (double *) malloc(bandWidth*3*sizeof(double)); // thread-safe on Windows
	for (i = 1; i < bandWidth*3; i++)
		m[i] = R_NegInf; // initialize to -Inf
	if (anchor == -1) {
		I = lkup_row[(unsigned char)p[l1]]; // last position in pattern
		if (I == NA_INTEGER) {
			free(m);
			return l1; // signal error
		}
		J = lkup_col[(unsigned char)s[l2]]; // last position in subject
		if (J == NA_INTEGER) {
			free(m);
			return -1*l2; // signal error
		}
	} else {
		I = lkup_row[(unsigned char)p[pos1]]; // first position in pattern
		if (I == NA_INTEGER) {
			free(m);
			return pos1; // signal error
		}
		J = lkup_col[(unsigned char)s[pos2]]; // first position in subject
		if (J == NA_INTEGER) {
			free(m);
			return -1*pos2; // signal error
		}
	}
	m[0] = subMatrix[I + J];
	
	// initialize traceback matrix [(bandWidth + 2 rows) x (diags columns)]
	int *o = (int *) malloc(size*diags*sizeof(int)); // thread-safe on Windows
	
	// initialize variables for looping
	int m1 = pos1; // position of max in first sequence
	int m2 = pos2; // position of max in second sequence
	int M1 = m1; // start of traceback in first sequence
	int M2 = m2; // start of traceback in second sequence
	int C1 = 0; // start of traceback in max diagonal
	int C2 = 0; // start of traceback in diagonal
	double max_score = R_NegInf; // highest observed score (initialize to -Inf)
	double diag_max_score = R_NegInf; // highest score on diagonal (initialize to -Inf)
	int diagonal = -1; // current diagonal
	int col0 = 0; // ultimate score column
	int col1 = 2; // penultimate score column
	int col2 = 1; // antipenultimate score column
	
	// loop through each diagonal
	while (m1 <= l1 && // within pattern
		m2 <= l2 && // within subject
		(anchor == 0 || // global unanchored
		anchor == 2 || // global anchored
		diag_max_score >= max_score + dropScore)) { // extending locally
		diagonal++; // next diagonal
		if (diagonal == diags) { // need to expand traceback matrix (o)
			diags *= 2; // double the number of possible diagonals
			if (diags > tot)
				diags = tot; // use maximum possible diagonals
			o = (int *) realloc(o, size*diags*sizeof(int)); // thread-safe on Windows
		}
		
		// initialize location
		c1 = m1 - half; // current position in first sequence
		c2 = m2 + half; // current position in second sequence
		if (c1 < pos1) { // out of bounds in pattern
			c2 = c2 - pos1 + c1;
			c1 = pos1;
		}
		if (c2 > l2) { // out of bounds in subject
			c1 = c1 - l2 + c2;
			c2 = l2;
		}
		o[bandWidth + size*diagonal] = c1; // record starting position
		
		// calculate relative positions in previous diagonals
		if (diagonal > 0) {
			p1 = c1 - o[bandWidth + size*(diagonal - 1)]; // one diagonal ago
			if (diagonal > 1) {
				p2 = c1 - o[bandWidth + size*(diagonal - 2)]; // two diagonals ago
			} else {
				p2 = -1*bandWidth - 1;
			}
		} else { // set positions out of bounds
			p1 = -1*bandWidth - 1;
			p2 = -1*bandWidth - 1;
		}
		
		// fill diagonal
		m1 = c1;
		m2 = c2;
		count = -1;
		max_count = 0;
		while (c1 <= l1 && // within bounds of pattern
			c2 >= pos2 && // within bounds of subject
			count < bandWidth - 1) { // within bandWidth
			count++; // next position
			o[count + size*diagonal] = 0; // initialize to across
			
			// compute score for aligning pattern and subject positions
			if (anchor == -1) {
				I = lkup_row[(unsigned char)p[l1 - c1 + 1]]; // position in pattern
				if (I == NA_INTEGER) {
					free(m);
					free(o);
					return l1 - c1 + 1; // signal error
				}
				J = lkup_col[(unsigned char)s[l2 - c2 + 1]]; // position in subject
				if (J == NA_INTEGER) {
					free(m);
					free(o);
					return -1*(l2 - c2 + 1); // signal error
				}
			} else {
				I = lkup_row[(unsigned char)p[c1]]; // first position in pattern
				if (I == NA_INTEGER) {
					free(m);
					free(o);
					return c1; // signal error
				}
				J = lkup_col[(unsigned char)s[c2]]; // first position in subject
				if (J == NA_INTEGER) {
					free(m);
					free(o);
					return -1*c2; // signal error
				}
			}
			subScore = subMatrix[I + J];
			
			// compute score for adding a gap in the subject
			if (p1 >= 1 && p1 <= o[bandWidth + 1 + size*(diagonal - 1)]) { // up
				if (anchor == 0 && c2 == pos2) { // at edge without terminal gap penalties
					m[count + bandWidth*col0] = subScore;
				} else {
					score = m[p1 - 1 + bandWidth*col1];
					if (o[p1 - 1 + size*(diagonal - 1)] < 1) {
						score += GO;
					} else { // existing gap in subject
						score += GE;
					}
					if (score > m[count + bandWidth*col0]) { // add gap in subject
						m[count + bandWidth*col0] = score;
						if (o[p1 - 1 + size*(diagonal - 1)] < 1) {
							o[count + size*diagonal] = 1;
						} else { // extend existing gap
							o[count + size*diagonal] = o[p1 - 1 + size*(diagonal - 1)] + 1;
						}
					}
				}
			}
			p1++;
			
			// compute score for adding a gap in the pattern
			if (p1 >= 1 && p1 <= o[bandWidth + 1 + size*(diagonal - 1)]) { // left
				if (anchor == 0 && c1 == pos1) { // at edge without terminal gap penalties
					m[count + bandWidth*col0] = subScore;
				} else {
					score = m[p1 - 1 + bandWidth*col1];
					if (o[p1 - 1 + size*(diagonal - 1)] > -1) {
						score += GO;
					} else { // existing gap in pattern
						score += GE;
					}
					if (score > m[count + bandWidth*col0]) { // add gap in pattern
						m[count + bandWidth*col0] = score;
						if (o[p1 - 1 + size*(diagonal - 1)] > -1) {
							o[count + size*diagonal] = -1;
						} else { // extend existing gap
							o[count + size*diagonal] = o[p1 - 1 + size*(diagonal - 1)] - 1;
						}
					}
				}
			}
			
			// compare with score for aligning positions
			if (p2 >= 1 && p2 <= o[bandWidth + 1 + size*(diagonal - 2)]) { // across
				subScore += m[p2 - 1 + bandWidth*col2];
				if (subScore > m[count + bandWidth*col0]) {
					m[count + bandWidth*col0] = subScore;
					o[count + size*diagonal] = 0;
				}
			}
			p2++;
			
			if (m[count + bandWidth*col0] >= m[max_count + bandWidth*col0]) {
				m1 = c1;
				m2 = c2;
				max_count = count;
				if (anchor == 1 || anchor == -1 || c1 == l1 || c2 == l2) {
					// check for new maximum score
					if (anchor < 2) {
						diag_max_score = 0;
					} else if (c1 == l1) {
						diag_max_score = l2 - c2;
					} else {
						diag_max_score = l1 - c1;
					}
					if (diag_max_score > 1) {
						diag_max_score = GO + GE*(diag_max_score - 1);
					} else if (diag_max_score > 0) {
						diag_max_score = GO;
					}
					diag_max_score += m[max_count + bandWidth*col0];
					if (diag_max_score >= max_score) {
						max_score = diag_max_score;
						M1 = c1;
						M2 = c2;
						C1 = count;
						C2 = diagonal;
					}
				} // else inherit previous diag_max_score
			}
			c1++;
			c2--;
		}
		o[bandWidth + 1 + size*diagonal] = count + 1; // record diagonal length
		
		// rotate columns in the score matrix
		temp = col0;
		col0 = col2;
		col2 = col1;
		col1 = temp;
		for (i = 0; i < bandWidth; i++)
			m[i + bandWidth*col0] = R_NegInf; // reset scores
		
		// shift max to next diagonal
		if (m2 == l2) {
			m1++;
		} else {
			m2++;
		}
	}
	free(m);
	
	// perform traceback
	m1 = M1; // position of max in sequence 1
	m2 = M2; // position of max in sequence 2
	c1 = C1 + 1; // position of max in diagonal
	c2 = C2 + 1; // max diagonal
	int count1 = 0; // number of indels
	int max_count1 = 1; // max number of indels (>= 1)
	int count2 = 0; // number of indels
	int max_count2 = 1; // max number of indels (>= 1)
	*indels1 = (int *) malloc(max_count1*sizeof(int)); // thread-safe on Windows
	*lengths1 = (int *) malloc(max_count1*sizeof(int)); // thread-safe on Windows
	*indels2 = (int *) malloc(max_count2*sizeof(int)); // thread-safe on Windows
	*lengths2 = (int *) malloc(max_count2*sizeof(int)); // thread-safe on Windows
	int *ind1 = indels1[0];
	int *len1 = lengths1[0];
	int *ind2 = indels2[0];
	int *len2 = lengths2[0];
	int end1, end2;
	if (anchor == 0 || anchor == 2) {
		end1 = l1;
		end2 = l2;
		if (m1 < l1) {
			ind2[count2] = l2 - pos2 + 2; // relative to pos2
			len2[count2++] = l1 - m1;
		}
		if (m2 < l2) {
			ind1[count1] = l1 - pos1 + 2; // relative to pos1
			len1[count1++] = l2 - m2;
		}
	} else { // only one anchor
		end1 = m1;
		end2 = m2;
	}
	while (m1 >= pos1 && m2 >= pos2) {
		if (o[c1 - 1 + size*(c2 - 1)] == 0) { // across
			c1 += o[bandWidth + size*(c2 - 1)];
			c2 -= 2; // shift two diagonals left
			if (c2 >= 1) // still within bounds
				c1 -= o[bandWidth + size*(c2 - 1)] + 1;
			m1--;
			m2--;
		} else if (o[c1 - 1 + size*(c2 - 1)] > 0) { // up
			if (count2 == max_count2) {
				max_count2 *= 2;
				*indels2 = (int *) realloc(*indels2, max_count2*sizeof(int)); // thread-safe on Windows
				*lengths2 = (int *) realloc(*lengths2, max_count2*sizeof(int)); // thread-safe on Windows
				ind2 = indels2[0];
				len2 = lengths2[0];
			}
			ind2[count2] = m2 + 2 - pos2;
			len2[count2++] = o[c1 - 1 + size*(c2 - 1)];
			m1 -= o[c1 - 1 + size*(c2 - 1)];
			temp = o[bandWidth + size*(c2 - 1)] + c1 - o[c1 - 1 + size*(c2 - 1)];
			c2 -= o[c1 - 1 + size*(c2 - 1)];
			c1 = temp;
			if (c2 >= 1) // still within bounds
				c1 -= o[bandWidth + size*(c2 - 1)];
		} else { // left
			if (count1 == max_count1) {
				max_count1 *= 2;
				*indels1 = (int *) realloc(*indels1, max_count1*sizeof(int)); // thread-safe on Windows
				*lengths1 = (int *) realloc(*lengths1, max_count1*sizeof(int)); // thread-safe on Windows
				ind1 = indels1[0];
				len1 = lengths1[0];
			}
			ind1[count1] = m1 + 2 - pos1;
			len1[count1++] = -1*o[c1 - 1 + size*(c2 - 1)];
			m2 += o[c1 - 1 + size*(c2 - 1)];
			temp = o[bandWidth + size*(c2 - 1)] + c1;
			c2 += o[c1 - 1 + size*(c2 - 1)];
			c1 = temp;
			if (c2 >= 1) // still within bounds
				c1 -= o[bandWidth + size*(c2 - 1)];
		}
	}
	free(o);
	if (m1 >= pos1) {
		if (count2 == max_count2) {
			max_count2++;
			*indels2 = (int *) realloc(*indels2, max_count2*sizeof(int)); // thread-safe on Windows
			*lengths2 = (int *) realloc(*lengths2, max_count2*sizeof(int)); // thread-safe on Windows
			ind2 = indels2[0];
			len2 = lengths2[0];
		}
		ind2[count2] = 1; // relative to pos2
		len2[count2++] = m1 - pos1 + 1;
	}
	if (m2 >= pos2) {
		if (count1 == max_count1) {
			max_count1++;
			*indels1 = (int *) realloc(*indels1, max_count1*sizeof(int)); // thread-safe on Windows
			*lengths1 = (int *) realloc(*lengths1, max_count1*sizeof(int)); // thread-safe on Windows
			ind1 = indels1[0];
			len1 = lengths1[0];
		}
		ind1[count1] = 1; // relative to pos1
		len1[count1++] = m2 - pos2 + 1;
	}
	if (anchor == -1) { // make positions relative to end
		for (i = 0; i < count1; i++)
			ind1[i] = end1 - pos1 - ind1[i] + 3;
		for (i = 0; i < count2; i++)
			ind2[i] = end2 - pos2 - ind2[i] + 3;
		pos1 = l1 - end1 + 1;
		pos2 = l2 - end2 + 1;
		end1 = l1;
		end2 = l2;
	}
	results[0] = count1;
	results[1] = count2;
	results[2] = pos1;
	results[3] = end1;
	results[4] = pos2;
	results[5] = end2;
	
	return 0; // signal success
}

SEXP alignPairs(SEXP pattern, SEXP subject, SEXP query, SEXP target, SEXP position, SEXP bandWidth, SEXP gapOpening, SEXP gapExtension, SEXP dropScore, SEXP subMatrix, SEXP matchMatrix, SEXP letters, SEXP verbose, SEXP pBar, SEXP nThreads)
{
	int i;
	int *q = INTEGER(query);
	int *t = INTEGER(target);
	int l = length(query);
	int bW = asInteger(bandWidth);
	double GO = asReal(gapOpening);
	double GE = asReal(gapExtension);
	double dS = asReal(dropScore);
	double *sM = REAL(subMatrix);
	int *mM = INTEGER(matchMatrix);
	int nthreads = asInteger(nThreads);
	
	int before, v, *rPercentComplete;
	double complete, soFar;
	SEXP percentComplete, utilsPackage;
	v = asLogical(verbose);
	if (v) { // percent complete variables
		complete = 0; // completed iterations
		soFar = 0;
		before = 0;
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	// set up a timer to minimize interrupt checks
	time_t start, end;
	double elapsed;
	time(&start);
	
	// build a vector of thread-safe pointers
	int **ptrs = (int **) malloc(l*sizeof(int *)); // thread-safe on Windows
	int *tot = (int *) malloc((l + 1)*sizeof(int)); // thread-safe on Windows
	tot[0] = 0;
	for (i = 0; i < l; i++) {
		ptrs[i] = INTEGER(VECTOR_ELT(position, i)); // anchor matrix
		tot[i + 1] = tot[i] + length(VECTOR_ELT(position, i))/4 + 1;
	}
	
	XStringSet_holder p_set, s_set, l_set;
	p_set = hold_XStringSet(pattern);
	s_set = hold_XStringSet(subject);
	l_set = hold_XStringSet(letters);
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
	
	// built output vectors for alignment results
	int **res1 = (int **) malloc(tot[l]*sizeof(int *)); // thread-safe on Windows
	int **res2 = (int **) malloc(tot[l]*sizeof(int *)); // thread-safe on Windows
	int **res3 = (int **) malloc(tot[l]*sizeof(int *)); // thread-safe on Windows
	int **res4 = (int **) malloc(tot[l]*sizeof(int *)); // thread-safe on Windows
	int **res5 = (int **) malloc(tot[l]*sizeof(int *)); // thread-safe on Windows
	
	// initialize empty results
	for (i = 0; i < tot[l]; i++) {
		int *res = (int *) malloc(6*sizeof(int)); // thread-safe on Windows
		res[0] = -1; // no memory allocated for indels in pattern
		res[1] = -1; // no memory allocated for indels in subject
		res1[i] = res;
	}
	
	// initialize output vectors
	SEXP ans1, ans2, ans3, ans4;
	PROTECT(ans1 = allocVector(INTSXP, l));
	int *starts1 = INTEGER(ans1);
	PROTECT(ans2 = allocVector(INTSXP, l));
	int *ends1 = INTEGER(ans2);
	PROTECT(ans3 = allocVector(INTSXP, l));
	int *starts2 = INTEGER(ans3);
	PROTECT(ans4 = allocVector(INTSXP, l));
	int *ends2 = INTEGER(ans4);
	
	int abort[3] = {0, 0, 0};
	#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(dynamic) num_threads(nthreads)
	#endif
	for (i = 0; i < l; i++) {
		if (abort[0] == 0) {
			int index1 = q[i];
			int index2 = t[i];
			int *anchor = ptrs[i];
			Chars_holder p_i = get_elt_from_XStringSet_holder(&p_set, index1 - 1);
			Chars_holder s_i = get_elt_from_XStringSet_holder(&s_set, index2 - 1);
			int N = tot[i + 1] - tot[i] - 1; // number of anchors
			int n = tot[i];
			int signal; // completion signal (success == 0)
			int *res, p1, p2, q1, q2, k;
			if (N == 0) { // no anchor positions
				starts1[i] = 1;
				starts2[i] = 1;
				ends1[i] = p_i.length;
				ends2[i] = s_i.length;
				signal = alignRegion(&p_i, &s_i, 1, 1, p_i.length, s_i.length, 0, bW, GO, GE, dS, sM, lkup_row, lkup_col, res1[n], &res2[n], &res3[n], &res4[n], &res5[n]);
				if (signal != 0) {
					abort[0] = 1; // unknown character flag
					abort[1] = signal > 0 ? index1 : index2; // sequence flag
					abort[2] = signal; // letter number
					continue;
				}
			} else {
				p1 = anchor[0] - 1; // right bound in pattern
				p2 = anchor[2] - 1; // right bound in subject
				if (p1 > 0 && p2 > 0) {
					if (p1 < 1 || p1 > p_i.length) {
						abort[0] = 3; // out-of-bounds anchors flag
						abort[1] = i + 1; // sequence flag
						abort[2] = 1; // anchor number
						continue;
					}
					if (p2 < 1 || p2 > s_i.length) {
						abort[0] = 3; // out-of-bounds anchors flag
						abort[1] = -1 - i; // sequence flag
						abort[2] = 1; // anchor number
						continue;
					}
					signal = alignRegion(&p_i, &s_i, 1, 1, p1, p2, -1, bW, GO, GE, dS, sM, lkup_row, lkup_col, res1[n], &res2[n], &res3[n], &res4[n], &res5[n]);
					if (signal != 0) {
						abort[0] = 1; // unknown character flag
						abort[1] = signal > 0 ? index1 : index2; // sequence flag
						abort[2] = signal; // letter number
						continue;
					}
					res = res1[n];
					starts1[i] = res[2];
					starts2[i] = res[4];
				} else if (p1 >= 0 && p2 >= 0) {
					starts1[i] = anchor[0];
					starts2[i] = anchor[2];
				} else { // virtual anchor out of bounds
					starts1[i] = 1;
					starts2[i] = 1;
				}
				n++;
				for (int j = 1; j < N; j++) {
					p1 = anchor[1 + 4*(j - 1)];
					q1 = anchor[4*j];
					if (p1 > q1) {
						abort[0] = 2; // overlapping anchor flag
						abort[1] = i + 1; // anchor flag
						abort[2] = j; // anchor number
						continue;
					}
					p2 = anchor[3 + 4*(j - 1)];
					q2 = anchor[2 + 4*j];
					if (p2 > q2) {
						abort[0] = 2; // overlapping anchor flag
						abort[1] = i + 1; // anchor flag
						abort[2] = -1*j; // anchor number
						continue;
					}
					p1++;
					p2++;
					q1--;
					q2--;
					if (q1 < p1) {
						if (q2 >= p2) {
							res = res1[n];
							res[0] = 1;
							res = (int *) malloc(1*sizeof(int)); // thread-safe on Windows
							res[0] = p1 - starts1[i] + 1;
							res2[n] = res;
							res = (int *) malloc(1*sizeof(int)); // thread-safe on Windows
							res[0] = q2 - p2 + 1;
							res3[n] = res;
						}
					} else if (q2 < p2) {
						if (q1 >= p1) {
							res = res1[n];
							res[1] = 1;
							res = (int *) malloc(1*sizeof(int)); // thread-safe on Windows
							res[0] = p2 - starts2[i] + 1;
							res4[n] = res;
							res = (int *) malloc(1*sizeof(int)); // thread-safe on Windows
							res[0] = q1 - p1 + 1;
							res5[n] = res;
						}
					} else {
						if (p1 < 1 || q1 > p_i.length) {
							abort[0] = 3; // out-of-bounds anchor flag
							abort[1] = i + 1; // sequence flag
							abort[2] = p1 < 1 ? j : j + 1; // anchor number
							continue;
						}
						if (p2 < 1 || q2 > s_i.length) {
							abort[0] = 3; // out-of-bounds anchor flag
							abort[1] = -1 - i; // sequence flag
							abort[2] = p2 < 1 ? j : j + 1; // anchor number
							continue;
						}
						signal = alignRegion(&p_i, &s_i, p1, p2, q1, q2, 2, bW, GO, GE, dS, sM, lkup_row, lkup_col, res1[n], &res2[n], &res3[n], &res4[n], &res5[n]);
						if (signal != 0) {
							abort[0] = 1; // unknown character flag
							abort[1] = signal > 0 ? index1 : index2; // sequence flag
							abort[2] = signal; // letter number
							continue;
						}
						res = res1[n];
						k = res[0];
						res = res2[n];
						while (k > 0) {
							k--;
							res[k] += p1 - starts1[i];
						}
						res = res1[n];
						k = res[1];
						res = res4[n];
						while (k > 0) {
							k--;
							res[k] += p2 - starts2[i];
						}
					}
					n++;
				}
				p1 = anchor[1 + 4*(N - 1)] + 1; // left bound in pattern
				p2 = anchor[3 + 4*(N - 1)] + 1; // left bound in subject
				q1 = p_i.length;
				q2 = s_i.length;
				if (p1 <= q1 && p2 <= q2) {
					if (p1 < 1 || q1 > p_i.length) {
						abort[0] = 3; // out-of-bounds anchor flag
						abort[1] = i + 1; // sequence flag
						abort[2] = N; // anchor number
						continue;
					}
					if (p2 < 1 || q2 > s_i.length) {
						abort[0] = 3; // out-of-bounds anchor flag
						abort[1] = -1 - i; // sequence flag
						abort[2] = N; // anchor number
						continue;
					}
					signal = alignRegion(&p_i, &s_i, p1, p2, q1, q2, 1, bW, GO, GE, dS, sM, lkup_row, lkup_col, res1[n], &res2[n], &res3[n], &res4[n], &res5[n]);
					if (signal != 0) {
						abort[0] = 1; // unknown character flag
						abort[1] = signal > 0 ? index1 : index2; // sequence flag
						abort[2] = signal; // letter number
						continue;
					}
					res = res1[n];
					ends1[i] = res[3];
					ends2[i] = res[5];
					res = res1[n];
					k = res[0];
					res = res2[n];
					while (k > 0) {
						k--;
						res[k] += p1 - starts1[i];
					}
					res = res1[n];
					k = res[1];
					res = res4[n];
					while (k > 0) {
						k--;
						res[k] += p2 - starts2[i];
					}
				} else if (anchor[1 + 4*(N - 1)] <= q1 && anchor[3 + 4*(N - 1)] <= q2)  {
					ends1[i] = anchor[1 + 4*(N - 1)];
					ends2[i] = anchor[3 + 4*(N - 1)];
				} else { // virtual anchor out of bounds
					ends1[i] = q1;
					ends2[i] = q2;
				}
			}
			
			if (v) {
				#ifdef _OPENMP
				#pragma omp critical
				{
					complete++;
				}
				#else
				complete++;
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
					if (checkInterrupt() != 0)
						abort[0] = -1;
					if (abort[0] == 0 && v) {
						soFar = complete/l;
						*rPercentComplete = floor(100*soFar);
						if (*rPercentComplete > before) { // when the percent has changed
							// tell the progress bar to update in the R console
							eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
							before = *rPercentComplete;
						}
					}
				}
			}
		}
	}
	free(ptrs);
	
	if (abort[0] != 0) {
		// release memory
		for (i = 0; i < tot[l]; i++) {
			int *res = res1[i];
			if (res[0] >= 0) {
				free(res2[i]);
				free(res3[i]);
			}
			if (res[1] >= 0) {
				free(res4[i]);
				free(res5[i]);
			}
			free(res);
		}
		free(res1);
		free(res2);
		free(res3);
		free(res4);
		free(res5);
		free(tot);
		if (v) {
			UNPROTECT(6);
		} else {
			UNPROTECT(4);
		}
		if (abort[0] < 0) {
			error("Received user interrupt.");
		} else if (abort[0] == 1) {
			if (abort[2] > 0) {
				error("Unexpected character in pattern[%d] position %d.", abort[1], abort[2]);
			} else {
				error("Unexpected character in subject[%d] position %d.", abort[1], -1*abort[2]);
			}
		} else if (abort[0] == 2) {
			if (abort[2] > 0) {
				error("Overlapping pattern anchor in Position[%d] column %d.", abort[1], abort[2]);
			} else {
				error("Overlapping subject anchor in Position[%d] column %d.", abort[1], -1*abort[2]);
			}
		} else if (abort[0] == 3) {
			if (abort[1] > 0) {
				error("Out-of-bounds pattern anchor in Position[%d] column %d.", abort[1], abort[2]);
			} else {
				error("Out-of-bounds subject anchor in Position[%d] column %d.", -1*abort[1], abort[2]);
			}
		} else {
			error("Unknown error.");
		}
	}
	
	// initialize output vectors
	SEXP ans5, ans6, ans7, ans8, ans9, ans10, ans11, ans12;
	PROTECT(ans5 = allocVector(INTSXP, l));
	int *matches = INTEGER(ans5);
	PROTECT(ans6 = allocVector(INTSXP, l));
	int *mismatches = INTEGER(ans6);
	PROTECT(ans7 = allocVector(INTSXP, l));
	int *counts = INTEGER(ans7);
	PROTECT(ans8 = allocVector(REALSXP, l));
	double *scores = REAL(ans8);
	PROTECT(ans9 = allocVector(VECSXP, l));
	PROTECT(ans10 = allocVector(VECSXP, l));
	PROTECT(ans11 = allocVector(VECSXP, l));
	PROTECT(ans12 = allocVector(VECSXP, l));
	
	int j, k, p, p1, p2, c1, c2, ms, mms, count, count1, count2, signal;
	int *res, *ind1, *len1, *ind2, *len2;
	double score;
	SEXP indels1, lengths1, indels2, lengths2;
	unsigned char s1, s2; // sequence positions
	for (i = 0; i < l; i++) {
		count1 = 0;
		count2 = 0;
		for (j = tot[i]; j < tot[i + 1]; j++) {
			res = res1[j];
			if (res[0] > 0)
				count1 += res[0];
			if (res[1] > 0)
				count2 += res[1];
		}
		PROTECT(indels1 = allocVector(INTSXP, count1));
		ind1 = INTEGER(indels1);
		PROTECT(lengths1 = allocVector(INTSXP, count1));
		len1 = INTEGER(lengths1);
		PROTECT(indels2 = allocVector(INTSXP, count2));
		ind2 = INTEGER(indels2);
		PROTECT(lengths2 = allocVector(INTSXP, count2));
		len2 = INTEGER(lengths2);
		count1 = 0;
		count2 = 0;
		int N = tot[i + 1] - tot[i] - 1; // number of anchors
		for (j = tot[i]; j < tot[i + 1]; j++) {
			res = res1[j];
			p = res[0];
			if (j == tot[i] && N > 0) { // record forwards
				for (k = 0; k < p; k++) {
					res = res2[j];
					ind1[count1] = res[k];
					res = res3[j];
					len1[count1++] = res[k];
				}
			} else { // record backwards
				for (k = p - 1; k >= 0; k--) {
					res = res2[j];
					ind1[count1] = res[k];
					res = res3[j];
					len1[count1++] = res[k];
				}
			}
			if (p >= 0) {
				res = res2[j];
				free(res);
				res = res3[j];
				free(res);
			}
			res = res1[j];
			p = res[1];
			if (j == tot[i] && N > 0) { // record forwards
				for (k = 0; k < p; k++) {
					res = res4[j];
					ind2[count2] = res[k];
					res = res5[j];
					len2[count2++] = res[k];
				}
			} else { // record backwards
				for (k = p - 1; k >= 0; k--) {
					res = res4[j];
					ind2[count2] = res[k];
					res = res5[j];
					len2[count2++] = res[k];
				}
			}
			if (p >= 0) {
				res = res4[j];
				free(res);
				res = res5[j];
				free(res);
			}
			res = res1[j];
			free(res);
		}
		
		// calculate matches, mismatches, score, and alignment length
		Chars_holder p_i = get_elt_from_XStringSet_holder(&p_set, q[i] - 1);
		Chars_holder s_i = get_elt_from_XStringSet_holder(&s_set, t[i] - 1);
		p1 = 1; // position in pattern
		p2 = 1; // position in subject
		c1 = 0; // position in indels1
		c2 = 0; // position in indels2
		count = 0; // alignment length
		score = 0; // alignment score
		mms = 0; // number of mismatches
		ms = 0; // number of matches
		signal = 0; // completion signal (success == 0)
		while (p1 + starts1[i] - 1 <= ends1[i] ||
			p2 + starts2[i] - 1 <= ends2[i]) {
			if (c1 < count1 && p1 == ind1[c1]) {
				count += len1[c1]; // add gap to count
				score += GO;
				if (len1[c1] > 1)
					score += GE*(len1[c1] - 1);
				p2 += len1[c1]; // skip subject positions
				c1++; // advance to next gap in pattern
			} else if (c2 < count2 && p2 == ind2[c2]) {
				count += len2[c2]; // add gap to count
				score += GO;
				if (len2[c2] > 1)
					score += GE*(len2[c2] - 1);
				p1 += len2[c2]; // skip pattern positions
				c2++; // advance to next gap in subject
			} else {
				count++; // increment count
				s1 = (unsigned char)p_i.ptr[p1 + starts1[i] - 2];
				s2 = (unsigned char)s_i.ptr[p2 + starts2[i] - 2];
				if (lkup_row[s1] == NA_INTEGER) {
					abort[0] = 1;
					abort[1] = q[i];
					abort[2] = p1 + starts1[i] - 1;
					break;
				}
				if (lkup_col[s2] == NA_INTEGER) {
					abort[0] = 1;
					abort[1] = t[i];
					abort[2] = 1 - p2 - starts2[i];
					break;
				}
				if (mM[lkup_row[s1] + lkup_col[s2]]) {
					ms++;
				} else {
					mms++;
				}
				score += sM[lkup_row[s1] + lkup_col[s2]];
				p1++;
				p2++;
			}
		}
		mismatches[i] = mms;
		matches[i] = ms;
		counts[i] = count;
		scores[i] = score;
		
		SET_VECTOR_ELT(ans9, i, indels1);
		SET_VECTOR_ELT(ans10, i, lengths1);
		SET_VECTOR_ELT(ans11, i, indels2);
		SET_VECTOR_ELT(ans12, i, lengths2);
		UNPROTECT(4);
	}
	free(res1);
	free(res2);
	free(res3);
	free(res4);
	free(res5);
	free(tot);
	free(lkup_row);
	free(lkup_col);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 12));
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	SET_VECTOR_ELT(ret_list, 2, ans3);
	SET_VECTOR_ELT(ret_list, 3, ans4);
	SET_VECTOR_ELT(ret_list, 4, ans5);
	SET_VECTOR_ELT(ret_list, 5, ans6);
	SET_VECTOR_ELT(ret_list, 6, ans7);
	SET_VECTOR_ELT(ret_list, 7, ans8);
	SET_VECTOR_ELT(ret_list, 8, ans9);
	SET_VECTOR_ELT(ret_list, 9, ans10);
	SET_VECTOR_ELT(ret_list, 10, ans11);
	SET_VECTOR_ELT(ret_list, 11, ans12);
	
	UNPROTECT(13);
	if (v)
		UNPROTECT(2);
	
	if (abort[0] != 0) {
		if (abort[2] > 0) {
			error("Unexpected character in pattern[%d] position %d.", abort[1], abort[2]);
		} else {
			error("Unexpected character in subject[%d] position %d.", abort[1], -1*abort[2]);
		}
	}
	
	return ret_list;
}

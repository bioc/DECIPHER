/****************************************************************************
 *                       Performs Pairwise Alignment                        *
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

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

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
	
	#pragma omp parallel for private(i,l1,l2,i1,j1,i2,j2,d,l,u,p1,p2,p3,p4) num_threads(nthreads)
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

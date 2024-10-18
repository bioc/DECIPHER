/****************************************************************************
 *                   Predicts Protein Secondary Structure                   *
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// for math functions
#include <math.h>

// DECIPHER header file
#include "DECIPHER.h"

SEXP predictHEC(SEXP x, SEXP windowSize, SEXP background, SEXP HEC_MI1, SEXP HEC_MI2, SEXP output, SEXP names)
{
	int i, j, k, l, k1, k2, p1, p2, c = 0;
	int wS1 = INTEGER(windowSize)[0];
	int wS2 = INTEGER(windowSize)[1];
	double f1 = (2*(double)wS1 - 1)/(2*(double)wS1 + 1); // fraction of each single
	double f2 = 2/(2*(double)wS2 + 1); // fraction of each double
	double *bg = REAL(background); // background scores
	int w = length(background);
	double *MI1 = REAL(HEC_MI1); // [AA, pos, state]
	double *MI2 = REAL(HEC_MI2); // [AA1, AA2, pos1, pos2, state]
	const char *nms = CHAR(STRING_ELT(names, 0)); // state names
	int N = length(STRING_ELT(names, 0)); // number of states
	int total1 = length(HEC_MI1)/(20*N); // maximum window
	int center1 = (total1 - 1)/2; // center of window
	int total2 = (int)sqrt(length(HEC_MI2)/(400*N)); // maximum window
	int center2 = (total2 - 1)/2; // center of window
	int o = asInteger(output);
	
	double HEC[N], temp, *rans;
	char *states;
	
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	int x_length = get_length_from_XStringSet_holder(&x_set);
	
	SEXP ret, ans;
	if (o == 1) { // return a character vector
		PROTECT(ret = allocVector(STRSXP, x_length));
	} else { // return a list of matrices (x_i.length x N)
		PROTECT(ret = allocVector(VECSXP, x_length));
	}
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		
		int *residues = R_Calloc(x_i.length, int);
		l = 0;
		for (j = 0; j < x_i.length; j++) {
			switch (x_i.ptr[j]) {
				case 65: // A
					residues[l] = 0;
					l++;
					break;
				case 82: // R
					residues[l] = 1;
					l++;
					break;
				case 78: // N
					residues[l] = 2;
					l++;
					break;
				case 68: // D
					residues[l] = 3;
					l++;
					break;
				case 67: // C
					residues[l] = 4;
					l++;
					break;
				case 81: // Q
					residues[l] = 5;
					l++;
					break;
				case 69: // E
					residues[l] = 6;
					l++;
					break;
				case 71: // G
					residues[l] = 7;
					l++;
					break;
				case 72: // H
					residues[l] = 8;
					l++;
					break;
				case 73: // I
					residues[l] = 9;
					l++;
					break;
				case 76: // L
					residues[l] = 10;
					l++;
					break;
				case 75: // K
					residues[l] = 11;
					l++;
					break;
				case 77: // M
					residues[l] = 12;
					l++;
					break;
				case 70: // F
					residues[l] = 13;
					l++;
					break;
				case 80: // P
					residues[l] = 14;
					l++;
					break;
				case 83: // S
					residues[l] = 15;
					l++;
					break;
				case 84: // T
					residues[l] = 16;
					l++;
					break;
				case 87: // W
					residues[l] = 17;
					l++;
					break;
				case 89: // Y
					residues[l] = 18;
					l++;
					break;
				case 86: // V
					residues[l] = 19;
					l++;
					break;
				case 85: // U
					residues[l] = 20; // unknown
					l++;
					break;
				case 79: // O
					residues[l] = 20; // unknown
					l++;
					break;
				case 66: // B = N or D
					residues[l] = 2; // treat as N
					l++;
					break;
				case 90: // Z = Q or E
					residues[l] = 5; // treat as Q
					l++;
					break;
				case 74: // J = I or L
					residues[l] = 9; // treat as I
					l++;
					break;
				case 88: // X = any letter
					residues[l] = 20; // unknown
					l++;
					break;
				case 42: // * (stop)
					residues[l] = 20; // unknown
					l++;
					break;
				case 45: // -
					break; // does not participate
				case 43: // +
					break; // does not participate
				case 46: // . treated as -
					break; // does not participate
				default:
					error("not AA!");
					break;
			}
		}
		
		if (o == 1) {
			states = R_Calloc(l + 1, char); // last position is for null terminating
		} else  {
			PROTECT(ans = allocMatrix(REALSXP, N, l)); // [state][pos]
			rans = REAL(ans);
		}
		
		for (j = 0; j < l; j++) {
			for (k = 0; k < N; k++)
				HEC[k] = *(bg + c++);
			if (c == w)
				c = 0; // recycle background
			
			for (k1 = -1*wS1; k1 <= wS1; k1++) {
				p1 = j + k1;
				if (p1 < 0)
					continue;
				if (p1 >= l)
					break;
				
				if (residues[p1] < 20) {
					// add mutual information from single residues
					for (k = 0; k < N; k++)
						HEC[k] -= f1 * *(MI1 + residues[p1] + 20*(center1 + k1) + 20*k*total1);
				}
			}
			
			for (k1 = -1*wS2; k1 <= wS2; k1++) {
				p1 = j + k1;
				if (p1 < 0)
					continue;
				if (p1 >= l)
					break;
				
				if (residues[p1] >= 20)
					continue;
				
				for (k2 = wS2; k2 > k1; k2--) {
					p2 = j + k2;
					if (p2 < 0)
						break;
					if (p2 >= l)
						continue;
					
					if (residues[p2] < 20) {
						// add mutual information from pairs of residues
						for (k = 0; k < N; k++)
							HEC[k] += f2 * *(MI2 + residues[p1] + 20*residues[p2] + 400*(center2 + k1) + 400*total2*(center2 + k2) + 400*k*total2*total2);
					}
				}
			}
			
			if (o == 1) { // states
				temp = HEC[0];
				states[j] = nms[0];
				for (k = 0; k < N; k++) {
					if (HEC[k] > temp) {
						temp = HEC[k];
						states[j] = nms[k];
					}
				}
			} else if (o == 2) { // log-odds scores
				for (k = 0; k < N; k++)
					*(rans + N*j + k) = HEC[k];
			} else { // normalized probabilities
				temp = 0;
				for (k = 0; k < N; k++) {
					HEC[k] = exp(HEC[k])/N;
					temp += HEC[k];
				}
				for (k = 0; k < N; k++)
					*(rans + N*j + k) = HEC[k]/temp;
			}
		}
		
		if (o == 1) {
			states[l] = '\0'; // end (null terminate) the string
			SET_STRING_ELT(ret, i, mkChar(states));
			R_Free(states);
		} else {
			SET_VECTOR_ELT(ret, i, ans);
			UNPROTECT(1); // ans
		}
		
		R_Free(residues);
	}
	
	UNPROTECT(1);
	
	return ret;
}

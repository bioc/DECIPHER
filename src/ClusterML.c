/****************************************************************************
 *                        Cluster Maximum Likelihood                        *
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

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

/* for Calloc/Free */
#include <R_ext/RS.h>

// for math functions
#include <math.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// for calloc/free
#include <stdlib.h>

// for floating point limits
#include <float.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

static void L_known(const char *p, double *Ls, int gaps)
{
	switch (*p) {
		case 1: // A
			*(Ls) = 1;
			break;
		case 2: // C
			*(Ls + 1) = 1;
			break;
		case 3: // M
			*(Ls) = 1/(double)2; *(Ls + 1) = 1/(double)2; // AC
			break;
		case 4: // G
			*(Ls + 2) = 1;
			break;
		case 5: // R
			*(Ls) = 1/(double)2; *(Ls + 2) = 1/(double)2; // AG
			break;
		case 6: // S
			*(Ls + 1) = 1/(double)2; *(Ls + 2) = 1/(double)2; // CG
			break;
		case 7: // V
			*(Ls) = 1/(double)3; *(Ls + 1) = 1/(double)3; *(Ls + 2) = 1/(double)3; // ACG
			break;
		case 8: // T
			*(Ls + 3) = 1;
			break;
		case 9: // W
			*(Ls) = 1/(double)2; *(Ls + 3) = 1/(double)2; // AT
			break;
		case 10: // Y
			*(Ls + 1) = 1/(double)2; *(Ls + 3) = 1/(double)2; // CT
			break;
		case 11: // H
			*(Ls) = 1/(double)3; *(Ls + 1) = 1/(double)3; *(Ls + 3) = 1/(double)3; // ACT
			break;
		case 12: // K
			*(Ls + 2) = 1/(double)2; *(Ls + 3) = 1/(double)2; // GT
			break;
		case 13: // D
			*(Ls) = 1/(double)3; *(Ls + 2) = 1/(double)3; *(Ls + 3) = 1/(double)3; // AGT
			break;
		case 14: // B
			*(Ls + 1) = 1/(double)3; *(Ls + 2) = 1/(double)3; *(Ls + 3) = 1/(double)3; // CGT
			break;
		case 15: // N
			// ACGT (unknown)
			break;
		case 16: // -
			if (gaps)
				*(Ls + 4) = 1;
			break;
		case 32: // +
			// mask
			break;
		case 64: // .
			if (gaps)
				*(Ls + 4) = 1;
			break;
		default:
			error("not nucleotides!");
			break;
	}
}

static void L_known_AA(const char *p, double *Ls, int gaps)
{
	switch (*p) {
		case 65: // A
			*(Ls) = 1;
			break;
		case 82: // R
			*(Ls + 1) = 1;
			break;
		case 78: // N
			*(Ls + 2) = 1;
			break;
		case 68: // D
			*(Ls + 3) = 1;
			break;
		case 67: // C
			*(Ls + 4) = 1;
			break;
		case 81: // Q
			*(Ls + 5) = 1;
			break;
		case 69: // E
			*(Ls + 6) = 1;
			break;
		case 71: // G
			*(Ls + 7) = 1;
			break;
		case 72: // H
			*(Ls + 8) = 1;
			break;
		case 73: // I
			*(Ls + 9) = 1;
			break;
		case 76: // L
			*(Ls + 10) = 1;
			break;
		case 75: // K
			*(Ls + 11) = 1;
			break;
		case 77: // M
			*(Ls + 12) = 1;
			break;
		case 70: // F
			*(Ls + 13) = 1;
			break;
		case 80: // P
			*(Ls + 14) = 1;
			break;
		case 83: // S
			*(Ls + 15) = 1;
			break;
		case 84: // T
			*(Ls + 16) = 1;
			break;
		case 87: // W
			*(Ls + 17) = 1;
			break;
		case 89: // Y
			*(Ls + 18) = 1;
			break;
		case 86: // V
			*(Ls + 19) = 1;
			break;
		case 85: // U
			// ignore
			break;
		case 79: // O
			// ignore
			break;
		case 66: // B = N or D
			*(Ls + 2) = 1/(double)2; *(Ls + 3) = 1/(double)2;
			break;
		case 90: // Z = Q or E
			*(Ls + 5) = 1/(double)2; *(Ls + 6) = 1/(double)2;
			break;
		case 74: // J = I or L
			*(Ls + 9) = 1/(double)2; *(Ls + 10) = 1/(double)2;
			break;
		case 88: // X = any letter
			// ignore
			break;
		case 42: // * (stop)
			// ignore
			break;
		case 45: // -
			if (gaps)
				*(Ls + 20) = 1;
			break;
		case 43: // +
			// mask
			break;
		case 46: // .
			if (gaps)
				*(Ls + 20) = 1;
			break;
		default:
			error("not AA!");
			break;
	}
}

static void L_unknown(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P, const int j1, const int j2, const int s2, const double epsilon, const double inv_epsilon, const int root)
{
	int i, j;
	const int s0 = 4;
	const int s1 = 5;
	double L1[s0], L2[s0];
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	const double *P1 = P + j1*s2;
	const double *P2 = P + j2*s2;
	
	int Z1 = 0, Z2 = 0;
	for (i = 0; i < s0; i++) {
		if (*(Ls1 + i) != 0)
			Z1 = 1;
		if (*(Ls2 + i) != 0)
			Z2 = 1;
		if (Z1 && Z2)
			break;
	}
	
	if (root == 0 && Z1) {
		if (Z2) {
			// neither branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				P1++;
				P2++;
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i]*L2[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			// second branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
			}
			for (j = 1; j < s0; j++) {
				P1++;
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		}
	} else {
		if (Z2) {
			// first branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				P2++;
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L2[i];
			}
			
			if (root && Z1) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= *(Ls1 + i);
				*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			} else {
				*(Ls3 + s1) = *(Ls2 + s1);
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon) {
					Z3 = 1;
					break;
				}
			}
			
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			for (i = 0; i < s0; i++)
				*(Ls3 + i) = *(Ls1 + i);
			*(Ls3 + s1) = *(Ls1 + s1);
		}
	}
	*(Ls3 + s0) = 0;
}

static void L_unknown_Indels(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P, const int j1, const int j2, const int s2, const double epsilon, const double inv_epsilon, const int root)
{
	int i, j;
	const int s0 = 5;
	const int s1 = 5;
	double L1[s0], L2[s0];
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	const double *P1 = P + j1*s2;
	const double *P2 = P + j2*s2;
	
	int Z1 = 0, Z2 = 0;
	for (i = 0; i < s0; i++) {
		if (*(Ls1 + i) != 0)
			Z1 = 1;
		if (*(Ls2 + i) != 0)
			Z2 = 1;
		if (Z1 && Z2)
			break;
	}
	
	if (root == 0 && Z1) {
		if (Z2) {
			// neither branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i]*L2[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			// second branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
			}
			for (j = 1; j < s0; j++) {
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		}
	} else {
		if (Z2) {
			// first branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L2[i];
			}
			
			if (root && Z1) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= *(Ls1 + i);
				*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			} else {
				*(Ls3 + s1) = *(Ls2 + s1);
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon) {
					Z3 = 1;
					break;
				}
			}
			
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			for (i = 0; i < s0; i++)
				*(Ls3 + i) = *(Ls1 + i);
			*(Ls3 + s1) = *(Ls1 + s1);
		}
	}
}

static void L_unknown_AA(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P, const int j1, const int j2, const int s2, const double epsilon, const double inv_epsilon, const int root)
{
	int i, j;
	const int s0 = 20;
	const int s1 = 21;
	double L1[s0], L2[s0];
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	const double *P1 = P + j1*s2;
	const double *P2 = P + j2*s2;
	
	int Z1 = 0, Z2 = 0;
	for (i = 0; i < s0; i++) {
		if (*(Ls1 + i) != 0)
			Z1 = 1;
		if (*(Ls2 + i) != 0)
			Z2 = 1;
		if (Z1 && Z2)
			break;
	}
	
	if (root == 0 && Z1) {
		if (Z2) {
			// neither branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				P1++;
				P2++;
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i]*L2[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			// second branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
			}
			for (j = 1; j < s0; j++) {
				P1++;
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		}
	} else {
		if (Z2) {
			// first branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				P2++;
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L2[i];
			}
			
			if (root && Z1) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= *(Ls1 + i);
				*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			} else {
				*(Ls3 + s1) = *(Ls2 + s1);
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon) {
					Z3 = 1;
					break;
				}
			}
			
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			for (i = 0; i < s0; i++)
				*(Ls3 + i) = *(Ls1 + i);
			*(Ls3 + s1) = *(Ls1 + s1);
		}
	}
	*(Ls3 + s0) = 0;
}

static void L_unknown_AA_Indels(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P, const int j1, const int j2, const int s2, const double epsilon, const double inv_epsilon, const int root)
{
	int i, j;
	const int s0 = 21;
	const int s1 = 21;
	double L1[s0], L2[s0];
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	const double *P1 = P + j1*s2;
	const double *P2 = P + j2*s2;
	
	int Z1 = 0, Z2 = 0;
	for (i = 0; i < s0; i++) {
		if (*(Ls1 + i) != 0)
			Z1 = 1;
		if (*(Ls2 + i) != 0)
			Z2 = 1;
		if (Z1 && Z2)
			break;
	}
	
	if (root == 0 && Z1) {
		if (Z2) {
			// neither branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i]*L2[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			// second branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L1[i] = *(P1++)*(*(Ls1));
			}
			for (j = 1; j < s0; j++) {
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L1[i] += *(P1++)*(*(Ls1 + j));
				}
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L1[i];
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon)
					Z3 = 1;
			}
			
			*(Ls3 + s1) = *(Ls1 + s1);
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		}
	} else {
		if (Z2) {
			// first branch can be disregarded
			
			#pragma omp simd
			for (i = 0; i < s0; i++) {
				L2[i] = *(P2++)*(*(Ls2));
			}
			for (j = 1; j < s0; j++) {
				#pragma omp simd
				for (i = 0; i < s0; i++) {
					L2[i] += *(P2++)*(*(Ls2 + j));
				}
			}
			
			for (i = 0; i < s0; i++) {
				*(Ls3 + i) = L2[i];
			}
			
			if (root && Z1) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= *(Ls1 + i);
				*(Ls3 + s1) = *(Ls1 + s1) + *(Ls2 + s1);
			} else {
				*(Ls3 + s1) = *(Ls2 + s1);
			}
			
			int Z3 = 0;
			for (i = 0; i < s0; i++) {
				if (*(Ls3 + i) > 0 && *(Ls3 + i) < inv_epsilon) {
					Z3 = 1;
					break;
				}
			}
			
			if (Z3) {
				#pragma omp simd
				for (i = 0; i < s0; i++)
					*(Ls3 + i) *= epsilon;
				*(Ls3 + s1) += 1;
			}
		} else {
			for (i = 0; i < s0; i++)
				*(Ls3 + i) = *(Ls1 + i);
			*(Ls3 + s1) = *(Ls1 + s1);
		}
	}
}

static void ProbChangeExp(double *m, double *E, double v)
{
	v = (v < 1e-6) ? 1e-6 : v;
	double A = *(m), C = *(m + 1), G = *(m + 2), T = *(m + 3), I = *(m + 4);
	double k1 = *(m + 5), k2 = *(m + 6), k3 = *(m + 7), k4 = *(m + 8), k5 = *(m + 9), k6 = *(m + 10);
	v = v/(2*(A*I*k6 + C*I*k6 + G*I*k6 + I*T*k6 + G*T + A*C*k3 + A*G*k1 + C*G*k5 + A*T*k4 + C*T*k2));
	
	double *Q = (double *) calloc(25, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *F = (double *) calloc(25, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *H = (double *) calloc(25, sizeof(double)); // initialized to zero (thread-safe on Windows)
	
	*(Q + 1) = C*k3*v;
	*(Q + 2) = G*k1*v;
	*(Q + 3) = T*k4*v;
	*(Q + 4) = I*k6*v;
	*(Q + 0) = -1*(*(Q + 1) + *(Q + 2) + *(Q + 3) + *(Q + 4));
	
	*(Q + 5) = A*k3*v;
	*(Q + 7) = G*k5*v;
	*(Q + 8) = T*k2*v;
	*(Q + 9) = I*k6*v;
	*(Q + 6) = -1*(*(Q + 5) + *(Q + 7) + *(Q + 8) + *(Q + 9));
	
	*(Q + 10) = A*k1*v;
	*(Q + 11) = C*k5*v;
	*(Q + 13) = T*v;
	*(Q + 14) = I*k6*v;
	*(Q + 12) = -1*(*(Q + 10) + *(Q + 11) + *(Q + 13) + *(Q + 14));
	
	*(Q + 15) = A*k4*v;
	*(Q + 16) = C*k2*v;
	*(Q + 17) = G*v;
	*(Q + 19) = I*k6*v;
	*(Q + 18) = -1*(*(Q + 15) + *(Q + 16) + *(Q + 17) + *(Q + 19));
	
	*(Q + 20) = A*k6*v;
	*(Q + 21) = C*k6*v;
	*(Q + 22) = G*k6*v;
	*(Q + 23) = T*k6*v;
	*(Q + 24) = -1*(*(Q + 20) + *(Q + 21) + *(Q + 22) + *(Q + 23));
	
	*(F + 0) = 1;
	*(F + 6) = 1;
	*(F + 12) = 1;
	*(F + 18) = 1;
	*(F + 24) = 1;
	
	// calculate infinity norm of the matrix Q
	int i, j, k, l;
	double r, x = 0, y;
	for (i = 0; i < 5; i++) {
		r = 0;
		for (j = i; j < 25; j += 5) {
			if (*(Q + j) >= 0) {
				r += *(Q + j);
			} else {
				r -= *(Q + j);
			}
		}
		if (r > x)
			x = r; // maximum absolute row sum
	}
	
	// E = exp(Q) = exp(Q/x)^x
	x = ceil(log2(x));
	if (x > 0) {
		double m = exp2(x);
		for (i = 0; i < 25; i++)
			*(Q + i) /= m;
	}
	
	k = 0;
	do {
		k++;
		
		// E = E + F
		for (i = 0; i < 25; i++)
			*(E + i) += *(F + i);
		
		// F = (A*F)/k
		for (i = 0; i < 5; i++) {
			for (j = 0; j < 5; j++) {
				*(H + 5*j + i) = 0;
				for (l = 0; l < 5; l++)
					*(H + 5*j + i) += *(Q + 5*l + i) * *(F + 5*j + l);
			}
		}
		for (i = 0; i < 25; i++)
			*(F + i) = *(H + i)/k;
		
		// G = E + F - E gives machine precision
		for (i = 0; i < 25; i++) {
			*(H + i) = *(E + i);
			*(H + i) += *(F + i);
			*(H + i) -= *(E + i);
		}
		
		// calculate the one norm of the matrix G
		y = 0;
		for (i = 0; i < 25;) {
			r = 0;
			for (j = i; j < i + 5; j++) {
				if (*(H + j) >= 0) {
					r += *(H + j);
				} else {
					r -= *(H + j);
				}
			}
			if (r > y)
				y = r; // maximum absolute column sum
			i = j;
		}
	} while (y > 0);
	
	if (x > 0) {
		for (k = 1; k <= x; k++) {
			for (i = 0; i < 5; i++) {
				for (j = 0; j < 5; j++) {
					*(H + 5*j + i) = 0;
					for (l = 0; l < 5; l++)
						*(H + 5*j + i) += *(E + 5*l + i) * *(E + 5*j + l);
				}
			}
			for (i = 0; i < 25; i++)
				*(E + i) = *(H + i);
		}
	}
	
	free(Q);
	free(F);
	free(H);
}

static void ProbChangeExpAA(double *m, double *E, double v)
{
	int i, j, k, l;
	double r, x = 0, y;
	
	double *Q = (double *) calloc(441, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *F = (double *) calloc(441, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *H = (double *) calloc(441, sizeof(double)); // initialized to zero (thread-safe on Windows)
	
	// fill Q matrix
	k = 0; // starting index in m
	for (j = 1; j < 20; j++) {
		for (i = 0; i < j; i++) {
			*(Q + 21*j + i) = *(m + k) * *(m + 190 + i);
			*(Q + 21*i + j) = *(m + k) * *(m + 190 + j);
			k++;
		}
	}
	for (i = 0; i < 20; i++) { // j = 20 (indels)
		*(Q + 21*j + i) = *(m + 211) * *(m + 190 + i);
		*(Q + 21*i + j) = *(m + 211) * *(m + 210);
	}
	
	r = 0; // sum of Q*freqs
	for (j = 0; j < 21; j++)
		for (i = 0; i < 21; i++)
			r += *(Q + 21*j + i) * *(m + 190 + j);
	
	v = (v < 1e-6) ? 1e-6 : v;
	v /= r;
	for (j = 0; j < 21; j++) {
		for (i = 0; i < 21; i++) {
			if (i != j) {
				*(Q + 21*j + i) *= v;
				*(Q + 21*j + j) -= *(Q + 21*j + i);
			}
		}
	}
	
	*(F + 0) = 1;
	*(F + 22) = 1;
	*(F + 44) = 1;
	*(F + 66) = 1;
	*(F + 88) = 1;
	*(F + 110) = 1;
	*(F + 132) = 1;
	*(F + 154) = 1;
	*(F + 176) = 1;
	*(F + 198) = 1;
	*(F + 220) = 1;
	*(F + 242) = 1;
	*(F + 264) = 1;
	*(F + 286) = 1;
	*(F + 308) = 1;
	*(F + 330) = 1;
	*(F + 352) = 1;
	*(F + 374) = 1;
	*(F + 396) = 1;
	*(F + 418) = 1;
	*(F + 440) = 1;
	
	// calculate infinity norm of the matrix Q
	for (i = 0; i < 21; i++) {
		r = 0;
		for (j = i; j < 441; j += 21) {
			if (*(Q + j) >= 0) {
				r += *(Q + j);
			} else {
				r -= *(Q + j);
			}
		}
		if (r > x)
			x = r; // maximum absolute row sum
	}
	
	// E = exp(Q) = exp(Q/x)^x
	x = ceil(log2(x));
	if (x > 0) {
		double m = exp2(x);
		for (i = 0; i < 441; i++)
			*(Q + i) /= m;
	}
	
	k = 0;
	do {
		k++;
		
		// E = E + F
		for (i = 0; i < 441; i++)
			*(E + i) += *(F + i);
		
		// F = (A*F)/k
		for (i = 0; i < 21; i++) {
			for (j = 0; j < 21; j++) {
				*(H + 21*j + i) = 0;
				for (l = 0; l < 21; l++)
					*(H + 21*j + i) += *(Q + 21*l + i) * *(F + 21*j + l);
			}
		}
		for (i = 0; i < 441; i++)
			*(F + i) = *(H + i)/k;
		
		// G = E + F - E gives machine precision
		for (i = 0; i < 441; i++) {
			*(H + i) = *(E + i);
			*(H + i) += *(F + i);
			*(H + i) -= *(E + i);
		}
		
		// calculate the one norm of the matrix G
		y = 0;
		for (i = 0; i < 441;) {
			r = 0;
			for (j = i; j < i + 21; j++) {
				if (*(H + j) >= 0) {
					r += *(H + j);
				} else {
					r -= *(H + j);
				}
			}
			if (r > y)
				y = r; // maximum absolute column sum
			i = j;
		}
	} while (y > 0);
	
	if (x > 0) {
		for (k = 1; k <= x; k++) {
			for (i = 0; i < 21; i++) {
				for (j = 0; j < 21; j++) {
					*(H + 21*j + i) = 0;
					for (l = 0; l < 21; l++)
						*(H + 21*j + i) += *(E + 21*l + i) * *(E + 21*j + l);
				}
			}
			for (i = 0; i < 441; i++)
				*(E + i) = *(H + i);
		}
	}
	
	free(Q);
	free(F);
	free(H);
}

static void Transpose(double *E, int n)
{
	double temp;
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			temp = *(E + n*i + j);
			*(E + n*i + j) = *(E + n*j + i);
			*(E + n*j + i) = temp;
		}
	}
}

SEXP clusterML(SEXP x, SEXP y, SEXP model, SEXP branches, SEXP lengths, SEXP states, SEXP type, SEXP weights, SEXP uLengths, SEXP nThreads)
{
	// initialize variables
	XStringSet_holder y_set;
	Chars_holder y_i;
	y_set = hold_XStringSet(y);
	int length = get_length_from_XStringSet_holder(&y_set);
	int i, j, k, o, p, numRates, indels, row, params;
	int *T = INTEGER(x); // Tree Topology
	double *m = REAL(model); // Substitution Model
	double s = asReal(states); // > 0 for reconstruct
	int t = asInteger(type);
	const int s1 = (t == 3) ? 21 : 5;
	const int width = s1 + 1;
	const int s2 = s1*s1;
	int *W = INTEGER(weights);
	int nthreads = asInteger(nThreads);
	
	SEXP dims;
	PROTECT(dims = GET_DIM(x));
	int l1 = INTEGER(dims)[0]; // number of rows in x (at most the length of x - 1)
	UNPROTECT(1);
	
	// alternative branch lengths
	int altL = length(lengths); // number of altered lengths
	int *ls = INTEGER(lengths); // index of new branch lengths
	int altB = length(branches); // number of altered branches
	int *bs = INTEGER(branches); // altered branches
	int lu = length(uLengths); // number of unique branch lengths
	double *ul = REAL(uLengths); // unique branch lengths
	// If altL != altB then NNI mode, in which case,
	// bs is the index of the central branch of the NNI,
	// and ls are the five associated branch lengths:
	// center, opposite, down-left, down-right, up
	// If bs > 0 then swap down-left with opposite
	// If bs < 0 then swap down-right with opposite
	if (s != 0 && altL != 0)
		error("branches or lengths specified when states is not zero.");
	
	// calculate a vector of sequence lengths
	y_i = get_elt_from_XStringSet_holder(&y_set, 0);
	int maxWidth = y_i.length; // assume all sequences are maxWidth (aligned)
	
	double *node;
	// matrix of base probabilities at each internal node (including the root)
	// [node base site]
	if (s > 0)
		node = (double *) calloc(l1*s1*maxWidth, sizeof(double)); // initialized to zero (thread-safe on Windows)
	
	// L_unknown is split into multiple versions to allow fixed loop iterations for compiler optimizations
	static void (*L_unknown_)(double *, const int, const int, const int, const double *, const int, const int, const int, const double, const double, const int);
	static void (*L_known_)(const char *, double *, int);
	if (t == 3) {
		params = 212;
		numRates = (length(model) - params)/2; // number of bins for the gamma distribution
		indels = (*(m + 210) == 0) ? 0 : 1;
		if (indels) {
			L_unknown_ = &L_unknown_AA_Indels;
		} else {
			L_unknown_ = &L_unknown_AA;
		}
		L_known_ = &L_known_AA;
	} else {
		params = 11;
		numRates = (length(model) - params)/2; // number of bins for the gamma distribution
		indels = (*(m + 4) == 0) ? 0 : 1;
		if (indels) {
			L_unknown_ = &L_unknown_Indels;
		} else {
			L_unknown_ = &L_unknown;
		}
		L_known_ = &L_known;
	}
	
	const double epsilon = pow(2, DBL_MAX_EXP/4);
	const double inv_epsilon = 1/epsilon;
	const double log_epsilon = log(epsilon);
	int size = lu*s2;
	
	int *Up;
	double *P;
	if (s > 0 || altB > 0) { // initialize pointers up the tree
		Up = Calloc(l1 - 1, int); // root node is untouched
		for (j = 1; j < l1; j++) { // for each node
			if (T[2*l1 + j] > 0) { // first branch is a node
				row = T[2*l1 + j] - 1;
				if (s > 0) {
					Up[row] = j;
				} else {
					Up[row] = -j;
				}
			}
			if (T[3*l1 + j] > 0) { // second branch is a node
				row = T[3*l1 + j] - 1;
				if (s > 0) {
					Up[row] = j;
				} else {
					Up[row] = -j;
				}
			}
		}
		for (j = 0; j < altB; j++) {
			p = *(bs + j);
			if (p < 0) // only applicable in NNI mode
				p *= -1;
			if (p > l1)
				p -= l1;
			p--;
			while (p < l1 - 1 && Up[p] < 0) {
				Up[p] *= -1;
				p = Up[p];
			}
		}
		P = Calloc(size*numRates + s2, double); // initialized to zero
		// last matrix is identity
		if (t == 3) {
			P[size*numRates + 0] = 1;
			P[size*numRates + 22] = 1;
			P[size*numRates + 44] = 1;
			P[size*numRates + 66] = 1;
			P[size*numRates + 88] = 1;
			P[size*numRates + 110] = 1;
			P[size*numRates + 132] = 1;
			P[size*numRates + 154] = 1;
			P[size*numRates + 176] = 1;
			P[size*numRates + 198] = 1;
			P[size*numRates + 220] = 1;
			P[size*numRates + 242] = 1;
			P[size*numRates + 264] = 1;
			P[size*numRates + 286] = 1;
			P[size*numRates + 308] = 1;
			P[size*numRates + 330] = 1;
			P[size*numRates + 352] = 1;
			P[size*numRates + 374] = 1;
			P[size*numRates + 396] = 1;
			P[size*numRates + 418] = 1;
			P[size*numRates + 440] = 1;
		} else {
			P[size*numRates + 0] = 1;
			P[size*numRates + 6] = 1;
			P[size*numRates + 12] = 1;
			P[size*numRates + 18] = 1;
			P[size*numRates + 24] = 1;
		}
	} else {
		P = Calloc(size*numRates, double); // initialized to zero
	}
	
	for (k = 0; k < numRates; k++) { // for each bin of the gamma distribution determined by alpha
		// P = expM(Q*v)
		// transpose for cache efficiency
		#pragma omp parallel for num_threads(nthreads)
		for (i = 0; i < lu; i++) {
			if (t == 3) {
				ProbChangeExpAA(m, (P + i*s2 + k*size), ul[i] * *(m + k + params));
			} else {
				ProbChangeExp(m, (P + i*s2 + k*size), ul[i] * *(m + k + params));
			}
			Transpose(P + i*s2 + k*size, s1);
		}
	}
	
	double *sumL = Calloc(maxWidth*(altB + 1), double);
	#pragma omp parallel for private(j,k,o,p,y_i,row) num_threads(nthreads)
	for (i = 0; i < maxWidth; i++) { // for each position
		int weight;
		if (s > 0) { // reconstruct ancestral states
			weight = 1;
		} else {
			weight = *(W + i);
		}
		if (weight > 0) {
			double *Ls;
			if (altL != altB) {
				Ls = (double *) calloc((length + 1)*3*width, sizeof(double)); // initialized to zero (thread-safe on Windows)
			} else {
				Ls = (double *) calloc(l1*3*width, sizeof(double)); // initialized to zero (thread-safe on Windows)
			}
			// Ls[row, 0] = likelihood left of node
			// Ls[row, width] = likelihood right of node
			// Ls[row, 2*width] = likelihood top of node
			double *mins = (double *) calloc(altB + 1, sizeof(double)); // initialized to zero (thread-safe on Windows)
			
			for (k = 0; k < numRates; k++) { // for each bin of the gamma distribution determined by alpha
				for (j = 0; j < l1; j++) { // for each node
					// if first branch is leaf then its base L is 1
					if (T[2*l1 + j] < 0) { // first branch is a leaf
						if (k == 0) {
							y_i = get_elt_from_XStringSet_holder(&y_set, (-1*T[2*l1 + j] - 1));
							(*L_known_)(&y_i.ptr[i], (Ls + 0 + j*3*width), indels);
						}
					} else  { // first branch is a node
						// L is probability(branch lengths) * L(previous nodes)
						row = T[2*l1 + j] - 1;
						(*L_unknown_)(Ls, j*3*width + 0, row*3*width + 0, row*3*width + width, P, T[0*l1 + row] - 1 + k*lu, T[1*l1 + row] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
					}
					// if second branch is leaf then its base L is 1
					if (T[3*l1 + j] < 0) { // second branch is a leaf
						if (k == 0) {
							y_i = get_elt_from_XStringSet_holder(&y_set, (-1*T[3*l1 + j] - 1));
							(*L_known_)(&y_i.ptr[i], (Ls + width + j*3*width), indels);
						}
					} else { // second branch is a node
						// L is probability(branch lengths) * L(previous nodes)
						row = T[3*l1 + j] - 1;
						(*L_unknown_)(Ls, j*3*width + width, row*3*width + 0, row*3*width + width, P, T[0*l1 + row] - 1 + k*lu, T[1*l1 + row] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
					}
				}
				
				if (k > 0) // clear Ls[root, 2*width]
					for (j = 2*width; j < 3*width; j++)
						*(Ls + (l1 - 1)*3*width + j) = 0;
				
				if (s > 0 || altL > 0) { // descend tree
					for (j = l1 - 2; j >= 0; j--) { // for each node below the root
						row = Up[j]; // calculate likelihood at node above
						if (row > 0) {
							int side = (T[2*l1 + row] == j + 1) ? width : 0; // edge opposite j
							int off1 = (side == 0) ? 0 : l1;
							if (row == (l1 - 1)) { // root is above
								// no node above so ignore by giving Ls[root, 2*width] = {0} currently
								(*L_unknown_)(Ls, j*3*width + 2*width, row*3*width + 2*width, row*3*width + side, P, 0, T[off1 + row] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
							} else {
								int off2 = (T[2*l1 + Up[row]] == row + 1) ? 0 : l1; // edge of row
								(*L_unknown_)(Ls, j*3*width + 2*width, row*3*width + 2*width, row*3*width + side, P, T[off2 + Up[row]] - 1 + k*lu, T[off1 + row] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
							}
						}
					}
				}
				
				// calculate Likelihood of each base at final node
				row = l1 - 1;
				(*L_unknown_)(Ls, row*3*width + 2*width, row*3*width + 0, row*3*width + width, P, T[0*l1 + row] - 1 + k*lu, T[1*l1 + row] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
				
				if (s > 0) { // final node
					if (t == 3) {
						if (!((*(Ls + 0 + row*3*width) == 0 &&
							*(Ls + 1 + row*3*width) == 0 &&
							*(Ls + 2 + row*3*width) == 0 &&
							*(Ls + 3 + row*3*width) == 0 &&
							*(Ls + 4 + row*3*width) == 0 &&
							*(Ls + 5 + row*3*width) == 0 &&
							*(Ls + 6 + row*3*width) == 0 &&
							*(Ls + 7 + row*3*width) == 0 &&
							*(Ls + 8 + row*3*width) == 0 &&
							*(Ls + 9 + row*3*width) == 0 &&
							*(Ls + 10 + row*3*width) == 0 &&
							*(Ls + 11 + row*3*width) == 0 &&
							*(Ls + 12 + row*3*width) == 0 &&
							*(Ls + 13 + row*3*width) == 0 &&
							*(Ls + 14 + row*3*width) == 0 &&
							*(Ls + 15 + row*3*width) == 0 &&
							*(Ls + 16 + row*3*width) == 0 &&
							*(Ls + 17 + row*3*width) == 0 &&
							*(Ls + 18 + row*3*width) == 0 &&
							*(Ls + 19 + row*3*width) == 0) ||
							(*(Ls + width + 0 + row*3*width) == 0 &&
							*(Ls + width + 1 + row*3*width) == 0 &&
							*(Ls + width + 2 + row*3*width) == 0 &&
							*(Ls + width + 3 + row*3*width) == 0 &&
							*(Ls + width + 4 + row*3*width) == 0 &&
							*(Ls + width + 5 + row*3*width) == 0 &&
							*(Ls + width + 6 + row*3*width) == 0 &&
							*(Ls + width + 7 + row*3*width) == 0 &&
							*(Ls + width + 8 + row*3*width) == 0 &&
							*(Ls + width + 9 + row*3*width) == 0 &&
							*(Ls + width + 10 + row*3*width) == 0 &&
							*(Ls + width + 11 + row*3*width) == 0 &&
							*(Ls + width + 12 + row*3*width) == 0 &&
							*(Ls + width + 13 + row*3*width) == 0 &&
							*(Ls + width + 14 + row*3*width) == 0 &&
							*(Ls + width + 15 + row*3*width) == 0 &&
							*(Ls + width + 16 + row*3*width) == 0 &&
							*(Ls + width + 17 + row*3*width) == 0 &&
							*(Ls + width + 18 + row*3*width) == 0 &&
							*(Ls + width + 19 + row*3*width) == 0))) { // neither branch is a gap
							for (j = 0; j < s1; j++)
								node[row*maxWidth*s1 + i*s1 + j] += *(Ls + 2*width + j + row*3*width) * *(m + numRates + k + params);
						}
					} else {
						if (!((*(Ls + 0 + row*3*width) == 0 &&
							*(Ls + 1 + row*3*width) == 0 &&
							*(Ls + 2 + row*3*width) == 0 &&
							*(Ls + 3 + row*3*width) == 0) ||
							(*(Ls + width + 0 + row*3*width) == 0 &&
							*(Ls + width + 1 + row*3*width) == 0 &&
							*(Ls + width + 2 + row*3*width) == 0 &&
							*(Ls + width + 3 + row*3*width) == 0))) { // neither branch is a gap
							for (j = 0; j < s1; j++)
								node[row*maxWidth*s1 + i*s1 + j] += *(Ls + 2*width + j + row*3*width) * *(m + numRates + k + params);
						}
					}
				}
				
				// calculate overall Likelihood
				int count = 0;
				for (o = 0; o <= altB; o++) {
					if (o > 0) {
						if (altL != altB) { // NNI mode
							int down, flip, side;
							p = *(bs + o - 1);
							if (p < 0) {
								flip = -1; // swap down-right with opposite
								p *= -1;
							} else {
								flip = 1; // swap down-left with opposite
							}
							if (p > l1) {
								side = 0;
								p -= l1;
								p--;
								down = T[3*l1 + p] - 1; // center is right branch
							} else {
								side = width;
								p--;
								down = T[2*l1 + p] - 1; // center is left branch
							}
							
							if (p == row) { // root node
								int opposite;
								if (side == 0) {
									opposite = T[2*l1 + p] - 1;
								} else {
									opposite = T[3*l1 + p] - 1;
								}
								if (flip > 0) {
									// opposite-left merging with down-right
									(*L_unknown_)(Ls, l1*3*width + 0, opposite*3*width + 0, down*3*width + width, P, *(ls + (o - 1)*5 + 2) - 1 + k*lu, *(ls + (o - 1)*5 + 3) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
									// opposite-right merging with down-left
									(*L_unknown_)(Ls, l1*3*width + width, opposite*3*width + width, down*3*width + 0, P, *(ls + (o - 1)*5 + 4) - 1 + k*lu, *(ls + (o - 1)*5 + 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								} else {
									// opposite-left merging with down-left
									(*L_unknown_)(Ls, l1*3*width + 0, opposite*3*width + 0, down*3*width + 0, P, *(ls + (o - 1)*5 + 3) - 1 + k*lu, *(ls + (o - 1)*5 + 2) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
									// opposite-right merging with down-right
									(*L_unknown_)(Ls, l1*3*width + width, opposite*3*width + width, down*3*width + width, P, *(ls + (o - 1)*5 + 4) - 1 + k*lu, *(ls + (o - 1)*5 + 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								}
								// merge both across center
								if (side == 0) { // center is right branch
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, P, T[0*l1 + row] - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								} else { // center is left branch
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, P, T[1*l1 + row] - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								}
							} else {
								if (count == 11) {
									count = 1;
								} else {
									count++;
								}
								if (count == 1) {
									if (flip > 0) {
										// opposite merging with down-right
										(*L_unknown_)(Ls, l1*3*width + 0, p*3*width + side, down*3*width + width, P, *(ls + (o - 1)*5 + 2) - 1 + k*lu, *(ls + (o - 1)*5 + 3) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
										// up merging with down-left
										(*L_unknown_)(Ls, l1*3*width + width, p*3*width + 2*width, down*3*width + 0, P, *(ls + (o - 1)*5 + 4) - 1 + k*lu, *(ls + (o - 1)*5 + 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
									} else {
										// opposite merging with down-left
										(*L_unknown_)(Ls, l1*3*width + 0, p*3*width + side, down*3*width + 0, P, *(ls + (o - 1)*5 + 3) - 1 + k*lu, *(ls + (o - 1)*5 + 2) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
										// up merging with down-right
										(*L_unknown_)(Ls, l1*3*width + width, p*3*width + 2*width, down*3*width + width, P, *(ls + (o - 1)*5 + 4) - 1 + k*lu, *(ls + (o - 1)*5 + 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
									}
									// merge both across center
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, P, numRates*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
									
									// up merging with below
									(*L_unknown_)(Ls, length*3*width + 0, p*3*width + 2*width, l1*3*width + 0, P, *(ls + (o - 1)*5 + 4) - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
									if (flip > 0) {
										// down-left merging with below
										(*L_unknown_)(Ls, l1*3*width + 2*width, down*3*width + 0, l1*3*width + 0, P, *(ls + (o - 1)*5 + 1) - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
										// down-right merging with above
										(*L_unknown_)(Ls, length*3*width + width, down*3*width + width, l1*3*width + width, P, *(ls + (o - 1)*5 + 3) - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
										// opposite merging with above
										(*L_unknown_)(Ls, length*3*width + 2*width, p*3*width + side, l1*3*width + width, P, *(ls + (o - 1)*5 + 2) - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
									} else {
										// down-right merging with below
										(*L_unknown_)(Ls, l1*3*width + 2*width, down*3*width + width, l1*3*width + 0, P, *(ls + (o - 1)*5 + 1) - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
										// down-left merging with above
										(*L_unknown_)(Ls, length*3*width + width, down*3*width + 0, l1*3*width + width, P, *(ls + (o - 1)*5 + 2) - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
										// opposite merging with above
										(*L_unknown_)(Ls, length*3*width + 2*width, p*3*width + side, l1*3*width + width, P, *(ls + (o - 1)*5 + 3) - 1 + k*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
									}
								} else if (count == 2 || count == 7) {
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, P, numRates*lu, *(ls + (o - 1)*5 + 0) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
								} else if (count == 3 || count == 8) {
									if (flip > 0) {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 0, down*3*width + 0, P, numRates*lu, *(ls + (o - 1)*5 + 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
									} else {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 0, down*3*width + width, P, numRates*lu, *(ls + (o - 1)*5 + 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
									}
								} else if (count == 4 || count == 9) {
									if (flip > 0) {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + width, p*3*width + side, P, numRates*lu, *(ls + (o - 1)*5 + 2) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
									} else {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 2*width, down*3*width + 0, P, numRates*lu, *(ls + (o - 1)*5 + 2) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
									}
								} else if (count == 5 || count == 10) {
									if (flip > 0) {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 2*width, down*3*width + width, P, numRates*lu, *(ls + (o - 1)*5 + 3) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
									} else {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + width, p*3*width + side, P, numRates*lu, *(ls + (o - 1)*5 + 3) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
									}
								} else { // count == 6 || count == 11
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 2*width, p*3*width + 2*width, P, numRates*lu, *(ls + (o - 1)*5 + 4) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
								}
							}
						} else {
							p = *(bs + o - 1);
							int side;
							if (p > l1) {
								side = width;
								p -= l1;
							} else {
								side = 0;
							}
							
							if (p == l1) { // root node
								if (side == width) {
									(*L_unknown_)(Ls, row*3*width + 2*width, row*3*width + 0, row*3*width + width, P, T[0*l1 + row] - 1 + k*lu, *(ls + o - 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								} else {
									(*L_unknown_)(Ls, row*3*width + 2*width, row*3*width + 0, row*3*width + width, P, *(ls + o - 1) - 1 + k*lu, T[1*l1 + row] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								}
							} else { // root at p - 1
								// calculate likelihood at node
								int off1 = (T[2*l1 + Up[p - 1]] == p) ? 0 : l1; // branch of row
								if (side == 0) {
									(*L_unknown_)(Ls, row*3*width + 2*width, (p - 1)*3*width + width, (p - 1)*3*width + 2*width, P, T[1*l1 + p - 1] - 1 + k*lu, T[off1 + Up[p - 1]] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								} else {
									(*L_unknown_)(Ls, row*3*width + 2*width, (p - 1)*3*width + 0, (p - 1)*3*width + 2*width, P, T[0*l1 + p - 1] - 1 + k*lu, T[off1 + Up[p - 1]] - 1 + k*lu, s2, epsilon, inv_epsilon, 0);
								}
								(*L_unknown_)(Ls, row*3*width + 2*width, (p - 1)*3*width + side, row*3*width + 2*width, P, numRates*lu, *(ls + o - 1) - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
							}
						}
					}
					
					if (k > 0) {
						double delta = *(Ls + 2*width + s1 + row*3*width) - mins[o];
						if (delta < 0) {
							for (j = 0; j > delta; j--)
								*(sumL + i + o*maxWidth) *= inv_epsilon;
							mins[o] = *(Ls + 2*width + s1 + row*3*width);
						} else if (delta > 0) {
							while (delta > 0) {
								for (j = 0; j < s1; j++)
									*(Ls + 2*width + j + row*3*width) *= inv_epsilon;
								delta--;
							}
						}
					} else {
						mins[o] = *(Ls + 2*width + s1 + row*3*width);
					}
					
					if (t == 3) {
						for (j = 0; j < s1; j++)
							*(sumL + i + o*maxWidth) += (*(m + j + 190) * *(Ls + 2*width + j + row*3*width)) * *(m + numRates + k + params);
					} else {
						for (j = 0; j < s1; j++)
							*(sumL + i + o*maxWidth) += (*(m + j) * *(Ls + 2*width + j + row*3*width)) * *(m + numRates + k + params);
					}
				}
				
				if (s > 0) {
					for (j = l1 - 2; j >= 0; j--) { // for each node below the root
						// apply three-way parsimony to resolve gaps
						int c1 = 0, c2 = 0;
						for (o = 0; o < s1 - 1; o++) {
							c1 |= *(Ls + o + j*3*width) != 0;
							c2 |= *(Ls + width + o + j*3*width) != 0;
						}
						int c = c1 + c2;
						if (c >= 1) { // at least one non-gap
							int side = (T[2*l1 + Up[j]] == j + 1) ? 0 : width; // edge of row
							if (c == 2 || c2) {
								int off1 = (side == 0) ? 0 : l1;
								(*L_unknown_)(Ls, row*3*width + 2*width, Up[j]*3*width + side, j*3*width + 2*width, P, numRates*lu, T[off1 + Up[j]] - 1 + k*lu, s2, epsilon, inv_epsilon, 1);
								for (o = 0; o < s1; o++)
									node[j*maxWidth*s1 + i*s1 + o] += *(Ls + 2*width + o + row*3*width) * *(m + numRates + k + params);
							}
						}
					}
				}
			}
			
			for (o = 0; o <= altB; o++)
				if (*(sumL + i + o*maxWidth) > 0)
					*(sumL + i + o*maxWidth) = log(*(sumL + i + o*maxWidth)) - mins[o]*log_epsilon;
			
			free(Ls);
			free(mins);
		}
	}
	Free(P);
	if (altL > 0)
		Free(Up);
	
	SEXP ans;
	if (s > 0) { // return character states
		SEXP ans1, ans2;
		PROTECT(ans = allocVector(VECSXP, 2));
		PROTECT(ans1 = allocVector(STRSXP, l1));
		s = (s >= 1) ? 0.9999999 : s; // machine precision
		if (t == 3) {
			for (j = 0; j < l1; j++) {
				char *rans1 = Calloc(maxWidth + 1, char);
				for (i = 0; i < maxWidth; i++) {
					double La = *(m) * node[j*maxWidth*s1 + i*s1 + 0];
					double Lr = *(m + 1) * node[j*maxWidth*s1 + i*s1 + 1];
					double Ln = *(m + 2) * node[j*maxWidth*s1 + i*s1 + 2];
					double Ld = *(m + 3) * node[j*maxWidth*s1 + i*s1 + 3];
					double Lc = *(m + 4) * node[j*maxWidth*s1 + i*s1 + 4];
					double Lq = *(m + 5) * node[j*maxWidth*s1 + i*s1 + 5];
					double Le = *(m + 6) * node[j*maxWidth*s1 + i*s1 + 6];
					double Lg = *(m + 7) * node[j*maxWidth*s1 + i*s1 + 7];
					double Lh = *(m + 8) * node[j*maxWidth*s1 + i*s1 + 8];
					double Li = *(m + 9) * node[j*maxWidth*s1 + i*s1 + 9];
					double Ll = *(m + 10) * node[j*maxWidth*s1 + i*s1 + 10];
					double Lk = *(m + 11) * node[j*maxWidth*s1 + i*s1 + 11];
					double Lm = *(m + 12) * node[j*maxWidth*s1 + i*s1 + 12];
					double Lf = *(m + 13) * node[j*maxWidth*s1 + i*s1 + 13];
					double Lp = *(m + 14) * node[j*maxWidth*s1 + i*s1 + 14];
					double Lz = *(m + 15) * node[j*maxWidth*s1 + i*s1 + 15];
					double Lt = *(m + 16) * node[j*maxWidth*s1 + i*s1 + 16];
					double Lw = *(m + 17) * node[j*maxWidth*s1 + i*s1 + 17];
					double Ly = *(m + 18) * node[j*maxWidth*s1 + i*s1 + 18];
					double Lv = *(m + 19) * node[j*maxWidth*s1 + i*s1 + 19];
					double L = *(m + 20) * node[j*maxWidth*s1 + i*s1 + 20];
					if (s*L >= La && s*L >= Lr && s*L >= Ln && s*L >= Ld && s*L >= Lc && s*L >= Lq && s*L >= Le && s*L >= Lg && s*L >= Lh && s*L >= Li && s*L >= Ll && s*L >= Lk && s*L >= Lm && s*L >= Lf && s*L >= Lp && s*L >= Lz && s*L >= Lt && s*L >= Lw && s*L >= Ly && s*L >= Lv) {
						rans1[i] = '-';
					} else if (s*La >= Lr && s*La >= Ln && s*La >= Ld && s*La >= Lc && s*La >= Lq && s*La >= Le && s*La >= Lg && s*La >= Lh && s*La >= Li && s*La >= Ll && s*La >= Lk && s*La >= Lm && s*La >= Lf && s*La >= Lp && s*La >= Lz && s*La >= Lt && s*La >= Lw && s*La >= Ly && s*La >= Lv) {
						rans1[i] = 'A';
					} else if (s*Lr >= La && s*Lr >= Ln && s*Lr >= Ld && s*Lr >= Lc && s*Lr >= Lq && s*Lr >= Le && s*Lr >= Lg && s*Lr >= Lh && s*Lr >= Li && s*Lr >= Ll && s*Lr >= Lk && s*Lr >= Lm && s*Lr >= Lf && s*Lr >= Lp && s*Lr >= Lz && s*Lr >= Lt && s*Lr >= Lw && s*Lr >= Ly && s*Lr >= Lv) {
						rans1[i] = 'R';
					} else if (s*Ln >= La && s*Ln >= Lr && s*Ln >= Ld && s*Ln >= Lc && s*Ln >= Lq && s*Ln >= Le && s*Ln >= Lg && s*Ln >= Lh && s*Ln >= Li && s*Ln >= Ll && s*Ln >= Lk && s*Ln >= Lm && s*Ln >= Lf && s*Ln >= Lp && s*Ln >= Lz && s*Ln >= Lt && s*Ln >= Lw && s*Ln >= Ly && s*Ln >= Lv) {
						rans1[i] = 'N';
					} else if (s*Ld >= La && s*Ld >= Lr && s*Ld >= Ln && s*Ld >= Lc && s*Ld >= Lq && s*Ld >= Le && s*Ld >= Lg && s*Ld >= Lh && s*Ld >= Li && s*Ld >= Ll && s*Ld >= Lk && s*Ld >= Lm && s*Ld >= Lf && s*Ld >= Lp && s*Ld >= Lz && s*Ld >= Lt && s*Ld >= Lw && s*Ld >= Ly && s*Ld >= Lv) {
						rans1[i] = 'D';
					} else if (s*Lc >= La && s*Lc >= Lr && s*Lc >= Ln && s*Lc >= Ld && s*Lc >= Lq && s*Lc >= Le && s*Lc >= Lg && s*Lc >= Lh && s*Lc >= Li && s*Lc >= Ll && s*Lc >= Lk && s*Lc >= Lm && s*Lc >= Lf && s*Lc >= Lp && s*Lc >= Lz && s*Lc >= Lt && s*Lc >= Lw && s*Lc >= Ly && s*Lc >= Lv) {
						rans1[i] = 'C';
					} else if (s*Lq >= La && s*Lq >= Lr && s*Lq >= Ln && s*Lq >= Ld && s*Lq >= Lc && s*Lq >= Le && s*Lq >= Lg && s*Lq >= Lh && s*Lq >= Li && s*Lq >= Ll && s*Lq >= Lk && s*Lq >= Lm && s*Lq >= Lf && s*Lq >= Lp && s*Lq >= Lz && s*Lq >= Lt && s*Lq >= Lw && s*Lq >= Ly && s*Lq >= Lv) {
						rans1[i] = 'Q';
					} else if (s*Le >= La && s*Le >= Lr && s*Le >= Ln && s*Le >= Ld && s*Le >= Lc && s*Le >= Lq && s*Le >= Lg && s*Le >= Lh && s*Le >= Li && s*Le >= Ll && s*Le >= Lk && s*Le >= Lm && s*Le >= Lf && s*Le >= Lp && s*Le >= Lz && s*Le >= Lt && s*Le >= Lw && s*Le >= Ly && s*Le >= Lv) {
						rans1[i] = 'E';
					} else if (s*Lg >= La && s*Lg >= Lr && s*Lg >= Ln && s*Lg >= Ld && s*Lg >= Lc && s*Lg >= Lq && s*Lg >= Le && s*Lg >= Lh && s*Lg >= Li && s*Lg >= Ll && s*Lg >= Lk && s*Lg >= Lm && s*Lg >= Lf && s*Lg >= Lp && s*Lg >= Lz && s*Lg >= Lt && s*Lg >= Lw && s*Lg >= Ly && s*Lg >= Lv) {
						rans1[i] = 'G';
					} else if (s*Lh >= La && s*Lh >= Lr && s*Lh >= Ln && s*Lh >= Ld && s*Lh >= Lc && s*Lh >= Lq && s*Lh >= Le && s*Lh >= Lg && s*Lh >= Li && s*Lh >= Ll && s*Lh >= Lk && s*Lh >= Lm && s*Lh >= Lf && s*Lh >= Lp && s*Lh >= Lz && s*Lh >= Lt && s*Lh >= Lw && s*Lh >= Ly && s*Lh >= Lv) {
						rans1[i] = 'H';
					} else if (s*Li >= La && s*Li >= Lr && s*Li >= Ln && s*Li >= Ld && s*Li >= Lc && s*Li >= Lq && s*Li >= Le && s*Li >= Lg && s*Li >= Lh && s*Li >= Ll && s*Li >= Lk && s*Li >= Lm && s*Li >= Lf && s*Li >= Lp && s*Li >= Lz && s*Li >= Lt && s*Li >= Lw && s*Li >= Ly && s*Li >= Lv) {
						rans1[i] = 'I';
					} else if (s*Ll >= La && s*Ll >= Lr && s*Ll >= Ln && s*Ll >= Ld && s*Ll >= Lc && s*Ll >= Lq && s*Ll >= Le && s*Ll >= Lg && s*Ll >= Lh && s*Ll >= Li && s*Ll >= Lk && s*Ll >= Lm && s*Ll >= Lf && s*Ll >= Lp && s*Ll >= Lz && s*Ll >= Lt && s*Ll >= Lw && s*Ll >= Ly && s*Ll >= Lv) {
						rans1[i] = 'L';
					} else if (s*Lk >= La && s*Lk >= Lr && s*Lk >= Ln && s*Lk >= Ld && s*Lk >= Lc && s*Lk >= Lq && s*Lk >= Le && s*Lk >= Lg && s*Lk >= Lh && s*Lk >= Li && s*Lk >= Ll && s*Lk >= Lm && s*Lk >= Lf && s*Lk >= Lp && s*Lk >= Lz && s*Lk >= Lt && s*Lk >= Lw && s*Lk >= Ly && s*Lk >= Lv) {
						rans1[i] = 'K';
					} else if (s*Lm >= La && s*Lm >= Lr && s*Lm >= Ln && s*Lm >= Ld && s*Lm >= Lc && s*Lm >= Lq && s*Lm >= Le && s*Lm >= Lg && s*Lm >= Lh && s*Lm >= Li && s*Lm >= Ll && s*Lm >= Lk && s*Lm >= Lf && s*Lm >= Lp && s*Lm >= Lz && s*Lm >= Lt && s*Lm >= Lw && s*Lm >= Ly && s*Lm >= Lv) {
						rans1[i] = 'M';
					} else if (s*Lf >= La && s*Lf >= Lr && s*Lf >= Ln && s*Lf >= Ld && s*Lf >= Lc && s*Lf >= Lq && s*Lf >= Le && s*Lf >= Lg && s*Lf >= Lh && s*Lf >= Li && s*Lf >= Ll && s*Lf >= Lk && s*Lf >= Lm && s*Lf >= Lp && s*Lf >= Lz && s*Lf >= Lt && s*Lf >= Lw && s*Lf >= Ly && s*Lf >= Lv) {
						rans1[i] = 'F';
					} else if (s*Lp >= La && s*Lp >= Lr && s*Lp >= Ln && s*Lp >= Ld && s*Lp >= Lc && s*Lp >= Lq && s*Lp >= Le && s*Lp >= Lg && s*Lp >= Lh && s*Lp >= Li && s*Lp >= Ll && s*Lp >= Lk && s*Lp >= Lm && s*Lp >= Lf && s*Lp >= Lz && s*Lp >= Lt && s*Lp >= Lw && s*Lp >= Ly && s*Lp >= Lv) {
						rans1[i] = 'P';
					} else if (s*Lz >= La && s*Lz >= Lr && s*Lz >= Ln && s*Lz >= Ld && s*Lz >= Lc && s*Lz >= Lq && s*Lz >= Le && s*Lz >= Lg && s*Lz >= Lh && s*Lz >= Li && s*Lz >= Ll && s*Lz >= Lk && s*Lz >= Lm && s*Lz >= Lf && s*Lz >= Lp && s*Lz >= Lt && s*Lz >= Lw && s*Lz >= Ly && s*Lz >= Lv) {
						rans1[i] = 'S';
					} else if (s*Lt >= La && s*Lt >= Lr && s*Lt >= Ln && s*Lt >= Ld && s*Lt >= Lc && s*Lt >= Lq && s*Lt >= Le && s*Lt >= Lg && s*Lt >= Lh && s*Lt >= Li && s*Lt >= Ll && s*Lt >= Lk && s*Lt >= Lm && s*Lt >= Lf && s*Lt >= Lp && s*Lt >= Lz && s*Lt >= Lw && s*Lt >= Ly && s*Lt >= Lv) {
						rans1[i] = 'T';
					} else if (s*Lw >= La && s*Lw >= Lr && s*Lw >= Ln && s*Lw >= Ld && s*Lw >= Lc && s*Lw >= Lq && s*Lw >= Le && s*Lw >= Lg && s*Lw >= Lh && s*Lw >= Li && s*Lw >= Ll && s*Lw >= Lk && s*Lw >= Lm && s*Lw >= Lf && s*Lw >= Lp && s*Lw >= Lz && s*Lw >= Lt && s*Lw >= Ly && s*Lw >= Lv) {
						rans1[i] = 'W';
					} else if (s*Ly >= La && s*Ly >= Lr && s*Ly >= Ln && s*Ly >= Ld && s*Ly >= Lc && s*Ly >= Lq && s*Ly >= Le && s*Ly >= Lg && s*Ly >= Lh && s*Ly >= Li && s*Ly >= Ll && s*Ly >= Lk && s*Ly >= Lm && s*Ly >= Lf && s*Ly >= Lp && s*Ly >= Lz && s*Ly >= Lt && s*Ly >= Lw && s*Ly >= Lv) {
						rans1[i] = 'Y';
					} else if (s*Lv >= La && s*Lv >= Lr && s*Lv >= Ln && s*Lv >= Ld && s*Lv >= Lc && s*Lv >= Lq && s*Lv >= Le && s*Lv >= Lg && s*Lv >= Lh && s*Lv >= Li && s*Lv >= Ll && s*Lv >= Lk && s*Lv >= Lm && s*Lv >= Lf && s*Lv >= Lp && s*Lv >= Lz && s*Lv >= Lt && s*Lv >= Lw && s*Lv >= Ly) {
						rans1[i] = 'V';
					} else if (s*Ln >= La && s*Ln >= Lr && s*Ln >= Lc && s*Ln >= Lq && s*Ln >= Le && s*Ln >= Lg && s*Ln >= Lh && s*Ln >= Li && s*Ln >= Ll && s*Ln >= Lk && s*Ln >= Lm && s*Ln >= Lf && s*Ln >= Lp && s*Ln >= Lz && s*Ln >= Lt && s*Ln >= Lw && s*Ln >= Ly && s*Ln >= Lv &&
						s*Ld >= La && s*Ld >= Lr && s*Ld >= Lc && s*Ld >= Lq && s*Ld >= Le && s*Ld >= Lg && s*Ld >= Lh && s*Ld >= Li && s*Ld >= Ll && s*Ld >= Lk && s*Ld >= Lm && s*Ld >= Lf && s*Ld >= Lp && s*Ld >= Lz && s*Ld >= Lt && s*Ld >= Lw && s*Ld >= Ly && s*Ld >= Lv) {
						rans1[i] = 'B';
					} else if (s*Lq >= La && s*Lq >= Lr && s*Lq >= Ln && s*Lq >= Ld && s*Lq >= Lc && s*Lq >= Lg && s*Lq >= Lh && s*Lq >= Li && s*Lq >= Ll && s*Lq >= Lk && s*Lq >= Lm && s*Lq >= Lf && s*Lq >= Lp && s*Lq >= Lz && s*Lq >= Lt && s*Lq >= Lw && s*Lq >= Ly && s*Lq >= Lv &&
						s*Le >= La && s*Le >= Lr && s*Le >= Ln && s*Le >= Ld && s*Le >= Lc && s*Le >= Lg && s*Le >= Lh && s*Le >= Li && s*Le >= Ll && s*Le >= Lk && s*Le >= Lm && s*Le >= Lf && s*Le >= Lp && s*Le >= Lz && s*Le >= Lt && s*Le >= Lw && s*Le >= Ly && s*Le >= Lv) {
						rans1[i] = 'Z';
					} else if (s*Li >= La && s*Li >= Lr && s*Li >= Ln && s*Li >= Ld && s*Li >= Lc && s*Li >= Lq && s*Li >= Le && s*Li >= Lg && s*Li >= Lh && s*Li >= Lk && s*Li >= Lm && s*Li >= Lf && s*Li >= Lp && s*Li >= Lz && s*Li >= Lt && s*Li >= Lw && s*Li >= Ly && s*Li >= Lv &&
						s*Ll >= La && s*Ll >= Lr && s*Ll >= Ln && s*Ll >= Ld && s*Ll >= Lc && s*Ll >= Lq && s*Ll >= Le && s*Ll >= Lg && s*Ll >= Lh && s*Ll >= Lk && s*Ll >= Lm && s*Ll >= Lf && s*Ll >= Lp && s*Ll >= Lz && s*Ll >= Lt && s*Ll >= Lw && s*Ll >= Ly && s*Ll >= Lv) {
						rans1[i] = 'J';
					} else {
						rans1[i] = 'X';
					}
				}
				rans1[i] = '\0'; // end (null terminate) the string
				SET_STRING_ELT(ans1, j, mkChar(rans1));
				Free(rans1);
			}
		} else {
			char TU = t == 1 ? 'T' : 'U';
			for (j = 0; j < l1; j++) {
				char *rans1 = Calloc(maxWidth + 1, char);
				for (i = 0; i < maxWidth; i++) {
					double La = *(m) * node[j*maxWidth*s1 + i*s1 + 0];
					double Lc = *(m + 1) * node[j*maxWidth*s1 + i*s1 + 1];
					double Lg = *(m + 2) * node[j*maxWidth*s1 + i*s1 + 2];
					double Lt = *(m + 3) * node[j*maxWidth*s1 + i*s1 + 3];
					double Li = *(m + 4) * node[j*maxWidth*s1 + i*s1 + 4];
					if (s*Li >= La && s*Li >= Lc && s*Li >= Lg && s*Li >= Lt) {
						rans1[i] = '-';
					} else if (s*La > Lc && s*La > Lg && s*La > Lt) {
						rans1[i] = 'A';
					} else if (s*Lc > La && s*Lc > Lg && s*Lc > Lt) {
						rans1[i] = 'C';
					} else if (s*Lg > La && s*Lg > Lc && s*Lg > Lt) {
						rans1[i] = 'G';
					} else if (s*Lt > La && s*Lt > Lc && s*Lt > Lg) {
						rans1[i] = TU;
					} else if (s*La > Lg && s*La > Lt && s*Lc > Lg && s*Lc > Lt) { // M = A or C
						rans1[i] = 'M';
					} else if (s*La > Lc && s*La > Lt && s*Lg > Lc && s*Lg > Lt) { // R = A or G
						rans1[i] = 'R';
					} else if (s*La > Lc && s*La > Lg && s*Lt > Lc && s*Lt > Lg) { // W = A or T
						rans1[i] = 'W';
					} else if (s*Lc > La && s*Lc > Lt && s*Lg > La && s*Lg > Lt) { // S = C or G
						rans1[i] = 'S';
					} else if (s*Lc > La && s*Lc > Lg && s*Lt > La && s*Lt > Lg) { // Y = C or T
						rans1[i] = 'Y';
					} else if (s*Lg > La && s*Lg > Lc && s*Lt > La && s*Lt > Lc) { // K = G or T
						rans1[i] = 'K';
					} else if (s*La > Lt && s*Lc > Lt && s*Lg > Lt) { // V = A or C or G
						rans1[i] = 'V';
					} else if (s*La > Lg && s*Lc > Lg && s*Lt > Lg) { // H = A or C or T
						rans1[i] = 'H';
					} else if (s*La > Lc && s*Lg > Lc && s*Lt > Lc) { // D = A or G or T
						rans1[i] = 'D';
					} else if (s*Lc > La && s*Lg > La && s*Lt > La) { // B = C or G or T
						rans1[i] = 'B';
					} else { // N = A or C or G or T
						rans1[i] = 'N';
					}
				}
				rans1[i] = '\0'; // end (null terminate) the string
				SET_STRING_ELT(ans1, j, mkChar(rans1));
				Free(rans1);
			}
		}
		free(node);
		PROTECT(ans2 = allocVector(REALSXP, maxWidth));
		double *rans2 = REAL(ans2);
		for (i = 0; i < maxWidth; i++) {
			if (*(sumL + i) < 0) {
				rans2[i] = *(sumL + i);
			} else {
				rans2[i] = 0;
			}
		}
		
		SET_VECTOR_ELT(ans, 0, ans1);
		SET_VECTOR_ELT(ans, 1, ans2);
		UNPROTECT(2);
	} else { // return -LnL
		PROTECT(ans = allocVector(REALSXP, altB + 1));
		double *rans = REAL(ans);
		for (o = 0; o <= altB; o++) {
			rans[o] = 0;
			for (i = 0; i < maxWidth; i++) {
				int weight = *(W + i);
				if (weight > 0)
					rans[o] -= ((double)weight) * *(sumL + i + o*maxWidth);
			}
		}
	}
	
	Free(sumL);
	UNPROTECT(1);
	
	return ans;
}

SEXP expM(SEXP x, SEXP model, SEXP type)
{
	// initialize variables
	double l = asReal(x);
	double *m = REAL(model); // Substitution Model
	int t = asInteger(type);
	int size = (t == 3) ? 21 : 5;
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size, size));
	double *rans = REAL(ans);
	size *= size;
	for (int i = 0; i < size; i++)
		rans[i] = 0;
	
	if (t == 3) {
		ProbChangeExpAA(m, rans, l);
	} else {
		ProbChangeExp(m, rans, l);
	}
	
	UNPROTECT(1);
	
	return ans;
}

/****************************************************************************
 *                        Cluster Maximum Parsimony                         *
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

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// for calloc/free
#include <stdlib.h>

// DECIPHER header file
#include "DECIPHER.h"

static void allStates(double *R, int *P, double *S, int c, int k0, int o0, int k1, int o1, int k2, int o2, int scoreOnly)
{
	int i, j, w1, w2;
	double t1, t2, min1, min2;
	double *R0 = R + k0*3*c + o0;
	double *R1 = R + k1*3*c + o1;
	double *R2 = R + k2*3*c + o2;
	
	if (scoreOnly) {
		for (i = 0; i < c; i++) {
			min1 = R_PosInf;
			min2 = R_PosInf;
			for (j = 0; j < c; j++) {
				t1 = *(R1 + j) + *(S + i*c + j);
				t2 = *(R2 + j) + *(S + i*c + j);
				if (t1 < min1)
					min1 = t1;
				if (t2 < min2)
					min2 = t2;
			}
			if (min1 != R_PosInf) {
				*(R0 + i) = min1;
				if (min2 != R_PosInf)
					*(R0 + i) += min2;
			} else if (min2 != R_PosInf) {
				*(R0 + i) = min2;
			} // else already R_PosInf
		}
	} else {
		for (i = 0; i < c; i++) {
			min1 = R_PosInf;
			min2 = R_PosInf;
			w1 = 0;
			w2 = 0;
			for (j = 0; j < c; j++) {
				t1 = *(R1 + j) + *(S + i*c + j);
				t2 = *(R2 + j) + *(S + i*c + j);
				if (t1 < min1) {
					w1 = j;
					min1 = t1;
				}
				if (t2 < min2) {
					w2 = j;
					min2 = t2;
				}
			}
			if (min1 != R_PosInf) {
				*(R0 + i) = min1;
				if (min2 != R_PosInf)
					*(R0 + i) += min2;
			} else if (min2 != R_PosInf) {
				*(R0 + i) = min2;
			} // else already R_PosInf
			*(P + k1*2*c + o1 + i) = w1 + 1;
			*(P + k2*2*c + o2 + i) = w2 + 1;
		}
	}
}

SEXP clusterMP(SEXP x, SEXP z, SEXP s, SEXP sizes, SEXP scoreOnly, SEXP add, SEXP weights, SEXP nThreads)
{
	// initialize variables
	int i, j, k, m, w;
	int *T = INTEGER(x); // Tree Topology
	int *Z = INTEGER(z); // Sequences
	double *S = REAL(s); // Substitution Matrix
	int c = INTEGER(sizes)[0]; // Dimensions of S
	int l = INTEGER(sizes)[1]; // Number of Sites
	int n = INTEGER(sizes)[2]; // Number of Nodes
	int N = INTEGER(sizes)[3]; // Number of Sequences
	int only = asInteger(scoreOnly);
	int a = asInteger(add); // Sequence To Add
	int *W = INTEGER(weights);
	int nthreads = asInteger(nThreads);
	
	double *lengths, *score;
	int *nodes, *subM;
	if (only == 0) {
		lengths = (double *) calloc(2*n, sizeof(double)); // initialized to zero (thread-safe on Windows)
		nodes = (int *) calloc(n*l, sizeof(int)); // initialized to zero (thread-safe on Windows)
		subM = (int *) calloc(c*c, sizeof(int)); // initialized to zero (thread-safe on Windows)
	}
	if (a > 0) { // insert leaf
		score = (double *) calloc(2*n + 2, sizeof(double)); // initialized to zero (thread-safe on Windows)
	} else if (a < 0) { // NNIs
		score = (double *) calloc(2*n - 1, sizeof(double)); // initialized to zero (thread-safe on Windows)
	} else {
		score = (double *) calloc(1, sizeof(double)); // initialized to zero (thread-safe on Windows)
	}
	
	int *Up;
	if (a != 0) {
		Up = (int *) calloc(n - 1, sizeof(int)); // initialized to zero (thread-safe on Windows)
		for (j = 0; j < n; j++) {
			k = *(T + j);
			if (k > 0)
				Up[k - 1] = j;
			k = *(T + n + j);
			if (k > 0)
				Up[k - 1] = j;
		}
	}
	
	#pragma omp parallel for private(i,j,k,m,w) num_threads(nthreads)
	for (i = 0; i < l; i++) {
		int weight;
		if (only == 0) { // reconstruct ancestral states
			weight = 1;
		} else {
			weight = *(W + i);
		}
		if (weight > 0) {
			double *R = (double *) calloc(3*c*(n + 1), sizeof(double)); // initialized to zero (thread-safe on Windows)
			for (j = 0; j < 3*c*(n + 1); j++)
				*(R + j) = R_PosInf;
			int *P;
			if (only == 0)
				P = (int *) calloc(2*c*n, sizeof(int)); // initialized to zero (thread-safe on Windows)
			
			// determine states going up the tree
			for (j = 0; j < n; j++) {
				k = *(T + j);
				if (k < 0) {
					m = *(Z + i*N - k - 1);
					if (m != NA_INTEGER)
						*(R + j*3*c + m - 1) = 0;
				} else {
					allStates(R, P, S, c, j, 0, k - 1, 0, k - 1, c, only);
				}
				k = *(T + n + j);
				if (k < 0) {
					m = *(Z + i*N - k - 1);
					if (m != NA_INTEGER)
						*(R + j*3*c + m - 1 + c) = 0;
				} else {
					allStates(R, P, S, c, j, c, k - 1, 0, k - 1, c, only);
				}
			}
			allStates(R, P, S, c, n - 1, 2*c, n - 1, 0, n - 1, c, only);
			
			w = 0;
			double temp[c];
			for (j = 0; j < c; j++) {
				temp[j] = *(R + 3*c*(n - 1) + 2*c + j);
				if (temp[j] < temp[w])
					w = j;
			}
			if (temp[w] != R_PosInf) {
				#pragma omp critical
				{
					score[0] += weight*temp[w];
				}
			}
			
			if (only == 0) {
				if (temp[w] != R_PosInf) {
					*(nodes + i*n + n - 1) = w + 1;
				} else {
					*(nodes + i*n + n - 1) = NA_INTEGER;
				}
				*(P + (n - 1)*2*c + w) *= -1;
				*(P + (n - 1)*2*c + c + w) *= -1;
			}
			
			// determine states going down the tree
			if (a < 0 ||
				(a > 0 &&
				*(Z + i*N + a - 1) != NA_INTEGER)) {
				j = n - 2;
				while (j >= 0) {
					if (Up[j] == n - 1) { // root is above
						// pass through opposite node
						if (*(T + Up[j]) == j + 1) {
							for (k = 0; k < c; k++)
								*(R + j*3*c + 2*c + k) = *(R + (n - 1)*3*c + c + k);
						} else {
							for (k = 0; k < c; k++)
								*(R + j*3*c + 2*c + k) = *(R + (n - 1)*3*c + k);
						}
					} else {
						if (*(T + Up[j]) == j + 1) {
							k = c;
						} else {
							k = 0;
						}
						allStates(R, P, S, c, j, 2*c, Up[j], 2*c, Up[j], k, 1);
					}
					j--;
				}
				
				if (a < 0) { // NNIs
					int count = 0;
					for (j = n - 1; j >= 0; j--) {
						k = *(T + j);
						if (k > 0) {
							// swap left-left with right
							count++;
							for (m = 0; m < 3*c; m++)
								*(R + 3*c*n + m) = R_PosInf;
							allStates(R, P, S, c, n, 0, k - 1, c, j, c, 1);
							allStates(R, P, S, c, n, c, k - 1, 0, n, 0, 1);
							if (j < n - 1) {
								allStates(R, P, S, c, n, 2*c, j, 2*c, n, c, 1);
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + 2*c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							} else {
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							}
							
							// swap left-right with right
							count++;
							for (m = 0; m < 3*c; m++)
								*(R + 3*c*n + m) = R_PosInf;
							allStates(R, P, S, c, n, 0, k - 1, 0, j, c, 1);
							allStates(R, P, S, c, n, c, k - 1, c, n, 0, 1);
							if (j < n - 1) {
								allStates(R, P, S, c, n, 2*c, j, 2*c, n, c, 1);
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + 2*c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							} else {
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							}
						}
						
						k = *(T + n + j);
						if (k > 0) {
							// swap right-left with left
							count++;
							for (m = 0; m < 3*c; m++)
								*(R + 3*c*n + m) = R_PosInf;
							allStates(R, P, S, c, n, 0, k - 1, c, j, 0, 1);
							allStates(R, P, S, c, n, c, k - 1, 0, n, 0, 1);
							if (j < n - 1) {
								allStates(R, P, S, c, n, 2*c, j, 2*c, n, c, 1);
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + 2*c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							} else {
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							}
							
							// swap right-right with left
							count++;
							for (m = 0; m < 3*c; m++)
								*(R + 3*c*n + m) = R_PosInf;
							allStates(R, P, S, c, n, 0, k - 1, 0, j, 0, 1);
							allStates(R, P, S, c, n, c, k - 1, c, n, 0, 1);
							if (j < n - 1) {
								allStates(R, P, S, c, n, 2*c, j, 2*c, n, c, 1);
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + 2*c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							} else {
								w = 0;
								for (m = 0; m < c; m++) {
									temp[m] = *(R + 3*c*n + c + m);
									if (temp[m] < temp[w])
										w = m;
								}
								if (temp[w] != R_PosInf) {
									#pragma omp critical
									{
										score[count] += weight*temp[w];
									}
								}
							}
						}
					}
				} else { // insert leaf
					m = *(Z + i*N + a - 1) - 1;
					
					// add to root
					*(R + 3*c*n + m) = 0;
					allStates(R, P, S, c, n, c, n - 1, 2*c, n, 0, 1);
					w = 0;
					for (j = 0; j < c; j++) {
						temp[j] = *(R + 3*c*n + c + j);
						if (temp[j] < temp[w])
							w = j;
					}
					if (*(R + 3*c*n + c + w) != R_PosInf) {
						#pragma omp critical
						{
							score[2*n + 1] += weight*temp[w];
						}
					}
					
					for (j = 0; j < c; j++)
						*(R + 3*c*(n - 1) + 2*c + j) = R_PosInf;
					*(R + 3*c*(n - 1) + 2*c + m) = 0;
					for (j = 0; j < n; j++) {
						// add to the first column of row j
						for (k = 0; k < 3*c; k++)
							*(R + 3*c*n + k) = R_PosInf;
						allStates(R, P, S, c, n, 0, j, 0, n - 1, 2*c, 1);
						allStates(R, P, S, c, n, c, j, c, n, 0, 1);
						if (j < n - 1) {
							allStates(R, P, S, c, n, 2*c, j, 2*c, n, c, 1);
							w = 0;
							for (k = 0; k < c; k++) {
								temp[k] = *(R + 3*c*n + 2*c + k);
								if (temp[k] < temp[w])
									w = k;
							}
						} else {
							w = 0;
							for (k = 0; k < c; k++) {
								temp[k] = *(R + 3*c*n + c + k);
								if (temp[k] < temp[w])
									w = k;
							}
						}
						if (temp[w] != R_PosInf) {
							#pragma omp critical
							{
								score[j + 1] += weight*temp[w];
							}
						}
						
						// add to the second column of row j
						for (k = 0; k < 3*c; k++)
							*(R + 3*c*n + k) = R_PosInf;
						allStates(R, P, S, c, n, 0, j, c, n - 1, 2*c, 1);
						allStates(R, P, S, c, n, c, j, 0, n, 0, 1);
						if (j < n - 1) {
							allStates(R, P, S, c, n, 2*c, j, 2*c, n, c, 1);
							w = 0;
							for (k = 0; k < c; k++) {
								temp[k] = *(R + 3*c*n + 2*c + k);
								if (temp[k] < temp[w])
									w = k;
							}
						} else {
							w = 0;
							for (k = 0; k < c; k++) {
								temp[k] = *(R + 3*c*n + c + k);
								if (temp[k] < temp[w])
									w = k;
							}
						}
						if (temp[w] != R_PosInf) {
							#pragma omp critical
							{
								score[n + j + 1] += weight*temp[w];
							}
						}
					}
				}
			}
			
			if (only != 0) {
				free(R);
				continue;
			}
			
			for (j = n - 1; j >= 0; j--) {
				for (w = 0; w < c; w++) {
					m = *(P + 2*c*j + w);
					if (m < 0)
						break;
				}
				m *= -1;
				m--;
				#pragma omp critical
				{
					*(lengths + j) += *(S + w*c + m);
				}
				k = *(T + j);
				if (k > 0) {
					k--;
					if (*(R + 3*c*j + m) != R_PosInf) {
						*(nodes + i*n + k) = m + 1;
						w = *(nodes + i*n + j);
						if (w != NA_INTEGER) {
							#pragma omp critical
							{
								*(subM + c*m + w - 1) += 1;
							}
						}
					} else {
						*(nodes + i*n + k) = NA_INTEGER;
					}
					*(P + 2*c*k + m) *= -1;
					*(P + 2*c*k + c + m) *= -1;
				} else {
					m = *(Z + i*N - k - 1);
					if (m != NA_INTEGER) {
						w = *(nodes + i*n + j);
						if (w != NA_INTEGER) {
							#pragma omp critical
							{
								*(subM + c*(m - 1) + w - 1) += 1;
							}
						}
					}
				}
				
				for (w = 0; w < c; w++) {
					m = *(P + 2*c*j + c + w);
					if (m < 0)
						break;
				}
				m *= -1;
				m--;
				#pragma omp critical
				{
					*(lengths + n + j) += *(S + w*c + m);
				}
				k = *(T + n + j);
				if (k > 0) {
					k--;
					if (*(R + 3*c*j + c + m) != R_PosInf) {
						*(nodes + i*n + k) = m + 1;
						w = *(nodes + i*n + j);
						if (w != NA_INTEGER) {
							#pragma omp critical
							{
								*(subM + c*m + w - 1) += 1;
							}
						}
					} else {
						*(nodes + i*n + k) = NA_INTEGER;
					}
					*(P + 2*c*k + m) *= -1;
					*(P + 2*c*k + c + m) *= -1;
				} else {
					m = *(Z + i*N - k - 1);
					if (m != NA_INTEGER) {
						w = *(nodes + i*n + j);
						if (w != NA_INTEGER) {
							#pragma omp critical
							{
								*(subM + c*(m - 1) + w - 1) += 1;
							}
						}
					}
				}
			}
			
			free(R);
			if (only == 0)
				free(P);
		}
	}
	
	SEXP ans1;
	if (a > 0) {
		j = 2*n + 2;
	} else if (a < 0) {
		j = 2*n - 1;
	} else {
		j = 1;
	}
	PROTECT(ans1 = allocVector(REALSXP, j));
	double *rans = REAL(ans1);
	for (i = 0; i < j; i++)
		rans[i] = score[i];
	
	if (a != 0)
		free(Up);
	free(score);
	
	if (only != 0) {
		UNPROTECT(1);
		
		return ans1;
	}
	
	SEXP ans2;
	PROTECT(ans2 = allocMatrix(INTSXP, n, l));
	int *rans_int = INTEGER(ans2);
	for (i = 0; i < n*l; i++)
		rans_int[i] = nodes[i];
	free(nodes);
	
	SEXP ans3;
	PROTECT(ans3 = allocMatrix(REALSXP, n, 2));
	rans = REAL(ans3);
	for (i = 0; i < 2*n; i++)
		rans[i] = lengths[i];
	free(lengths);
	
	SEXP ans4;
	PROTECT(ans4 = allocMatrix(REALSXP, c, c));
	rans = REAL(ans4);
	for (i = 0; i < c*c; i++)
		rans[i] = subM[i];
	free(subM);
	
	SEXP ret_list;
	PROTECT(ret_list = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(ret_list, 0, ans1);
	SET_VECTOR_ELT(ret_list, 1, ans2);
	SET_VECTOR_ELT(ret_list, 2, ans3);
	SET_VECTOR_ELT(ret_list, 3, ans4);
	
	UNPROTECT(5);
	
	return ret_list;
}

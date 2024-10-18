/****************************************************************************
 *                        Cluster Maximum Parsimony                         *
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

// DECIPHER header file
#include "DECIPHER.h"

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

static void allStates(double *R, int *P, double *S, int c, int k0, int o0, int k1, int o1, int k2, int o2, int scoreOnly)
{
	int i, j, w1, w2;
	double t1, t2, min1, min2;
	double *R0 = R + k0*3*c + o0;
	double *R1 = R + k1*3*c + o1;
	double *R2 = R + k2*3*c + o2;
	
	if (scoreOnly == 1) {
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
			} // else R0 already R_PosInf
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
				*(P + k1*2*c + o1 + i) = w1 + 1;
				if (min2 != R_PosInf) {
					*(R0 + i) += min2;
					*(P + k2*2*c + o2 + i) = w2 + 1;
				}
			} else if (min2 != R_PosInf) {
				*(R0 + i) = min2;
				*(P + k2*2*c + o2 + i) = w2 + 1;
			} // else R0 already R_PosInf and P already zero
		}
	}
}

SEXP clusterMP(SEXP z, SEXP x, SEXP s, SEXP letters, SEXP scoreOnly, SEXP add, SEXP weights, SEXP nThreads)
{
	// initialize variables
	int i, j, k, m, w;
	int *T = INTEGER(z); // Tree Topology
	int n = length(z)/2; // number of nodes
	double *S = REAL(s); // Substitution Matrix
	int only = asInteger(scoreOnly);
	int a = asInteger(add); // Sequence To Add
	int *W = INTEGER(weights);
	int nthreads = asInteger(nThreads);
	
	XStringSet_holder x_set, l_set;
	Chars_holder x_i, l_i;
	x_set = hold_XStringSet(x);
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	int l = x_i.length; // number of sites
	l_set = hold_XStringSet(letters);
	l_i = get_elt_from_XStringSet_holder(&l_set, 0);
	int c = l_i.length; // number of states
	int *lkup = (int *) malloc(256*sizeof(int)); // thread-safe on Windows
	for (i = 0; i < 256; i++)
		lkup[i] = NA_INTEGER;
	for (i = 0; i < c; i++)
		lkup[(unsigned char)l_i.ptr[i]] = i;
	
	double *lengths, *score;
	int size, *nodes, *subM;
	if (only != 1) {
		lengths = (double *) calloc(2*n, sizeof(double)); // initialized to zero (thread-safe on Windows)
		subM = (int *) calloc(c*c, sizeof(int)); // initialized to zero (thread-safe on Windows)
		if (only != 0)
			nodes = (int *) malloc(n*l*sizeof(int)); // thread-safe on Windows
	}
	if (a > 0) { // insert leaf
		size = 2*n + 2;
	} else if (a < 0) { // NNIs
		size = 2*n - 1;
	} else {
		size = 1;
	}
	score = (double *) calloc(size, sizeof(double)); // initialized to zero (thread-safe on Windows)
	
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
	
	#ifdef _OPENMP
	#pragma omp parallel num_threads(nthreads)
	{
	#endif
		double *temp_score, *temp_lengths;
		int *temp_subM, *temp_nodes;
		temp_score = (double *) calloc(size, sizeof(double)); // initialized to zero (thread-safe on Windows)
		if (only != 1) {
			temp_lengths = (double *) calloc(2*n, sizeof(double)); // initialized to zero (thread-safe on Windows)
			temp_subM = (int *) calloc(c*c, sizeof(int)); // initialized to zero (thread-safe on Windows)
			temp_nodes = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
		}
		
		#ifdef _OPENMP
		#pragma omp for private(i,j,k,m,w)
		#endif
		for (i = 0; i < l; i++) {
			int weight;
			if (only == 0 || only == 1) {
				weight = *(W + i);
			} else { // reconstruct ancestral states
				weight = 1;
			}
			if (weight > 0) {
				double *R = (double *) calloc(3*c*(n + 1), sizeof(double)); // initialized to zero (thread-safe on Windows)
				for (j = 0; j < 3*c*(n + 1); j++)
					*(R + j) = R_PosInf;
				int *P;
				if (only != 1) {
					P = (int *) calloc(2*c*n, sizeof(int)); // initialized to zero (thread-safe on Windows)
					if (only != 0)
						for (j = 0; j < n; j++)
							temp_nodes[j] = 0;
				}
				
				// determine states going up the tree
				for (j = 0; j < n; j++) {
					k = *(T + j);
					if (k < 0) {
						x_i = get_elt_from_XStringSet_holder(&x_set, -k - 1);
						m = lkup[(unsigned char)x_i.ptr[i]];
						if (m != NA_INTEGER)
							*(R + j*3*c + m) = 0;
					} else {
						allStates(R, P, S, c, j, 0, k - 1, 0, k - 1, c, only);
					}
					k = *(T + n + j);
					if (k < 0) {
						x_i = get_elt_from_XStringSet_holder(&x_set, -k - 1);
						m = lkup[(unsigned char)x_i.ptr[i]];
						if (m != NA_INTEGER)
							*(R + j*3*c + m + c) = 0;
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
				if (temp[w] != R_PosInf)
					temp_score[0] += weight*temp[w];
				
				if (only != 1) {
					if (temp[w] != R_PosInf) {
						*(temp_nodes + n - 1) = w + 1;
					} else {
						*(temp_nodes + n - 1) = NA_INTEGER;
					}
					*(P + (n - 1)*2*c + w) *= -1;
					*(P + (n - 1)*2*c + c + w) *= -1;
				}
				
				// determine states going down the tree
				if (a > 0) {
					x_i = get_elt_from_XStringSet_holder(&x_set, a - 1);
					m = lkup[(unsigned char)x_i.ptr[i]];
				}
				if (a < 0 ||
					(a > 0 &&
					m != NA_INTEGER)) {
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
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
								} else {
									w = 0;
									for (m = 0; m < c; m++) {
										temp[m] = *(R + 3*c*n + c + m);
										if (temp[m] < temp[w])
											w = m;
									}
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
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
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
								} else {
									w = 0;
									for (m = 0; m < c; m++) {
										temp[m] = *(R + 3*c*n + c + m);
										if (temp[m] < temp[w])
											w = m;
									}
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
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
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
								} else {
									w = 0;
									for (m = 0; m < c; m++) {
										temp[m] = *(R + 3*c*n + c + m);
										if (temp[m] < temp[w])
											w = m;
									}
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
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
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
								} else {
									w = 0;
									for (m = 0; m < c; m++) {
										temp[m] = *(R + 3*c*n + c + m);
										if (temp[m] < temp[w])
											w = m;
									}
									if (temp[w] != R_PosInf)
										temp_score[count] += weight*temp[w];
								}
							}
						}
					} else { // insert leaf
						// add to root
						*(R + 3*c*n + m) = 0;
						allStates(R, P, S, c, n, c, n - 1, 2*c, n, 0, 1);
						w = 0;
						for (j = 0; j < c; j++) {
							temp[j] = *(R + 3*c*n + c + j);
							if (temp[j] < temp[w])
								w = j;
						}
						if (*(R + 3*c*n + c + w) != R_PosInf)
							temp_score[2*n + 1] += weight*temp[w];
						
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
							if (temp[w] != R_PosInf)
								temp_score[j + 1] += weight*temp[w];
							
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
							if (temp[w] != R_PosInf)
								temp_score[n + j + 1] += weight*temp[w];
						}
					}
				}
				
				if (only == 1) {
					free(R);
					continue;
				}
				
				for (j = n - 1; j >= 0; j--) {
					for (w = 0; w < c; w++) {
						m = *(P + 2*c*j + w);
						if (m < 0)
							break;
					}
					if (m < 0) {
						m *= -1;
						m--;
						*(temp_lengths + j) += *(S + w*c + m)*weight;
						k = *(T + j);
						if (k > 0) {
							k--;
							if (*(R + 3*c*j + m) != R_PosInf) {
								*(temp_nodes + k) = m + 1;
								w = *(temp_nodes + j);
								if (w != NA_INTEGER && w - 1 != m)
									*(temp_subM + c*m + w - 1) += weight;
							} else {
								*(temp_nodes + k) = NA_INTEGER;
							}
							*(P + 2*c*k + m) *= -1;
							*(P + 2*c*k + c + m) *= -1;
						} else {
							x_i = get_elt_from_XStringSet_holder(&x_set, -k - 1);
							m = lkup[(unsigned char)x_i.ptr[i]];
							if (m != NA_INTEGER) {
								w = *(temp_nodes + j);
								if (w != NA_INTEGER && w - 1 != m)
									*(temp_subM + c*m + w - 1) += weight;
							}
						}
					} else {
						k = *(T + j);
						if (k > 0)
							*(temp_nodes + k - 1) = NA_INTEGER;
					}
					
					for (w = 0; w < c; w++) {
						m = *(P + 2*c*j + c + w);
						if (m < 0)
							break;
					}
					if (m < 0) {
						m *= -1;
						m--;
						*(temp_lengths + n + j) += *(S + w*c + m)*weight;
						k = *(T + n + j);
						if (k > 0) {
							k--;
							if (*(R + 3*c*j + c + m) != R_PosInf) {
								*(temp_nodes + k) = m + 1;
								w = *(temp_nodes + j);
								if (w != NA_INTEGER && w - 1 != m)
									*(temp_subM + c*m + w - 1) += weight;
							} else {
								*(temp_nodes + k) = NA_INTEGER;
							}
							*(P + 2*c*k + m) *= -1;
							*(P + 2*c*k + c + m) *= -1;
						} else {
							x_i = get_elt_from_XStringSet_holder(&x_set, -k - 1);
							m = lkup[(unsigned char)x_i.ptr[i]];
							if (m != NA_INTEGER) {
								w = *(temp_nodes + j);
								if (w != NA_INTEGER && w - 1 != m)
									*(temp_subM + c*m + w - 1) += weight;
							}
						}
					} else {
						k = *(T + n + j);
						if (k > 0)
							*(temp_nodes + k - 1) = NA_INTEGER;
					}
				}
				
				free(R);
				free(P);
				
				// perform array reduction
				#ifdef _OPENMP
				#pragma omp critical
				{
				#endif
					if (only != 0)
						for (j = 0; j < n; j++)
							nodes[i*n + j] = temp_nodes[j];
				#ifdef _OPENMP
				}
				#endif
			}
		}
		
		// perform array reductions
		#ifdef _OPENMP
		#pragma omp critical
		{
		#endif
			for (i = 0; i < size; i++)
				score[i] += temp_score[i];
			if (only != 1) {
				for (i = 0; i < 2*n; i++)
					lengths[i] += temp_lengths[i];
				for (i = 0; i < c*c; i++)
					subM[i] += temp_subM[i];
			}
		#ifdef _OPENMP
		}
		#endif
		free(temp_score);
		if (only != 1) {
			free(temp_lengths);
			free(temp_subM);
			free(temp_nodes);
		}
	#ifdef _OPENMP
	}
	#endif
	free(lkup);
	
	SEXP ans1;
	PROTECT(ans1 = allocVector(REALSXP, size));
	double *rans = REAL(ans1);
	for (i = 0; i < size; i++)
		rans[i] = score[i];
	
	if (a != 0)
		free(Up);
	free(score);
	
	if (only == 1) {
		UNPROTECT(1);
		
		return ans1;
	}
	
	SEXP ans2;
	if (only == 0) {
		PROTECT(ans2 = allocMatrix(INTSXP, 0, 0));
	} else {
		PROTECT(ans2 = allocMatrix(INTSXP, n, l));
		int *rans_int = INTEGER(ans2);
		for (i = 0; i < n*l; i++)
			rans_int[i] = nodes[i];
		free(nodes);
	}
	
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

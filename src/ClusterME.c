/****************************************************************************
 *                        Cluster Minimum Evolution                         *
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

// for math functions
#include <math.h>

// for calloc/free
#include <stdlib.h>

// DECIPHER header file
#include "DECIPHER.h"

static R_len_t index1D(const int r, const int c, const int N)
{
	if (r > c) {
		return (2*(R_len_t)N - (R_len_t)c - 1)*(R_len_t)c/2 + (R_len_t)r - (R_len_t)c - 1;
	} else {
		return (2*(R_len_t)N - (R_len_t)r - 1)*(R_len_t)r/2 + (R_len_t)c - (R_len_t)r - 1;
	}
}

static R_len_t index2D(const int r, const int c, const int N)
{
	return (R_len_t)r + (R_len_t)c*(R_len_t)N;
}

// row sums of a matrix or 'dist' object
SEXP rowSums(SEXP dist, SEXP n)
{
	int i, j;
	double *D = REAL(dist);
	int N = asInteger(n);
	
	R_len_t (*ind)(const int, const int, const int);
	if (N < 0) {
		ind = &index1D;
		N *= -1;
	} else {
		ind = &index2D;
	}
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, N));
	double *rans = REAL(ans);
	for (i = 0; i < N; i++)
		rans[i] = 0;
	
	double d;
	for (i = 1; i < N; i++) {
		for (j = 0; j < i; j++) {
			d = D[ind(i, j, N)];
			rans[i] += d;
			rans[j] += d;
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

// patristic distances from clusters
// multiplies distance between partitions by elements of z
SEXP patristic(SEXP x, SEXP y, SEXP z)
{
	int i, j, k, posE, count = 0;
	int *C = INTEGER(x); // node indices
	double *H = REAL(y); // node heights
	double *mult = REAL(z); // standard deviations
	int num = length(z);
	int n = length(x)/2; // number of rows in x
	int N = n + 1;
	R_len_t l = (R_len_t)N*(R_len_t)n/2;
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, l));
	double *rans = REAL(ans);
	
	// record subtrees
	int *sides; // active side
	sides = (int *) calloc(n, sizeof(int)); // initialized to zero (thread-safe on Windows)
	int *index; // index of leaves
	index = (int *) malloc(N*sizeof(int)); // thread-safe on Windows
	int *ups; // pointers up tree
	ups = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	int *posL; // start of left branch indices
	posL = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	int *posR; // start of right branch indices
	posR = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	double *node_depth; // distance from node to root
	node_depth = (double *) malloc(n*sizeof(double)); // thread-safe on Windows
	double *leaf_depth; // distance from leaf to root
	leaf_depth = (double *) malloc(N*sizeof(double)); // thread-safe on Windows
	j = 0; // position in index
	i = n - 1; // row in clusters
	node_depth[i] = 0;
	ups[i] = n;
	do {
		if (sides[i] == 2) { // both sides visited
			posE = j;
			node_depth[i] *= 2;
			for (k = posL[i]; k < posR[i]; k++)
				for (j = posR[i]; j < posE; j++)
					rans[index1D(index[j], index[k], N)] = mult[count]*(leaf_depth[index[k]] + leaf_depth[index[j]] - node_depth[i]);
			count++;
			if (count == num)
				count = 0;
			i = ups[i]; // ascend tree
		} else {
			k = C[i + n*sides[i]];
			if (sides[i]) { // right branch
				posR[i] = j;
			} else { // left branch
				posL[i] = j;
			}
			if (k < 0) { // add leaf to index
				index[j] = -k - 1;
				leaf_depth[index[j]] = node_depth[i] + H[i + n*sides[i]];
				sides[i]++; // switch side
				j++;
			} else { // descend tree
				k--;
				ups[k] = i;
				node_depth[k] = node_depth[i] + H[i + n*sides[i]];
				sides[i]++; // switch side
				i = k;
			}
		}
	} while (i < n);
	
	free(sides);
	free(index);
	free(ups);
	free(posL);
	free(posR);
	free(leaf_depth);
	free(node_depth);
	
	UNPROTECT(1);
	
	return ans;
}

SEXP clusterME(SEXP x, SEXP y, SEXP l, SEXP f)
{
	// initialize variables
	int i, j, k, r, c;
	int *C = INTEGER(x); // tree topology
	double *D = REAL(y); // distance matrix
	int n = asInteger(l); // dimension of x
	int flag = asInteger(f); // 0 = total length, 1 = branch lengths, 2 = NNIs
	
	// index is split into multiple versions to allow full or half distance matrices
	R_len_t (*ind)(const int, const int, const int);
	if (n < 0) {
		ind = &index1D;
		n *= -1;
	} else {
		ind = &index2D;
	}
	int N = n + 1;
	
	SEXP ans;
	double *rans;
	int *B;
	if (flag == 0L) { // total length
		PROTECT(ans = allocVector(REALSXP, 1));
		rans = REAL(ans);
		rans[0] = 0;
	} else {
		PROTECT(ans = allocMatrix(REALSXP, n, 2));
		rans = REAL(ans);
		for (i = 0; i < 2*n; i++)
			rans[i] = 0;
		
		B = (int *) malloc((R_len_t)N*(R_len_t)n/2*sizeof(int)); // thread-safe on Windows
	}
	
	// record depth of each split
	int *leaf_depth = (int *) malloc(N*sizeof(int)); // thread-safe on Windows
	int *depth = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	i = n - 1;
	depth[i] = 0;
	while (i >= 0) {
		if (C[i] > 0) {
			depth[C[i] - 1] = depth[i] + 1;
		} else {
			leaf_depth[-C[i] - 1] = depth[i] + 1;
		}
		if (C[i + n] > 0) {
			depth[C[i + n] - 1] = depth[i] + 1;
		} else {
			leaf_depth[-C[i + n] - 1] = depth[i] + 1;
		}
		depth[i] *= -2;
		i--;
	}
	depth[n - 1] = -1; // root
	
	// record subtrees
	int *sides; // active side
	sides = (int *) calloc(n, sizeof(int)); // initialized to zero (thread-safe on Windows)
	int *index; // index of leaves
	index = (int *) malloc(N*sizeof(int)); // thread-safe on Windows
	int *ups; // pointers up tree
	ups = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	int *posL; // start of left branch indices
	posL = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	int *posR; // start of right branch indices
	posR = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	int *posE; // start of right branch indices
	posE = (int *) malloc(n*sizeof(int)); // thread-safe on Windows
	j = 0; // position in index
	i = n - 1; // row in C
	ups[i] = n;
	do {
		if (sides[i] == 2) { // both sides visited
			posE[i] = j;
			if (flag == 0) { // total length
				for (k = posL[i]; k < posR[i]; k++)
					for (j = posR[i]; j < posE[i]; j++)
						rans[0] += D[ind(index[j], index[k], N)]/pow(2, leaf_depth[index[j]] + leaf_depth[index[k]] + depth[i]);
			} else {
				for (k = posL[i]; k < posR[i]; k++)
					for (j = posR[i]; j < posE[i]; j++)
						B[index1D(index[j], index[k], N)] = leaf_depth[index[j]] + leaf_depth[index[k]] + depth[i];
			}
			i = ups[i]; // ascend tree
		} else {
			k = C[i + n*sides[i]];
			if (sides[i]) { // right branch
				posR[i] = j;
			} else { // left branch
				posL[i] = j;
			}
			sides[i]++; // switch side
			if (k < 0) { // add leaf to index
				index[j++] = -k - 1;
			} else { // descend tree
				k--;
				ups[k] = i;
				i = k;
			}
		}
	} while (i < n);
	
	free(leaf_depth);
	free(depth);
	free(sides);
	free(ups);
	
	if (flag == 0) { // total length
		rans[0] *= 2;
	} else {
		// calculate row sums
		double d;
		double *sumDist = (double *) calloc(N, sizeof(double)); // initialized to zero (thread-safe on Windows)
		for (i = 1; i < N; i++) {
			for (j = 0; j < i; j++) {
				d = D[ind(i, j, N)]/pow(2, B[index1D(i, j, N)]);
				sumDist[i] += d;
				sumDist[j] += d;
			}
		}
		
		double LLRL, LRRL, LLRR, LRRR;
		double *subDists; // a, b, ab, ac, bc, ac + bc
		subDists = (double *) calloc((n - 1)*6, sizeof(double)); // initialized to zero (thread-safe on Windows)
		for (i = 0; i < n - 1; i++) {
			if (C[i] > 0 && C[i + n] > 0) { // two branches
				LLRL = 0; // left-left right-left
				for (r = posL[C[i] - 1]; r < posR[C[i] - 1]; r++)
					for (c = posL[C[i + n] - 1]; c < posR[C[i + n] - 1]; c++)
						LLRL += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
				LRRL = 0; // left-right right-left
				for (r = posR[C[i] - 1]; r < posE[C[i] - 1]; r++)
					for (c = posL[C[i + n] - 1]; c < posR[C[i + n] - 1]; c++)
						LRRL += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
				LLRR = 0; // left-left right-right
				for (r = posL[C[i] - 1]; r < posR[C[i] - 1]; r++)
					for (c = posR[C[i + n] - 1]; c < posE[C[i + n] - 1]; c++)
						LLRR += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
				LRRR = 0; // left-right right-right
				for (r = posR[C[i] - 1]; r < posE[C[i] - 1]; r++)
					for (c = posR[C[i + n] - 1]; c < posE[C[i + n] - 1]; c++)
						LRRR += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
			}
			for (j = 0; j < 2; j++) {
				if (C[i + j*n] > 0) { // branch
					subDists[i + j*(n - 1)] = subDists[C[i + j*n] + 0*(n - 1) - 1] + subDists[C[i + j*n] + 1*(n - 1) - 1];
					double ab = subDists[C[i + j*n] + 5*(n - 1) - 1]; // next ab = previous ac + previous bc
					if (ab == 0) {
						for (r = posL[C[i + j*n] - 1]; r < posR[C[i + j*n] - 1]; r++)
							for (c = posR[C[i + j*n] - 1]; c < posE[C[i + j*n] - 1]; c++)
								ab += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
					}
					subDists[i + j*(n - 1)] -= 2*ab;
					subDists[C[i + j*n] + 2*(n - 1) - 1] = ab; // ab
					if (C[i] > 0 && C[i + n] > 0) { // two branches
						if (j == 0) {
							subDists[C[i + j*n] + 3*(n - 1) - 1] = LLRL + LLRR; // ac
							subDists[C[i + j*n] + 4*(n - 1) - 1] = LRRL + LRRR; // bc
						} else {
							subDists[C[i + j*n] + 3*(n - 1) - 1] = LLRL + LRRL; // ac
							subDists[C[i + j*n] + 4*(n - 1) - 1] = LLRR + LRRR; // bc
						}
					} else { // one branch
						if (j == 0) { // c is right
							for (r = posL[C[i + j*n] - 1]; r < posR[C[i + j*n] - 1]; r++)
								for (c = posR[i]; c < posE[i]; c++)
									subDists[C[i + j*n] + 3*(n - 1) - 1] += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]); // ac
							for (r = posR[C[i + j*n] - 1]; r < posE[C[i + j*n] - 1]; r++)
								for (c = posR[i]; c < posE[i]; c++)
									subDists[C[i + j*n] + 4*(n - 1) - 1] += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]); // bc
						} else { // c is left
							for (r = posL[C[i + j*n] - 1]; r < posR[C[i + j*n] - 1]; r++)
								for (c = posL[i]; c < posR[i]; c++)
									subDists[C[i + j*n] + 3*(n - 1) - 1] += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]); // ac
							for (r = posR[C[i + j*n] - 1]; r < posE[C[i + j*n] - 1]; r++)
								for (c = posL[i]; c < posR[i]; c++)
									subDists[C[i + j*n] + 4*(n - 1) - 1] += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]); // bc
						}
					}
					subDists[i + 5*(n - 1)] = subDists[C[i + j*n] + 3*(n - 1) - 1] + subDists[C[i + j*n] + 4*(n - 1) - 1]; // ac + bc (equal in both directions)
				} else {
					subDists[i + j*(n - 1)] = sumDist[-C[i + j*n] - 1];
				}
			}
			
			if (flag == 1) { // branch lengths
				for (j = 0; j < 2; j++) {
					if (C[i + j*n] > 0) { // branch
						double ab = subDists[C[i + j*n] + 2*(n - 1) - 1];
						double ac = subDists[C[i + j*n] + 3*(n - 1) - 1];
						double bc = subDists[C[i + j*n] + 4*(n - 1) - 1];
						
						double ad = subDists[C[i + j*n] + 0*(n - 1) - 1] - ab - ac;
						double bd = subDists[C[i + j*n] + 1*(n - 1) - 1] - ab - bc;
						double cd = subDists[i + (1 - j)*(n - 1)] - ac - bc;
						
						rans[i + j*n] = 2*(ac + bd + ad + bc - ab - cd);
					}
				}
			} else { // nearest neighbor interchanges (NNIs)
				for (j = 0; j < 2; j++) {
					if (C[i + j*n] > 0) { // branch
						double ab = subDists[C[i + j*n] + 2*(n - 1) - 1];
						double ac = subDists[C[i + j*n] + 3*(n - 1) - 1];
						double bc = subDists[C[i + j*n] + 4*(n - 1) - 1];
						
						double ad = subDists[C[i + j*n] + 0*(n - 1) - 1] - ab - ac;
						double bd = subDists[C[i + j*n] + 1*(n - 1) - 1] - ab - bc;
						double cd = subDists[i + (1 - j)*(n - 1)] - ac - bc;
						
						double one = ab + cd;
						
						// swap down-left (a) with right (c)
						double two = ad + bc;
						double delta = one - 2*two;
						if (delta > 0)
							rans[i + j*n] = delta;
						
						// swap down-right (b) with right (c)
						two = ac + bd;
						delta = 2*two - one;
						if (delta < -rans[i + j*n])
							rans[i + j*n] = delta;
					}
				}
			}
		}
		
		// consider (virtual) root branch
		if (C[n - 1] > 0 && C[n + n - 1] > 0) { // left and right are branches
			double ab = subDists[C[n - 1] + 5*(n - 1) - 1];
			if (ab == 0) {
				for (r = posL[C[n - 1] - 1]; r < posR[C[n - 1] - 1]; r++)
					for (c = posR[C[n - 1] - 1]; c < posE[C[n - 1] - 1]; c++)
						ab += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
			}
			double cd = subDists[C[n + n - 1] + 5*(n - 1) - 1];
			if (cd == 0) {
				for (r = posL[C[n + n - 1] - 1]; r < posR[C[n + n - 1] - 1]; r++)
					for (c = posR[C[n + n - 1] - 1]; c < posE[C[n + n - 1] - 1]; c++)
						cd += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
			}
			double ac = 0; // left-left right-left
			for (r = posL[C[n - 1] - 1]; r < posR[C[n - 1] - 1]; r++)
				for (c = posL[C[n + n - 1] - 1]; c < posR[C[n + n - 1] - 1]; c++)
					ac += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
			double ad = 0; // left-left right-right
			for (r = posL[C[n - 1] - 1]; r < posR[C[n - 1] - 1]; r++)
				for (c = posR[C[n + n - 1] - 1]; c < posE[C[n + n - 1] - 1]; c++)
					ad += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
			double bd = 0; // left-right right-right
			for (r = posR[C[n - 1] - 1]; r < posE[C[n - 1] - 1]; r++)
				for (c = posR[C[n + n - 1] - 1]; c < posE[C[n + n - 1] - 1]; c++)
					bd += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
			double bc = 0; // left-right right-left
			for (r = posR[C[n - 1] - 1]; r < posE[C[n - 1] - 1]; r++)
				for (c = posL[C[n + n - 1] - 1]; c < posR[C[n + n - 1] - 1]; c++)
					bc += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
			
			if (flag == 1) { // branch lengths
				subDists[C[n - 1] + 2*(n - 1) - 1] = ab;
				subDists[C[n + n - 1] + 2*(n - 1) - 1] = cd;
				
				rans[n - 1] = ac + bd + ad + bc - ab - cd;
				if (rans[n - 1] < 0)
					rans[n - 1] = 0; // enforce strictly positive
				rans[n + n - 1] = rans[n - 1];
			} else { // nearest neighbor interchanges (NNIs)
				double one = ab + cd;
				
				// swap left-left (a) with right-left (c)
				double two = ad + bc;
				double delta = one - 2*two;
				if (delta > 0)
					rans[n - 1] = delta;
				
				// swap left-right (b) with right-left (c)
				two = ac + bd;
				delta = 2*two - one;
				if (delta < -rans[n - 1])
					rans[n - 1] = delta;
			}
		} else if (flag == 1) { // branch lengths
			if (C[n - 1] > 0) { // left is branch (therefore right is leaf)
				r = -C[n + n - 1] - 1; // a
				
				double ab = 0; // right left-left
				for (c = posL[C[n - 1] - 1]; c < posR[C[n - 1] - 1]; c++)
					ab += D[ind(r, index[c], N)]/pow(2, B[index1D(r, index[c], N)]);
				double ac = 0; // right left-right
				for (c = posR[C[n - 1] - 1]; c < posE[C[n - 1] - 1]; c++)
					ac += D[ind(r, index[c], N)]/pow(2, B[index1D(r, index[c], N)]);
				double bc = subDists[C[n - 1] + 5*(n - 1) - 1];
				if (bc == 0) {
					for (r = posL[C[n - 1] - 1]; r < posR[C[n - 1] - 1]; r++)
						for (c = posR[C[n - 1] - 1]; c < posE[C[n - 1] - 1]; c++)
							bc += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
				}
				subDists[C[n - 1] + 2*(n - 1) - 1] = bc;
				
				rans[n - 1] = ab + ac - bc;
			} else if (C[n + n - 1] > 0) { // right is branch (therefore left is leaf)
				r = -C[n - 1] - 1; // a
				
				double ab = 0; // right left-left
				for (c = posL[C[n + n - 1] - 1]; c < posR[C[n + n - 1] - 1]; c++)
					ab += D[ind(r, index[c], N)]/pow(2, B[index1D(r, index[c], N)]);
				double ac = 0; // right left-right
				for (c = posR[C[n + n - 1] - 1]; c < posE[C[n + n - 1] - 1]; c++)
					ac += D[ind(r, index[c], N)]/pow(2, B[index1D(r, index[c], N)]);
				double bc = subDists[C[n + n - 1] + 5*(n - 1) - 1];
				if (bc == 0) {
					for (r = posL[C[n + n - 1] - 1]; r < posR[C[n + n - 1] - 1]; r++)
						for (c = posR[C[n + n - 1] - 1]; c < posE[C[n + n - 1] - 1]; c++)
							bc += D[ind(index[r], index[c], N)]/pow(2, B[index1D(index[r], index[c], N)]);
				}
				subDists[C[n + n - 1] + 2*(n - 1) - 1] = bc;
				
				rans[n - 1] = ab + ac - bc;
			} else { // neither is a branch (therefore left and right are leaves)
				r = -C[n - 1] - 1;
				c = -C[n + n - 1] - 1;
				rans[n - 1] = D[ind(r, c, N)]/pow(2, B[index1D(r, c, N)])/2;
			}
			
			if (rans[n - 1] < 0)
				rans[n - 1] = 0; // enforce strictly positive
			rans[n + n - 1] = rans[n - 1];
		}
		
		if (flag == 1) { // branch lengths
			// set leaf lengths
			for (i = 0; i < n - 1; i++) {
				for (j = 0; j < 2; j++) {
					if (C[i + j*n] < 0) { // leaf
						double ab = subDists[i + 2*(n - 1)];
						double ac = subDists[i + j*(n - 1)] - ab;
						double bc = subDists[i + (1 - j)*(n - 1)] - ab;
						
						rans[i + j*n] = 2*(ab + ac - bc);
					}
					if (rans[i + j*n] < 0)
						rans[i + j*n] = 0; // enforce strictly positive
				}
			}
		}
		
		free(subDists);
		free(sumDist);
		free(B);
	}
	
	free(index);
	free(posL);
	free(posR);
	free(posE);
	
	UNPROTECT(1);
	
	return ans;
}

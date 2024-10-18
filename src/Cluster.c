/****************************************************************************
 *                         Cluster a Distance Matrix                        *
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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

void binUPGMA(double *rans, int i, int clusterNumber, double maxHeight, int length) {
	// NOTE:  cluster number is going to be the index i + 1
	if (rans[8*(length - 1) + i] == 0 || rans[9*(length - 1) + i] == 0) { // cluster number unassigned
		// then assign the node a cluster number
		if (rans[8*(length - 1) + i] == 0) { // the first fork is unassigned
			if (rans[6*(length - 1) + i] < 0) { // the first fork is a leaf
				rans[8*(length - 1) + i] = clusterNumber;
			} else { // the first fork is a branch
				rans[8*(length - 1) + i] = -1;
			}
		}
		if (rans[9*(length - 1) + i] == 0) { // the second fork is unassigned
			if (rans[7*(length - 1) + i] < 0) { // the second fork is a leaf
				rans[9*(length - 1) + i] = clusterNumber;
			} else { // the second fork is a branch
				rans[9*(length - 1) + i] = -1;
			}
		}
		
		// keep going up the tree
		for (int j = i + 1; j < length - 1; j++) {
			if (rans[6*(length - 1) + j] == (i + 1) || rans[7*(length - 1) + j] == (i + 1)) { // node is merged again
				if (rans[5*(length - 1) + j] <= maxHeight) { // if the node is within reach
					// then assign the same clusterNumber to this node
					binUPGMA(rans, j, clusterNumber, maxHeight, length);
					break;
				}
			}
		}
	}
	
	// follow the branches down the tree
	if (rans[6*(length - 1) + i] > 0) { // if first fork is a branch
		binUPGMA(rans, (int)(rans[6*(length - 1) + i] - 1), clusterNumber, maxHeight, length);
	}
	if (rans[7*(length - 1) + i] > 0) { // if second fork is a branch
		binUPGMA(rans, (int)(rans[7*(length - 1) + i] - 1), clusterNumber, maxHeight, length);
	}
}

void FollowBranch(double *rans, int i, double *branchLength, int length) {
	// add the longest branch length
	if (rans[8*(length - 1) + i] == 0) { // cluster number unassigned
		double alternative;
		if (rans[6*(length - 1) + i] < 0 && // first fork is a leaf and
			rans[7*(length - 1) + i] < 0) { // second fork is a leaf
			// add the longest one to branch length
			if (rans[3*(length - 1) + i] < rans[4*(length - 1) + i] && // second leaf is longest and
				rans[9*(length - 1) + i] == 0) { // second leaf is not assigned to a cluster
				// add second leaf to branch length
				*branchLength += rans[4*(length - 1) + i];
			} else {
				// add first leaf to branch length
				*branchLength += rans[3*(length - 1) + i];
			}
		} else if (rans[6*(length - 1) + i] > 0) { // first fork is a branch
			alternative = *branchLength + rans[4*(length - 1) + i];
			*branchLength += rans[3*(length - 1) + i];
			FollowBranch(rans, (int)(rans[6*(length - 1) + i] - 1), branchLength, length);
			if (*branchLength < alternative) {
				*branchLength = alternative;
			}
		} else if (rans[7*(length - 1) + i] > 0) { // second fork is a branch
			alternative = *branchLength + rans[3*(length - 1) + i];
			*branchLength += rans[4*(length - 1) + i];
			FollowBranch(rans, (int)(rans[7*(length - 1) + i] - 1), branchLength, length);
			if (*branchLength < alternative) {
				*branchLength = alternative;
			}
		}
	}
}

void assignNumber(double *rans, int i, int clusterNumber, double maxHeight, double floorHeight, int length) {
	// NOTE:  cluster number is going to be the index i + 1
	if (rans[8*(length - 1) + i] == 0 || rans[9*(length - 1) + i] == 0) { // cluster number unassigned
		// then assign the node a cluster number
		if (rans[8*(length - 1) + i] == 0) { // the first fork is unassigned
			if (rans[6*(length - 1) + i] < 0) { // the first fork is a leaf
				rans[8*(length - 1) + i] = clusterNumber;
			} else { // the first fork is a branch
				rans[8*(length - 1) + i] = -1;
			}
		}
		if (rans[9*(length - 1) + i] == 0) { // the second fork is unassigned
			if (rans[7*(length - 1) + i] < 0) { // the second fork is a leaf
				rans[9*(length - 1) + i] = clusterNumber;
			} else { // the second fork is a branch
				rans[9*(length - 1) + i] = -1;
			}
		}
		
		// keep going up the tree
		double branchLength;
		for (int j = i + 1; j < length - 1; j++) {
			if (rans[6*(length - 1) + j] == (i + 1) || rans[7*(length - 1) + j] == (i + 1)) { // node is merged again
				branchLength = 0;
				FollowBranch(rans, j, &branchLength, length);
				if ((rans[5*(length - 1) + j] + branchLength) <= maxHeight) { // if the node is within reach
					// then assign the same clusterNumber to this node
					assignNumber(rans, j, clusterNumber, maxHeight, floorHeight, length);
					break;
				}
			}
		}
	}
	// follow the branches down the tree
	if (rans[6*(length - 1) + i] > 0) { // if first fork is a branch
		if (rans[5*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] >= floorHeight) { // condition to stop following branches
			// then recursively number the branch
			if (((rans[5*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] -
				  rans[3*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] >= floorHeight) || // first fork is within reach or
				 rans[8*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] != 0) && // the fork is already assigned and
				((rans[5*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] -
				  rans[4*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] >= floorHeight) || // second fork is within reach or
				 rans[9*(length - 1) + (int)(rans[6*(length - 1) + i] - 1)] != 0)) { // the fork is already assigned
				assignNumber(rans, (int)(rans[6*(length - 1) + i] - 1), clusterNumber, maxHeight, floorHeight, length);
			}
		}
	}
	if (rans[7*(length - 1) + i] > 0) { // if second fork is a branch
		if (rans[5*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] >= floorHeight) { // condition to stop following branches
			// then recursively number the branch
			if (((rans[5*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] -
				  rans[3*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] >= floorHeight) || // first fork is within reach or
				 rans[8*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] != 0) && // the fork is already assigned and
				((rans[5*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] -
				  rans[4*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] >= floorHeight) || // second fork is within reach or
				 rans[9*(length - 1) + (int)(rans[7*(length - 1) + i] - 1)] != 0)) { // the fork is already assigned
					assignNumber(rans, (int)(rans[7*(length - 1) + i] - 1), clusterNumber, maxHeight, floorHeight, length);
			}
		}
	}
}

void Offset(int i, double *rans, double *offset, int length) {
	for (int j = i + 1; j < length - 1; j++) {
		if (rans[6*(length - 1) + j] == (i + 1)) {
			*offset = *offset + rans[5*(length - 1) + j] - rans[5*(length - 1) + i] - rans[3*(length - 1) + j];
			Offset(j, rans, offset, length);
			break;
		}
		if (rans[7*(length - 1) + j] == (i + 1)) {
			*offset = *offset + rans[5*(length - 1) + j] - rans[5*(length - 1) + i] - rans[4*(length - 1) + j];
			Offset(j, rans, offset, length);
			break;
		}
	}
}

double LHM(int w1, double lp1, int w2, double lp2) {
	double m; // weighted harmonic mean
	double lkup[5] = {0, 0, 0.6931471805599452862268, 1.098612288668109782108, 1.386294361119890572454}; // ln(index)
	
	if (lp1 > lp2) {
		m = lkup[w2] + lp1 + log1p(exp(lp2 - lp1)*w1/w2);
	} else {
		m = lkup[w1] + lp2 + log1p(exp(lp1 - lp2)*w2/w1);
	}
	
	m = lkup[w1 + w2] + lp1 + lp2 - m;
	if (m > 0) // finite precision
		m = 0;
	
	return m;
}

SEXP cluster(SEXP x, SEXP cutoff, SEXP method, SEXP l, SEXP verbose, SEXP pBar, SEXP nThreads)
{	
	/*
	 * **** Input (x) ****
	 *
	 *    Dist Matrix
	 *   A B C D E F G
	 * A 0
	 * B x 0
	 * C x x 0
	 * D x x x 0
	 * E x x x x 0
	 * F x x x x x 0
	 * G x x x x x x 0
	 */
	
	/*
	 * **********   Output (ans)   **********
	 * rans[0*(length - 1) + k] // 1st col = row merged
	 * rans[1*(length - 1) + k] // 2nd col = col merged
	 * rans[2*(length - 1) + k] // 3rd col = cluster number
	 * rans[3*(length - 1) + k] // 4th col = row branch length
	 * rans[4*(length - 1) + k] // 5th col = col branch length
	 * rans[5*(length - 1) + k] // 6th col = height of merger
	 * rans[6*(length - 1) + k] // 7th col = alternative numbering
	 * rans[7*(length - 1) + k] // 8th col = alternative numbering
	 * rans[8*(length - 1) + k] // 9th col = cutoff numbering
	 * rans[9*(length - 1) + k] // 10th col = cutoff numbering
	 */
	
	// distanceMatrix is a pointer to x
	// clusterNum hold the number of clusters
	// dMatrix2 is a manipulatable copy of distance matrix
	// *rans is a pointer to the output ans
	// nDiv is the net divergence
	// single leaves are negative, clusters are positive
	// Cluster numbering notes:
	// starting at lowest leaf, climb the tree until over maxHeight
	// maxHeight is merge height plus cutoff minus longest branch
	// do not go down the tree past floorHeight
	
	int k, dobj, clusterNum, minRow, minCol, index, minC, met;
	R_xlen_t i, j, length, size;
	int before, v, *rPercentComplete;
	double soFar, total, minHeight, *cut, *rans, *distanceMatrix, minH;
	SEXP ans, percentComplete, utilsPackage;
	int nthreads = asInteger(nThreads);
	
	// initialize variables
	clusterNum = 0; // increments with each new cluster
	size = asInteger(l); // square distance matrix dimension
	if (size < 0) { // x is of class "dist"
		dobj = 1;
		size *= -1;
	} else { // x is of class "matrix"
		dobj = 0;
	}
	length = size; // doesn't decrease in size
	PROTECT(ans = allocMatrix(REALSXP, (size - 1), 10));
	rans = REAL(ans);
	distanceMatrix = REAL(x);
	double *dMatrix2 = (double *) R_alloc(size*(size - 1)/2, sizeof(double));
	double dist1, dist2;
	double maxval = 1;
	int foundNA = 0;
	double *dTemp = (double *) R_alloc(size - 1, sizeof(double));
	double *cumHeight = (double *) R_alloc(size - 1, sizeof(double));
	int *clusterNums = (int *) R_alloc(size - 1, sizeof(int));
	int *rowNums = (int *) R_alloc(size - 1, sizeof(int));
	int *colNums = (int *) R_alloc(size - 1, sizeof(int));
	R_xlen_t *rowIndices = (R_xlen_t *) R_alloc(size - 1, sizeof(R_xlen_t));
	R_xlen_t *colIndices = (R_xlen_t *) R_alloc(size - 1, sizeof(R_xlen_t));
	int *minCols;
	met = asInteger(method);
	if (met != 1) // non-NJ method
		minCols = (int *) R_alloc(size - 1, sizeof(int));
	cut = REAL(cutoff);
	v = asLogical(verbose);
	double *nDiv;
	
	if (v) { // initialize progress variables
		soFar = 0;
		before = 0;
		total = length*(length - 1);
		PROTECT(percentComplete = NEW_INTEGER(1));
		rPercentComplete = INTEGER(percentComplete);
		// make it possible to access R functions from the utils package for the progress bar
		PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
	}
	
	if (dobj) { // x is of class "dist"
		for (i = 0; i < size*(size - 1)/2; i++) {
			dMatrix2[i] = *(distanceMatrix + i);
			if (!R_FINITE(dMatrix2[i])) {
				foundNA = 1;
			} else if (dMatrix2[i] > maxval) {
				maxval = dMatrix2[i];
			}
		}
	} else { // x is of class "matrix"
		// copy the lower triangle of the distance matrix (x)
		for (i = 1; i < length; i++) {
			for (j = 0; j < i; j++) {
				dMatrix2[length*j - j*(j + 1)/2 + i - j - 1] = *(distanceMatrix + j*length + i);
				if (!R_FINITE(dMatrix2[length*j - j*(j + 1)/2 + i - j - 1])) {
					foundNA = 1;
				} else if (dMatrix2[length*j - j*(j + 1)/2 + i - j - 1] > maxval) {
					maxval = dMatrix2[length*j - j*(j + 1)/2 + i - j - 1];
				}
			}
		}
	}
	
	// replace all non-finite values with maxval >= 1
	if (foundNA) {
		warning("Substituting %1.2f for non-finite values in myDistMatrix.", maxval);
		for (i = 0; i < size*(size - 1)/2; i++) {
			if (!R_FINITE(dMatrix2[i])) {
				dMatrix2[i] = maxval;
			}
		}
	}
	
	// initialize the order of groupings in the matrix
	for (i = 0; i < (size - 1); i++) {
		*(rowNums + i) = -(i + 2);
		*(colNums + i) = -(i + 1);
		*(rowIndices + i) = i;
		*(colIndices + i) = i;
		if (met != 1) // non-NJ method
			*(minCols + i) = -1;
	}
	
	// initialize rans to zero representing incomplete answer
	for (i = 0; i < (size - 1); i++) {
		for (j = 0; j < 10; j++) {
			*(rans + j*(length - 1) + i) = 0;
		}
	}
	
	if (met <= 1) { // NJ method
		nDiv = (double *) R_alloc(size, sizeof(double));
		
		for (i = 0; i < size; i++)
			nDiv[i] = 0;
		
		// calculate the net divergence
		for (i = 0; i < (size - 1); i++) {
			for (j = 0; j <= i; j++) {
				nDiv[j] += dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]]; // col sums
				nDiv[i + 1] += dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]]; // row sums
			}
		}
	}
	
	// start the loop that goes from tree leaf to root
	for (k = 0; k < (length - 1); k++) {
		// calculate the Q matrix & find the smallest element in the Q matrix
		minRow = 0;
		minCol = 0;
		minHeight = 1e50;
		if (nthreads == 1) {
			if (met < 1) { // heuristic NJ method
				for (i = (size - 2); i >= 0; i--) {
					if (minCols[rowIndices[i]] < 0) { // need to find minimum column
						minH = 1e50;
						for (j = i; j >= 0; j--) {
							if (dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2) < minH) {
								minH = dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2);
								minC = j;
							}
						}
						minCols[rowIndices[i]] = minC;
					} else { // minimum column is unchanged
						minC = minCols[rowIndices[i]];
						minH = dMatrix2[length*colIndices[minC] - colIndices[minC]*(colIndices[minC] + 1)/2 + rowIndices[i] - colIndices[minC]] - (nDiv[i + 1] + nDiv[minC])/(size - 2);
					}
					if (minH < minHeight) {
						minHeight = minH;
						minRow = i;
						minCol = minC;
					}
				}
			} else if (met == 1) { // NJ method
				for (i = (size - 2); i >= 0; i--) {
					minH = 1e50;
					for (j = i; j >= 0; j--) {
						if (dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2) < minH) {
							minH = dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2);
							minC = j;
						}
					}
					if (minH < minHeight) {
						minHeight = minH;
						minRow = i;
						minCol = minC;
					}
				}
			} else {
				for (i = (size - 2); i >= 0; i--) {
					if (minCols[rowIndices[i]] < 0) { // need to find minimum column
						minH = 1e50;
						for (j = i; j >= 0; j--) {
							if (dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] < minH) {
								minH = dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]];
								minC = j;
							}
						}
						minCols[rowIndices[i]] = minC;
					} else { // minimum column is unchanged
						minC = minCols[rowIndices[i]];
						minH = dMatrix2[length*colIndices[minC] - colIndices[minC]*(colIndices[minC] + 1)/2 + rowIndices[i] - colIndices[minC]];
					}
					if (minH < minHeight) {
						minHeight = minH;
						minRow = i;
						minCol = minC;
					}
				}
			}
		} else { // nthreads > 1
			if (met < 1) { // heuristic NJ method
				double *minHs = R_Calloc(size - 1, double); // initialized to zero
				int *minCs = R_Calloc(size - 1, int); // initialized to zero
				#pragma omp parallel for private(i,j,minC,minH) schedule(guided) num_threads(nthreads)
				for (i = (size - 2); i >= 0; i--) {
					if (minCols[rowIndices[i]] < 0) { // need to find minimum column
						minH = 1e50;
						for (j = i; j >= 0; j--) {
							if (dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2) < minH) {
								minH = dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2);
								minC = j;
							}
						}
						minCols[rowIndices[i]] = minC;
					} else { // minimum column is unchanged
						minC = minCols[rowIndices[i]];
						minH = dMatrix2[length*colIndices[minC] - colIndices[minC]*(colIndices[minC] + 1)/2 + rowIndices[i] - colIndices[minC]] - (nDiv[i + 1] + nDiv[minC])/(size - 2);
					}
					minHs[i] = minH;
					minCs[i] = minC;
				}
				for (i = (size - 2); i >= 0; i--) {
					if (minHs[i] < minHeight) {
						minHeight = minHs[i];
						minRow = i;
						minCol = minCs[i];
					}
				}
				R_Free(minHs);
				R_Free(minCs);
			} else if (met == 1) { // NJ method
				double *minHs = R_Calloc(size - 1, double); // initialized to zero
				int *minCs = R_Calloc(size - 1, int); // initialized to zero
				#pragma omp parallel for private(i,j,minC,minH) schedule(guided) num_threads(nthreads)
				for (i = (size - 2); i >= 0; i--) {
					minH = 1e50;
					for (j = i; j >= 0; j--) {
						if (dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2) < minH) {
							minH = dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] - (nDiv[i + 1] + nDiv[j])/(size - 2);
							minC = j;
						}
					}
					minHs[i] = minH;
					minCs[i] = minC;
				}
				for (i = (size - 2); i >= 0; i--) {
					if (minHs[i] < minHeight) {
						minHeight = minHs[i];
						minRow = i;
						minCol = minCs[i];
					}
				}
				R_Free(minHs);
				R_Free(minCs);
			} else {
				double *minHs = R_Calloc(size - 1, double); // initialized to zero
				int *minCs = R_Calloc(size - 1, int); // initialized to zero
				#pragma omp parallel for private(i,j,minC,minH) schedule(guided) num_threads(nthreads)
				for (i = (size - 2); i >= 0; i--) {
					if (minCols[rowIndices[i]] < 0) { // need to find minimum column
						minH = 1e50;
						for (j = i; j >= 0; j--) {
							if (dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]] < minH) {
								minH = dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]];
								minC = j;
							}
						}
						minCols[rowIndices[i]] = minC;
					} else { // minimum column is unchanged
						minC = minCols[rowIndices[i]];
						minH = dMatrix2[length*colIndices[minC] - colIndices[minC]*(colIndices[minC] + 1)/2 + rowIndices[i] - colIndices[minC]];
					}
					minHs[i] = minH;
					minCs[i] = minC;
				}
				for (i = (size - 2); i >= 0; i--) {
					if (minHs[i] < minHeight) {
						minHeight = minHs[i];
						minRow = i;
						minCol = minCs[i];
					}
				}
				R_Free(minHs);
				R_Free(minCs);
			}
		}
		
		// merge into a cluster
		rans[0*(length - 1) + k] = *(rowNums + rowIndices[minRow]); // row merged
		rans[1*(length - 1) + k] = *(colNums + colIndices[minCol]); // column merged
		
		// cluster
		if (((rans[0*(length - 1) + k] < 0) && (rans[1*(length - 1) + k] < 0)) ||
			((rans[0*(length - 1) + k] > 0) && (rans[1*(length - 1) + k] > 0))) {
			// merge into a new cluster
			clusterNum++;
			rans[2*(length - 1) + k] = clusterNum; // cluster formed
			// calculate both branch lengths
			if (met <= 1) { // NJ method
				if ((size - 2) == 0) { // case of (0/0 == NaN)
					rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2; // col
					rans[3*(length - 1) + k] = rans[4*(length - 1) + k]; // row
				} else {
					rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2 + (nDiv[minCol] - nDiv[minRow + 1])/(2*(size-2)); // col
					rans[3*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]] - rans[4*(length - 1) + k]; // row
				}
				
				// transfer negative branch lengths to positive branch lengths
				if (rans[3*(length - 1) + k] < 0 && rans[4*(length - 1) + k] < 0) {
					rans[3*(length - 1) + k] = -1*rans[3*(length - 1) + k];
					rans[4*(length - 1) + k] = -1*rans[4*(length - 1) + k];
				} else if (rans[4*(length - 1) + k] < 0) {
					rans[3*(length - 1) + k] += -1*rans[4*(length - 1) + k]; // add difference to other branch
					rans[4*(length - 1) + k] = 0;
				} else if (rans[3*(length - 1) + k] < 0) {
					rans[4*(length - 1) + k] += -1*rans[3*(length - 1) + k]; // add difference to other branch
					rans[3*(length - 1) + k] = 0;
				}
			} else if (met == 8) { // UPGMH
				rans[4*(length - 1) + k] = exp(dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]])/2; // col
				rans[3*(length - 1) + k] = rans[4*(length - 1) + k]; // row
			} else {
				rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2; // col
				rans[3*(length - 1) + k] = rans[4*(length - 1) + k]; // row
			}
			
			// add longest path to cumulative height
			if ((rans[0*(length - 1) + k] < 0) && (rans[1*(length - 1) + k] < 0)) {
				if (met <= 1) { // NJ method
					if (rans[4*(length - 1) + k] > rans[3*(length - 1) + k]) {
						cumHeight[clusterNum - 1] = rans[4*(length - 1) + k];
					} else {
						cumHeight[clusterNum - 1] = rans[3*(length - 1) + k];
					}
				} else {
					cumHeight[clusterNum - 1] = rans[4*(length - 1) + k];
				}
				
				// alternative cluster numbering
				rans[6*(length - 1) + k] = rans[0*(length - 1) + k];
				rans[7*(length - 1) + k] = rans[1*(length - 1) + k];
			} else {
				if (met <= 1) { // NJ method
					if ((cumHeight[(int)rans[0*(length - 1) + k] - 1] + rans[3*(length - 1) + k]) >
						(cumHeight[(int)rans[1*(length - 1) + k] - 1] + rans[4*(length - 1) + k])) {
						cumHeight[clusterNum - 1] = cumHeight[(int)rans[0*(length - 1) + k] - 1] + rans[3*(length - 1) + k];
					} else {
						cumHeight[clusterNum - 1] = cumHeight[(int)rans[1*(length - 1) + k] - 1] + rans[4*(length - 1) + k];
					}
				} else {
					cumHeight[clusterNum - 1] = rans[3*(length - 1) + k];
					rans[4*(length - 1) + k] -= cumHeight[(int)rans[0*(length - 1) + k] - 1]; // col
					rans[3*(length - 1) + k] -= cumHeight[(int)rans[1*(length - 1) + k] - 1]; // row
				}
				
				// alternative cluster numbering
				rans[6*(length - 1) + k] = clusterNums[(int)rans[0*(length - 1) + k] - 1];
				rans[7*(length - 1) + k] = clusterNums[(int)rans[1*(length - 1) + k] - 1];
			}
			clusterNums[clusterNum - 1] = k + 1;
			rans[5*(length - 1) + k] = cumHeight[clusterNum - 1];
		} else if (rans[0*(length - 1) + k] > 0) {
			// row is a cluster from before
			rans[2*(length - 1) + k] = rans[0*(length - 1) + k]; // merge with previous cluster
			// calculate both branch lengths
			if (met <= 1) { // NJ method
				if ((size - 2) == 0) { // case of (0/0 == NaN)
					rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]; // col
					rans[3*(length - 1) + k] = 0; // row
				} else {
					rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2 + (nDiv[minCol] - nDiv[minRow + 1])/(2*(size-2)); // col
					rans[3*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]] - rans[4*(length - 1) + k]; // row
				}
				
				// zero negative branch lengths
				if (rans[3*(length - 1) + k] < 0 && rans[4*(length - 1) + k] < 0) {
					rans[3*(length - 1) + k] = 0;
					rans[4*(length - 1) + k] = 0;
				} else if (rans[4*(length - 1) + k] < 0) {
					rans[3*(length - 1) + k] += -1*rans[4*(length - 1) + k]; // add difference to other branch
					rans[4*(length - 1) + k] = 0;
				} else if (rans[3*(length - 1) + k] < 0) {
					rans[4*(length - 1) + k] += -1*rans[3*(length - 1) + k]; // add difference to other branch
					rans[3*(length - 1) + k] = 0;
				}
			} else if (met == 8) { // UPGMH
				rans[4*(length - 1) + k] = exp(dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]])/2; // col
				rans[3*(length - 1) + k] = exp(dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]])/2 - cumHeight[(int)rans[0*(length - 1) + k] - 1]; // row
			} else {
				rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2; // col
				rans[3*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2 - cumHeight[(int)rans[0*(length - 1) + k] - 1]; // row
			}
			
			// add longest path to cumulative height
			cumHeight[(int)rans[0*(length - 1) + k] - 1] += rans[3*(length - 1) + k];
			if (rans[4*(length - 1) + k] > cumHeight[(int)rans[0*(length - 1) + k] - 1]) {
				cumHeight[(int)rans[0*(length - 1) + k] - 1] = rans[4*(length - 1) + k];
			}
			rans[5*(length - 1) + k] = cumHeight[(int)rans[0*(length - 1) + k] - 1];
			
			// alternative cluster numbering
			rans[6*(length - 1) + k] = clusterNums[(int)rans[0*(length - 1) + k] - 1];
			rans[7*(length - 1) + k] = rans[1*(length - 1) + k];
			clusterNums[(int)rans[0*(length - 1) + k] - 1] = k + 1;
		} else if (rans[1*(length - 1) + k] > 0) {
			rans[2*(length - 1) + k] = rans[1*(length - 1) + k]; // merge with previous cluster
			// calculate both branch lengths
			if (met <= 1) { // NJ method
				if ((size - 2) == 0) { // case of (0/0 == NaN)
					// unclear what to do in this case - splitting branch (edge) lengths as compromise
					rans[4*(length - 1) + k] = 0; // col
					rans[3*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]; // row
				} else {
					rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2 + (nDiv[minCol] - nDiv[minRow + 1])/(2*(size-2)); // col
					rans[3*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]] - rans[4*(length - 1) + k]; // row
				}
				
				// zero negative branch lengths
				if (rans[3*(length - 1) + k] < 0 && rans[4*(length - 1) + k] < 0) {
					rans[3*(length - 1) + k] = 0;
					rans[4*(length - 1) + k] = 0;
				} else if (rans[4*(length - 1) + k] < 0) {
					rans[3*(length - 1) + k] += -1*rans[4*(length - 1) + k]; // add difference to other branch
					rans[4*(length - 1) + k] = 0;
				} else if (rans[3*(length - 1) + k] < 0) {
					rans[4*(length - 1) + k] += -1*rans[3*(length - 1) + k]; // add difference to other branch
					rans[3*(length - 1) + k] = 0;
				}
			} else if (met == 8) { // UPGMH
				rans[4*(length - 1) + k] = exp(dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]])/2 - cumHeight[(int)rans[1*(length - 1) + k] - 1]; // col
				rans[3*(length - 1) + k] = exp(dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]])/2; // row
			} else {
				rans[4*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2 - cumHeight[(int)rans[1*(length - 1) + k] - 1]; // col
				rans[3*(length - 1) + k] = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]]/2; // row
			}
			
			// add longest path to cumulative height
			cumHeight[(int)rans[1*(length - 1) + k] - 1] += rans[4*(length - 1) + k];
			if (rans[3*(length - 1) + k] > cumHeight[(int)rans[1*(length - 1) + k] - 1]) {
				cumHeight[(int)rans[1*(length - 1) + k] - 1] = rans[3*(length - 1) + k];
			}
			rans[5*(length - 1) + k] = cumHeight[(int)rans[1*(length - 1) + k] - 1];
			
			// alternative cluster numbering
			rans[7*(length - 1) + k] = clusterNums[(int)rans[1*(length - 1) + k] - 1];
			rans[6*(length - 1) + k] = rans[0*(length - 1) + k];
			clusterNums[(int)rans[1*(length - 1) + k] - 1] = k + 1;
		}
		
		// calculate distances to the new node/cluster
		index = 0;
		if (met == 4) { // complete
			for (i = -1; i < (size - 1); i++) {
				if (!(i == minRow) && !(i == minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dist1 = dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[minRow] - colIndices[(i + 1)]];
					} else {
						dist1 = dMatrix2[length*colIndices[(minRow + 1)] - colIndices[(minRow + 1)]*(colIndices[(minRow + 1)] + 1)/2 + rowIndices[i] - colIndices[(minRow + 1)]];
					}
					if (i >= minCol) {
						dist2 = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[i] - colIndices[minCol]];
					} else {
						dist2 = dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[(i + 1)]];
					}
					if (dist1 > dist2) { // pick max distance
						dTemp[index] = dist1;
					} else {
						dTemp[index] = dist2;
					}
					index++;
				}
			}
		} else if (met == 5) { // single
			for (i = -1; i < (size - 1); i++) {
				if (!(i == minRow) && !(i == minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dist1 = dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[minRow] - colIndices[(i + 1)]];
					} else {
						dist1 = dMatrix2[length*colIndices[(minRow + 1)] - colIndices[(minRow + 1)]*(colIndices[(minRow + 1)] + 1)/2 + rowIndices[i] - colIndices[(minRow + 1)]];
					}
					if (i >= minCol) {
						dist2 = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[i] - colIndices[minCol]];
					} else {
						dist2 = dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[(i + 1)]];
					}
					if (dist1 < dist2) { // pick min distance
						dTemp[index] = dist1;
					} else {
						dTemp[index] = dist2;
					}
					index++;
				}
			}
		} else if (met == 2) { // UPGMA
			int weight1, weight2;
			if (*(rowNums + rowIndices[minRow]) < 0) {
				weight1 = 1;
			} else {
				weight1 = 2;
			}
			if (*(colNums + colIndices[minCol]) < 0) {
				weight2 = 1;
			} else {
				weight2 = 2;
			}
			for (i = -1; i < (size - 1); i++) {
				if (!(i == minRow) && !(i == minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dTemp[index] += dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[minRow] - colIndices[(i + 1)]]*weight1;
					} else {
						dTemp[index] += dMatrix2[length*colIndices[(minRow + 1)] - colIndices[(minRow + 1)]*(colIndices[(minRow + 1)] + 1)/2 + rowIndices[i] - colIndices[(minRow + 1)]]*weight1;
					}
					if (i >= minCol) {
						dTemp[index] += dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[i] - colIndices[minCol]]*weight2;
					} else {
						dTemp[index] += dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[(i + 1)]]*weight2;
					}
					dTemp[index] /= (weight1 + weight2); // average distance
					index++;
				}
			}
		} else if (met == 6) { // WPGMA
			for (i = -1; i < (size - 1); i++) {
				if (!(i == minRow) && !(i == minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dTemp[index] += dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[minRow] - colIndices[(i + 1)]];
					} else {
						dTemp[index] += dMatrix2[length*colIndices[(minRow + 1)] - colIndices[(minRow + 1)]*(colIndices[(minRow + 1)] + 1)/2 + rowIndices[i] - colIndices[(minRow + 1)]];
					}
					if (i >= minCol) {
						dTemp[index] += dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[i] - colIndices[minCol]];
					} else {
						dTemp[index] += dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[(i + 1)]];
					}
					dTemp[index] /= 2; // average distance
					index++;
				}
			}
		} else if (met == 8) { // UPGMH
			int weight1, weight2;
			if (*(rowNums + rowIndices[minRow]) < 0) {
				weight1 = 1;
			} else {
				weight1 = 2;
			}
			if (*(colNums + colIndices[minCol]) < 0) {
				weight2 = 1;
			} else {
				weight2 = 2;
			}
			for (i = -1; i < (size - 1); i++) {
				if (!(i == minRow) && !(i == minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dist1 = dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[minRow] - colIndices[(i + 1)]];
					} else {
						dist1 = dMatrix2[length*colIndices[(minRow + 1)] - colIndices[(minRow + 1)]*(colIndices[(minRow + 1)] + 1)/2 + rowIndices[i] - colIndices[(minRow + 1)]];
					}
					if (i >= minCol) {
						dist2 = dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[i] - colIndices[minCol]];
					} else {
						dist2 = dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[(i + 1)]];
					}
					dTemp[index] = LHM(weight1, dist1, weight2, dist2); // weighted harmonic mean of distance
					index++;
				}
			}
		} else { // NJ (met <= 1)
			for (i = -1; i < (size - 1); i++) {
				if (!(i == minRow) && !(i == minCol - 1)) {
					dTemp[index] = 0;
					// calculate distance from the new node
					if (minRow >= i) {
						dTemp[index] += dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[minRow] - colIndices[(i + 1)]];
					} else {
						dTemp[index] += dMatrix2[length*colIndices[(minRow + 1)] - colIndices[(minRow + 1)]*(colIndices[(minRow + 1)] + 1)/2 + rowIndices[i] - colIndices[(minRow + 1)]];
					}
					if (i >= minCol) {
						dTemp[index] += dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[i] - colIndices[minCol]];
					} else {
						dTemp[index] += dMatrix2[length*colIndices[(i + 1)] - colIndices[(i + 1)]*(colIndices[(i + 1)] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[(i + 1)]];
					}
					dTemp[index] -= dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[minRow] - colIndices[minCol]];
					dTemp[index] /= 2;
					index++;
				}
			}
		}
		
		// subtract net divergence of removed row
		if (met <= 1) { // NJ method
			// i = minRow
			for (j = 0; j <= minRow; j++) {
				nDiv[j] -= dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[minRow] - colIndices[j]]; // col sums
				nDiv[minRow + 1] -= dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[minRow] - colIndices[j]]; // row sums
			}
			for (i = j; i < (size - 1); i++) {
				nDiv[j] -= dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]]; // col sums
				nDiv[i + 1] -= dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[i] - colIndices[j]]; // row sums
			}
		}
		
		// move each index left after minRow
		for (i = minRow + 1; i < (size - 1); i++) {
			rowIndices[i - 1] = rowIndices[i];
			if (i > minRow + 1)
				colIndices[i - 1] = colIndices[i];
			
			if (met <= 1) { // NJ method
				nDiv[i] = nDiv[i + 1];
			}
			if (met != 1) { // non-NJ method
				// shift pointers left after minRow + 1
				if (minCols[rowIndices[i]] == (minRow + 1)) {
					minCols[rowIndices[i]] = -1; // eliminate
				} else if (minCols[rowIndices[i]] > (minRow + 1)) {
					minCols[rowIndices[i]]--; // shift left
				}
			}
		}
		
		// give the cluster its new number
		*(colNums + colIndices[minCol]) = rans[2*(length - 1) + k];
		if ((minCol - 1) >= 0) {
			*(rowNums + rowIndices[minCol - 1]) = rans[2*(length - 1) + k];
			if (met != 1) // non-NJ method
				minCols[rowIndices[minCol - 1]] = -1;
		}
		
		// decrement size of the matricies; length remains constant
		size--;
		
		// put new distances into the new cluster at minCol
		for (j = 0; j < (size - 1); j++) {
			if (j < minCol) {
				if (met <= 1) { // NJ method
					// remove original value
					nDiv[j] -= dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[j]]; // col sums
					nDiv[minCol] -= dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[j]]; // row sums
				}
				dMatrix2[length*colIndices[j] - colIndices[j]*(colIndices[j] + 1)/2 + rowIndices[(minCol - 1)] - colIndices[j]] = dTemp[j];
				if (met <= 1) { // NJ method
					// add new value
					nDiv[j] += dTemp[j]; // col sums
					nDiv[minCol] += dTemp[j]; // row sums
				}
			} else {
				if (met <= 1) { // NJ method
					// remove original value
					nDiv[minCol] -= dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[j] - colIndices[minCol]]; // col sums
					nDiv[j + 1] -= dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[j] - colIndices[minCol]]; // row sums
				}
				dMatrix2[length*colIndices[minCol] - colIndices[minCol]*(colIndices[minCol] + 1)/2 + rowIndices[j] - colIndices[minCol]] = dTemp[j];
				if (met <= 1) { // NJ method
					// add new value
					nDiv[minCol] += dTemp[j]; // col sums
					nDiv[j + 1] += dTemp[j]; // row sums
				}
				if (met < 1) { // heuristic NJ method
					if (minCols[rowIndices[j]] >= 0) {
						if (minCols[rowIndices[j]] == minCol) { // replacing self value
							minCols[rowIndices[j]] = -1; // re-evaluate row
						} else if (dTemp[j] < dMatrix2[length*colIndices[minCols[rowIndices[j]]] - colIndices[minCols[rowIndices[j]]]*(colIndices[minCols[rowIndices[j]]] + 1)/2 + rowIndices[j] - colIndices[minCols[rowIndices[j]]]] - (nDiv[j + 1] + nDiv[minCols[rowIndices[j]]])/(size - 2)) {
							minCols[rowIndices[j]] = minCol;
						}
					}
				} else if (met > 1) { // non-NJ method
					if (minCols[rowIndices[j]] >= 0) {
						if (minCols[rowIndices[j]] == minCol) { // replacing self value
							minCols[rowIndices[j]] = -1; // re-evaluate row
						} else if (dTemp[j] < dMatrix2[length*colIndices[minCols[rowIndices[j]]] - colIndices[minCols[rowIndices[j]]]*(colIndices[minCols[rowIndices[j]]] + 1)/2 + rowIndices[j] - colIndices[minCols[rowIndices[j]]]]) {
							minCols[rowIndices[j]] = minCol;
						}
					}
				}
			}
		}
		
		if (v) {
			// print the percent completed so far
			soFar = (2*length - 2 - k)*(k + 1);
			*rPercentComplete = floor(100*soFar/total);
			
			if (*rPercentComplete > before) { // when the percent has changed
				// tell the progress bar to update in the R console
				eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
				before = *rPercentComplete;
			}
		} else {
			R_CheckUserInterrupt();
		}
	}
	
	// bin sequences
	// assign cluster numbers to nodes within cutoff percent different
	int clusterNumber = 1;
	if (met <= 1) { // NJ method
		// insure that nodes are at the correct heights
		double offset;
		for (i = 0; i < length - 1; i++) {
			offset = 0;
			Offset(i, rans, &offset, length);
			rans[5*(length - 1) + i] = rans[5*(length - 1) + i] + offset;
		}
		
		double maxHeight, lowHeight, tempHeight;
		double longLeaf, longestLeaf, floorHeight;
		double tempSwap; // for swapping order of leaves
		int lowestBranch;
		for (i = 0; i < length; i++) {
			// find the lowest unassigned leaf
			lowestBranch = -1;
			lowHeight = 1e50;
			for (j = 0; j < length - 1; j++) {
				if (rans[8*(length - 1) + j] == 0 && // cluster number unassigned
					(rans[6*(length - 1) + j] < 0 || // first fork is a leaf or
					rans[7*(length - 1) + j] < 0)) { // second fork is a leaf
					if (rans[6*(length - 1) + j] < 0 && // first fork is a leaf and
						rans[7*(length - 1) + j] < 0) { // second fork is a leaf
						if (rans[3*(length - 1) + j] < rans[4*(length - 1) + j] && // second leaf is longest and
							rans[9*(length - 1) + j] == 0) { // second cluster number is unassigned
							tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
							longLeaf = rans[4*(length - 1) + j];
						} else {
							tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
							longLeaf = rans[3*(length - 1) + j];
						}
					} else if (rans[6*(length - 1) + j] < 0) { // first fork is a leaf
						tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
						longLeaf = rans[3*(length - 1) + j];
					} else if (rans[7*(length - 1) + j] < 0) { // second fork is a leaf
						tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
						longLeaf = rans[4*(length - 1) + j];
					}
					if (tempHeight < lowHeight) { // found a lower leaf
						lowestBranch = j;
						lowHeight = tempHeight;
						longestLeaf = longLeaf;
					}
				}
			}
			// assign a number to the lowest leaf
			if (lowestBranch >= 0) {
				// handle the special case where the two leaves are far apart
				if (rans[6*(length - 1) + lowestBranch] < 0 && // first fork is a leaf and
					rans[7*(length - 1) + lowestBranch] < 0 && // second fork is a leaf and
					rans[9*(length - 1) + lowestBranch] == 0 && // second cluster number is unassigned and
					((rans[3*(length - 1) + lowestBranch] + rans[4*(length - 1) + lowestBranch]) > *cut)) { // leaves farther apart then cutoff
					// then assign the longest leaf its own cluster number
					rans[9*(length - 1) + lowestBranch] = clusterNumber;
					if (rans[3*(length - 1) + lowestBranch] > rans[4*(length - 1) + lowestBranch]) { // first leaf is longest
						// then swap leaves
						tempSwap = rans[4*(length - 1) + lowestBranch];
						rans[4*(length - 1) + lowestBranch] = rans[3*(length - 1) + lowestBranch];
						rans[3*(length - 1) + lowestBranch] = tempSwap;
						tempSwap = rans[1*(length - 1) + lowestBranch];
						rans[1*(length - 1) + lowestBranch] = rans[0*(length - 1) + lowestBranch];
						rans[0*(length - 1) + lowestBranch] = tempSwap;
						tempSwap = rans[7*(length - 1) + lowestBranch];
						rans[7*(length - 1) + lowestBranch] = rans[6*(length - 1) + lowestBranch];
						rans[6*(length - 1) + lowestBranch] = tempSwap;
					}
				} else {
					// assign clusters going down the tree
					floorHeight = lowHeight - *cut + 2*longestLeaf;
					maxHeight = lowHeight + *cut;
					assignNumber(rans, lowestBranch, clusterNumber, maxHeight, floorHeight, length);
				}
				clusterNumber++;
			} else { // no leaves left to assign
				break;
			}
		}
	} else {
		for (i = 0; i < length - 1; i++) {
			if (rans[5*(length - 1) + i] > *cut/2 &&
				rans[8*(length - 1) + i] == 0 && // first fork is unassigned
				rans[6*(length - 1) + i] < 0 && // first fork is a leaf
				rans[9*(length - 1) + i] == 0 && // second fork is unassigned
				rans[7*(length - 1) + i] < 0) { // second fork is a leaf
				rans[8*(length - 1) + i] = clusterNumber;
				clusterNumber++;
				rans[9*(length - 1) + i] = clusterNumber;
				clusterNumber++;
			} else {
				if (rans[8*(length - 1) + i] == 0 && // first fork is unassigned
					rans[6*(length - 1) + i] < 0) { // first fork is a leaf
					binUPGMA(rans, i, clusterNumber, *cut/2, length);
					clusterNumber++;
				}
				if (rans[9*(length - 1) + i] == 0 && // second fork is unassigned
					rans[7*(length - 1) + i] < 0) { // second fork is a leaf
					binUPGMA(rans, i, clusterNumber, *cut/2, length);
					clusterNumber++;
				}
			}
		}
	}
	
	if (v) {
		UNPROTECT(3);
	} else {
		UNPROTECT(1);
	}
	
	return ans;
}

SEXP reclusterUPGMA(SEXP x, SEXP cutoff)
{
	int i;
	double *cut, *rans;
	cut = REAL(cutoff);
	SEXP ans;
	PROTECT(ans = duplicate(x));
	rans = REAL(ans);
	const int length = length(ans)/10 + 1; // number of rows
	
	// zero out previous clusters
	for (i = 0; i < length - 1; i++) {
		rans[8*(length - 1) + i] = 0;
		rans[9*(length - 1) + i] = 0;
	}
	
	// rebin sequences
	// assign cluster numbers to nodes within cutoff percent different
	int clusterNumber = 1;
	for (i = 0; i < length - 1; i++) {
		if (rans[5*(length - 1) + i] > *cut/2 &&
			rans[8*(length - 1) + i] == 0 && // first fork is unassigned
			rans[6*(length - 1) + i] < 0 && // first fork is a leaf
			rans[9*(length - 1) + i] == 0 && // second fork is unassigned
			rans[7*(length - 1) + i] < 0) { // second fork is a leaf
			rans[8*(length - 1) + i] = clusterNumber;
			clusterNumber++;
			rans[9*(length - 1) + i] = clusterNumber;
			clusterNumber++;
		} else {
			if (rans[8*(length - 1) + i] == 0 && // first fork is unassigned
				rans[6*(length - 1) + i] < 0) { // first fork is a leaf
				binUPGMA(rans, i, clusterNumber, *cut/2, length);
				clusterNumber++;
			}
			if (rans[9*(length - 1) + i] == 0 && // second fork is unassigned
				rans[7*(length - 1) + i] < 0) { // second fork is a leaf
				binUPGMA(rans, i, clusterNumber, *cut/2, length);
				clusterNumber++;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP reclusterNJ(SEXP x, SEXP cutoff)
{
	// bin sequences
	// assign cluster numbers to nodes within cutoff percent different
	int clusterNumber = 1;
	double maxHeight, lowHeight, tempHeight;
	double longLeaf, longestLeaf, floorHeight;
	double tempSwap; // for swapping order of leaves
	int lowestBranch;
	
	int i, j;
	double *cut, *rans;
	cut = REAL(cutoff);
	SEXP ans;
	PROTECT(ans = duplicate(x));
	rans = REAL(ans);
	const int length = length(ans)/10 + 1; // number of rows
	
	// zero out previous clusters
	for (i = 0; i < length - 1; i++) {
		rans[8*(length - 1) + i] = 0;
		rans[9*(length - 1) + i] = 0;
	}
	
	for (i = 0; i < length; i++) {
		// find the lowest unassigned leaf
		lowestBranch = -1;
		lowHeight = 1e50;
		for (j = 0; j < length - 1; j++) {
			if (rans[8*(length - 1) + j] == 0 && // cluster number unassigned
				(rans[6*(length - 1) + j] < 0 || // first fork is a leaf or
				 rans[7*(length - 1) + j] < 0)) { // second fork is a leaf
					if (rans[6*(length - 1) + j] < 0 && // first fork is a leaf and
						rans[7*(length - 1) + j] < 0) { // second fork is a leaf
						if (rans[3*(length - 1) + j] < rans[4*(length - 1) + j] && // second leaf is longest and
							rans[9*(length - 1) + j] == 0) { // second cluster number is unassigned
							tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
							longLeaf = rans[4*(length - 1) + j];
						} else {
							tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
							longLeaf = rans[3*(length - 1) + j];
						}
					} else if (rans[6*(length - 1) + j] < 0) { // first fork is a leaf
						tempHeight = rans[5*(length - 1) + j] - rans[3*(length - 1) + j];
						longLeaf = rans[3*(length - 1) + j];
					} else if (rans[7*(length - 1) + j] < 0) { // second fork is a leaf
						tempHeight = rans[5*(length - 1) + j] - rans[4*(length - 1) + j];
						longLeaf = rans[4*(length - 1) + j];
					}
					if (tempHeight < lowHeight) { // found a lower leaf
						lowestBranch = j;
						lowHeight = tempHeight;
						longestLeaf = longLeaf;
					}
				}
		}
		// assign a number to the lowest leaf
		if (lowestBranch >= 0) {
			// handle the special case where the two leaves are far apart
			if (rans[6*(length - 1) + lowestBranch] < 0 && // first fork is a leaf and
				rans[7*(length - 1) + lowestBranch] < 0 && // second fork is a leaf and
				rans[9*(length - 1) + lowestBranch] == 0 && // second cluster number is unassigned and
				((rans[3*(length - 1) + lowestBranch] + rans[4*(length - 1) + lowestBranch]) > *cut)) { // leaves farther apart then cutoff
				// then assign the longest leaf its own cluster number
				rans[9*(length - 1) + lowestBranch] = clusterNumber;
				if (rans[3*(length - 1) + lowestBranch] > rans[4*(length - 1) + lowestBranch]) { // first leaf is longest
					// then swap leaves
					tempSwap = rans[4*(length - 1) + lowestBranch];
					rans[4*(length - 1) + lowestBranch] = rans[3*(length - 1) + lowestBranch];
					rans[3*(length - 1) + lowestBranch] = tempSwap;
					tempSwap = rans[1*(length - 1) + lowestBranch];
					rans[1*(length - 1) + lowestBranch] = rans[0*(length - 1) + lowestBranch];
					rans[0*(length - 1) + lowestBranch] = tempSwap;
					tempSwap = rans[7*(length - 1) + lowestBranch];
					rans[7*(length - 1) + lowestBranch] = rans[6*(length - 1) + lowestBranch];
					rans[6*(length - 1) + lowestBranch] = tempSwap;
				}
			} else {
				// assign clusters going down the tree
				floorHeight = lowHeight - *cut + 2*longestLeaf;
				maxHeight = lowHeight + *cut;
				assignNumber(rans, lowestBranch, clusterNumber, maxHeight, floorHeight, length);
			}
			clusterNumber++;
		} else { // no leaves left to assign
			break;
		}
		
	}
	
	UNPROTECT(1);
	
	return ans;
}

SEXP adjustHeights(SEXP x)
{
	// insure that nodes are at the correct heights
	double offset, *rans;
	int length = length(x)/10 + 1;
	rans = REAL(x);
	
	for (int i = 0; i < length - 1; i++) {
		offset = 0;
		Offset(i, rans, &offset, length);
		rans[5*(length - 1) + i] = rans[5*(length - 1) + i] + offset;
	}
	
	return x;
}

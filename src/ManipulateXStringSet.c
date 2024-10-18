/****************************************************************************
 *                         Manipulate XStringSets                           *
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

/* for Calloc/Free */
#include <R_ext/RS.h>

// for math functions
#include <math.h>

// strcpy
#include <string.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

// insert gaps into aligned XStringSets
SEXP insertGaps(SEXP x, SEXP positions, SEXP lengths, SEXP type, SEXP nThreads)
{
	int i, j, x_length, *width, sum, start;
	SEXP ans_width, ans;
	int t = asInteger(type);
	int *p = INTEGER(positions);
	int *l = INTEGER(lengths);
	int n = length(positions);
	int nthreads = asInteger(nThreads);
	
	// determine the cumulative width of insertions
	sum = 0;
	for (i = 0; i < n; i++)
		sum += l[i];
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the length of the XStringSet
	XStringSet_holder x_set, ans_holder;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	// determine the widths of the aligned (equal width) XStringSet
	PROTECT(ans_width = NEW_INTEGER(x_length));
	Chars_holder x_s;
	x_s = get_elt_from_XStringSet_holder(&x_set, 0);
	for (i = 0, width = INTEGER(ans_width); i < x_length; i++, width++)
		*width = x_s.length + sum;
	
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
	
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,ans_elt_holder,x_s,sum,start) schedule(guided) num_threads(nthreads)
	#endif
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		//ans_elt_holder.length = 0;
		x_s = get_elt_from_XStringSet_holder(&x_set, i);
		sum = 0; // position in ans_elt_holder.ptr
		start = 0; // position in x_s.ptr
		for (j = 0; j < n; j++) {
			if ((p[j] - 1) > start) { // copy over sequence
				memcpy((char *) ans_elt_holder.ptr + sum, x_s.ptr + start, (p[j] - 1 - start) * sizeof(char));
				sum += (p[j] - 1 - start);
				start += (p[j] - 1 - start);
			}
			if (l[j] > 0) { // insert gaps
				if (t == 3) { // AAStringSet
					memset((char *) ans_elt_holder.ptr + sum, 45, l[j] * sizeof(char));
				} else { // DNAStringSet or RNAStringSet
					memset((char *) ans_elt_holder.ptr + sum, 16, l[j] * sizeof(char));
				}
				sum += l[j];
			}
		}
		if (sum < ans_elt_holder.length) {
			memcpy((char *) ans_elt_holder.ptr + sum, x_s.ptr + start, (ans_elt_holder.length - sum) * sizeof(char));
		}
	}
	
	UNPROTECT(2);
	return ans;
}

// remove common gaps in an XStringSet
SEXP commonGaps(SEXP x)
{
	int i, j, l, p, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	// initialize vector of positions
	int bits[longest];
	for (i = 0; i < longest; i++)
		bits[i] = 0;
	
	// determine which positions are not gaps
	p = 0; // number of non-gap positions
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		for (j = 0; j < l; j++) {
			if (bits[j] == 0) {
				if (!(seq[j] == '-' || seq[j] == '.')) {
					bits[j] = 1;
					p++;
				}
			}
		}
	}
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char s[p + 1]; // each sequence
	
	// write new character vector
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		count = 0;
		for (j = 0; j < l; j++) {
			if (bits[j] == 1) {
				s[count] = seq[j];
				count ++;
			}
		}
		s[count] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(s));
	}
	
	UNPROTECT(1);
	
	return seqs;	
}

// consolidate gaps in an XStringSet
SEXP consolidateGaps(SEXP x, SEXP type)
{
	int i, j, k, l, mp, x_length, count, c;
	int t = asInteger(type);
	char p;
	
	// determine the length of the XStringSet
	XStringSet_holder x_set;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	Chars_holder x_i;
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	l = x_i.length; // assume sequences are aligned
	
	// pos keeps track of last (non-gap) position
	int *pos = (int *) R_alloc(x_length, sizeof(int));
	for (i = 0; i < x_length; i++) {
		pos[i] = NA_INTEGER;
	}
	
	// w keeps track of which sequences are non-gap
	int *w = (int *) R_alloc(x_length, sizeof(int));
	
	for (i = 0; i < l; i++) {
		count = 0;
		mp = 0; // max pos
		for (j = 0; j < x_length; j++) {
			x_i = get_elt_from_XStringSet_holder(&x_set, j);
			
			if ((t == 3 && x_i.ptr[i] ^ 0x2D && x_i.ptr[i] ^ 0x2E) ||
				(t != 3 && !(x_i.ptr[i] & 0x10 || x_i.ptr[i] & 0x40))) {
				// not a gap ("-") or unknown (".") in this position
				w[count] = j;
				count++;
				
				if (pos[j] == NA_INTEGER) // pos is NA
					pos[j] = i;
				
				if (pos[j] > mp)
					mp = pos[j];
			}
		}
		
		if (count == 0) // column is all gaps
			continue;
		
		if ((i - mp) > 1) { // swap columns
			c = 0;
			for (j = 0; j < x_length; j++) {
				x_i = get_elt_from_XStringSet_holder(&x_set, j);
				
				if (c < count && j == w[c]) { // shift to left
					c++;
					p = x_i.ptr[mp + 1];
					*((char *)x_i.ptr + mp + 1) = *((char *)x_i.ptr + i);
					*((char *)x_i.ptr + i) = p;
					pos[j] = mp + 1;
				} else { // shift to right
					p = x_i.ptr[i];
					for (k = i; k >= (mp + 2); k--) {
						*((char *)x_i.ptr + k) = *((char *)x_i.ptr + k - 1);
					}
					*((char *)x_i.ptr + mp + 1) = p;
					
					if (pos[j] >= mp + 1)
						pos[j]++;
				}
			}
		} else {
			for (j = 0; j < count; j++)
				pos[w[j]] = i;
		}
	}
	
	return x;
}

// quickly replace characters
SEXP replaceChars(SEXP x, SEXP r, SEXP t)
{
	int i, j, l, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	const char *repChar = CHAR(STRING_ELT(r, 0));
	int fail = (STRING_ELT(r, 0) == NA_STRING);
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char *s = R_Calloc(longest + 1, char); // each sequence
	
	// write new character vector
	if (asInteger(t) == 1) {
		for (i = 0; i < n; i++) {
			l = length(STRING_ELT(x, i));
			seq = CHAR(STRING_ELT(x, i));
			count = 0;
			for (j = 0; j < l; j++) {
				if (seq[j] != 'U' && seq[j] != 'u') {
					switch (seq[j]) {
						case '-':
						case 'A':
						case 'a':
						case 'C':
						case 'c':
						case 'G':
						case 'g':
						case 'T':
						case 't':
						case 'N':
						case 'n':
						case 'M':
						case 'm':
						case 'R':
						case 'r':
						case 'W':
						case 'w':
						case 'S':
						case 's':
						case 'Y':
						case 'y':
						case 'K':
						case 'k':
						case 'V':
						case 'v':
						case 'H':
						case 'h':
						case 'D':
						case 'd':
						case 'B':
						case 'b':
						case '+':
						case '.':
							s[count] = seq[j];
							count++;
							break;
						default:
							if (fail) {
								error("Incompatible character ('%c') found in DNAStringSet when replaceChar = NA.", seq[j]);
							} else if (repChar[0] != '\0') {
								s[count] = repChar[0];
								count++;
							}
							break;
					}
				} else {
					s[count] = 'T';
					count++;
				}
			}
			s[count] = '\0'; // null-terminate
			SET_STRING_ELT(seqs, i, mkChar(s));
		}
	} else if (asInteger(t) == 2) {
		for (i = 0; i < n; i++) {
			l = length(STRING_ELT(x, i));
			seq = CHAR(STRING_ELT(x, i));
			count = 0;
			for (j = 0; j < l; j++) {
				if (seq[j] != 'T' && seq[j] != 't') {
					switch (seq[j]) {
						case '-':
						case 'A':
						case 'a':
						case 'C':
						case 'c':
						case 'G':
						case 'g':
						case 'U':
						case 'u':
						case 'N':
						case 'n':
						case 'M':
						case 'm':
						case 'R':
						case 'r':
						case 'W':
						case 'w':
						case 'S':
						case 's':
						case 'Y':
						case 'y':
						case 'K':
						case 'k':
						case 'V':
						case 'v':
						case 'H':
						case 'h':
						case 'D':
						case 'd':
						case 'B':
						case 'b':
						case '+':
						case '.':
							s[count] = seq[j];
							count++;
							break;
						default:
							if (fail) {
								error("Incompatible character ('%c') in RNAStringSet found when replaceChar = NA.", seq[j]);
							} else if (repChar[0] != '\0') {
								s[count] = repChar[0];
								count++;
							}
							break;
					}
				} else {
					s[count] = 'U';
					count++;
				}
			}
			s[count] = '\0'; // null-terminate
			SET_STRING_ELT(seqs, i, mkChar(s));
		}
	} else {
		for (i = 0; i < n; i++) {
			l = length(STRING_ELT(x, i));
			seq = CHAR(STRING_ELT(x, i));
			count = 0;
			for (j = 0; j < l; j++) {
				switch (seq[j]) {
					case '-':
					case 'A':
					case 'R':
					case 'N':
					case 'D':
					case 'C':
					case 'Q':
					case 'E':
					case 'G':
					case 'H':
					case 'I':
					case 'L':
					case 'K':
					case 'M':
					case 'F':
					case 'P':
					case 'S':
					case 'T':
					case 'W':
					case 'Y':
					case 'V':
					case 'U':
					case 'O':
					case 'B':
					case 'Z':
					case 'X':
					case '*':
					case '+':
					case '.':
						s[count] = seq[j];
						count++;
						break;
					case 'a':
						s[count] = 'A';
						count++;
						break;
					case 'r':
						s[count] = 'R';
						count++;
						break;
					case 'n':
						s[count] = 'N';
						count++;
						break;
					case 'd':
						s[count] = 'D';
						count++;
						break;
					case 'c':
						s[count] = 'C';
						count++;
						break;
					case 'q':
						s[count] = 'Q';
						count++;
						break;
					case 'e':
						s[count] = 'E';
						count++;
						break;
					case 'g':
						s[count] = 'G';
						count++;
						break;
					case 'h':
						s[count] = 'H';
						count++;
						break;
					case 'i':
						s[count] = 'I';
						count++;
						break;
					case 'l':
						s[count] = 'L';
						count++;
						break;
					case 'k':
						s[count] = 'K';
						count++;
						break;
					case 'm':
						s[count] = 'M';
						count++;
						break;
					case 'f':
						s[count] = 'F';
						count++;
						break;
					case 'p':
						s[count] = 'P';
						count++;
						break;
					case 's':
						s[count] = 'S';
						count++;
						break;
					case 't':
						s[count] = 'T';
						count++;
						break;
					case 'w':
						s[count] = 'W';
						count++;
						break;
					case 'y':
						s[count] = 'Y';
						count++;
						break;
					case 'v':
						s[count] = 'V';
						count++;
						break;
					case 'u':
						s[count] = 'U';
						count++;
						break;
					case 'o':
						s[count] = 'O';
						count++;
						break;
					case 'b':
						s[count] = 'B';
						count++;
						break;
					case 'j':
						s[count] = 'J';
						count++;
						break;
					case 'z':
						s[count] = 'Z';
						count++;
						break;
					case 'x':
						s[count] = 'X';
						count++;
						break;
					default:
						if (fail) {
							error("Incompatible character ('%c') in AAStringSet found when replaceChar = NA.", seq[j]);
						} else if (repChar[0] != '\0') {
							s[count] = repChar[0];
							count++;
						}
						break;
				}
			}
			s[count] = '\0'; // null-terminate
			SET_STRING_ELT(seqs, i, mkChar(s));
		}
	}
	
	R_Free(s);
	
	UNPROTECT(1);
	
	return seqs;
}

// replace characters in an XStringSet
SEXP replaceChar(SEXP x, SEXP c, SEXP r)
{
	int i, j, l, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	const char *repChar = CHAR(STRING_ELT(r, 0));
	const char *charRep = CHAR(STRING_ELT(c, 0));
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char *s = R_Calloc(longest + 1, char); // each sequence
	
	// write new character vector
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		count = 0;
		for (j = 0; j < l; j++) {
			if (seq[j] == charRep[0]) {
				if (repChar[0] != '\0') {
					s[count] = repChar[0];
					count++;
				}
			} else {
				s[count] = seq[j];
				count++;
			}
		}
		s[count] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(s));
	}
	
	R_Free(s);
	
	UNPROTECT(1);
	
	return seqs;
}

// replace Gaps in an XStringSet
SEXP replaceGaps(SEXP x, SEXP y, SEXP start, SEXP type)
{
	int i, j, x_length;
	SEXP ans_width, ans;
	int s = asInteger(start) - 1; // starting position in y
	int t = asInteger(type);
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the length of the XStringSet
	Chars_holder x_i, y_holder, ans_elt_holder;
	XStringSet_holder x_set, ans_holder;
	y_holder = hold_XRaw(y);
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	// count the sequence lengths
	PROTECT(ans_width = NEW_INTEGER(x_length));
	int *width = INTEGER(ans_width);
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		width[i] = x_i.length;
	}
	
	// set the class of the XStringSet
	char ans_classname[40];
	if (t == 1) {
		strcpy(ans_classname, "DNAStringSet");
	} else if (t == 2) {
		strcpy(ans_classname, "RNAStringSet");
	} else { // t == 3
		strcpy(ans_classname, "AAStringSet");
	}
	
	// initialize a new XStringSet
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (t == 3) { // AAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E)) { // position is a gap
					memcpy((char *) ans_elt_holder.ptr + j, x_i.ptr + j, 1 * sizeof(char));
				} else {
					memcpy((char *) ans_elt_holder.ptr + j, y_holder.ptr + s, 1 * sizeof(char));
					s++;
				}
			}
		} else { // DNAStringSet or RNAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40) { // position is a gap
					memcpy((char *) ans_elt_holder.ptr + j, x_i.ptr + j, 1 * sizeof(char));
				} else {
					memcpy((char *) ans_elt_holder.ptr + j, y_holder.ptr + s, 1 * sizeof(char));
					s++;
				}
			}
		}
	}
	
	UNPROTECT(2);
	return ans;
}

// remove common gaps in an XStringSet
SEXP removeCommonGaps(SEXP x, SEXP type, SEXP mask, SEXP nThreads)
{
	int i, j, k, l, w, x_length, *width, sum, start, delta;
	SEXP ans_width, ans;
	int t = asInteger(type);
	int m = asInteger(mask);
	int nthreads = asInteger(nThreads);
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the length of the XStringSet
	XStringSet_holder x_set, ans_holder;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	Chars_holder x_i;
	x_i = get_elt_from_XStringSet_holder(&x_set, 0);
	l = x_i.length;
	
	// initialize a vector of columns that are 100% gaps
	int *remove = (int *) R_alloc(l, sizeof(int));
	sum = 0; // number of columns to remove
	if (t == 3) { // AAStringSet
		for (i = 0; i < l; i++) {
			if (!(x_i.ptr[i] ^ 0x2D) || !(x_i.ptr[i] ^ 0x2E) || (m && !(x_i.ptr[i] ^ 0x2B))) { // position is a gap
				remove[sum] = i;
				sum++;
			}
		}
		
		// find columns to remove
		for (i = 1; i < x_length; i++) { // each sequence
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			if (x_i.length != l)
				error("Sequences are not equal length.");
			
			for (j = 0; j < sum; j++) {
				if (x_i.ptr[remove[j]] ^ 0x2D && x_i.ptr[remove[j]] ^ 0x2E && (!m || x_i.ptr[remove[j]] ^ 0x2B)) { // non-gap in this position
					// stop the position from being removed
					sum--;
					// shift all positions over by one
					for (k = j; k < sum; k++)
						remove[k] = remove[k + 1];
					j--;
				}
			}
		}
	} else { // DNAStringSet or RNAStringSet
		for (i = 0; i < l; i++) {
			if (x_i.ptr[i] & 0x10 || x_i.ptr[i] & 0x40 || (m && x_i.ptr[i] & 0x20)) { // position is a gap
				remove[sum] = i;
				sum++;
			}
		}
		
		// find columns to remove
		for (i = 1; i < x_length; i++) { // each sequence
			x_i = get_elt_from_XStringSet_holder(&x_set, i);
			if (x_i.length != l)
				error("Sequences are not equal length.");
			
			for (j = 0; j < sum; j++) {
				if (!(x_i.ptr[remove[j]] & 0x10 || x_i.ptr[remove[j]] & 0x40 || (m && x_i.ptr[remove[j]] & 0x20))) { // non-gap in this position
					// stop the position from being removed
					sum--;
					// shift all positions over by one
					for (k = j; k < sum; k++)
						remove[k] = remove[k + 1];
					j--;
				}
			}
		}
	}
	
	// determine the widths of the aligned (equal width) XStringSet
	PROTECT(ans_width = NEW_INTEGER(x_length));
	w = l - sum;
	for (i = 0, width = INTEGER(ans_width); i < x_length; i++, width++)
		*width = w;
	
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
	
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,ans_elt_holder,x_i,start,delta) schedule(guided) num_threads(nthreads)
	#endif
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		start = 0; // position in ans_elt_holder.ptr
		for (j = 0; j < sum; j++) {
			if (j > 0) {
				delta = remove[j] - remove[j - 1] - 1;
				if (delta > 0) {
					memcpy((char *) ans_elt_holder.ptr + start, x_i.ptr + remove[j - 1] + 1, delta * sizeof(char));
					start += delta;
				}
			} else {
				memcpy((char *) ans_elt_holder.ptr, x_i.ptr, remove[j] * sizeof(char));
				start = remove[j];
			}
		}
		if (start < w) {
			delta = w - start;
			memcpy((char *) ans_elt_holder.ptr + start, x_i.ptr + l - delta, delta * sizeof(char));
		}
	}
	
	UNPROTECT(2);
	return ans;
}

// remove all gaps in an XStringSet
SEXP removeGaps(SEXP x, SEXP type, SEXP mask, SEXP nThreads)
{
	int i, j, x_length, p, sum;
	SEXP ans_width, ans;
	int t = asInteger(type);
	int m = asInteger(mask);
	int nthreads = asInteger(nThreads);
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the length of the XStringSet
	Chars_holder x_i, ans_elt_holder;
	XStringSet_holder x_set, ans_holder;
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	// count the sequence lengths
	PROTECT(ans_width = NEW_INTEGER(x_length));
	int *width = INTEGER(ans_width);
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,x_i) schedule(guided) num_threads(nthreads)
	#endif
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		width[i] = x_i.length;
		if (t == 3) { // AAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E) || (m && !(x_i.ptr[j] ^ 0x2B))) { // position is a gap
					width[i]--;
				}
			}
		} else { // DNAStringSet or RNAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40 || (m && x_i.ptr[j] & 0x20)) { // position is a gap
					width[i]--;
				}
			}
		}
	}
	
	// set the class of the XStringSet
	char ans_classname[40];
	if (t == 1) {
		strcpy(ans_classname, "DNAStringSet");
	} else if (t == 2) {
		strcpy(ans_classname, "RNAStringSet");
	} else { // t == 3
		strcpy(ans_classname, "AAStringSet");
	}
	
	// initialize a new XStringSet
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	
	#ifdef _OPENMP
	#pragma omp parallel for private(i,j,ans_elt_holder,x_i,p,sum) schedule(guided) num_threads(nthreads)
	#endif
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		sum = 0; // number of positions to copy
		p = 0; // position in ans_elt_holder
		if (t == 3) { // AAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E) || (m && !(x_i.ptr[j] ^ 0x2B))) { // position is a gap
					if (sum > 0) {
						memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
						p += sum;
						sum = 0;
					}
				} else {
					sum++;
				}
			}
			if (sum > 0)
				memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
		} else { // DNAStringSet or RNAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40 || (m && x_i.ptr[j] & 0x20)) { // position is a gap
					if (sum > 0) {
						memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
						p += sum;
						sum = 0;
					}
				} else {
					sum++;
				}
			}
			if (sum > 0)
				memcpy((char *) ans_elt_holder.ptr + p, x_i.ptr + j - sum, sum * sizeof(char));
		}
	}
	
	UNPROTECT(2);
	return ans;
}

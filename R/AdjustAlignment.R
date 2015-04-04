AdjustAlignment <- function(myXStringSet,
	perfectMatch=6,
	misMatch=0,
	gapLetter=-5,
	gapOpening=-0.5,
	gapExtension=0,
	substitutionMatrix=NULL,
	shiftPenalty=-0.2,
	threshold=0.1,
	weight=1,
	processors=NULL) {
	
	# error checking
	type <- switch(class(myXStringSet),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet."))
	if (length(myXStringSet) < 2)
		return(myXStringSet)
	u <- unique(width(myXStringSet))
	if (length(u)!=1)
		stop("Sequences in myXStringSet must be the same width (aligned).")
	if (u < 4) # no changes can be made
		return(myXStringSet)
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapLetter))
		stop("gapLetter must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold < 0)
		stop("threshold must be at least zero.")
	if (!is.numeric(shiftPenalty))
		stop("shiftPenalty must be a numeric.")
	if (shiftPenalty > 0)
		stop("shiftPenalty must be less than or equal to zero.")
	if (!is.numeric(weight))
		stop("weight must be a numeric.")
	if (length(weight)!=1 && length(weight)!=length(myXStringSet))
		stop("Length of weight must equal one or the length of the myXStringSet.")
	if (length(weight)==1) {
		weight <- rep(1, length(myXStringSet))
	} else {
		if (!isTRUE(all.equal(1, mean(weight))))
			stop("The mean of weight must be 1.")
	}
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors)!=processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	if (type==3) { # AAStringSet
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.null(substitutionMatrix)) {
			# use MIQS
			substitutionMatrix <- matrix(c(3, -1, 0, 0, 2, 0, 0, 0, -1, -1, -1, -1, -1, -2, 0, 1, 1, -4, -2, 0, -6, -1, 6, 0, -1, -3, 2, -1, -2, 1, -2, -3, 3, -1, -3, -1, 0, -1, -4, -2, -2, -6, 0, 0, 5, 3, -2, 1, 1, 0, 1, -4, -4, 1, -2, -3, -1, 1, 0, -5, -1, -3, -6, 0, -1, 3, 6, -4, 1, 3, -1, 0, -5, -5, 0, -3, -6, 0, 0, 0, -5, -4, -3, -6, 2, -3, -2, -4, 12, -3, -3, -2, -1, 0, -2, -3, 0, -3, -3, 1, 0, -6, -1, 2, -6, 0, 2, 1, 1, -3, 4, 2, -2, 1, -2, -2, 2, 0, -2, 0, 0, 0, -5, -3, -2, -6, 0, -1, 1, 3, -3, 2, 4, -1, 0, -3, -3, 1, -2, -4, 0, 0, 0, -6, -2, -2, -6, 0, -2, 0, -1, -2, -2, -1, 8, -2, -5, -5, -2, -4, -5, -2, 0, -2, -5, -4, -4, -6, -1, 1, 1, 0, -1, 1, 0, -2, 7, -2, -2, 0, -2, 0, -2, 0, 0, 0, 2, -2, -6, -1, -2, -4, -5, 0, -2, -3, -5, -2, 5, 3, -2, 2, 1, -4, -3, -1, -1, -1, 3, -6, -1, -3, -4, -5, -2, -2, -3, -5, -2, 3, 5, -2, 3, 2, -3, -3, -2, 0, 0, 2, -6, -1, 3, 1, 0, -3, 2, 1, -2, 0, -2, -2, 4, -1, -4, 0, 0, 0, -4, -2, -2, -6, -1, -1, -2, -3, 0, 0, -2, -4, -2, 2, 3, -1, 5, 1, -3, -2, -1, -2, -1, 1, -6, -2, -3, -3, -6, -3, -2, -4, -5, 0, 1, 2, -4, 1, 7, -4, -3, -2, 4, 5, 0, -6, 0, -1, -1, 0, -3, 0, 0, -2, -2, -4, -3, 0, -3, -4, 8, 0, 0, -4, -5, -3, -6, 1, 0, 1, 0, 1, 0, 0, 0, 0, -3, -3, 0, -2, -3, 0, 3, 2, -4, -2, -1, -6, 1, -1, 0, 0, 0, 0, 0, -2, 0, -1, -2, 0, -1, -2, 0, 2, 4, -5, -2, 0, -6, -4, -4, -5, -5, -6, -5, -6, -5, 0, -1, 0, -4, -2, 4, -4, -4, -5, 15, 5, -3, -6, -2, -2, -1, -4, -1, -3, -2, -4, 2, -1, 0, -2, -1, 5, -5, -2, -2, 5, 8, -1, -6, 0, -2, -3, -3, 2, -2, -2, -4, -2, 3, 2, -2, 1, 0, -3, -1, 0, -3, -1, 4, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (is.character(substitutionMatrix)) {
			if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
		"PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")))
				stop("Invalid substitutionMatrix.")
		}
		if (is.matrix(substitutionMatrix)) {
			if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix, envir=environment())))
		}
		subMatrix <- subMatrix[AAs, AAs]
		subMatrix <- as.numeric(subMatrix)
		
		functionCall <- "shiftGapsAA"
	} else { # DNAStringSet or RNAStringSet
		if (is.matrix(substitutionMatrix)) {
			bases <- c("A", "C", "G", "T")
			if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix[bases, bases]
		} else if (is.null(substitutionMatrix)) {
			bases <- c("A", "C", "G", "T")
			subMatrix <- matrix(misMatch,
				nrow=4, ncol=4)
			diag(subMatrix) <- perfectMatch
		} else {
			stop("substitutionMatrix must be NULL or a matrix.")
		}
		
		functionCall <- "shiftGaps"
	}
	
	# add gaps to both ends of the sequences
	ns <- names(myXStringSet)
	myXStringSet <- .Call("insertGaps",
		myXStringSet,
		c(1L, u + 1L),
		c(1L, 1L),
		type,
		processors,
		PACKAGE="DECIPHER")
	
	# adjust the alignment
	changes <- .Call(functionCall,
		myXStringSet,
		as.numeric(subMatrix),
		gapOpening,
		gapExtension,
		gapLetter,
		shiftPenalty,
		threshold,
		weight,
		PACKAGE="DECIPHER")
	
	# remove all 100% gap columns
	myXStringSet <- .Call("removeCommonGaps",
		myXStringSet,
		type,
		processors,
		PACKAGE="DECIPHER")
	
	names(myXStringSet) <- ns
	
	return(myXStringSet)
}
RemoveGaps <- function(myXStringSet,
	removeGaps="all",
	includeMask=FALSE,
	processors=1) {
	
	# error checking
	if (is(myXStringSet, "DNAStringSet")) {
		type <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		type <- 2L
	} else if (is(myXStringSet, "AAStringSet")) {
		type <- 3L
	} else {
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	GAPS <- c("none", "all", "common")
	removeGaps <- pmatch(removeGaps[1], GAPS)
	if (is.na(removeGaps))
		stop("Invalid removeGaps method.")
	if (removeGaps == -1)
		stop("Ambiguous removeGaps method.")
	if (!is.logical(includeMask))
		stop("includeMask must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors) != processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- .detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	if (removeGaps == 2L) { # all gaps
		ns <- names(myXStringSet)
		myXStringSet <- .Call("removeGaps",
			myXStringSet,
			type,
			includeMask,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- ns
	} else if (removeGaps == 3L) { # common gaps
		ns <- names(myXStringSet)
		myXStringSet <- .Call("removeCommonGaps",
			myXStringSet,
			type,
			includeMask,
			processors,
			PACKAGE="DECIPHER")
		names(myXStringSet) <- ns
	}
	
	return(myXStringSet)
}

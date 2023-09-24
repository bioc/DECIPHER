DistanceMatrix <- function(myXStringSet,
	method="overlap",
	type="matrix",
	includeTerminalGaps=FALSE,
	penalizeGapLetterMatches=TRUE,
	minCoverage=0,
	correction="none",
	processors=1,
	verbose=TRUE) {
	
	# error checking
	TYPES <- c("matrix", "dist")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	METHODS <- c("overlap", "shortest", "longest")
	method <- pmatch(method[1], METHODS)
	if (is.na(method))
		stop("Invalid method.")
	if (method == -1)
		stop("Ambiguous method.")
	CORRECTIONS <- c("none", "Jukes-Cantor", "JC", "JC69", "F81")
	correction <- pmatch(correction, CORRECTIONS)
	if (is.na(correction))
		stop("Invalid distance correction method.")
	if (correction == -1)
		stop("Ambiguous distance correction method.")
	if (correction == 3 || correction == 4)
		correction <- 2
	if (!is.logical(includeTerminalGaps))
		stop("includeTerminalGaps must be a logical.")
	if (!is.logical(penalizeGapLetterMatches))
		stop("penalizeGapLetterMatches must be a logical.")
	if (!is.numeric(minCoverage))
		stop("maxCoverage must be a numeric.")
	if (minCoverage < -1)
		stop("minCoverage must be at least -1.")
	if (minCoverage > 1)
		stop("minCoverage can be at most 1.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is(myXStringSet, "XStringSet"))
		stop("myXStringSet must be an XStringSet.")
	if (is(myXStringSet, "BStringSet"))
		stop("myXStringSet cannot be a BStringSet.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors) != processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	maxW <- unique(width(myXStringSet))
	if (length(maxW) != 1) {
		warning("\n",
			length(maxW),
			" different sequence lengths.\n",
			"Using shorter length in each comparison.\n")
	}
	numF <- length(myXStringSet)
	if (numF < 2) {
		stop("At least two sequences are required.")
	}
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
	}
	
	# calculate distance correction
	if (correction != 1L) {
		if (is(myXStringSet, "DNAStringSet")) {
			if (penalizeGapLetterMatches) {
				alphabet <- c(DNA_BASES, "-")
			} else {
				alphabet <- DNA_BASES
			}
		} else if (is(myXStringSet, "RNAStringSet")) {
			if (penalizeGapLetterMatches) {
				alphabet <- c(RNA_BASES, "-")
			} else {
				alphabet <- RNA_BASES
			}
		} else if (is(myXStringSet, "AAStringSet")) {
			if (penalizeGapLetterMatches) {
				alphabet <- c(AA_STANDARD, "-")
			} else {
				alphabet <- AA_STANDARD
			}
		}
		if (correction == 5L) {
			E <- letterFrequency(myXStringSet,
				alphabet)
			E <- colSums(E)
			E <- E/sum(E)
			E <- 1 - sum(E^2)
		} else { # JC69
			E <- 1/length(alphabet) # assumes even frequencies
			E <- 1 - length(alphabet)*sum(E^2)
		}
	} else {
		E <- 0
	}
	
	# calculate the distance matrix
	distMatrix <- .Call("distMatrix",
		myXStringSet,
		ifelse(is(myXStringSet, "AAStringSet"), 3L, 1L),
		includeTerminalGaps,
		penalizeGapLetterMatches,
		TRUE, # full matrix
		type,
		E,
		minCoverage,
		method,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	if (type == 1) { # matrix
		dimnames(distMatrix) <- list(names(myXStringSet),
			names(myXStringSet))
	} else { # dist
		attr(distMatrix, "Size") <- length(myXStringSet)
		if (!is.null(names(myXStringSet)))
			attr(distMatrix, "Labels") <- names(myXStringSet)
		attr(distMatrix, "Diag") <- TRUE
		attr(distMatrix, "Upper") <- TRUE
		class(distMatrix) <- "dist"
	}
	attr(distMatrix, "correction") <- CORRECTIONS[correction]
	
	if (verbose) {
		close(pBar)
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	return(distMatrix)
}

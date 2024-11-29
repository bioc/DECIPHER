DistanceMatrix <- function(myXStringSet,
	method="overlap",
	type="matrix",
	includeTerminalGaps=FALSE,
	penalizeGapLetterMatches=FALSE,
	minCoverage=0,
	correction=NA,
	substitutionMatrix=NULL,
	frequencies=NULL,
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
	if (length(correction) != 1L)
		stop("correction must be a single value.")
	if (is.na(correction)) { # Hamming distance
		empirical <- FALSE
		correction <- 0L
	} else {
		if (!is.character(correction))
			stop("correction must be a character string.")
		empirical <- endsWith(correction, "+F")
		if (empirical)
			correction <- substring(correction, 1L, nchar(correction) - 2L)
		CORRECTIONS <- c("JC69", "F81", "Poisson", "K80", "HKY85", "T92", "TN93", colnames(.ProtModels))
		correction <- pmatch(correction, CORRECTIONS)
		if (is.na(correction))
			stop("Invalid distance correction method.")
		if (correction == -1)
			stop("Ambiguous distance correction method.")
		if (empirical) {
			if (correction == 1L)
				stop("correction cannot be 'JC69+F'.")
			if (correction == 3L)
				stop("correction cannot be 'Poisson+F'.")
			if (correction == 4L)
				stop("correction cannot be 'K80+F'.")
		} else {
			if (correction == 2L)
				stop("correction cannot be 'F81'.")
			if (correction == 5L)
				stop("correction cannot be 'HKY85'.")
			if (correction == 6L)
				stop("correction cannot be 'TN92'.")
			if (correction == 7L)
				stop("correction cannot be 'TN93'.")
		}
	}
	if (correction != 0L && !is.null(substitutionMatrix)) {
		stop("substitutionMatrix cannot be specified when correction is not NA.")
	} else if (correction > 7L) {
		correction <- correction - 7L
		substitutionMatrix <- matrix(0, 20L, 20L)
		substitutionMatrix[upper.tri(substitutionMatrix)] <- .ProtModels[1:190, correction]
		substitutionMatrix <- substitutionMatrix + t(substitutionMatrix)
		if (is.null(frequencies)) {
			if (empirical) {
				frequencies <- letterFrequency(myXStringSet,
					AA_STANDARD, # same order as substitutionMatrix
					collapse=TRUE)
			} else if (is.null(frequencies)) {
				frequencies <- .ProtModels[191:210, correction]
			}
			frequencies <- frequencies/sum(frequencies)
		} else if (!empirical) {
			stop("frequencies can only be specified when correction includes '+F'.")
		} else {
			if (is.null(names(frequencies)))
				stop("frequencies must be a named vector.")
			if (length(frequencies) != 20L)
				stop("frequencies must be a length 20 vector.")
			if (!all(AA_STANDARD %in% names(frequencies)))
				stop("All amino acids in AA_STANDARD must be present in frequencies.")
			frequences <- frequencies[AA_STANDARD]
		}
		correction <- 0L
	}
	if (!is.logical(penalizeGapLetterMatches))
		stop("penalizeGapLetterMatches must be a logical.")
	if (is.null(frequencies)) {
		lkup <- NULL
	} else if (correction == 0L && is.null(substitutionMatrix)) {
		stop("substitutionMatrix must be specified with frequencies.")
	} else {
		if (correction == 1L)
			stop("frequencies cannot be specified when correction is 'JC69'.")
		if (correction == 3L)
			stop("frequencies cannot be specified when correction is 'Poisson'.")
		if (correction == 4L)
			stop("frequencies cannot be specified when correction is 'K80'.")
		if (length(frequencies) <= 1L)
			stop("frequencies must contain multiple values.")
		if (!is.numeric(frequencies))
			stop("frequencies must be a numeric vector.")
		if (any(is.na(frequencies)))
			stop("frequencies cannot contain NA values.")
		if (!isTRUE(all.equal(1, sum(frequencies))))
			stop("frequencies must sum to 1.")
		lkup <- names(frequencies)
	}
	if (!is.null(substitutionMatrix)) {
		if (is.character(substitutionMatrix))
			substitutionMatrix <- -.getSubMatrix(substitutionMatrix)
		if (!is.numeric(substitutionMatrix))
			stop("substitutionMatrix must be a numeric matrix.")
		if (is(substitutionMatrix, "matrix")) {
			if (!isSymmetric(substitutionMatrix))
				stop("substitutionMatrix must be a symmetric matrix.")
			if (!is.null(frequencies)) {
				if (any(diag(substitutionMatrix) != 0))
					stop("The diagonal of substitutionMatrix must contain zeros.")
				if (any(substitutionMatrix[lower.tri(substitutionMatrix)] <= 0))
					stop("substitutionMatrix must contain positive values.")
			}
		} else if (is(substitutionMatrix, "dist")) {
			if (!is.null(frequencies) && any(substitutionMatrix <= 0))
				stop("substitutionMatrix must contain positive values.")
			substitutionMatrix <- as.matrix(substitutionMatrix)
		} else {
			stop("substitutionMatrix must be a dist or matrix object.")
		}
		if (any(is.na(substitutionMatrix)))
			stop("substitutionMatrix cannot contain NA values.")
		if (nrow(substitutionMatrix) < 4L)
			stop("substitutionMatrix must have at least four rows/columns.")
		if (!is.null(frequencies) && nrow(substitutionMatrix) != length(frequencies))
			stop("The dimensions of substitutionMatrix do not match the number of frequencies.")
		if (is.null(lkup)) {
			lkup <- rownames(substitutionMatrix)
			if (is.null(lkup))
				lkup <- colnames(substitutionMatrix)
			if (is.null(lkup))
				stop("substitutionMatrix must have row or column names.")
		}
		if ("-" %in% lkup)
			stop("frequencies cannot contain '-'.")
		if ("." %in% lkup)
			stop("frequencies cannot contain '.'.")
		if ("+" %in% lkup)
			stop("frequencies cannot contain '+'.")
		if (!isFALSE(penalizeGapLetterMatches) && is.null(frequencies))
			stop("penalizeGapLetterMatches must be FALSE when substitutionMatrix is provided without frequencies.")
	}
	if (!is.null(lkup)) {
		if (any(duplicated(lkup))) {
			if (!is.null(frequencies)) {
				stop("Names of frequencies are duplicated.")
			} else {
				stop("Names of substitutionMatrix are duplicated.")
			}
		}
		if (!all(lkup %in% alphabet(myXStringSet))) {
			if (!is.null(frequencies)) {
				stop("Names of frequencies missing from the alphabet of myXStringSet.")
			} else {
				stop("Names of substitutionMatrix missing from the alphabet of myXStringSet.")
			}
		}
		lkup <- as(paste(lkup, collapse=""),
			class(myXStringSet))
	}
	if (!is.logical(includeTerminalGaps))
		stop("includeTerminalGaps must be a logical.")
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
		processors <- .detectCores()
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
	
	# calculate distance correction
	if (correction == 3L) { # Poisson
		E <- 1
	} else if (correction != 0L) {
		if (is(myXStringSet, "DNAStringSet")) {
			alphabet <- DNA_BASES
		} else if (is(myXStringSet, "RNAStringSet")) {
			alphabet <- RNA_BASES
		} else if (is(myXStringSet, "AAStringSet")) {
			if (correction > 3L)
				stop("Inapplicable correction when myXStringSet is an AAStringSet.")
			alphabet <- AA_STANDARD
		}
		if (correction == 1L) { # JC69
			E <- 1/length(alphabet) # assumes even frequencies
			E <- 1 - length(alphabet)*sum(E^2)
		} else if (correction == 4L) { # K80
			E <- numeric()
		} else {
			if (is.null(frequencies)) {
				E <- letterFrequency(myXStringSet,
					alphabet,
					collapse=TRUE)
				E <- unname(E)
				E <- E/sum(E)
			} else {
				if (length(frequencies) != length(alphabet))
					stop("Incorrect number of frequencies.")
				if (is.null(names(frequencies)))
					stop("frequencies must be a named vector.")
				if (!all(alphabet %in% names(frequencies)))
					stop("Names of frequencies missing from the alphabet of myXStringSet.")
				E <- frequencies[alphabet]
			}
			if (correction == 2L) { # F81+F
				E <- 1 - sum(E^2)
			} else if (correction == 6L) { # T92+F
				E <- E[2L] + E[3L]
				E <- 2*E*(1 - E)
				E <- c(E, 1 - E)
			} else {
				E <- c(E[1L]*E[3L], E[1L] + E[3L], E[2L]*E[4L], E[2L] + E[4L])
				if (correction == 5L) { # HKY85+F
					E <- c(E[1L]/E[2L] + E[3L]/E[4L], E[1L] + E[3L], E[2L]*E[4L])
				} # else TN93+F
			}
		}
	} else if (!is.null(frequencies)) {
		substitutionMatrix <- frequencies*substitutionMatrix
		norm <- sum(t(substitutionMatrix)*frequencies)
		diag(substitutionMatrix) <- -colSums(substitutionMatrix)
		substitutionMatrix <- substitutionMatrix/norm
		E <- eigen(substitutionMatrix, FALSE)
		E <- c(E$vectors, E$values, solve.default(E$vectors))
	} else if (!is.null(substitutionMatrix)) {
		E <- substitutionMatrix
		if (is(E, "integer"))
			mode(E) <- "double"
	} else {
		E <- 0
	}
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
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
		lkup,
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

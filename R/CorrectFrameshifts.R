CorrectFrameshifts <- function(myXStringSet,
	myAAStringSet,
	type="indels",
	acceptDistance=0.01,
	rejectDistance=0.60,
	maxComparisons=10,
	gapOpening=-13,
	gapExtension=-1,
	frameShift=-15,
	geneticCode=GENETIC_CODE,
	substitutionMatrix="PFASUM50",
	verbose=TRUE,
	processors=1) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	if (!is(myAAStringSet, "AAStringSet"))
		stop("myAAStringSet must be an AAStringSet.")
	if (length(myXStringSet) == 0)
		stop("At least one sequence is required in myXStringSet.")
	if (length(myAAStringSet) == 0)
		stop("At least one sequence is required in myAAStringSet.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("+", myXStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') in myXStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("-", myAAStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') in myAAStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern("+", myAAStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') in myAAStringSet must be removed before correcting frameshifts.")
	a <- vcountPattern(".", myAAStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') in myAAStringSet must be removed before correcting frameshifts.")
	org_index <- which(!duplicated(myAAStringSet))
	myAAStringSet <- myAAStringSet[org_index]
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
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	TYPES <- c("indels", "sequences", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(frameShift))
		stop("frameShift must be a numeric.")
	if (!is.numeric(rejectDistance))
		stop("rejectDistance must be a numeric.")
	if (rejectDistance > 1)
		stop("rejectDistance can be at most 1.")
	if (rejectDistance <= 0)
		stop("rejectDistance must be greater than zero.")
	if (!is.numeric(acceptDistance))
		stop("acceptDistance must be a numeric.")
	if (acceptDistance > rejectDistance)
		stop("acceptDistance can be at most rejectDistance.")
	if (acceptDistance < 0)
		stop("acceptDistance must be at least zero.")
	if (!is.numeric(maxComparisons))
		stop("maxComparisons must be a numeric.")
	if (maxComparisons < 1)
		stop("maxComparisons must be at least one.")
	if (floor(maxComparisons) != maxComparisons)
		stop("maxComparisons must be a whole number.")
	AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
	if (is.character(substitutionMatrix)) {
		if (substitutionMatrix == "PFASUM50") {
			subMatrix <- matrix(c(4.1181,-1.1516,-1.3187,-1.4135,0.4271,-0.5467,-0.6527,0.1777,-1.6582,-1.1243,-1.1843,-1.0235,-0.5685,-1.9515,-0.6072,0.8284,0.0361,-2.5368,-2.1701,0.0661,-11,-1.1516,6.341,0.0543,-0.6628,-3.2085,1.6006,0.5067,-1.961,0.7706,-3.5053,-3.0357,2.938,-1.9894,-3.7846,-1.3455,-0.4194,-0.5594,-2.1629,-1.7957,-2.9403,-11,-1.3187,0.0543,6.4672,2.3024,-2.5179,0.8192,0.5566,0.1585,1.104,-4.1629,-4.0977,0.8743,-2.6216,-3.805,-1.0904,1.1291,0.3253,-3.7763,-1.874,-3.6076,-11,-1.4135,-0.6628,2.3024,6.8156,-4.358,0.6705,2.582,-0.5667,-0.196,-5.475,-5.1661,0.226,-3.9595,-5.3456,-0.5662,0.4273,-0.5218,-4.7691,-3.4644,-4.5477,-11,0.4271,-3.2085,-2.5179,-4.358,13.5349,-3.3641,-4.3086,-2.1614,-1.8945,-0.7546,-0.9453,-3.8239,-0.5923,-0.8182,-3.6019,-0.3927,-0.801,-1.9317,-1.1607,0.0673,-11,-0.5467,1.6006,0.8192,0.6705,-3.3641,5.5795,2.1372,-1.5923,1.0862,-3.3001,-2.7545,1.872,-1.1216,-3.6631,-1.0426,0.1982,-0.0434,-3.061,-1.9214,-2.6993,-11,-0.6527,0.5067,0.5566,2.582,-4.3086,2.1372,5.5684,-1.6462,-0.2488,-4.1849,-4.0275,1.4821,-2.7964,-4.8311,-0.7028,0.0283,-0.312,-4.1969,-2.9489,-3.281,-11,0.1777,-1.961,0.1585,-0.5667,-2.1614,-1.5923,-1.6462,7.6508,-1.8185,-4.7058,-4.4215,-1.5991,-3.2786,-3.9992,-1.4409,0.184,-1.4823,-3.8328,-3.7343,-3.7264,-11,-1.6582,0.7706,1.104,-0.196,-1.8945,1.0862,-0.2488,-1.8185,9.7543,-3.3812,-2.8685,0.1425,-1.8724,-1.2545,-1.5333,-0.4285,-0.8896,-0.9385,1.6476,-2.8729,-11,-1.1243,-3.5053,-4.1629,-5.475,-0.7546,-3.3001,-4.1849,-4.7058,-3.3812,5.1229,2.5319,-3.5454,1.8309,0.9346,-3.4603,-3.0985,-1.2543,-1.5006,-1.117,3.3961,-11,-1.1843,-3.0357,-4.0977,-5.1661,-0.9453,-2.7545,-4.0275,-4.4215,-2.8685,2.5319,4.7049,-3.4581,2.5303,1.687,-3.365,-3.1578,-1.8626,-0.5308,-0.6881,1.4829,-11,-1.0235,2.938,0.8743,0.226,-3.8239,1.872,1.4821,-1.5991,0.1425,-3.5454,-3.4581,5.5476,-2.164,-4.3516,-0.7583,0.0275,-0.1516,-3.5889,-2.4422,-3.0453,-11,-0.5685,-1.9894,-2.6216,-3.9595,-0.5923,-1.1216,-2.7964,-3.2786,-1.8724,1.8309,2.5303,-2.164,7.0856,1.2339,-3.0823,-1.7587,-0.7402,-0.5841,-0.3946,0.9477,-11,-1.9515,-3.7846,-3.805,-5.3456,-0.8182,-3.6631,-4.8311,-3.9992,-1.2545,0.9346,1.687,-4.3516,1.2339,7.4322,-3.6222,-3.0316,-2.2851,2.6305,3.8302,0.1942,-11,-0.6072,-1.3455,-1.0904,-0.5662,-3.6019,-1.0426,-0.7028,-1.4409,-1.5333,-3.4603,-3.365,-0.7583,-3.0823,-3.6222,9.1796,-0.0652,-0.8587,-3.3634,-3.3006,-2.5443,-11,0.8284,-0.4194,1.1291,0.4273,-0.3927,0.1982,0.0283,0.184,-0.4285,-3.0985,-3.1578,0.0275,-1.7587,-3.0316,-0.0652,4.2366,1.8491,-3.1454,-2.1838,-2.1839,-11,0.0361,-0.5594,0.3253,-0.5218,-0.801,-0.0434,-0.312,-1.4823,-0.8896,-1.2543,-1.8626,-0.1516,-0.7402,-2.2851,-0.8587,1.8491,4.8833,-2.8511,-1.8993,-0.2699,-11,-2.5368,-2.1629,-3.7763,-4.7691,-1.9317,-3.061,-4.1969,-3.8328,-0.9385,-1.5006,-0.5308,-3.5889,-0.5841,2.6305,-3.3634,-3.1454,-2.8511,13.6485,3.3017,-1.851,-11,-2.1701,-1.7957,-1.874,-3.4644,-1.1607,-1.9214,-2.9489,-3.7343,1.6476,-1.117,-0.6881,-2.4422,-0.3946,3.8302,-3.3006,-2.1838,-1.8993,3.3017,8.7568,-1.2438,-11,0.0661,-2.9403,-3.6076,-4.5477,0.0673,-2.6993,-3.281,-3.7264,-2.8729,3.3961,1.4829,-3.0453,0.9477,0.1942,-2.5443,-2.1839,-0.2699,-1.851,-1.2438,4.6928,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,14),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else {
			subMatrix <- .getSubMatrix(substitutionMatrix)
		}
	} else if (is.matrix(substitutionMatrix)) {
		if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
			any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
			stop("substitutionMatrix is incomplete.")
		subMatrix <- substitutionMatrix
	} else {
		stop("Invalid substitutionMatrix must be NULL, a character string, or a matrix.")
	}
	subMatrix <- subMatrix[AAs, AAs]
	mode(subMatrix) <- "numeric"
	
	# de-replicate
	w <- which(!duplicated(myXStringSet))
	l <- length(myXStringSet)
	ns <- names(myXStringSet)
	if (length(w) < l) {
		dupes <- TRUE
		derep <- numeric(l)
		derep[w] <- seq_along(w)
		derep[-w] <- match(myXStringSet[-w],
			myXStringSet[w])
		myXStringSet <- myXStringSet[w]
		l <- length(w)
	} else {
		dupes <- FALSE
	}
	
	ends <- width(myXStringSet)
	if (min(ends) < 2)
		stop("All sequences in myXStringSet must be at least two nucleotides long.")
	
	frames <- AAStringSet()
	for (i in 1:3) {
		end <- ends
		offset <- end - i + 1
		end <- end - offset %% 3
		end <- ifelse(end < i - 1,
			i - 1,
			end)
		
		AA <- translate(subseq(myXStringSet,
				i,
				end),
			genetic.code=geneticCode,
			if.fuzzy.codon="solve")
		frames <- c(frames,
			AA)
	}
	
	if (length(myAAStringSet) > 1) {
		if (verbose) {
			time.1 <- Sys.time()
			cat("Finding the closest reference amino acid sequences:\n",
				sep="")
			flush.console()
			pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
		} else {
			pBar <- NULL
		}
		
		v1 <- .Call("enumerateSequenceAA",
			frames,
			7L,
			PACKAGE="DECIPHER")
		v2 <- .Call("enumerateSequenceAA",
			myAAStringSet,
			7L,
			PACKAGE="DECIPHER")
		
		v1 <- lapply(v1,
			sort)
		v2 <- lapply(v2,
			sort)
		
		d <- .Call("matchListsDual",
			v1,
			v2,
			verbose,
			pBar,
			processors,
			PACKAGE="DECIPHER")
		d <- rowsum(d,
			rep(seq_len(l),
				3),
			na.rm=TRUE)
		
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
		}
	} else {
		d <- matrix(0, nrow=l)
	}
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Assessing frameshifts in nucleotide sequences:\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
	}
	
	widths <- width(myAAStringSet)
	
	index <- matrix(nrow=nrow(d), ncol=ncol(d))
	for (i in seq_len(nrow(d)))
		index[i,] <- order(d[i,], widths, decreasing=TRUE)
	if (ncol(index) > maxComparisons) {
		index <- index[, seq_len(maxComparisons), drop=FALSE]
	} else {
		maxComparisons <- ncol(index)
	}
	
	X <- .Call("findFrameshifts",
		myAAStringSet,
		ends,
		frames,
		index,
		org_index,
		maxComparisons,
		gapOpening,
		gapExtension,
		frameShift,
		acceptDistance,
		rejectDistance,
		subMatrix,
		verbose,
		pBar,
		PACKAGE="DECIPHER")
	if (type != 2) {
		X <- lapply(X,
			setNames,
			c("insertions", "deletions", "distance", "index"))
	}
	
	if (type > 1) {
		pos <- lapply(X,
			function(x) {
				ins <- rle(x[[1]])
				dels <- rle(x[[2]])
				IRanges(start=c(ins[[2]],
						dels[[2]]),
					width=c(ins[[1]],
						rep(0,
							length(dels[[2]]))))
			})
		pos <- unname(pos)
		
		val <- lapply(X,
			function(x) {
				ins <- rep("",
					length(unique(x[[1]])))
				dels <- rle(x[[2]])
				types <- c("N", "NN")
				dels <- types[dels[[1]]]
				c(ins, dels)
			})
		
		if (type == 2L) { # sequences
			X <- replaceAt(myXStringSet,
				as(pos, "IRangesList"),
				val)
		} else { # both
			X <- list(indels=X,
				sequences=replaceAt(myXStringSet,
					as(pos, "IRangesList"),
					val))
		}
	}
	
	if (dupes) { # re-replicate
		if (type == 1 || type == 2) {
			X <- X[derep]
		} else {
			X[[1]] <- X[[1]][derep]
			X[[2]] <- X[[2]][derep]
		}
	}
	
	if (type == 1 || type == 2) {
		names(X) <- ns
	} else {
		names(X[[1]]) <- names(X[[2]]) <- ns
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
	}
	
	return(X)
}

SearchIndex <- function(myXStringSet,
	invertedIndex,
	type="one",
	minScore=NA,
	scoreOnly=FALSE,
	sepCost=-0.4,
	gapCost=-2.5,
	maskRepeats=TRUE,
	maskLCRs=TRUE,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	TYPES <- c("all", "one", "top")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (length(minScore) != 1L)
		stop("minScore must be a single numeric.")
	if (!is.na(minScore) && !is.numeric(minScore))
		stop("minScore must be a numeric.")
	if (!isTRUEorFALSE(scoreOnly))
		stop("scoreOnly must be TRUE or FALSE.")
	if (length(sepCost) != 1L)
		stop("sepCost must be a single numeric.")
	if (is.na(sepCost) || !is.numeric(sepCost))
		stop("sepCost must be a numeric.")
	if (length(gapCost) != 1L)
		stop("gapCost must be a single numeric.")
	if (is.na(gapCost) || !is.numeric(gapCost))
		stop("gapCost must be a numeric.")
	if (!isTRUEorFALSE(maskRepeats))
		stop("maskRepeats must be TRUE or FALSE.")
	if (!isTRUEorFALSE(maskLCRs))
		stop("maskLCRs must be TRUE or FALSE.")
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be TRUE or FALSE.")
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
	
	K <- invertedIndex$k
	step <- invertedIndex$step
	alphabet <- invertedIndex$alphabet
	freqs <- invertedIndex$frequency
	num <- invertedIndex$count
	len <- invertedIndex$length
	loc <- invertedIndex$location
	ind <- invertedIndex$index
	
	if (length(alphabet) == 20L) {
		if (!is(myXStringSet, "AAStringSet"))
			stop("myXStringSet must be an AAStringSet.")
		xtype <- 3L
	} else if (length(alphabet) == 4L) {
		if (!(is(myXStringSet, "DNAStringSet") || is(myXStringSet, "RNAStringSet")))
			stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
		xtype <- 1L
	} else {
		stop("Corrupted invertedIndex.")
	}
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
	}
	
	L <- length(num)
	maxSep <- as.integer(sqrt(L))
	
	offset <- numeric(L)
	for (i in seq_len(L - 1L))
		offset[i + 1L] <- offset[i] + num[i]
	total <- sum(len) + 1 # size of the target database
	
	if (xtype == 3L) {
		kmers <- .Call("enumerateSequenceReducedAA",
			myXStringSet,
			K,
			alphabet,
			maskRepeats,
			maskLCRs,
			processors,
			PACKAGE="DECIPHER")
	} else {
		kmers <- .Call("enumerateSequence",
			myXStringSet,
			K,
			maskRepeats,
			maskLCRs,
			processors,
			PACKAGE="DECIPHER")
	}
	
	ans <- .Call("searchIndex",
		kmers, # query k-mers [0 to length(num)]
		K, # wordSize
		step, # separation between k-mers (>= 1 and <= K)
		-log(freqs), # -log of normalized letter frequencies
		num, # count
		loc, # location
		ind, # index
		len, # positions
		sepCost, # sepC
		gapCost, # gapC
		type, # output
		total, # size of target database
		minScore, # minimum score or NA to calculate
		scoreOnly, # FALSE to output anchor positions
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (!scoreOnly) {
		pos <- ans[[4L]]
		ans <- ans[1:3]
	}
	names(ans) <- c("Pattern", "Subject", "Score")
	ans <- data.frame(ans)
	if (!scoreOnly)
		ans$Position <- pos
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(ans)
}

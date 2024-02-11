IndexSeqs <- function(subject,
	K,
	sensitivity,
	percentIdentity,
	patternLength,
	step=1,
	alphabet=AA_REDUCED[[171]],
	maskRepeats=TRUE,
	maskLCRs=TRUE,
	batchSize=1e+07,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (is(subject, "DNAStringSet")) {
		xtype <- 1L
		
		alphabet <- setNames(0L:3L, DNA_BASES)
		size <- 4 # alphabet size
		maxK <- 15L
	} else if (is(subject, "RNAStringSet")) {
		xtype <- 2L
		
		alphabet <- setNames(0L:3L, RNA_BASES)
		size <- 4 # alphabet size
		maxK <- 15L
	} else if (is(subject, "AAStringSet")) {
		xtype <- 3L
		
		if (!is.character(alphabet))
			stop("alphabet must be a character vector.")
		if (any(alphabet == ""))
			stop("No elements of alphabet can be empty.")
		r <- strsplit(alphabet, "", fixed=TRUE)
		alphabet <- setNames(rep(0L, 20),
			AA_STANDARD)
		for (i in seq_along(r)) {
			w <- which(!(r[[i]] %in% AA_STANDARD))
			if (length(w) > 0)
				stop("Unrecognized letter(s) found in alphabet:  ",
					paste(r[[i]][w], collapse=", "),
					".")
			w <- which(alphabet[r[[i]]] != 0L)
			if (length(w) > 0)
				stop("Repeated amino acids found in alphabet:  ",
					paste(r[[i]][w], collapse=", "),
					".")
			alphabet[r[[i]]] <- i
		}
		w <- which(alphabet == 0L)
		if (length(w) > 0)
			stop("Standard amino acids missing from alphabet:  ",
				paste(names(w), collapse=", "),
				".")
		size <- max(alphabet) # alphabet size
		alphabet <- alphabet - 1L
		maxK <- as.integer(log(2147483647L, size)) # 2147483647L == 2^31 - 1
	} else {
		stop("subject must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	l <- length(subject)
	if (l == 0L)
		stop("subject must contain at least one sequence.")
	if (l > 2147483647L) # 2^31 - 1 == 2147483647
		stop("subject can contain at most 2,147,483,647 sequences.")
	if (step < 1)
		stop("step must be at least 1.")
	if (step != floor(step))
		stop("step must be a whole number.")
	if (!isTRUEorFALSE(maskRepeats))
		stop("maskRepeats must be TRUE or FALSE.")
	if (!isTRUEorFALSE(maskLCRs))
		stop("maskLCRs must be TRUE or FALSE.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (floor(batchSize) != batchSize)
		stop("batchSize must be a whole number.")
	if (batchSize <= 0)
		stop("batchSize must be greater than zero.")
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
	
	if (missing(K)) {
		if (missing(patternLength))
			stop("patternLength must be specified if K is not specified.")
		if (missing(sensitivity))
			stop("sensitivity must be specified if K is not specified.")
		if (missing(percentIdentity))
			stop("percentIdentity must be specified if K is not specified.")
		if (length(patternLength) != 1L)
			stop("patternLength must be a single number.")
		if (length(sensitivity) != 1L)
			stop("sensitivity must be a single number.")
		if (length(percentIdentity) != 1L)
			stop("percentIdentity must be a single number.")
		if (!is.numeric(patternLength))
			stop("patternLength must be a numeric.")
		if (!is.numeric(sensitivity))
			stop("sensitivity must be a numeric.")
		if (!is.numeric(percentIdentity))
			stop("percentIdentity must be a numeric.")
		if (is.na(patternLength))
			stop("patternLength cannot be NA.")
		if (is.na(sensitivity))
			stop("sensitivity cannot be NA.")
		if (is.na(percentIdentity))
			stop("percentIdentity cannot be NA.")
		if (patternLength < 15L)
			stop("patternLength must be at least 15.")
		if (sensitivity < 0.5)
			stop("sensitivity must be at least 0.5.")
		if (sensitivity >= 1)
			stop("sensitivity must be less than 1.")
		if (percentIdentity <= 1)
			stop("percentIdentity must be greater than 1.")
		if (percentIdentity > 100)
			stop("percentIdentity must be less than 100.")
		K <- ceiling(patternLength/step*log(1 - exp(log(percentIdentity/100)*log(patternLength*max(width(subject))/step)/log(size)))/log(1 - sensitivity))
		if (K > maxK) {
			K <- maxK
		} else if (K < 2) {
			stop("The desired sensitivity is unachievable at the specified percentIdentity for queries of patternLength.")
		} else if (K < step) {
			stop("The desired sensitivity is unachievable with the specified step.")
		}
	} else {
		if (length(K) != 1L)
			stop("K must be a single number.")
		if (!is.numeric(K))
			stop("K must be a numeric.")
		if (is.na(K))
			stop("K cannot be NA.")
		if (K != floor(K))
			stop("K must be a whole number.")
		if (K < 2L)
			stop("K must be at least 2.")
		if (K > maxK)
			stop("K can be at most ", maxK, ".")
	}
	if (step > K)
		stop("step can be at most K.")
	
	K <- as.integer(K)
	step <- as.integer(step)
	L <- as.integer(size^K) # number of possible k-mers
	num <- integer(L) # assumes fewer than 2^31 hits per k-mer
	
	# process k-mers in batches to conserve memory
	cum_pos <- cumsum(as.double(width(subject)))
	if (verbose) {
		pBar <- txtProgressBar(max=ifelse(cum_pos[l] <= batchSize, cum_pos[l], 2*cum_pos[l]), style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	}
	last <- 0 # number of positions processed
	prev <- 0L # previous index
	batches <- 0L
	while (prev < l) {
		batches <- batches + 1L
		# determine number of sequences to process
		prev <- prev + 1L
		N <- prev # number of sequences processed
		while(N < l && cum_pos[N + 1L] - last <= batchSize)
			N <- N + 1L
		last <- cum_pos[N]
		
		if (xtype == 3L) {
			kmers <- .Call("enumerateSequenceReducedAA",
				subject[prev:N],
				K,
				alphabet,
				maskRepeats,
				maskLCRs,
				processors,
				PACKAGE="DECIPHER")
		} else {
			kmers <- .Call("enumerateSequence",
				subject[prev:N],
				K,
				maskRepeats,
				maskLCRs,
				processors,
				PACKAGE="DECIPHER")
		}
		
		.Call("countIndex",
			num, # in-place change of num
			kmers, # query k-mers
			step, # separation between k-mers (>= 1 and <= K)
			PACKAGE="DECIPHER")
		
		prev <- N
		if (verbose)
			setTxtProgressBar(pBar, last)
	}
	
	offset <- numeric(L) # starting index of each k-mer
	freqs <- numeric(size) # background frequencies
	.Call("approxFreqs", # returns NULL
		offset, # in-place change of offset
		freqs, # in-place change of freqs
		num, # count per k-mer
		PACKAGE="DECIPHER")
	total <- offset[L] + num[L]
	freqs[freqs == 0] <- 1 # add pseudocount
	freqs <- freqs/sum(freqs) # normalize frequencies
	
	loc <- integer(total) # location of k-mer in the sequence
	ind <- integer(total) # sequence index of the k-mer
	len <- integer(N) # number of unmasked positions per sequence
	
	if (batches == 1L) {
		.Call("updateIndex", # returns NULL
			offset, # in-place change of offset
			kmers, # query k-mers
			K, # wordSize
			step, # separation between k-mers (>= 1 and <= K)
			loc, # in-place change of loc
			ind, # in-place change of ind
			len, # in-place change of len
			0L, # previous count
			PACKAGE="DECIPHER")
	} else { # batch processing requires re-computing k-mers
		last <- 0 # number of positions processed
		prev <- 0L # previous index
		while (prev < l) {
			# determine number of sequences to process
			prev <- prev + 1L
			N <- prev # number of sequences processed
			while(N < l && cum_pos[N + 1L] - last <= batchSize)
				N <- N + 1L
			last <- cum_pos[N]
			
			if (xtype == 3L) {
				kmers <- .Call("enumerateSequenceReducedAA",
					subject[prev:N],
					K,
					alphabet,
					maskRepeats,
					maskLCRs,
					processors,
					PACKAGE="DECIPHER")
			} else {
				kmers <- .Call("enumerateSequence",
					subject[prev:N],
					K,
					maskRepeats,
					maskLCRs,
					processors,
					PACKAGE="DECIPHER")
			}
			
			.Call("updateIndex", # returns NULL
				offset, # in-place change of offset
				kmers, # query k-mers
				K, # wordSize
				step, # separation between k-mers (>= 1 and <= K)
				loc, # in-place change of loc
				ind, # in-place change of ind
				len, # in-place change of len
				prev - 1L, # previous count
				PACKAGE="DECIPHER")
			
			prev <- N
			if (verbose)
				setTxtProgressBar(pBar, last + cum_pos[l])
		}
	}
	
	ans <- list(k=K,
		step=step,
		alphabet=alphabet,
		frequency=freqs,
		count=num,
		length=len,
		location=loc,
		index=ind)
	class(ans) <- "InvertedIndex"
	
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
	
	return(ans)
}

Clusterize <- function(myXStringSet,
	cutoff=0,
	method="overlap",
	includeTerminalGaps=FALSE,
	penalizeGapLetterMatches=NA,
	minCoverage=0.5,
	maxPhase1=2e4,
	maxPhase2=2e3,
	maxPhase3=2e3,
	maxAlignments=2e2,
	rareKmers=50,
	probability=0.99,
	invertCenters=FALSE,
	singleLinkage=FALSE,
	alphabet=AA_REDUCED[[186]],
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.numeric(cutoff))
		stop("cutoff must be a numeric.")
	if (is.integer(cutoff))
		cutoff <- as.numeric(cutoff)
	if (any(is.na(cutoff)))
		stop("cutoff must not contain NA values.")
	if (any(cutoff < 0))
		stop("cutoff must be at least zero.")
	if (any(cutoff >= 1))
		stop("cutoff must be less than one.")
	if (any(duplicated(cutoff)))
		stop("cutoff cannot contain duplicated values.")
	METHODS <- c("overlap", "shortest", "longest")
	method <- pmatch(method[1], METHODS)
	if (is.na(method))
		stop("Invalid method.")
	if (method == -1)
		stop("Ambiguous method.")
	if (!is.logical(includeTerminalGaps))
		stop("includeTerminalGaps must be a logical.")
	if (!is.logical(penalizeGapLetterMatches))
		stop("penalizeGapLetterMatches must be a logical.")
	ASC <- TRUE
	if (is.unsorted(cutoff)) {
		if (is.unsorted(rev(cutoff))) {
			stop("cutoff must be sorted.")
		} else {
			ASC <- FALSE
		}
	}
	if (!is.numeric(minCoverage))
		stop("maxCoverage must be a numeric.")
	if (minCoverage < -1)
		stop("minCoverage must be at least -1.")
	if (minCoverage > 1)
		stop("minCoverage can be at most 1.")
	if (!is.numeric(maxPhase1))
		stop("maxPhase1 must be a numeric.")
	if (length(maxPhase1) != 1)
		stop("maxPhase1 must only be a single number.")
	if (maxPhase1 < 10000)
		stop("maxPhase1 must be at least 10000.")
	if (floor(maxPhase1) != maxPhase1)
		stop("maxPhase1 must be a whole number.")
	if (!is.numeric(maxPhase2))
		stop("maxPhase2 must be a numeric.")
	if (length(maxPhase2) != 1)
		stop("maxPhase2 must only be a single number.")
	if (maxPhase2 < 1)
		stop("maxPhase2 must be at least 1.")
	if (floor(maxPhase2) != maxPhase2)
		stop("maxPhase2 must be a whole number.")
	if (!is.numeric(maxPhase3))
		stop("maxPhase3 must be a numeric.")
	if (length(maxPhase3) != 1)
		stop("maxPhase3 must only be a single number.")
	if (maxPhase3 < 2)
		stop("maxPhase3 must be at least 2.")
	if (floor(maxPhase3) != maxPhase3)
		stop("maxPhase3 must be a whole number.")
	if (!is.numeric(maxAlignments))
		stop("maxAlignments must be a numeric.")
	if (length(maxAlignments) != 1)
		stop("maxAlignments must only be a single number.")
	if (maxAlignments > maxPhase3)
		stop("maxAlignments must be less than or equal to maxPhase3.")
	if (floor(maxAlignments) != maxAlignments)
		stop("maxAlignments must be a whole number.")
	if (!is.numeric(rareKmers))
		stop("rareKmers must be a numeric.")
	if (length(rareKmers) != 1)
		stop("rareKmers must only be a single number.")
	if (rareKmers < 1)
		stop("rareKmers must be at least 1.")
	if (floor(rareKmers) != rareKmers)
		stop("rareKmers must be a whole number.")
	rareKmers <- as.integer(rareKmers)
	if (!is.numeric(probability))
		stop("probability must be a numeric.")
	if (length(probability) != 1)
		stop("probability must only be a single number.")
	if (probability <= 0 || probability >= 1)
		stop("probability must be between zero and one (exclusive).")
	if (!is.logical(singleLinkage))
		stop("singleLinkage must be a logical.")
	if (!is.logical(invertCenters))
		stop("invertCenters must be a logical.")
	if (invertCenters && singleLinkage)
		stop("invertCenters must be FALSE if singleLinkage is TRUE.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
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
	if (is(myXStringSet, "DNAStringSet")) {
		typeX <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		typeX <- 2L
	} else if (is(myXStringSet, "AAStringSet")) {
		typeX <- 3L
	} else {
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') are not allowed in myXStringSet.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') are not allowed in myXStringSet.")
	rm(a)
	widths <- width(myXStringSet)
	L <- length(myXStringSet)
	if (L == 0L)
		stop("myXStringSet contains no sequences.")
	if (L > 2147483647L)
		stop("myXStringSet can have at most 2,147,483,647 sequences.")
	if (all(widths == 0L))
		stop("All sequences in myXStringSet are zero width.")
	
	# initialize parameters
	minCount <- 10L # minimum number of replicate timings to stop optimizing processors (>= 1)
	quantile <- 0.99 # width quantile determining k-mer size (>> 0 and <= 1)
	fracRandom <- 0.1 # target fraction of random occurrences per rare k-mer (> 0 and << 1)
	pow <- 0.5 # power (> 0 and <= 1) applied to the number of k-mers to randomly project into a smaller space
	N <- 100L # approximate frequency of random k-mers (1 per every N on average)
	threshold <- 0.9 # inverse of false discovery rate for partly disjoint set inclusion in phase 1 (>> 0 and <= 1)
	maxSample <- 100L # number of k-mers to sample in phase 3 (>> 1 and ideally <= cache size)
	alpha1 <- 0.05 # weight of exponential moving average to smooth changes in variance of rank order (>= 0 and <= 1)
	split <- maxPhase3/2 # decrease in rank order variance required to split partitions (> 0 and < maxPhase3)
	attempts <- 20L # alignment attempts before possibly skipping (>> 0)
	binSize <- 0.05 # size of bins for distributing k-mer similarities (> 0 and << 1)
	minSimilarities <- 1000L # minimum similarities per cutoff to stop collecting data (>> 0)
	needAlignments <- 10L # minimum alignments within binSize needed to set k-mer similarity limit (>= 0 and <= minSimilarities)
	minAlignments <- 1L # minimum alignments to perform per sequence in phase 3 (> 0)
	batchSize <- maxPhase3 # number of k-mer similarities to compute per batch (ideally >> processors)
	boundComparisons <- 10 # fold-limit on variation beyond maxPhase3 (average comparisons per sequence)
	numRandom <- as.integer(maxPhase1/rareKmers*fracRandom) # expected number of random occurrences of rare k-mers
	minFraction <- 0.1 # minimum fraction of hits to consider further comparison in phase 3 (>= 0 and <= 1)
	stdDevs <- 1.96 # standard deviations above the mean to continue looking for a cluster (> 0)
	alpha2 <- 0.01 # weight of exponential moving average to smooth proportions in phase 3 (>= 0 and <= 1)
	updateRate <- 100L # update subsetRate at least every updateRate iterations in phase 3 (>> 1)
	
	if (maxAlignments < minAlignments)
		stop("maxAlignments must be at least ", minAlignments, ".")
	if (maxPhase2*maxAlignments < minSimilarities)
		stop("maxPhase2 must be at least", ceiling(minSimilarities/maxAlignments), ".")
	interval <- log((1 - probability)/probability) # logit(probability of error) when setting k-mer similarity limit
	if (verbose)
		time.1 <- Sys.time()
	
	if (typeX == 3L) { # AAStringSet
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
		words <- max(alphabet)
		alphabet <- alphabet - 1L
		
		# use rounded PFASUM
		sM <- matrix(c(4L, -1L, -1L, -1L, 0L, -1L, -1L, 0L, -2L, -1L, -1L, -1L, -1L, -2L, -1L, 1L, 0L, -3L, -2L, 0L, -1L, 6L, 0L, -1L, -3L, 2L, 1L, -2L, 1L, -4L, -3L, 3L, -2L, -4L, -1L, 0L, -1L, -2L, -2L, -3L, -1L, 0L, 6L, 2L, -3L, 1L, 1L, 0L, 1L, -4L, -4L, 1L, -3L, -4L, -1L, 1L, 0L, -4L, -2L, -4L, -1L, -1L, 2L, 7L, -4L, 1L, 3L, -1L, 0L, -5L, -5L, 0L, -4L, -5L, -1L, 0L, -1L, -5L, -3L, -5L, 0L, -3L, -3L, -4L, 14L, -3L, -4L, -2L, -2L, -1L, -1L, -4L, -1L, -1L, -4L, 0L, -1L, -2L, -1L, 0L, -1L, 2L, 1L, 1L, -3L, 6L, 2L, -2L, 1L, -3L, -3L, 2L, -1L, -4L, -1L, 0L, 0L, -3L, -2L, -3L, -1L, 1L, 1L, 3L, -4L, 2L, 6L, -2L, 0L, -4L, -4L, 1L, -3L, -5L, -1L, 0L, 0L, -4L, -3L, -3L, 0L, -2L, 0L, -1L, -2L, -2L, -2L, 8L, -2L, -5L, -4L, -2L, -3L, -4L, -1L, 0L, -1L, -4L, -4L, -4L, -2L, 1L, 1L, 0L, -2L, 1L, 0L, -2L, 10L, -3L, -3L, 0L, -2L, -1L, -2L, 0L, -1L, -1L, 2L, -3L, -1L, -4L, -4L, -5L, -1L, -3L, -4L, -5L, -3L, 5L, 3L, -4L, 2L, 1L, -3L, -3L, -1L, -2L, -1L, 3L, -1L, -3L, -4L, -5L, -1L, -3L, -4L, -4L, -3L, 3L, 5L, -3L, 3L, 2L, -3L, -3L, -2L, -1L, -1L, 1L, -1L, 3L, 1L, 0L, -4L, 2L, 1L, -2L, 0L, -4L, -3L, 6L, -2L, -4L, -1L, 0L, 0L, -4L, -2L, -3L, -1L, -2L, -3L, -4L, -1L, -1L, -3L, -3L, -2L, 2L, 3L, -2L, 7L, 1L, -3L, -2L, -1L, -1L, 0L, 1L, -2L, -4L, -4L, -5L, -1L, -4L, -5L, -4L, -1L, 1L, 2L, -4L, 1L, 7L, -4L, -3L, -2L, 3L, 4L, 0L, -1L, -1L, -1L, -1L, -4L, -1L, -1L, -1L, -2L, -3L, -3L, -1L, -3L, -4L, 9L, 0L, -1L, -3L, -3L, -3L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, -3L, -3L, 0L, -2L, -3L, 0L, 4L, 2L, -3L, -2L, -2L, 0L, -1L, 0L, -1L, -1L, 0L, 0L, -1L, -1L, -1L, -2L, 0L, -1L, -2L, -1L, 2L, 5L, -3L, -2L, 0L, -3L, -2L, -4L, -5L, -2L, -3L, -4L, -4L, -1L, -2L, -1L, -4L, -1L, 3L, -3L, -3L, -3L, 14L, 3L, -2L, -2L, -2L, -2L, -3L, -1L, -2L, -3L, -4L, 2L, -1L, -1L, -2L, 0L, 4L, -3L, -2L, -2L, 3L, 9L, -1L, 0L, -3L, -4L, -5L, 0L, -3L, -3L, -4L, -3L, 3L, 1L, -3L, 1L, 0L, -3L, -2L, 0L, -2L, -1L, 5L),
			nrow=20,
			ncol=20)
		entropy <- .Call("alphabetSizeReducedAA",
			myXStringSet,
			alphabet,
			PACKAGE="DECIPHER")
		if (words == 1L)
			stop("More than one grouping of amino acids is required in the alphabet.")
		maxK <- as.integer(log(2147483647L, words)) # 2147483647L == 2^31 - 1
	} else { # DNAStringSet or RNAStringSet
		# penalize transitions less than transversions
		sM <- matrix(c(3L, -3L, 0L, -3L, -3L, 3L, -3L, 0L, 0L, -3L, 3L, -3L, -3L, 0L, -3L, 3L),
			nrow=4,
			ncol=4)
		entropy <- .Call("alphabetSize",
			myXStringSet,
			PACKAGE="DECIPHER")
		words <- 4L
		maxK <- 15L
	}
	
	# identify duplicated sequences
	x <- selfmatch(myXStringSet)
	u <- unique(x)
	l <- length(u)
	t <- tabulate(x, L)[u]
	wu <- widths[u]
	
	wordSizeQuantile <- quantile(wu, quantile, names=FALSE)
	wordSize <- ceiling(log(N*wordSizeQuantile, entropy))
	if (wordSize < 1)
		wordSize <- 1L
	if (wordSize > maxK)
		wordSize <- maxK
	sumL <- sum(pmin.int(widths, rareKmers))
	kmerSize <- ceiling(log(sumL/numRandom, entropy))
	if (kmerSize < 1)
		kmerSize <- 1L
	if (kmerSize > maxK)
		kmerSize <- maxK
	
	minPhase1 <- as.integer(rareKmers^2/fracRandom/2147483647*l)
	if (maxPhase1 < minPhase1)
		stop("maxPhase1 must be at least ", minPhase1, ".")
	rm(minPhase1)
	if (rareKmers > 2147483647/l)
		stop("rareKmers can be at most ", floor(2147483647/l), ".")
	
	if (verbose) {
		cat("Partitioning sequences by ",
			kmerSize,
			"-mer similarity:\n",
			sep="")
		flush.console()
	}
	
	lc <- length(cutoff)
	C <- vector("list", lc) # cluster numbers
	if (lc > 1) {
		names(C) <- paste("cluster",
			gsub(".", "_", prettyNum(cutoff), fixed=TRUE),
			sep="_")
	} else {
		names(C) <- "cluster"
	}
	bins <- seq(0, 1, binSize) # similarity ranges
	counts <- rep(1L, length(bins) - 1L)
	incBins <- .bincode(c(cutoff - binSize, cutoff, cutoff + binSize),
		bins,
		include.lowest=TRUE)
	incBins <- unique(which(tabulate(incBins, length(counts)) > 0))
	cutoff <- 1 - cutoff
	
	overlap <- function(x, y, processors=1L) {
		.Call("matchOverlap",
			x,
			y,
			v,
			wordSize,
			processors,
			PACKAGE="DECIPHER")
	}
	
	similarity <- function(res, w1, w2, processors=1L, coverage=minCoverage) {
		.Call("similarities",
			res,
			w1,
			w2,
			includeTerminalGaps,
			penalizeGapLetterMatches,
			minCoverage,
			method,
			processors,
			PACKAGE="DECIPHER")
	}
	
	countHits <- function(x, v, processors=1L) {
		.Call("countHits",
			x,
			v,
			processors,
			PACKAGE="DECIPHER")
	}
	
	dist <- function(ali) {
		.Call("distMatrix",
			ali,
			typeX,
			includeTerminalGaps,
			penalizeGapLetterMatches,
			TRUE, # full matrix
			2L, # type = "dist"
			0, # correction (none)
			minCoverage,
			method,
			FALSE, # verbose
			NULL, # progress bar
			1L, # processors
			PACKAGE="DECIPHER")
	}
	
	align <- function(pair, # pair of indices in myXStringSet
		anchor=NULL, # optional anchors
		GO=-10, # gap opening
		GE=-2, # gap extension
		TG=0, # terminal gap
		maxLength=5, # maximum length to skip alignment in equal length regions
		perProcessor=5L, # minimum number of regions to align per processor
		processors=1L) {
		n <- as.integer(length(anchor)/4) # ncol(anchor) but works when anchor is NULL
		start1 <- start2 <- end1 <- end2 <- integer(n + 1L)
		
		l1 <- 0L
		l2 <- 0L
		for (i in seq_len(n)) {
			start1[i] <- l1 + 1L
			end1[i] <- anchor[1L, i] - 1L
			l1 <- anchor[2L, i]
			start2[i] <- l2 + 1L
			end2[i] <- anchor[3L, i] - 1L
			l2 <- anchor[4L, i]
		}
		start1[n + 1L] <- l1 + 1L
		end1[n + 1L] <- widths[pair[1L]]
		start2[n + 1L] <- l2 + 1L
		end2[n + 1L] <- widths[pair[2L]]
		
		if (processors > 1L) {
			l <- as.integer(length(start1)/perProcessor)
			if (l < processors)
				processors <- max(1L, l)
		}
		.Call("alignPair",
			myXStringSet,
			pair,
			start1,
			end1,
			start2,
			end2,
			GO,
			GE,
			TG,
			maxLength,
			typeX,
			sM,
			processors,
			PACKAGE="DECIPHER")
	}
	
	sample2 <- function(x, prob, seed) {
		# generate a uniformly distributed 8 bit number
		seed <- as.raw(seed %% 256L)
		seed <- xor(seed, rawShift(seed, 7L))
		seed <- xor(seed, rawShift(seed, -5L))
		seed <- xor(seed, rawShift(seed, 3L))
		seed <- as.numeric(seed)/255 # [0, 1]
		
		# project onto the probability vector
		prob <- cumsum(prob)
		prob <- prob/prob[length(prob)]
		prob <- c(0, prob) # [0, 1]
		x[.bincode(seed, prob, include.lowest=TRUE)]
	}
	
	if (typeX == 3L) { # AAStringSet
		v <- .Call("enumerateSequenceReducedAA",
			.subset(myXStringSet, u),
			kmerSize,
			alphabet,
			FALSE, # mask repeats
			FALSE, # mask low complexity regions
			processors,
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		v <- .Call("enumerateSequence",
			.subset(myXStringSet, u),
			kmerSize,
			FALSE, # mask repeats
			FALSE, # mask low complexity regions
			processors,
			PACKAGE="DECIPHER")
	}
	sizes <- lengths(v)
	for (i in seq_along(v)) {
		o <- .Call("radixOrder", v[[i]], 1L, 0L, PACKAGE="DECIPHER")
		v[[i]] <- list(v[[i]][o], seq_along(v[[i]])[o])
	}
	KMERS <- as.integer(words^kmerSize)
	words <- as.integer(words^wordSize)
	lf <- as.integer(max(KMERS^pow, mean(sizes)))
	freqs <- .Call("sumBins", v, lf, PACKAGE="DECIPHER")
	
	# record groups sharing the rarest k-mers
	select <- as.integer(min(rareKmers,
		lf, # bounded by number of possible mapped k-mers
		max(sizes))) # bounded by number of possible k-mers
	j <- 0L
	z <- seq_len(select)
	y <- rep(KMERS, select) # use max k-mer + 1 as placeholder
	kmers <- integer(l*select)
	for (i in seq_len(l)) {
		s <- v[[i]][[1L]]
		s <- .Call("sortedUnique", s, PACKAGE="DECIPHER") # fast unique(s)
		if (length(s) == 0) {
			s <- y
		} else {
			if (length(s) > select) {
				o <- .Call("xorShift", s, lf, PACKAGE="DECIPHER")
				o <- .Call("radixOrder", freqs[o], 1L, 1L, PACKAGE="DECIPHER")
				s <- s[o[z]]
			} else if (length(s) < select) {
				s <- c(s, rep(KMERS, select - length(s)))
			} # else length(s) == select
		}
		kmers[j + z] <- s
		j <- j + select
	}
	rm(y, z)
	
	o <- order(kmers, method="radix")
	ini <- integer(length(kmers)) # index of initial k-mer
	prev <- -1L
	j <- 0L
	for (i in seq_along(o)) {
		curr <- kmers[o[i]]
		if (curr != prev) {
			if (prev == KMERS) { # placeholder group
				for (k in j:(i - 1L))
					kmers[o[k]] <- ini[o[k]] # length 1 group
			} else {
				prev <- i - 1L
				for (k in j:prev)
					kmers[o[k]] <- prev # reuse for fin
			}
			prev <- curr
			j <- i
		}
		ini[o[i]] <- j
	}
	if (kmers[o[j]] == KMERS) { # placeholder group
		for (k in j:length(o))
			kmers[o[k]] <- ini[o[k]] # length 1 group
	} else {
		for (k in j:length(o))
			kmers[o[k]] <- length(o)
	}
	for (i in seq_along(o))
		o[i] <- (o[i] - 1L) %/% select + 1L # convert to sequence index
	
	# Phase 1: Partition into groups
	partition <- integer(l)
	stack <- integer(l)
	pos <- 0L
	group <- 0L
	j <- 1L
	if (verbose) {
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
		count <- 0L
		lastValue <- 0
	}
	repeat {
		if (verbose) {
			count <- count + 1L
			value <- round(count/l, 2)
			if (value > lastValue) {
				lastValue <- value
				setTxtProgressBar(pBar, value)
			}
		}
		if (pos == 0L) {
			group <- group + 1L
			while (j <= l) {
				if (partition[j] == 0L)
					break # new partition
				j <- j + 1L
			}
			if (j > l)
				break # done
			i <- j
			partition[i] <- group
			first <- TRUE
		} else {
			i <- stack[pos]
			pos <- pos - 1L
			first <- FALSE
		}
		
		w <- (i - 1L)*select + seq_len(select)
		w <- w[.Call("radixOrder", kmers[w] - ini[w], 1L, 1L, PACKAGE="DECIPHER")]
		groups <- .Call("selectGroups",
			o,
			ini[w],
			kmers[w],
			maxPhase1,
			PACKAGE="DECIPHER")
		
		d <- .Call("radixOrder", groups, 1L, 1L, PACKAGE="DECIPHER")
		d <- .Call("dereplicate", groups, d, PACKAGE="DECIPHER")
		
		# fit exponential decrease in probability due to chance
		weights <- tabulate(d[[2L]], select)
		k1 <- which.max(weights)
		w <- weights[k1]
		if (k1 == select) # only one datapoint
			next
		for (k2 in (k1 + 1L):select) {
			if (weights[k2] > w)
				break
			w <- weights[k2]
			if (w == 0L)
				break
		}
		k2 <- k2 - 1L
		if (k1 == k2) # only one datapoint
			next
		w <- k1:k2 # monotonically decreasing from the max
		weights <- weights[w]
		log_weights <- log(weights)
		w_mean <- sum(w*weights)/sum(weights)
		log_weights_mean <- sum(log_weights*weights)/sum(weights)
		m <- sum(weights*(w - w_mean)*(log_weights - log_weights_mean))/sum(weights*(w - w_mean)^2)
		b <- log_weights_mean - m*w_mean
		probs <- 1 - exp(b + m*seq_len(select)) # can be negative
		
		# remove insufficient probabilities
		sufficient <- which.max(probs >= threshold) # if 1 then likely all < threshold
		w <- which(d[[2L]] >= sufficient)
		d[[1L]] <- d[[1L]][w]
		d[[2L]] <- d[[2L]][w]
		
		# apply multiple testing correction
		w <- .Call("radixOrder", d[[2L]], 0L, 1L, PACKAGE="DECIPHER")
		d[[1L]] <- d[[1L]][w]
		d[[2L]] <- d[[2L]][w]
		d[[2L]] <- probs[d[[2L]]] # convert to probabilities
		d[[2L]] <- cumprod(d[[2L]]) # multiple testing correction
		w <- which(d[[2L]] >= threshold)
		groups <- groups[d[[1L]][w]]
		
		# add unassigned groups to stack
		w <- which(partition[groups] == 0L)
		if (length(w) == 0L) {
			if (first && # singleton
				length(groups) > 1L) { # related groups already exist
				# assign singleton to highest weighted existing group
				w <- which(groups != i)
				w <- w[.Call("radixOrder", partition[groups[w]], 0L, 1L, PACKAGE="DECIPHER")]
				groups <- partition[groups[w]]
				partition[i] <- .Call("maxGroup", groups, d[[2L]][w], PACKAGE="DECIPHER")
				group <- group - 1L
			}
		} else {
			groups <- groups[w]
			partition[groups] <- group
			stack[pos + seq_along(groups)] <- groups
			pos <- pos + length(groups)
		}
	}
	rm(o, ini, kmers, stack, probs, sufficient)
	
	keep <- sizes > 0L
	O <- order(partition, decreasing=TRUE, method="radix")
	G <- .Call("indexPartitions", O, partition, keep, PACKAGE="DECIPHER")
	inPlay <- G[[2L]]
	G <- G[[1L]]
	ls <- lengths(G)
	
	if (verbose) {
		close(pBar)
		
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
		
		time.1 <- time.2
		cat("Sorting by relatedness within",
			length(G),
			ifelse(length(G) == 1L, "group:\n", "groups:\n"))
		flush.console()
	}
	
	# re-select rare k-mers after weighting partitions equally
	freqs <- numeric(lf)
	for (i in seq_along(ls)) {
		temp <- .Call("sumBins", v[G[[i]]], lf, PACKAGE="DECIPHER")
		w <- which(temp > 0) # randomly projected k-mer exists in group
		if (length(w) > 0L) {
			temp <- temp[w]
			temp <- temp/max(temp)
			freqs[w] <- freqs[w] + temp
		}
	}
	rm(temp)
	freqs <- as.integer(freqs)
	
	j <- 0L
	z <- seq_len(select)
	y <- rep(KMERS, select) # use max k-mer + 1 as placeholder
	kmers <- integer(l*select)
	for (i in seq_len(l)) {
		s <- v[[i]][[1L]]
		s <- .Call("sortedUnique", s, PACKAGE="DECIPHER")
		if (length(s) == 0) {
			s <- y
		} else {
			if (length(s) > select) {
				o <- .Call("xorShift", s, lf, PACKAGE="DECIPHER")
				o <- .Call("radixOrder", freqs[o], 1L, 1L, PACKAGE="DECIPHER")
				s <- s[o[z]]
			} else if (length(s) < select) {
				s <- c(s, rep(KMERS, select - length(s)))
			} # else length(s) == select
		}
		kmers[j + z] <- s
		j <- j + select
	}
	rm(y, z)
	
	if (wordSize != kmerSize ||
		(typeX == 3L && any(alphabet != 0:19))) {
		rm(v, sizes, freqs, lf)
		if (typeX == 3L) { # AAStringSet
			alphabet <- 0L:19L # switch to the standard alphabet
			entropy <- .Call("alphabetSizeReducedAA",
				myXStringSet,
				alphabet,
				PACKAGE="DECIPHER")
			maxK <- 7L
			wordSize <- ceiling(log(N*wordSizeQuantile, entropy))
			if (wordSize < 1)
				wordSize <- 1L
			if (wordSize > maxK)
				wordSize <- maxK
			words <- as.integer(20L^wordSize)
			
			v <- .Call("enumerateSequenceReducedAA",
				.subset(myXStringSet, u),
				wordSize,
				alphabet,
				FALSE, # mask repeats
				FALSE, # mask low complexity regions
				processors,
				PACKAGE="DECIPHER")
		} else { # DNAStringSet or RNAStringSet
			words <- as.integer(4L^wordSize)
			
			v <- .Call("enumerateSequence",
				.subset(myXStringSet, u),
				wordSize,
				FALSE, # mask repeats
				FALSE, # mask low complexity regions
				processors,
				PACKAGE="DECIPHER")
		}
		sizes <- lengths(v)
		for (i in seq_along(v)) {
			o <- .Call("radixOrder", v[[i]], 1L, 0L, PACKAGE="DECIPHER")
			v[[i]] <- list(v[[i]][o], seq_along(v[[i]])[o])
		}
		lf <- as.integer(words^pow)
		freqs <- numeric(lf)
		for (i in seq_along(ls)) {
			temp <- .Call("sumBins", v[G[[i]]], lf, PACKAGE="DECIPHER")
			w <- which(temp > 0) # randomly projected k-mer exist in group
			if (length(w) > 0L) {
				temp <- temp[w]
				temp <- temp/max(temp)
				freqs[w] <- freqs[w] + temp
			}
		}
		rm(temp)
		freqs <- as.integer(freqs)
	}
	
	# Phase 2: Relatedness sorting
	if (!any(inPlay)) {
		maxPhase2 <- 0L
	} else {
		maxPhase2 <- as.integer(min(maxPhase2, max(ls/2)))
	}
	ksims <- psims <- vector("list", maxPhase2)
	var <- pmin(maxPhase3, ls[partition]/2) # doubles to avoid overflow
	avg_cor <- sum(var < maxPhase3/2)/l # initialize to expected value by chance
	resTime1 <- resTime2 <- rep(NA_real_, processors)
	count1 <- count2 <- integer(processors)
	optProcessors1 <- processors
	optProcessors2 <- 1L # assume 1 is most efficient for alignment
	alignments <- integer(length(cutoff))
	for (i in seq_len(maxPhase2)) {
		if (verbose && interactive()) {
			cat("\riteration ",
				i,
				" of up to ",
				maxPhase2,
				" (",
				formatC(100*avg_cor, digits=1, format="f"),
				"% stability) ",
				sep="")
		}
		
		batches <- unique(ls)
		batches <- lapply(batches,
			function(l) {
				batches <- seq.int(0L, l, batchSize)
				if (batches[length(batches)] < l)
					batches <- c(batches, l)
				batches
			})[match(ls, batches)]
		
		if (processors > 1L) {
			if (i == 1L) {
				# do nothing (optProcessors1 == processors)
			} else if (i == 2L) {
				optProcessors1 <- processors - 1L
			} else if (max(count1) < minCount) {
				w <- which(!is.na(resTime1))
				optProcessors1 <- sample2(w, 1/resTime1[w], i)
				w <- w[1L]
				if (w > 1L && optProcessors1 == w)
					optProcessors1 <- w - 1L
			}
		}
		
		m1 <- m2 <- rep(NA_real_, l)
		rand <- matrix(NA_integer_, length(ls), 2L)
		time.3 <- Sys.time()
		for (k in which(inPlay)) {
			rand[k, 1L] <- G[[k]][sample(ls[k], 1L, replace=TRUE, prob=var[G[[k]]]*sizes[G[[k]]]*keep[G[[k]]])]
			keep[rand[k, 1L]] <- FALSE
			
			# process sequences in batches to reduce memory
			ov <- integer(ls[k])
			for (j in seq_len(length(batches[[k]]) - 1L)) {
				b <- (batches[[k]][j] + 1L):batches[[k]][j + 1L]
				s <- G[[k]][b]
				res1 <- overlap(rand[k, 1L], s, optProcessors1)
				ov[b] <- .Call("overlap",
					res1,
					wu[rand[k, 1L]],
					wu[s],
					PACKAGE="DECIPHER")
				m1[s] <- similarity(res1, wu[rand[k, 1L]], wu[s], optProcessors1, coverage=0)
			}
			
			rand[k, 2L] <- G[[k]][sample(ls[k], 1L, replace=TRUE, prob=(1.001 - m1[G[[k]]])*keep[G[[k]]]/ov)] # maximize distance and overlap to other random sequence
			keep[rand[k, 2L]] <- FALSE
			
			# process sequences in batches to reduce memory
			for (j in seq_len(length(batches[[k]]) - 1L)) {
				s <- G[[k]][(batches[[k]][j] + 1L):batches[[k]][j + 1L]]
				res2 <- overlap(rand[k, 2L], s, optProcessors1)
				m2[s] <- similarity(res2, wu[rand[k, 2L]], wu[s], optProcessors1, coverage=0)
			}
		}
		time.4 <- Sys.time()
		if (is.na(resTime1[optProcessors1])) {
			resTime1[optProcessors1] <- difftime(time.4, time.3, units='secs')/sum(ls[inPlay])
		} else {
			resTime1[optProcessors1] <- (count1[optProcessors1]*resTime1[optProcessors1] + difftime(time.4, time.3, units='secs')/sum(ls[inPlay]))/(count1[optProcessors1] + 1L)
		}
		count1[optProcessors1] <- count1[optProcessors1] + 1L
		optProcessors1 <- which.min(resTime1)
		
		if (any(alignments < minSimilarities)) {
			b1 <- .bincode(m1, bins, include.lowest=TRUE)
			w <- which(m1 < 1)
			n <- min(length(w), minSimilarities/2)
			if (n > 0) {
				s1 <- sample(length(w),
					n,
					replace=TRUE, # for speed
					prob=(counts/(1L + tabulate(b1, length(counts))))[b1[w]])
				s1 <- w[s1]
				s1 <- s1[!duplicated(s1)]
			} else {
				s1 <- integer()
			}
			b2 <- .bincode(m2, bins, include.lowest=TRUE)
			w <- which(m2 < 1)
			n <- min(length(w), minSimilarities/2)
			if (n > 0) {
				s2 <- sample(length(w),
					n,
					replace=TRUE, # for speed
					prob=(counts/(1L + tabulate(b2, length(counts))))[b2[w]])
				s2 <- w[s2]
				s2 <- s2[!duplicated(s2)]
			} else {
				s2 <- integer()
			}
			
			if (processors > 1L) {
				if (i == 1L) {
					optProcessors2 <- processors
				} else if (i == 2L) {
					optProcessors2 <- 1L
				} else {
					w <- which(is.na(resTime2))
					if (length(w) == 0L) {
						optProcessors2 <- which.min(resTime2)
					} else if (length(w) == 1L) {
						optProcessors2 <- w
					} else {
						o <- order(resTime2)
						optProcessors2 <- w[which.min(abs(mean(o[1:2]) - w))]
					}
				}
			}
			
			ksim1 <- m1[s1]
			pdist1 <- numeric(length(s1))
			time.5 <- Sys.time()
			for (j in seq_along(s1)) {
				res <- overlap(rand[partition[s1[j]], 1L], s1[j], optProcessors1) # recompute
				if (minCoverage != 0) # need to replace k-mer similarity
					ksim1[j] <- similarity(res, wu[rand[partition[s1[j]], 1L]], wu[s1[j]], optProcessors1)
				ali <- align(c(u[rand[partition[s1[j]], 1L]], u[s1[j]]),
					anchor=res[[1L]],
					processors=optProcessors2)
				pdist1[j] <- dist(ali)
			}
			time.6 <- Sys.time()
			
			ksim2 <- m2[s2]
			pdist2 <- numeric(length(s2))
			time.7 <- Sys.time()
			for (j in seq_along(s2)) {
				res <- overlap(rand[partition[s2[j]], 2L], s2[j], optProcessors1) # recompute
				if (minCoverage != 0) # need to replace k-mer similarity
					ksim2[j] <- similarity(res, wu[rand[partition[s2[j]], 2L]], wu[s2[j]], optProcessors1)
				ali <- align(c(u[rand[partition[s2[j]], 2L]], u[s2[j]]),
					anchor=res[[1L]],
					processors=optProcessors2)
				pdist2[j] <- dist(ali)
			}
			time.8 <- Sys.time()
			if (is.na(resTime2[optProcessors2])) {
				resTime2[optProcessors2] <- difftime(time.8, time.7, units='secs') + difftime(time.6, time.5, units='secs')
			} else {
				resTime2[optProcessors2] <- (count2[optProcessors2]*resTime2[optProcessors2] + difftime(time.8, time.7, units='secs') + difftime(time.6, time.5, units='secs'))/(count2[optProcessors2] + 1L)
			}
			count2[optProcessors2] <- count2[optProcessors2] + 1L
			
			w <- .bincode(pdist1, bins, include.lowest=TRUE)
			w <- which(w %in% incBins)
			counts <- counts + (tabulate(b1[s1[w]], length(counts)) > 0)
			w <- .bincode(pdist2, bins, include.lowest=TRUE)
			w <- which(w %in% incBins)
			counts <- counts + (tabulate(b2[s2[w]], length(counts)) > 0)
			ksims[[i]] <- c(ksim1, ksim2)
			psims[[i]] <- 1 - c(pdist1, pdist2)
			
			for (k in seq_along(cutoff)) {
				w <- which(psims[[i]] >= cutoff[k] - binSize & psims[[i]] <= cutoff[k] + binSize)
				alignments[k] <- alignments[k] + length(w)
			}
		}
		
		d <- m1 - m2
		if (i > 1L) {
			for (k in which(inPlay)) {
				if (diff(range(rS[G[[k]]])) > 0) {
					if (diff(range(d[G[[k]]])) > 0) {
						temp <- suppressWarnings(cov(rS[G[[k]]], d[G[[k]]]))
					} else {
						temp <- NA
					}
					if (!is.na(temp)) {
						# calculate first principle component
						temp <- matrix(c(var(rS[G[[k]]]), temp, temp, var(d[G[[k]]])), nrow=2L)
						temp <- eigen(temp)$vectors
						temp <- temp[1L]*(rS[G[[k]]] - mean(rS[G[[k]]])) + temp[2L]*(d[G[[k]]] - mean(d[G[[k]]]))
						if (cor(temp, rS[G[[k]]]) < 0)
							temp <- -temp # flip sign of principle component
						rS[G[[k]]] <- temp
					}
				} else { # no variation in relative similarity
					rS[G[[k]]] <- d[G[[k]]]
				}
			}
			
			o <- order(partition,
				rS,
				wu,
				t, # frequency
				decreasing=TRUE,
				method="radix")
			r <- o
			r[r] <- seq_along(r) # rank
			
			w <- which(d != 0) # only update positions that could have moved
			var[w] <- alpha1*pmin(maxPhase3, abs(R[w] - r[w])) + (1 - alpha1)*var[w] # exponential moving average of change in rank order
			avg_cor <- sum(var*keep < maxPhase3/2)/l # exponential moving average of fraction stabilized in rank order
			
			O <- o
			R <- r
		} else {
			rS <- d # relative similarity
			
			# initialize rank order
			O <- order(partition,
				rS,
				wu,
				t, # frequency
				decreasing=TRUE,
				method="radix")
			R <- O
			R[R] <- seq_along(R) # rank
		}
		
		if (i < maxPhase2) { # possibly split partitions
			partition <- .Call("splitPartitions",
				O, # ordering partition
				partition, # maybe in-place change of `partition`
				var, # split when var <= split
				maxPhase3, # minimum split partition size
				split, # split threshold
				PACKAGE="DECIPHER")
			
			G <- .Call("indexPartitions", O, partition, keep, PACKAGE="DECIPHER")
			inPlay <- G[[2L]]
			G <- G[[1L]]
			ls <- lengths(G)
		}
		
		for (k in which(inPlay)) {
			if (all(alignments >= minSimilarities) &&
				all(var[G[[k]]] < maxPhase3/2)) {
				inPlay[k] <- FALSE # 100% stability
			}
		}
		
		if (!any(inPlay) ||
			(avg_cor >= 0.9995 && # rounds to 100%
			(all(alignments >= minSimilarities) || # sufficient alignments
			all(alignments*maxPhase2/i < minSimilarities)))) # insufficient sequences
			break
	}
	if (maxPhase2 > 0L) {
		optProcessors2 <- which.min(resTime2)
	} else {
		O <- R <- seq_along(u)
	}
	
	# fit similarity limit based on k-mer similarity
	ksims <- unlist(ksims)
	limits <- rep(1e-6, length(cutoff)) # initialize above zero
	if (length(ksims) > 1) {
		psims <- unlist(psims)
		psims[is.na(psims)] <- 0
		z <- c(0, sqrt(ksims), 1) # perform logistic regression in sqrt-space
		w <- .bincode(z, bins, include.lowest=TRUE)
		w <- (1/tabulate(w, length(counts)))[w] # weight bins equally
		fit <- function(ksim) {
			p <- predict(g,
				newdata=data.frame(z=ksim))
			abs(p - interval)
		}
		for (j in which(alignments >= needAlignments)) {
			y <- psims >= cutoff[j]
			if (sum(y) >= 1) {
				y <- c(FALSE, y, TRUE) # add fixed endpoints
				
				g <- suppressWarnings(glm(y ~ z,
					binomial(link="logit"),
					weights=w,
					control=glm.control(1e-6)))
				limits[j] <- optimize(fit, c(0, 1))$minimum # exclusive of bounds
				limits[j] <- limits[j]^2
			} else { # no data above above cutoff
				limits[j] <- max(ksims)
			}
		}
	}
	if (lc > 1) {
		if (ASC) {
			limits <- rev(cummax(rev(limits)))
		} else {
			limits <- cummax(limits)
		}
	}
	
	maxComparisons <- boundComparisons*maxPhase3
	minComparisons <- ceiling(maxPhase3/boundComparisons)
	var[var > maxComparisons] <- maxComparisons
	var[var < minComparisons] <- minComparisons
	var <- var[O]/mean(var) # normalized variability
	bL <- seq_len(l) - as.integer(var*(maxPhase3/2))
	bR <- seq_len(l) + as.integer(var*(maxPhase3/2))
	bL[bL < 1L] <- 1L
	bR[bR > l] <- l
	
	if (verbose) {
		if (maxPhase2 > 0L) {
			cat("\riteration ",
				i,
				" of up to ",
				maxPhase2,
				" (",
				formatC(100*avg_cor, digits=1, format="f"),
				"% stability) ",
				sep="")
			
			time.2 <- Sys.time()
			cat("\n\n")
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\n")
			
			time.1 <- time.2
		}
		lastValue <- 0
		cat("Clustering sequences by ",
			wordSize,
			"-mer similarity:\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	}
	
	if (maxPhase2 > 0) {
		rm(b1, b2, bins, d, inPlay, keep, ksims, psims, ls, m1, m2, ov, ksim1, ksim2, pdist1, pdist2, res1, res2, s1, s2, var, b, rS)
	} else {
		rm(bins, inPlay, keep, ksims, psims, ls, var)
	}
	
	partition <- partition[O]
	sizes <- sizes[O]
	t <- t[O]
	P <- order(sizes,
		t,
		partition,
		decreasing=TRUE,
		method="radix")
	P <- O[P]
	Q <- R[P]
	
	# order rare k-mers into groups
	E <- seq_len(l)
	E[P] <- E
	E <- rep(E, each=select)
	o <- order(kmers, E, method="radix")
	rm(E)
	ini <- integer(length(kmers)) # index of initial k-mer
	prev <- -1L
	j <- 0L
	for (i in seq_along(o)) {
		curr <- kmers[o[i]]
		if (curr != prev) {
			prev <- curr
			j <- i
		}
		ini[o[i]] <- j
	}
	prev <- -1L
	for (i in seq_along(o)) {
		curr <- ini[o[i]]
		if (curr != prev) {
			j <- 0L
			prev <- curr
		} else {
			j <- j + 1L
		}
		if (kmers[o[i]] == KMERS) { # missing k-mer
			kmers[o[i]] <- 0L # zero length group
		} else{
			kmers[o[i]] <- j # reuse as fin
		}
	}
	for (i in seq_len(l)) { # prioritize smaller k-mer groups
		j <- seq((i - 1L)*select + 1L, length.out=select)
		k <- j[.Call("radixOrder", kmers[j], 1L, 1L, PACKAGE="DECIPHER")]
		ini[j] <- ini[k]
		kmers[j] <- kmers[k]
	}
	for (i in seq_along(o))
		o[i] <- (o[i] - 1L) %/% select + 1L # convert to sequence index
	
	# Phase 3: Cluster sequences
	for (i in seq_len(lc))
		C[[i]] <- integer(L)
	V <- v # original ordering of `v`
	origin <- c(1, 1) # origin of existing clusters
	subsetRate <- 0 # rate of subsetting
	if (verbose) {
		matches <- integer(maxPhase3) # relative location of cluster matches
		totals <- integer(4L) # count of match origins
	}
	for (i in seq_len(lc)) {
		c <- C[[i]]
		C[[i]] <- list(NULL) # release memory
		if (!ASC && i > 1) {
			if (invertCenters) {
				P <- order(abs(C[[i - 1L]][u][O]),
					sizes,
					t,
					partition,
					decreasing=TRUE,
					method="radix")
			} else {
				P <- order(C[[i - 1L]][u][O],
					sizes,
					t,
					partition,
					decreasing=TRUE,
					method="radix")
			}
			P <- O[P]
			Q <- R[P]
		}
		
		v <- V[P] # reorder k-mers
		p <- u[P] # index of sequence
		
		cluster_num <- 1L
		offset <- 0L
		if (invertCenters) {
			c[p[1L]] <- -1L
		} else {
			c[p[1L]] <- 1L
		}
		seeds.index <- integer(l)
		seeds.index[Q[1L]] <- 1L
		recent <- integer(maxPhase3)
		count <- 1L
		recent[count] <- 1L
		
		j <- 1L
		while (j < l) {
			if (verbose) {
				value <- round(((i - 1)*l + j)/lc/l, 2)
				if (value > lastValue) {
					lastValue <- value
					setTxtProgressBar(pBar, value)
				}
			}
			j <- j + 1L
			
			if (!ASC && i > 1L) {
				if ((invertCenters &&
					abs(C[[i - 1L]][p[j]]) != abs(C[[i - 1L]][p[j - 1L]])) ||
					(!invertCenters &&
					C[[i - 1L]][p[j]] != C[[i - 1L]][p[j - 1L]])) {
					# different clusters in last cutoff
					offset <- offset + cluster_num
					cluster_num <- 1L
					c[p[j]] <- cluster_num + offset
					if (invertCenters)
						c[p[j]] <- -c[p[j]]
					recent <- integer(maxPhase3)
					count <- 1L
					recent[count] <- j
					seeds.index[Q[j]] <- j
					next
				} else if (c[p[j]] > 0L) { # cluster pre-assigned
					if (singleLinkage) {
						seeds.index[Q[j]] <- j
					} else {
						seeds.index[Q[j]] <- seeds.index[R[c[p[j]]]]
					}
					if (invertCenters) {
						c[p[j]] <- abs(c[u[c[p[j]]]])
					} else {
						c[p[j]] <- c[u[c[p[j]]]]
					}
					next
				}
			}
			
			compare <-  c(bL[Q[j]]:bR[Q[j]], # surrounding centers
				Q[recent]) # most recent centers
			compare <- compare[.Call("radixOrder", abs(Q[j] - compare), 1L, 1L, PACKAGE="DECIPHER")] # order by proximity
			compare <- seeds.index[compare]
			compare <- compare[compare != 0L] # remove unvisited seeds
			compare <- compare[!duplicated(compare)]
			
			# choose the most frequent sequences sharing rare k-mers
			groups <- .Call("selectIndices",
				P[j], # focal sequence
				select, # number of k-mers possible
				ini, # initial index of each k-mer group
				kmers, # number of k-mers to record in each group
				o, # ordering of k-mers
				maxPhase1, # maximum number of sequences to select
				PACKAGE="DECIPHER")
			d <- .Call("radixOrder", groups, 1L, 1L, PACKAGE="DECIPHER")
			d <- .Call("dereplicate", groups, d, PACKAGE="DECIPHER")
			d <- d[[1L]][.Call("radixOrder", d[[2L]], 0L, 1L, PACKAGE="DECIPHER")]
			if (length(d) > maxPhase3)
				length(d) <- maxPhase3
			groups <- groups[d]
			groups <- seeds.index[R[groups]]
			groups <- groups[!duplicated(groups)]
			
			if (!ASC && i > 1) { # check bounds
				groups <- groups[groups != 0L]
				if (invertCenters) {
					compare <- compare[abs(C[[i - 1L]][p[compare]]) == abs(C[[i - 1L]][p[j]])]
					groups <- groups[abs(C[[i - 1L]][p[groups]]) == abs(C[[i - 1L]][p[j]])]
				} else {
					compare <- compare[C[[i - 1L]][p[compare]] == C[[i - 1L]][p[j]]]
					groups <- groups[C[[i - 1L]][p[groups]] == C[[i - 1L]][p[j]]]
				}
			}
			
			if (minCoverage < 0) { # eliminate representatives that cannot meet minCoverage
				compare <- compare[wu[P[j]]/wu[P[compare]] >= -minCoverage]
				groups <- groups[wu[P[j]]/wu[P[groups]] >= -minCoverage]
			}
			
			if (length(compare) == 0L && length(groups) == 0L) {
				cluster_num <- cluster_num + 1L
				c[p[j]] <- cluster_num + offset
				if (invertCenters)
					c[p[j]] <- -c[p[j]]
				seeds.index[Q[j]] <- j
				count <- count + 1L
				if (count > maxPhase3)
					count <- 1L
				recent[count] <- j
				next
			}
			
			num <- bR[Q[j]] - bL[Q[j]] + 1L # total sequences to consider
			prop <- origin/sum(origin)
			num <- as.integer(ceiling(num*prop))
			if (num[2L] > length(groups)) {
				num[1L] <- num[1L] + num[2L] - length(groups)
			} else if (num[1L] > length(compare)) {
				num[2L] <- num[2L] + num[1L] - length(compare)
			}
			compare <- head(compare, num[1L])
			groups <- head(groups, num[2L])
			combined <- c(compare, groups)
			combined <- combined[!duplicated(combined)]
			
			# narrow to comparisons with most hits
			s <- unique(v[[j]][[1L]])
			if (length(s) > maxSample/(1 - subsetRate) || # worthwhile
				(j %% updateRate == 0 && length(s) > maxSample)) { # need to update
				m <- .Call("xorShift", s, lf, PACKAGE="DECIPHER")
				m <- .Call("radixOrder", freqs[m], 1L, 1L, PACKAGE="DECIPHER")
				length(m) <- maxSample
				s <- s[m]
				s <- s[.Call("radixOrder", s, 1L, 1L, PACKAGE="DECIPHER")]
				counts <- countHits(s, v[combined], processors)
				counts <- counts/sizes[Q[combined]]
				w <- which(!is.na(counts))
				if (length(w) > 0) {
					w <- which(counts >= max(counts[w])*minFraction)
					subsetRate <- (1 - alpha2)*subsetRate + alpha2*length(w)/length(combined)
					combined <- combined[w]
				}
			}
			
			res <- overlap(j, combined, optProcessors1)
			m <- similarity(res, wu[P[j]], wu[P[combined]], optProcessors1)
			w <- sort.list(m, method="radix", decreasing=TRUE)
			if (length(w) > maxAlignments)
				length(w) <- maxAlignments
			W <- which(m[w] >= limits[i])
			if (length(W) >= minAlignments) {
				w <- w[W]
			} else {
				w <- head(w, minAlignments)
			}
			
			if (m[w[1L]] < cutoff[i]) {
				for (k in seq_along(w)) {
					if (k > attempts) {
						s <- seq(to=k - 1L, length.out=attempts)
						if (max(m[w[s]]) + stdDevs*sd(m[w[s]]) < cutoff[i]) {
							w <- w[seq_len(k - 1L)] # prevent matching clusters without alignment
							break
						}
					}
					ali <- align(c(p[j], p[combined[w[k]]]),
						anchor=res[[w[k]]],
						processors=optProcessors2)
					m[w[k]] <- dist(ali)
					if (is.na(m[w[k]])) {
						m[w[k]] <- 0
					} else {
						m[w[k]] <- 1 - m[w[k]]
					}
					if (m[w[k]] >= cutoff[i])
						break
				}
				w <- w[which.max(m[w])]
			} else {
				w <- w[1L] # max similarity
			}
			
			if (length(w) == 0L ||
				m[w] < cutoff[i]) { # form a new group
				cluster_num <- cluster_num + 1L
				c[p[j]] <- cluster_num + offset
				if (invertCenters)
					c[p[j]] <- -c[p[j]]
				count <- count + 1L
				if (count > maxPhase3)
					count <- 1L
				recent[count] <- j
				seeds.index[Q[j]] <- j
			} else { # part of an existing group
				num <- c(match(combined[w], compare),
					match(combined[w], groups))
				isna <- is.na(num)
				if (verbose) {
					if (!isna[1L]) {
						totals[1L] <- totals[1L] + 1L
						if (isna[2L])
							totals[3L] <- totals[3L] + 1L
						if (num[1L] <= maxPhase3 &&
							(num[1L] < num[2L] || isna[2L]))
							matches[num[1L]] <- matches[num[1L]] + 1L
					}
					if (!isna[2L]) {
						totals[2L] <- totals[2L] + 1L
						if (isna[1L])
							totals[4L] <- totals[4L] + 1L
						if (num[2L] <= maxPhase3 &&
							(num[2L] <= num[1L] || isna[1L]))
							matches[num[2L]] <- matches[num[2L]] + 1L
					}
				}
				if (isna[1L] || num[1L] > maxPhase3)
					num[1L] <- maxPhase3
				if (isna[2L] || num[2L] > maxPhase3)
					num[2L] <- maxPhase3
				origin <- (1 - alpha2)*origin + alpha2/num
				if (invertCenters) {
					c[p[j]] <- -c[p[combined[w]]]
				} else {
					c[p[j]] <- c[p[combined[w]]]
				}
				if (!ASC && i < lc) {
					cols <- (i + 1):lc
					cols <- cols[which(m[w] >= cutoff[cols])]
					for (k in cols) # assign forward
						C[[k]][p[j]] <- P[combined[w]]
				}
				if (singleLinkage) {
					seeds.index[Q[j]] <- j
				} else {
					seeds.index[Q[j]] <- combined[w]
				}
			}
		}
		
		C[[i]] <- c[x]
	}
	rm(P, O, R, Q, ini, kmers, partition)
	
	C <- as.data.frame(C)
	if (!is.null(names(myXStringSet))) {
		dNames <- names(myXStringSet)
		w <- which(duplicated(dNames))
		if (length(w) > 0) {
			warning("Duplicated names of myXStringSet appended with index.")
			dNames[w] <- paste(dNames[w],
				w,
				sep="_")
		}
		rownames(C) <- dNames
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		
		num <- totals[1L] + totals[4L]
		if (num > 0 && l > maxPhase3) {
			cat("\nClusters via relatedness sorting: ",
				100*round(totals[1L]/num, 3),
				"% (",
				100*round(totals[3L]/num, 3),
				"% exclusively)",
				sep="")
			cat("\nClusters via rare ",
				kmerSize,
				"-mers: ",
				100*round(totals[2L]/num, 3),
				"% (",
				100*round(totals[4L]/num, 3),
				"% exclusively)",
				sep="")
			
			X <- seq_along(matches)
			Y <- cumsum(matches)
			Y <- Y/Y[length(Y)]
			Zipf <- function(s) {
				H <- cumsum(X^s)
				y <- H/H[length(H)]
				sum(X^3*abs(Y - y)) # up weight tail
			}
			s <- optimize(Zipf, c(-10, 0))$minimum
			
			Hy <- 0
			for (i in seq_len(maxPhase3))
				Hy <- Hy + i^s
			Hl <- Hy
			for (i in seq(maxPhase3 + 1L, l))
				Hl <- Hl + i^s
			
			cat("\nEstimated clustering effectiveness: ",
				round(100*Hy/Hl, 1),
				"%\n\n",
				sep="")
		} else {
			cat("\nClusters via relatedness sorting: 100% (0% exclusively)",
				"\nClusters via rare ",
				kmerSize,
				"-mers: 100% (0% exclusively)",
				"\nEstimated clustering effectiveness: 100%\n\n",
				sep="")
		}
	}
	
	return(C)
}

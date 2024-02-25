FindNonCoding <- function(x,
	myXStringSet,
	minScore=16,
	allScores=FALSE,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	if (length(myXStringSet) == 0L)
		stop("myXStringSet must contain sequences.")
	if (sum(width(myXStringSet)) > 1073741823L)
		stop("Total length of sequences in myXStringSet is too long.")
	if (is(x, "list")) {
		if (any(unlist(lapply(x, class)) != "NonCoding"))
			stop("x must be an object of class 'NonCoding'.")
	} else if (is(x, "NonCoding")) {
		x <- list(x)
	} else {
		stop("x must be an object of class 'NonCoding'.")
	}
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') are not allowed in myXStringSet.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') are not allowed in myXStringSet.")
	if (!is.numeric(minScore))
		stop("minScore must be a numeric.")
	if (any(minScore < 0))
		stop("minScore must be at least zero.")
	if (length(minScore) == 1L) {
		minScore <- rep(minScore, length(x))
	} else if (length(minScore) != length(x)) {
		stop("minScore must be a single numeric.")
	}
	if (!is.logical(allScores))
		stop("allScores must be a logcial.")
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
		stop("verbose must be a logcial.")
	
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	}
	
	# initialize parameters
	pad <- 10 # nucleotides to left and right of matches to search for hairpins
	minWindowSize <- 1e4L # size of k-mer background sliding window
	mult1 <- 0.6 # fraction of minScore required without hairpins
	mult2 <- 0.4 # fraction of minScore required from motifs alone
	mult3 <- 1 # fraction of minScore required from motifs and hairpins
	minConsider <- -10 # score required for consideration of begin or end
	myXString <- unlist(myXStringSet)
	myXString <- DNAString(myXString)
	myXString <- c(myXString,
		reverseComplement(myXString))
	l <- length(myXString)
	
	ions <- 1
	temp <- 37
	data("deltaHrulesRNA", envir=environment(), package="DECIPHER")
	data("deltaSrulesRNA", envir=environment(), package="DECIPHER")
	deltaSrulesRNA <- deltaSrulesRNA + 0.368*log(ions)/1000
	deltaGrulesRNA <- deltaHrulesRNA - (273.15 + temp)*deltaSrulesRNA
	dG_ini <- 4.065225 # 3.6 - (273.15 + temp)*(-1.5/1000 + 0.368*log(ions)/1000)
	max_dG <- 0
	
	K <- sapply(x, attr, which="K")
	annotations <- names(x)
	O <- order(K)
	x <- x[O]
	K <- K[O]
	minScore <- minScore[O]
	
	results <- matrix(0,
		nrow=0L,
		ncol=6L,
		dimnames=list(NULL,
			c("Index", "Strand", "Begin", "End", "TotalScore", "Gene")))
	
	.getIndex <- function(s,
		offset1,
		offset2,
		score) {
		
		s1 <- s + as.integer(offset1)
		s2 <- s + as.integer(offset2)
		
		.Call("getIndex",
			s1,
			s2,
			l,
			score,
			PACKAGE="DECIPHER")
	}
	
	for (k in seq_along(x)) {
		minLength <- attr(x[[k]], "minLength")
		maxLength <- attr(x[[k]], "maxLength")
		lenScores <- x[[k]]$lengthScores
		depend <- x[[k]]$dependence
		params <- attr(x[[k]], "background")
		mean <- params[1]
		if (any(is.na(mean))) {
			minS <- minScore[k]
		} else {
			sd <- params[2]
			minS <- qlnorm(exp(-minScore[k]),
				mean,
				sd,
				lower.tail=FALSE)
			if (minS < minScore[k])
				minS <- minScore[k]
		}
		
		windowSize <- maxLength*20L # maximally 5% of region
		if (windowSize < minWindowSize)
			windowSize <- minWindowSize
		
		if (k == 1L || K[k - 1L] != K[k]) {
			ints <- .Call("enumerateSequence",
				DNAStringSet(myXString),
				K[[k]],
				TRUE, # mask repeats
				FALSE, # mask low complexity regions
				1L, # left is fast moving side
				processors,
				PACKAGE="DECIPHER")[[1L]]
			
			ints <- ints + 1L
		}
		
		oligos <- x[[k]][[3L]]
		oligos <- oligos/sum(oligos)
		
		kmers <- .Call("kmerScores",
			oligos,
			ints,
			windowSize,
			K[k])
		
		starts <- numeric(l)
		ends <- numeric(l)
		
		# find motifs
		motifs <- x[[k]]$motifs
		n1 <- nrow(motifs)
		if (!is.null(depend))
			indices <- vector("list", n1)
		for (i in seq_len(n1)) {
			bins <- motifs[i, "minscore"][[1L]]
			m <- log(motifs[i, "pwm"][[1L]]/0.25)
			p <- .Call("scorePWM",
				m,
				myXString,
				bins[1L],
				processors,
				PACKAGE="DECIPHER")
			
			index_s <- .getIndex(p[[1L]],
				-motifs[i, "begin_high"],
				-motifs[i, "begin_low"],
				.bincode(p[[2L]], bins))
			index_e <- .getIndex(p[[1L]],
				motifs[i, "end_low"] + ncol(m) - 1L,
				motifs[i, "end_high"] + ncol(m) - 1L,
				.bincode(p[[2L]], bins))
			
			if (!is.null(depend))
				indices[[i]] <- list(index_s, index_e)
			
			delta_starts <- motifs[i, "begin_high"] - motifs[i, "begin_low"]
			delta_ends <- motifs[i, "end_high"] - motifs[i, "end_low"]
			if (delta_starts < delta_ends) {
				bg <- tabulate(index_s, length(bins) - 1L)
			} else {
				bg <- tabulate(index_e, length(bins) - 1L)
			}
			
			bg[bg == 0L] <- 1L # pseudocount
			bg <- bg/l
			bg <- c(1 - sum(bg), bg)
			
			scores <- motifs[i, "prevalence"][[1L]]/bg
			scores <- log(scores)
			
			if (scores[1L] < 0) {
				weights <- c(delta_starts + 1L,
					delta_ends + 1L)
				weights <- 1 - weights/sum(weights)
				
				# in-place change of starts
				starts <- .Call("addIfElse",
					starts,
					index_s,
					scores*weights[1L],
					PACKAGE="DECIPHER")
				# in-place change of ends
				ends <- .Call("addIfElse",
					ends,
					index_e,
					scores*weights[2L],
					PACKAGE="DECIPHER")
			}
		}
		
		# check whether there are any potential matches
		result <- .Call("getBounds",
			width(myXStringSet),
			starts,
			ends,
			minLength,
			maxLength,
			lenScores,
			kmers,
			K[[k]],
			-O[k],
			minS*mult1,
			minS*mult2,
			minConsider,
			PACKAGE="DECIPHER")
		
		if (nrow(result) == 0L) {
			if (verbose)
				setTxtProgressBar(pBar, k/length(x))
			next
		}
		
		if (!is.null(depend))
			indices <- lapply(indices,
				function(x)
					pmax(x[[1L]][result[, 3]],
						x[[2L]][result[, 4]]) + 1L)
		
		# find hairpins
		hairpins <- x[[k]]$hairpins
		n2 <- nrow(hairpins)
		if (n2 > 0) {
			if (!is.null(depend))
				indices <- c(indices,
					vector("list", n2))
			
			# select positions of potential matches
			r <- result
			r[, 3] <- r[, 3] - pad
			w <- which(r[, 3] < 1)
			if (length(w) > 0)
				r[w, 3] <- 1
			r[, 4] <- r[, 4] + pad
			w <- which(r[, 4] > l)
			if (length(w) > 0)
				r[w, 4] <- l
			r <- IRanges(r[, 3], r[, 4])
			r <- reduce(r) # union of all ranges
			myXSubseq <- unlist(extractAt(myXString, r))
			ce <- cumsum(width(r)) # ends in myXString
			cs <- c(1L, ce[-length(ce)] + 1L) # starts in myXString
			offset <- start(r) - cs # offset from start in myXString
			
			maxLoopLength <- max(result[, 8] - result[, 7] + 5)
			maxLoopLength <- min(maxLoopLength,
				attr(x[[k]], "maxLoopLength"))
			
			pals <- findPalindromes(myXSubseq,
				min.armlength=4,
				max.looplength=maxLoopLength,
				min.looplength=3,
				max.mismatch=1,
				allow.wobble=TRUE)
			p_starts <- start(pals)
			p_widths <- width(pals)
			p_ends <- p_starts + p_widths - 1L
			
			o <- order(p_starts)
			p_starts <- p_starts[o]
			p_widths <- p_widths[o]
			p_ends <- p_ends[o]
			pals <- pals[o]
			
			# shift positions from myXSubseq to myXString
			w <- integer(length(p_starts))
			current <- 1L
			for (j in seq_along(p_starts)) {
				while (p_starts[j] > ce[current])
					current <- current + 1L
				if (p_ends[j] <= ce[current])
					w[j] <- current
			}
			keep <- which(w != 0)
			
			if (length(keep) > 0) {
				w <- w[keep]
				pals <- pals[keep]
				p_starts <- p_starts[keep] + offset[w]
				p_ends <- p_ends[keep] + offset[w]
				p_widths <- p_widths[keep]
				
				dna <- DNAStringSet(pals)
				p_arms <- palindromeArmLength(dna,
					max.mismatch=1,
					allow.wobble=TRUE)
				max_arms <- as.integer((p_widths - 3)/2)
				w <- which(p_arms > max_arms)
				if (length(w) > 0)
					p_arms[w] <- max_arms[w]
				
				p_dG <- .Call("calculateHairpinDeltaG",
					dna,
					p_arms + 1L, # include end bases
					deltaGrulesRNA,
					PACKAGE="DECIPHER") + dG_ini
				
				w <- which(p_dG <= max_dG)
				
				if (length(w) > 0) {
					p_starts <- p_starts[w]
					p_ends <- p_ends[w]
					p_dG <- p_dG[w]
					
					o <- order(p_starts)
					p_starts <- p_starts[o]
					p_ends <- p_ends[o]
					p_dG <- p_dG[o]
					p_widths <- p_ends - p_starts + 1L
					
					o <- order(p_ends)
					p_ends2 <- p_ends[o]
					p_dG2 <- p_dG[o]
					p_widths2 <- p_widths[o]
					
					s <- as.integer(result[, 3])
					e <- as.integer(result[, 4])
					for (i in seq_len(n2)) {
						delta_start <- hairpins[i, "begin_high"] - hairpins[i, "begin_low"]
						delta_end <- hairpins[i, "end_high"] - hairpins[i, "end_low"]
						delta_width <- hairpins[i, "width_high"] - hairpins[i, "width_low"]
						if (delta_start < delta_end) {
							if (delta_end < delta_width) {
								index <- .Call("getHits",
									p_starts,
									p_ends,
									s + as.integer(hairpins[i, "begin_low"]),
									s + as.integer(hairpins[i, "begin_high"]),
									e - as.integer(hairpins[i, "end_high"]),
									e - as.integer(hairpins[i, "end_low"]),
									p_dG,
									PACKAGE="DECIPHER")
								index[index != 0] <- p_dG[index[index != 0]]
							} else {
								index <- .Call("getHits",
									p_starts,
									p_widths,
									s + as.integer(hairpins[i, "begin_low"]),
									s + as.integer(hairpins[i, "begin_high"]),
									rep(as.integer(hairpins[i, "width_low"]),
										length(s)),
									rep(as.integer(hairpins[i, "width_high"]),
										length(s)),
									p_dG,
									PACKAGE="DECIPHER")
								index[index != 0] <- p_dG[index[index != 0]]
							}
						} else {
							if (delta_start < delta_width) {
								index <- .Call("getHits",
									p_starts,
									p_ends,
									s + as.integer(hairpins[i, "begin_low"]),
									s + as.integer(hairpins[i, "begin_high"]),
									e - as.integer(hairpins[i, "end_high"]),
									e - as.integer(hairpins[i, "end_low"]),
									p_dG,
									PACKAGE="DECIPHER")
								index[index != 0] <- p_dG[index[index != 0]]
							} else {
								index <- .Call("getHits",
									p_ends2,
									p_widths2,
									e - as.integer(hairpins[i, "end_high"]),
									e - as.integer(hairpins[i, "end_low"]),
									rep(as.integer(hairpins[i, "width_low"]),
										length(e)),
									rep(as.integer(hairpins[i, "width_high"]),
										length(e)),
									p_dG2,
									PACKAGE="DECIPHER")
								index[index != 0] <- p_dG2[index[index != 0]]
							}
						}
						
						index <- .bincode(index,
							hairpins[i, "dG"][[1L]])
						if (!is.null(depend))
							indices[[i + n1]] <- index
						
						scores <- hairpins[i, "prevalence"][[1L]]/hairpins[i, "background"][[1L]]
						scores <- log(scores)
						result[, 5] <- result[, 5] + scores[index] # TotalScore
						result[, 9] <- result[, 9] + scores[index] # score from patterns alone
					}
				}
			}
		}
		
		if (!is.null(depend)) {
			# add score for dependencies among patterns
			count <- 0L
			for (i in seq_len(n1 + n2 - 1L)) {
				for (j in (i + 1L):(n1 + n2)) {
					count <- count + 1L
					if (!is.null(depend[[count]])) {
						index <- cbind(indices[[i]], indices[[j]])
						result[, 5] <- result[, 5] + depend[[count]][index]
					}
				}
			}
		}
		
		result <- result[result[, 9] >= minS*mult3,, drop=FALSE]
		result <- result[result[, 5] >= minS,, drop=FALSE]
		
		# replace coordinates in myXString with myXStringSet
		result[, 3:4] <- result[, 7:8]
		result <- result[, -7:-9, drop=FALSE]
		
		if (!is.na(mean)) # transform log-odds scores
			result[, 5] <- -plnorm(result[, 5],
				mean,
				sd,
				lower.tail=FALSE,
				log.p=TRUE)
		
		results <- rbind(results, result)
		
		if (verbose)
			setTxtProgressBar(pBar, k/length(x))
	}
	
	o <- order(results[, "Index"], results[, "Begin"], results[, "End"])
	results <- results[o,, drop=FALSE]
	bounds <- matrix(c(as.integer(results[, "Index"]),
			as.integer(results[, "Strand"]),
			as.integer(results[, "Begin"]),
			as.integer(results[, "End"])),
		ncol=4)
	
	if (!allScores) {
		indices <- .Call("chainGenes",
			bounds,
			results[, "TotalScore"],
			bounds[, 4] - bounds[, 3] + 1L,
			FALSE, # scoreIntergenic
			20L, # maxOverlapSame
			30L, # maxOverlapOpposite
			0.05, # maxFracOverlap
			numeric(), # sameScores
			numeric(), # oppoScores
			PACKAGE="DECIPHER")
		
		results <- results[indices,, drop=FALSE]
	}
	
	class(results) <- "Genes"
	attr(results, "widths") <- setNames(width(myXStringSet),
		names(myXStringSet))
	if (!is.null(annotations)) {
		lkup <- sort(unique(results[, "Gene"]),
			decreasing=TRUE)
		attr(results, "annotations") <- setNames(lkup,
			annotations[-lkup])
	}
	
	if (verbose) {
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
	}
	
	return(results)
}

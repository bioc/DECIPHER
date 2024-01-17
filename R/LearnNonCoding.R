LearnNonCoding <- function(myXStringSet,
	threshold=0.3,
	weight=NA,
	maxLoopLength=500,
	maxPatterns=20,
	scoreDependence=FALSE,
	structure=NULL,
	processors=1) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	l <- length(myXStringSet)
	if (l == 0)
		stop("myXStringSet must contain sequences.")
	uw <- unique(width(myXStringSet))
	if (length(uw) != 1)
		stop("Sequences in myXStringSet must be the same width (aligned).")
	if (uw <= 7)
		stop("Sequences in myXStringSet contain too few nucleotides.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (length(threshold) != 1L)
		stop("threshold must be a single numeric.")
	if (threshold <= 0)
		stop("threshold must be greater than zero.")
	if (threshold > 1)
		stop("threshold can be at most one.")
	if (!is.numeric(maxLoopLength))
		stop("maxLoopLength must be a numeric.")
	if (length(maxLoopLength) != 1L)
		stop("maxLoopLength must be a single numeric.")
	if (maxLoopLength <= 11)
		stop("maxLoopLength must be at least 12.")
	if (!is.numeric(maxPatterns))
		stop("maxPatterns must be a numeric.")
	if (length(maxPatterns) == 1L) {
		maxPatterns <- rep(maxPatterns, 2L)
	} else if (length(maxPatterns) != 2L) {
		stop("maxPatterns must be one or two numerics.")
	}
	if (!is.logical(scoreDependence))
		stop("scoreDependence must be a logical.")
	if (!is.null(structure)) {
		if (is.character(structure)) {
			if (length(structure) > 1)
				stop("structure must be a character vector of length one.")
			if (nchar(structure) != uw)
				stop("structure must be a character string with the same number of characters as the width of myXStringSet.")
			structure <- .parseDBN(structure)
		} else if (is.numeric(structure)) {
			if (!is.matrix(structure))
				stop("structure must be a matrix.")
			if (ncol(structure) < 2)
				stop("structure must be a matrix with at least two columns.")
			structure <- structure[, 1:2]
			if (any(floor(structure) != structure))
				stop("structure must be a matrix of whole numbers.")
			if (any(structure < 1))
				stop("All values in structure must be at least 1.")
		} else {
			stop("structure must be a character string or matrix.")
		}
	}
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
	
	weight <- .getWeights(weight,
		myXStringSet,
		FALSE, # verbose
		processors)
	
	myXStringSet <- DNAStringSet(myXStringSet)
	noGaps <- RemoveGaps(myXStringSet)
	wS <- width(noGaps)
	avgWidth <- mean(wS)
	if (avgWidth > 5e4)
		stop("Sequences in myXStringSet are too long.")
	# determine expected k-mer frequencies
	for (K in 1:4) { # optimize K
		o <- oligonucleotideFrequency(noGaps,
			K,
			as.prob=FALSE,
			fast.moving.side="left",
			with.labels=FALSE)
		o <- colSums(o*weight)
		if (any(o < 10L)) { # inaccurate frequencies
			K <- K - 1L
			break
		}
		o <- as.integer(o)
		oligos <- o
	}
	if (K == 0L) {
		if (any(o == 0)) {
			stop("Not enough diversity among sequences in myXStringSet.")
		} else {
			oligos <- o <- as.integer(o)
		}
	}
	if (avgWidth < K)
		stop("Sequences in myXStringSet are too short.")
	
	# initialize parameters
	maxFractionGaps <- 0.5 # maximum indels in considered positions
	maxEntropy <- 1.8 # maximum entropy of kept positions (in bits)
	quants_bounds <- 0.1 # quantile for boundaries
	quants_bounds <- c(quants_bounds, 1 - quants_bounds)
	maxIndels <- 0.05 # maximum cumulative indels in motifs
	length_params <- c(-2, mean(wS), 1) # initial parameters for fitting a distribution of sequence lengths
	N <- 1e6 # number of nucleotides in random sequences
	AT <- seq(0.25, 0.75, 0.05) # AT-content of random sequences
	maxNucs <- ceiling(log(1000*avgWidth, 4)) # maximum number of equivalent nucleotides in motifs
	windowSize <- 1 # points to left and right of center for moving average
	nBins <- 10 # maximum number of bins for scoring
	minBinSize <- 1 # minimum bin size for scoring free energy
	minBinCount <- 10 # minimum count per bin to score dependency among patterns
	alpha <- 0.01 # p-value threshold for scoring dependency among patterns
	maxDiscerningPower <- 100 # adequate total discerning power from motifs
	
	ions <- 1
	temp <- 37
	data("deltaHrulesRNA", envir=environment(), package="DECIPHER")
	data("deltaSrulesRNA", envir=environment(), package="DECIPHER")
	deltaSrulesRNA <- deltaSrulesRNA + 0.368*log(ions)/1000
	deltaGrulesRNA <- deltaHrulesRNA - (273.15 + temp)*deltaSrulesRNA
	dG_ini <- 4.065225 # 3.6 - (273.15 + temp)*(-1.5/1000 + 0.368*log(ions)/1000)
	max_dG <- 0 # maximum free energy of hairpins
	
	pos <- as.matrix(myXStringSet)
	pos <- lapply(seq_len(nrow(pos)),
		function(x)
			which(!(pos[x,] %in% c("-", "."))))
	
	a <- oligonucleotideFrequency(noGaps,
		1L,
		as.prob=TRUE)
	a <- colSums(weight*a)
	a <- a[DNA_BASES]
	a <- a/sum(a)
	random <- sapply(seq_len(N/avgWidth),
		function(x)
			paste(sample(DNA_BASES,
					avgWidth,
					replace=TRUE,
					prob=a),
				collapse=""))
	random <- DNAStringSet(random)
	ce <- cumsum(width(random))
	cs <- c(1L, ce[-length(ce)] + 1L)
	random <- unlist(random)
	
	# find conserved sequence patterns
	profile <- .Call("consensusProfile",
		myXStringSet,
		weight,
		NULL,
		PACKAGE="DECIPHER")
	ungapped <- profile[8,]*(profile[5,] - 1) + 1
	ungapped <- which(ungapped <= maxFractionGaps)
	
	ACGT <- profile[1:4, ungapped, drop=FALSE]
	rownames(ACGT) <- DNA_BASES
	ACGT <- t(t(ACGT)/colSums(ACGT))
	ACGT[ACGT == 0] <- min(weight)/l
	ACGT <- t(t(ACGT)/colSums(ACGT))
	entropy <- colSums(-ACGT*log2(ACGT), na.rm=TRUE)
	if (windowSize > 0) # apply moving average
		entropy <- .centerPoint(entropy, windowSize)
	keep <- entropy <= maxEntropy
	
	blocks <- rle(keep)
	w <- which(blocks$values)
	
	if (length(w) > 0) {
		end <- cumsum(blocks$lengths)
		start <- end - blocks$lengths + 1L
		end <- end[w]
		start <- start[w]
		
		# split motifs at indels
		starts <- ends <- vector("list", length(start))
		for (i in seq_along(start)) {
			starts[[i]] <- start[i]
			if (start[i] + 1L < end[i]) {
				for (j in (start[i] + 1L):(end[i] - 1L)) {
					s <- tail(starts[[i]], 1)
					indels <- subseq(myXStringSet,
						ungapped[s],
						ungapped[j])
					indels <- RemoveGaps(indels)
					fracIndels <- sum(weight[width(indels) != j - s + 1L])/l
					if (fracIndels > maxIndels) {
						starts[[i]] <- c(starts[[i]], j)
						ends[[i]] <- c(ends[[i]], j - 1L)
					}
				}
			}
			ends[[i]] <- c(ends[[i]], end[i])
		}
		start <- unlist(starts)
		end <- unlist(ends)
		
		# split motifs that are too long
		starts <- ends <- vector("list", length(start))
		for (i in seq_along(start)) {
			s <- start[i]
			e <- end[i]
			len <- 1 - entropy[s:e]/2 # equivalent length
			bins <- ceiling(sum(len)/maxNucs)
			starts[[i]] <- s
			if (bins > 1) {
				len <- cumsum(len)
				j <- seq(0,
					len[length(len)],
					length.out=bins + 1)
				j <- j[-c(1, length(j))]
				j <- sapply(j,
					function(x)
						which.min(abs(len - x)))
				j <- j + s - 1L
				starts[[i]] <- c(starts[[i]], j)
				ends[[i]] <- c(ends[[i]], j - 1L)
			}
			ends[[i]] <- c(ends[[i]], end[i])
		}
		start <- unlist(starts)
		end <- unlist(ends)
		
		consensus <- lapply(seq_along(start),
			function(i)
				ACGT[, start[i]:end[i], drop=FALSE])
	} else {
		stop("No conserved motifs found in myXStringSet.")
	}
	
	.vPWM <- function(pwm,
		subject,
		min.score) {
		
		ws <- width(subject)
		cs <- cumsum(ws)
		subject <- unlist(subject)
		
		suppressWarnings({
			v <- matchPWM(pwm,
				subject,
				min.score,
				with.score=TRUE)
			})
		
		s <- start(v)
		e <- end(v)
		
		index <- integer(length(s))
		j <- 1L
		for (i in seq_along(index)) {
			while (s[i] > cs[j])
				j <- j + 1L
			if (e[i] <= cs[j])
				index[i] <- j
		}
		
		score <- mcols(v)$score
		ws <- cs - ws
		starts <- scores <- vector("list", length(ws))
		j <- 1L
		for (i in seq_along(starts)) {
			w <- .Call("multiMatch",
				index,
				i,
				j,
				PACKAGE="DECIPHER")
			starts[[i]] <- s[w] - ws[i]
			scores[[i]] <- score[w]
		}
		
		return(list(starts, scores))
	}
	
	motifs <- integer(length(consensus))
	motifs <- data.frame(begin_low=motifs,
		begin_high=motifs,
		end_low=motifs,
		end_high=motifs,
		motif=character(length(consensus)),
		pwm=I(consensus),
		minscore=I(vector("list", length(consensus))),
		prevalence=I(vector("list", length(consensus))),
		background=I(vector("list", length(consensus))),
		hits=I(vector("list", length(consensus))))
	for (i in seq_along(consensus)) {
		motif <- log(motifs[i, "pwm"][[1]]/0.25)
		
		v <- .vPWM(motif,
			noGaps,
			min.score=0)
		v[[3L]] <- mapply(function(a, b) b[a],
			v[[1L]],
			pos,
			SIMPLIFY=FALSE)
		# v[[1]]: positions of hits in noGaps
		# v[[2]]: scores of hits in noGaps
		# v[[3]]: positions of hits in the alignment
		m <- sapply(v[[3]], match, x=ungapped[start[i]]) # index of hits to expected position
		hits <- mapply(`[`, v[[2]], m) # score of hits to expected position
		begin <- mapply(`[`, v[[1]], m) # position in noGaps of hits to expected position
		w <- which(!is.na(begin)) # indices of hits in noGaps
		
		len <- ncol(motif)
		motifs[i, "begin_low"] <- floor(quantile(begin[w] - 1L, quants_bounds[1]))
		motifs[i, "begin_high"] <- ceiling(quantile(begin[w] - 1L, quants_bounds[2]))
		motifs[i, "end_low"] <- floor(quantile(wS[w] - begin[w] - len + 1, quants_bounds[1]))
		motifs[i, "end_high"] <- ceiling(quantile(wS[w] - begin[w] - len + 1, quants_bounds[2]))
		delta_starts <- motifs[i, "begin_high"] - motifs[i, "begin_low"]
		delta_ends <- motifs[i, "end_high"] - motifs[i, "end_low"]
		
		# separate scores into appoximately even bins
		o <- order(hits, na.last=FALSE)
		bins <- hits[o[seq(1, length(o), length.out=nBins)]]
		bins <- bins[!is.na(bins)]
		bins <- unique(bins)
		bins[length(bins)] <- Inf
		bins <- c(0, bins)
		
		bg <- matchPWM(motif,
			random,
			min.score=0,
			with.score=TRUE)
		bg <- mcols(bg)$score
		bg <- .bincode(bg, bins)
		bg <- tabulate(bg, length(bins) - 1)
		bg[bg == 0L] <- 1L # pseudocount
		bg <- bg/length(random)
		delta <- min(delta_starts, delta_ends) + 1
		bg <- 1 - (1 - bg)^delta
		
		h <- .bincode(hits, bins)
		h <- tabulate(h)/l
		w <- which(h > bg)
		if (length(w) == 0L)
			next
		bins <- bins[w[1]:length(bins)]
		motifs[i, "minscore"][[1]] <- list(bins)
		
		bg <- bg[w[1]:length(bg)]
		bg <- c(1 - sum(bg), bg)
		names(bg) <- seq_along(bg) - 1L
		motifs[i, "background"][[1L]] <- list(bg)
		
		if (scoreDependence) {
			hits <- .bincode(hits, bins)
			hits[is.na(hits)] <- 0L
			motifs[i, "hits"][[1L]] <- list(hits)
		} else {
			vec <- NULL
		}
		
		found <- integer(length(v[[1]]))
		for (j in seq_along(v[[1]])) {
			s <- v[[1L]][[j]] - 1L
			if (delta_starts < delta_ends) {
				w <- which(s >= motifs[i, "begin_low"] &
					s <= motifs[i, "begin_high"])
			} else {
				e <- wS[j] - s - len
				w <- which(e >= motifs[i, "end_low"] &
					e <= motifs[i, "end_high"])
			}
			if (length(w) > 0) {
				w <- .bincode(v[[2L]][[j]][w], bins)
				if (any(!is.na(w)))
					found[j] <- max(w, na.rm=TRUE)
			}
		}
		t <- tapply(weight, found, sum)/l
		hits <- rep(1/l, length(bins))
		names(hits) <- seq_along(hits) - 1L
		hits[names(t)] <- t
		hits <- hits/sum(hits)
		if (hits[1L] >= bg[1L])
			next
		motifs[i, "prevalence"][[1L]] <- list(hits)
		
		motifs[i, "motif"] <- paste(apply(motif,
				2L,
				function(x) {
					w <- which(x >= 0)
					y <- DNA_BASES[w]
					y <- paste(y, collapse="")
					y <- match(y, IUPAC_CODE_MAP)
					y <- names(IUPAC_CODE_MAP)[y]
					if (all(x[w] < 1))
						y <- tolower(y)
					y
				}),
			collapse="")
	}
	motifs <- motifs[lengths(motifs[["prevalence"]]) > 0,]
	motifs <- unique(motifs)
	
	# find conserved structure patterns
	hairpins <- vector("list", l)
	patterns <- matrix(0L, uw, uw)
	for (i in seq_along(hairpins)) {
		pals <- findPalindromes(noGaps[[i]],
			min.armlength=4,
			max.looplength=maxLoopLength,
			min.looplength=3,
			max.mismatch=1,
			allow.wobble=TRUE)
		
		dna <- DNAStringSet(pals)
		arms <- palindromeArmLength(dna,
			max.mismatch=1,
			allow.wobble=TRUE)
		max_arms <- as.integer((width(pals) - 3)/2)
		w <- which(arms > max_arms)
		if (length(w) > 0)
			arms[w] <- max_arms[w]
		
		dG <- .Call("calculateHairpinDeltaG",
			dna,
			arms + 1L, # include end bases
			deltaGrulesRNA,
			PACKAGE="DECIPHER") + dG_ini
		o <- order(dG)
		o <- o[dG[o] <= max_dG]
		
		if (length(o) > 0) {
			s1 <- start(pals)
			s2 <- s1 + arms - 1L
			e2 <- end(pals)
			e1 <- e2 - arms + 1L
			
			hairpins[[i]] <- data.frame(index=i,
				begin=start(pals)[o] - 1L,
				end=wS[i] - end(pals)[o],
				width=width(pals)[o],
				length=arms[o],
				start1=pos[[i]][s1][o],
				start2=pos[[i]][s2][o],
				finish1=pos[[i]][e1][o],
				finish2=pos[[i]][e2][o],
				dG=dG[o])
		} else {
			hairpins[[i]] <- data.frame(index=integer(),
				begin=integer(),
				end=integer(),
				width=integer(),
				length=integer(),
				start1=integer(),
				start2=integer(),
				finish1=integer(),
				finish2=integer(),
				dG=numeric())
		}
	}
	
	if (is.null(structure)) {
		p <- PredictDBN(myXStringSet,
			type="pairs",
			minOccupancy=1 - maxFractionGaps,
			weight=weight,
			verbose=FALSE,
			processors=processors)
	} else {
		p <- structure
	}
	p <- p[, 1:2, drop=FALSE]
	p[] <- match(p, ungapped)
	p <- p[rowSums(is.na(p)) == 0L,, drop=FALSE]
	
	if (nrow(p) > 0) {
		# keep hairpins spanning arms in the alignment
		hairpins <- do.call(rbind, hairpins)
		keep <- mapply(function(a, b, c, d) {
				w <- which(hairpins[, "start1"] <= a)
				w <- w[hairpins[w, "start2"] >= b]
				w <- w[hairpins[w, "finish1"] <= c]
				w <- w[hairpins[w, "finish2"] >= d]
				w[!duplicated(hairpins[w, "index"])]
			},
			ungapped[p[, "("]],
			ungapped[p[, "("]],
			ungapped[p[, ")"]],
			ungapped[p[, ")"]],
			SIMPLIFY=FALSE)
		
		prevalence <- sapply(keep,
			function(x)
				sum(weight[hairpins[x, "index"]])/l)
		keep <- keep[prevalence >= threshold]
	} else {
		keep <- list()
	}
	
	if (length(keep) > 0) {
		pals <- findPalindromes(random,
			min.armlength=4,
			max.looplength=maxLoopLength,
			min.looplength=3,
			max.mismatch=1,
			allow.wobble=TRUE)
		
		dna <- DNAStringSet(pals)
		arms <- palindromeArmLength(dna,
			max.mismatch=1,
			allow.wobble=TRUE)
		max_arms <- as.integer((width(pals) - 3)/2)
		w <- which(arms > max_arms)
		if (length(w) > 0)
			arms[w] <- max_arms[w]
		
		p_dG <- .Call("calculateHairpinDeltaG",
			dna,
			arms + 1L, # include end bases
			deltaGrulesRNA,
			PACKAGE="DECIPHER") + dG_ini
		
		w <- which(p_dG <= max_dG)
		pals <- pals[w]
		p_dG <- p_dG[w]
		
		p_starts <- start(pals)
		p_ends <- end(pals)
		o <- order(p_starts)
		p_starts <- p_starts[o]
		p_ends <- p_ends[o]
		p_dG <- p_dG[o]
		p_widths <- p_ends - p_starts + 1L
		o <- order(p_ends)
		p_ends2 <- p_ends[o]
		p_widths2 <- p_widths[o]
		p_dG2 <- p_dG[o]
		
		hairpins <- sapply(seq_along(keep),
			function(x) {
				h <- hairpins[keep[[x]],, drop=FALSE]
				
				begin_low <- floor(quantile(h[, "begin"], quants_bounds[1]))
				begin_high <- ceiling(quantile(h[, "begin"], quants_bounds[2]))
				end_low <- floor(quantile(h[, "end"], quants_bounds[1]))
				end_high <- ceiling(quantile(h[, "end"], quants_bounds[2]))
				width_low <- floor(quantile(h[, "width"], quants_bounds[1]))
				width_high <- ceiling(quantile(h[, "width"], quants_bounds[2]))
				length_low <- floor(quantile(h[, "length"], quants_bounds[1]))
				length_high <- ceiling(quantile(h[, "length"], quants_bounds[2]))
				
				# correct for hairpins running off ends
				if (begin_low == 0L ||
					end_low == 0L) {
					# extend four nucleotides on each end
					begin_low <- begin_low - 4L
					end_low <- end_low - 4L
					width_high <- width_high + 8L
					length_high <- length_high + 4L
				}
				
				delta_start <- begin_high - begin_low
				delta_end <- end_high - end_low
				delta_width <- width_high - width_low
				if (delta_start < delta_end) {
					w <- which(hairpins[, "begin"] >= begin_low)
					w <- w[hairpins[w, "begin"] <= begin_high]
					if (delta_end < delta_width) {
						w <- w[hairpins[w, "end"] >= end_low]
						w <- w[hairpins[w, "end"] <= end_high]
					} else {
						w <- w[hairpins[w, "width"] >= width_low]
						w <- w[hairpins[w, "width"] <= width_high]
					}
				} else {
					w <- which(hairpins[, "end"] >= end_low)
					w <- w[hairpins[w, "end"] <= end_high]
					if (delta_start < delta_width) {
						w <- w[hairpins[w, "begin"] >= begin_low]
						w <- w[hairpins[w, "begin"] <= begin_high]
					} else {
						w <- w[hairpins[w, "width"] >= width_low]
						w <- w[hairpins[w, "width"] <= width_high]
					}
				}
				w <- w[!duplicated(hairpins[w, "index"])]
				
				dG <- hairpins[w, "dG"]
				dG <- dG[order(dG)]
				dG <- c(dG, rep(0, l - length(dG)))
				bins <- dG[seq(1, length(dG), length.out=nBins)]
				bins <- c(bins, 0)
				bins <- unique(bins)
				bins[1L] <- -Inf
				repeat {
					eliminate <- which.max(diff(bins) < minBinSize)
					eliminate <- eliminate + 1L
					if (eliminate == 2L)
						break
					if (eliminate == length(bins)) {
						bins <- bins[-eliminate + 1L]
						break
					} else {
						bins <- bins[-eliminate]
					}
				}
				
				if (scoreDependence) {
					vec <- numeric(l)
					vec[hairpins[w, "index"]] <- hairpins[w, "dG"]
					vec <- .bincode(vec, bins)
				}
				
				t <- .bincode(hairpins[w, "dG"], bins)
				t <- tapply(weight[hairpins[w, "index"]],
					t,
					sum)/l
				hits <- numeric(length(bins) - 1L)
				names(hits) <- seq_along(hits)
				hits[names(t)] <- t
				hits[length(hits)] <- hits[length(hits)] + 1 - sum(hits)
				if (isTRUE(all.equal(hits[length(hits)], 0, check.attributes=FALSE))) {
					hits[length(hits)] <- 1/l # pseudocount
					hits <- hits/sum(hits)
				}
				
				# estimate background frequency of hairpins
				if (delta_start < delta_end) {
					if (delta_end < delta_width) {
						index <- .Call("getHits",
							p_starts,
							p_ends,
							cs + as.integer(begin_low),
							cs + as.integer(begin_high),
							ce - as.integer(end_high),
							ce - as.integer(end_low),
							p_dG,
							PACKAGE="DECIPHER")
						bg <- p_dG[index]
					} else {
						index <- .Call("getHits",
							p_starts,
							p_widths,
							cs + as.integer(begin_low),
							cs + as.integer(begin_high),
							rep(as.integer(width_low),
								length(cs)),
							rep(as.integer(width_high),
								length(cs)),
							p_dG,
							PACKAGE="DECIPHER")
						bg <- p_dG[index]
					}
				} else {
					if (delta_start < delta_width) {
						index <- .Call("getHits",
							p_starts,
							p_ends,
							cs + as.integer(begin_low),
							cs + as.integer(begin_high),
							ce - as.integer(end_high),
							ce - as.integer(end_low),
							p_dG,
							PACKAGE="DECIPHER")
						bg <- p_dG[index]
					} else {
						index <- .Call("getHits",
							p_ends2,
							p_widths2,
							ce - as.integer(end_high),
							ce - as.integer(end_low),
							rep(as.integer(width_low),
								length(ce)),
							rep(as.integer(width_high),
								length(ce)),
							p_dG2,
							PACKAGE="DECIPHER")
						bg <- p_dG2[index]
					}
				}
				bg <- c(bg, rep(0, length(index) - length(bg)))
				bg <- .bincode(bg, bins)
				bg <- tabulate(bg, length(bins) - 1)
				bg[bg == 0L] <- 1L # pseudocount
				bg <- bg/length(ce)
				
				if (bg[1L] > hits[1L]) {
					data.frame(begin_low=numeric(),
						begin_high=numeric(),
						end_low=numeric(),
						end_high=numeric(),
						width_low=numeric(),
						width_high=numeric(),
						length_low=numeric(),
						length_high=numeric(),
						dG=I(list()),
						prevalence=I(list()),
						background=I(list()),
						hits=I(list()))
				} else {
					data.frame(begin_low=begin_low,
						begin_high=begin_high,
						end_low=end_low,
						end_high=end_high,
						width_low=width_low,
						width_high=width_high,
						length_low=length_low,
						length_high=length_high,
						dG=I(list(bins)),
						prevalence=I(list(hits)),
						background=I(list(bg)),
						hits=I(list(vec)))
				}
			},
			simplify=FALSE)
		hairpins <- do.call(rbind, hairpins)
	} else {
		hairpins <- data.frame(begin_low=numeric(),
			begin_high=numeric(),
			end_low=numeric(),
			end_high=numeric(),
			width_low=numeric(),
			width_high=numeric(),
			length_low=numeric(),
			length_high=numeric(),
			dG=I(list()),
			prevalence=I(list()),
			background=I(list()))
	}
	hairpins <- unique(hairpins)
	
	# filter patterns by discerning power
	f <- function(a, b) {
		# add pseudocounts
		a[a <= 0] <- 0.0001
		b[b <= 0] <- 0.0001
		sum(abs(log(a/b))*a)
	}
	score_motifs <- mapply(f,
		motifs[["prevalence"]],
		motifs[["background"]])
	score_hairpins <- mapply(f,
		hairpins[["prevalence"]],
		hairpins[["background"]])
	prob <- unlist(c(score_motifs, score_hairpins))
	keep <- prob > 0
	
	# remove overlapping motifs
	n1 <- nrow(motifs)
	n2 <- nrow(hairpins)
	w <- which(keep)
	o <- w[order(prob[w], decreasing=TRUE)]
	used_s1 <- used_e1 <- used_s2 <- used_e2 <- logical(uw)
	for (i in seq_along(o)) {
		if (o[i] > n1) {
			s <- hairpins[o[i] - n1, "begin_low"] + 1L
			if (s < 1L)
				s <- 1L
			s <- s:(s + hairpins[o[i] - n1, "length_low"] - 1L)
			e <- hairpins[o[i] - n1, "end_low"] + 1L
			if (e < 1L)
				e <- 1L
			e <- e:(e + hairpins[o[i] - n1, "length_low"] - 1L)
			
			if (any(used_s2[s]) || any(used_e2[e])) {
				keep[o[i]] <- FALSE
			} else {
				used_s2[s] <- TRUE
				used_e2[e] <- TRUE
			}
		} else {
			delta_starts <- motifs[o[i], "begin_high"] - motifs[o[i], "begin_low"]
			delta_ends <- motifs[o[i], "end_high"] - motifs[o[i], "end_low"]
			if (delta_starts < delta_ends) {
				s <- motifs[o[i], "begin_low"] + 1L
				s <- s:(s + ncol(motifs[o[i], "pwm"][[1]]) - 1L)
				if (all(used_s1[s])) {
					keep[o[i]] <- FALSE
				} else {
					used_s1[s] <- TRUE
				}
			} else {
				e <- motifs[o[i], "end_low"] + 1L
				e <- e:(e + ncol(motifs[o[i], "pwm"][[1]]) - 1L)
				if (all(used_e1[e])) {
					keep[o[i]] <- FALSE
				} else {
					used_e1[e] <- TRUE
				}
			}
		}
	}
	
	motifs <- motifs[keep[seq_len(n1)],]
	if (nrow(motifs) > maxPatterns[1L]) {
		w <- which(keep[seq_len(n1)])
		o <- order(prob[w],
			decreasing=TRUE)
		tot <- which(cumsum(prob[w[o]]) < maxDiscerningPower)
		tot <- tot[length(tot)]
		if (is.na(tot))
			tot <- maxPatterns[1L]
		motifs <- motifs[o,]
		motifs <- motifs[seq_len(min(maxPatterns[1L], tot)),]
	} else if (nrow(motifs) == 0L) {
		stop("No conserved motifs found in myXStringSet.")
	}
	motifs <- motifs[order(motifs[, "begin_low"], -motifs[, "end_low"]),]
	hairpins <- hairpins[keep[n1 + seq_len(n2)],]
	if (nrow(hairpins) > maxPatterns[2L]) {
		w <- which(keep[n1 + seq_len(n2)])
		o <- order(prob[w],
			decreasing=TRUE)
		hairpins <- hairpins[o,]
		hairpins <- hairpins[seq_len(maxPatterns[2L]),]
	}
	hairpins <- hairpins[order(hairpins[, "begin_low"]),]
	
	rownames(motifs) <- NULL
	rownames(hairpins) <- NULL
	
	if (scoreDependence) {
		# calculate dependence among motifs and hairpins
		n1 <- nrow(motifs)
		n2 <- nrow(hairpins)
		scores <- vector("list", (n1 + n2)*(n1 + n2 - 1)/2)
		alpha <- alpha/length(scores) # multiple testing correction
		k <- 0L
		for (i in seq_len(n1 + n2 - 1L)) {
			for (j in (i + 1L):(n1 + n2)) {
				k <- k + 1L
				if (i > n1) {
					t1 <- hairpins[[i - n1, "hits"]]
					p1 <- paste("hairpin", i - n1, sep="")
				} else {
					t1 <- motifs[[i, "hits"]] + 1L # offset because missing hits recorded as zero
					p1 <- paste("motif", i, sep="")
				}
				if (j > n1) {
					t2 <- hairpins[[j - n1, "hits"]]
					p2 <- paste("hairpin", j - n1, sep="")
				} else {
					t2 <- motifs[[j, "hits"]] + 1L # offset because missing hits recorded as zero
					p2 <- paste("motif", j, sep="")
				}
				names(scores)[k] <- paste(p1, p2, sep=" / ")
				
				m1 <- max(t1)
				m2 <- max(t2)
				observed <- matrix(0, m1, m2)
				for (i1 in seq_len(m1))
					for (i2 in seq_len(m2))
						observed[i1, i2] <- sum(weight[t1 == i1 & t2 == i2])
				
				expected <- outer(rowSums(observed), colSums(observed))/sum(observed)
				suppressWarnings(p <- chisq.test(observed, p=expected, rescale.p=TRUE)$p.value)
				if (is.na(p) || p >= alpha)
					next
				
				eliminate <- expected < minBinCount & observed < minBinCount
				eliminate <- eliminate | observed == 0
				eliminate <- eliminate |
					(observed < qbinom(0.975, l, expected/l) &
					observed > qbinom(0.025, l, expected/l))
				if (!all(eliminate)) {
					scores[[k]] <- log(observed/expected)
					scores[[k]][eliminate] <- 0
				}
			}
		}
	} else {
		scores <- NULL
	}
	motifs <- motifs[, -ncol(motifs)] # drop "hits"
	hairpins <- hairpins[, -ncol(hairpins)] # drop "hits"
	
	# compute log-odds length scores
	.PDF <- function(x, p)
		-p[1]*p[3]*exp(p[1]*log(x/p[2]))*(exp(p[1]*log(x/p[2])) + 1)^(-1 - p[3])/x
	.CDF <- function(x, p)
		sig <- (1 + exp(p[1]*log(x/p[2])))^-p[3]
	.fitSigmoid <- function(x, y, ini) {
		.SSE <- function(p) {
			# avoid conditions with numeric errors
			pdf <- .PDF(x, p)
			if (any(is.nan(pdf) |
				is.infinite(pdf) |
				pdf < 0))
				return(sum(y))
			
			expected <- .CDF(x, p)
			
			sum((y - expected)^2)
		}
		o <- suppressWarnings(nlminb(ini,
			.SSE))
		o$par
	}
	
	observed <- tabulate(wS, 2*max(wS))
	observed <- observed/sum(observed)
	observed <- cumsum(observed)
	include <- logical(length(observed))
	last <- 1
	for (i in which(observed > 0)) {
		if (observed[last] + 0.001 <= observed[i]) {
			last <- i
			include[i] <- TRUE
		}
	}
	include[1L] <- TRUE
	include[tail(which(observed == 0), 1)] <- TRUE
	include[length(include)] <- TRUE
	include <- which(include)
	observed <- observed[include]
	
	o <- .fitSigmoid(include,
		observed,
		length_params)
	
	minLength <- as.integer(min(wS)/2)
	maxLength <- as.integer(max(wS)*2)
	lenScores <- .PDF(minLength:maxLength, o) # fitted distribution of lengths
	w <- which(lenScores > 1e-100)
	minLength <- minLength + w[1] - 1L
	maxLength <- maxLength - length(lenScores) + w[length(w)]
	lenScores <- lenScores[w]
	lenScores <- lenScores/sum(lenScores)
	lenScores <- lenScores*(maxLength - minLength + 1) # assume a uniform distribution of background lengths
	lenScores <- log(lenScores)
	
	result <- list(motifs=motifs,
		hairpins=hairpins,
		kmers=oligos,
		lengthScores=lenScores,
		dependence=scores)
	class(result) <- "NonCoding"
	attr(result, "K") <- K
	attr(result, "minLength") <- minLength
	attr(result, "maxLength") <- maxLength
	attr(result, "maxLoopLength") <- maxLoopLength
	attr(result, "background") <- c(meanlog=NA_real_,
		sdlog=NA_real_)
	
	# determine inclusion threshold for reporting
	# based on the assumption that the scores of
	# false discoveries can be approximated by
	# the right tail of a normal distribution
	random <- paste(sample(DNA_BASES,
			N,
			replace=TRUE,
			prob=a),
		collapse="")
	bg <- FindNonCoding(result,
		DNAStringSet(random),
		minScore=0,
		processors=processors,
		verbose=FALSE)[, "TotalScore"]
	AT <- AT/2 # relative frequency of A or T
	for (i in seq_along(AT)) {
		p <- c(AT[i],
			rep(0.5 - AT[i], 2),
			AT[i])
		random <- paste(sample(DNA_BASES,
				N,
				replace=TRUE,
				prob=p),
			collapse="")
		temp <- FindNonCoding(result,
			DNAStringSet(random),
			minScore=0,
			processors=processors,
			verbose=FALSE)[, "TotalScore"]
		if (length(bg) == 0L ||
			(length(temp) > 0L &&
			max(temp) > max(bg))) {
			a <- p
			bg <- temp
		}
	}
	if (length(bg) > 0 &&
		max(bg) > log(N)/4) {
		# calibrate log-odds scores
		count <- 1L
		while (length(bg) < 1e4 &&
			count < 20) {
			random <- paste(sample(DNA_BASES,
					N,
					replace=TRUE,
					prob=a),
				collapse="")
			bg <- c(bg,
				FindNonCoding(result,
					DNAStringSet(random),
					minScore=0,
					processors=processors,
					verbose=FALSE)[, "TotalScore"])
			count <- count + 1L
		}
		
		# select tail of background distribution
		cutoff <- quantile(bg, 0.9)
		bg <- bg[bg > cutoff]
		
		# fit censored log-normal distribution
		missing <- count*N*2 - length(bg) # 2 for both strands
		.MLE <- function(params) {
			observed <- sum(-log(dlnorm(bg, params[1], exp(params[2]))))
			censored <- missing*plnorm(cutoff, params[1], exp(params[2]), log.p=TRUE)
			observed - censored
		}
		o <- optim(c(-5, 0.25),
			.MLE)
		mean <- o$par[[1]]
		sd <- exp(o$par[[2]])
		
		if (qlnorm(1.125352e-07, mean, sd, lower.tail=FALSE) > 16)
			attr(result, "background") <- c(meanlog=mean,
				sdlog=sd)
	}
	
	return(result)
}

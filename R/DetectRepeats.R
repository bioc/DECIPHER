DetectRepeats <- function(myXStringSet,
	type="tandem",
	minScore=15,
	allScores=FALSE,
	maxPeriod=1000,
	maxFailures=2,
	maxShifts=5,
	alphabet=AA_REDUCED[[125]],
	useEmpirical=TRUE,
	processors=1,
	verbose=TRUE,
	...) {
	
	# error checking
	if (length(type) == 0)
		stop("No type specified.")
	if (length(type) > 1)
		stop("Only one type may be specified.")
	TYPES <- c("tandem", "interspersed", "both")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (length(myXStringSet) == 0)
		stop("myXStringSet must contain at least one sequence.")
	if (is(myXStringSet, "DNAStringSet")) {
		xtype <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		xtype <- 2L
	} else if (is(myXStringSet, "AAStringSet")) {
		if (type > 1L)
			stop("type must be 'tandem' when myXStringSet is an AAStringSet.")
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
		sizeAA <- max(alphabet)
		if (sizeAA == 1L)
			stop("More than one grouping of amino acids is required in the alphabet.")
		n <- as.integer(floor(log(4294967295, sizeAA)))
		alphabet <- alphabet - 1L
	} else {
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	if (!is.numeric(minScore) || length(minScore) > 1)
		stop("minScore must be a numeric.")
	if (minScore < 0)
		stop("minScore must be at least zero.")
	if (!is.logical(allScores))
		stop("allScores must be a logical.")
	if (!is.numeric(maxPeriod) || length(maxPeriod) > 1)
		stop("maxPeriod must be a numeric.")
	if (maxPeriod < 1)
		stop("maxPeriod must be at least one.")
	maxPeriod <- min(maxPeriod, max(2000, width(myXStringSet))/2)
	if (!is.numeric(maxFailures))
		stop("maxFailures must be a numeric.")
	if (maxFailures != floor(maxFailures))
		stop("maxFailures must be a whole number.")
	if (maxFailures < 0)
		stop("maxFailures must be at least zero.")
	if (!is.numeric(maxShifts))
		stop("maxShifts must be a numeric.")
	if (maxShifts != floor(maxShifts))
		stop("maxShifts must be a whole number.")
	if (maxShifts < 0)
		stop("maxShifts must be at least zero.")
	if (maxPeriod < 1)
		stop("maxPeriod must be at least one.")
	if (!is.logical(useEmpirical))
		stop("useEmpirical must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') are not allowed in myXStringSet.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') are not allowed in myXStringSet.")
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
	
	if (verbose)
		time.1 <- Sys.time()
	
	# default parameters
	N <- 10 # find k-mers on average N times per maxPeriod
	gapCost <- c(-2.6, -0.1) # gap open and extension coefficients
	mult <- 2L # multiplier on K for lookahead length
	maxVisits <- 1L # maximum attempts to detect repeats in each position
	if (xtype == 3L) {
		residues <- c(A=0.05, R=-0.02, N=0.05, D=0.06, C=0.29, Q=-0.07, E=-0.09, G=0.01, H=-0.79, I=0.27, L=0.1, K=-0.12, M=-0.31, F=0.21, P=-0.15, S=-0.19, T=0.03, W=0.28, Y=0.2, V=0.22) # log-odds in tandem repeat
		periods <- c(0.6281454, -4.807124, 2.698869, 1.994083, 0.9163926, 10.58584) # sigmoid fit for log-odds of periodicity relative to background
	} else {
		periods <- c(0.6217107, -3.638397, 0.2537369, 0.1148134, 0.3390507, 11.63033) # sigmoid fit for log-odds of periodicity relative to background
	}
	lens <- c(NaN, 0.2, -2, 0, -0.5, 0.6, 0.9, 2.9, 1.1) # log-odds of repeat copy number relative to background
	
	.score <- function(x,
		gapCost,
		POSL=NA_integer_,
		POSR=NA_integer_) {
		if (length(POSL) > 1L) { # not NA
			period <- POSR - POSL + 1L
			period <- range(period)
			period <- period[1L]/period[2L]
			if (period < 0.3)
				return(-Inf)
		}
		
		if (xtype == 3L) {
			structures <- mapply(function(a, b) struct[, a:b, drop=FALSE],
				POSL,
				POSR,
				SIMPLIFY=FALSE)
			if (useEmpirical) {
				AAs <- colSums(alphabetFrequency(x)[, 1:20])
				addScore <- sum(residues*AAs) # score for amino acid enrichment
			} else {
				addScore <- 0
			}
		} else {
			addScore <- 0
		}
		
		if (useEmpirical &&
			length(POSL) > 1L) { # not NA
			# score periodicity of repeats relative to random background
			addScore <- addScore + periods[1] + (periods[2] - periods[1])/(periods[3] + periods[4]*exp(periods[5]*mean(POSR - POSL + 1L, na.rm=TRUE)))^(1/periods[6])
			
			# score copy number relative to background
			if (length(x) > length(lens)) {
				addScore <- addScore + lens[length(lens)] # use longest possible length
			} else {
				addScore <- addScore + lens[length(x)]
			}
		}
		
		addScore + ScoreAlignment(x,
			method="adjacent",
			gapOpening=gapCost[1L],
			gapExtension=gapCost[2L],
			substitutionMatrix=sM,
			structures=structures,
			structureMatrix=structureMatrix,
			includeTerminalGaps=TRUE)
	}
	
	if (xtype == 3L) {
		# convert substitution matrix to units of log-odds (rather than two-thirds bits)
		SM <- matrix(c(0.951, -0.266, -0.305, -0.327, 0.099, -0.126, -0.151, 0.041, -0.383, -0.26, -0.274, -0.236, -0.131, -0.451, -0.14, 0.191, 0.008, -0.586, -0.501, 0.015, -2.542, -0.266, 1.465, 0.013, -0.153, -0.741, 0.37, 0.117, -0.453, 0.178, -0.81, -0.701, 0.679, -0.46, -0.874, -0.311, -0.097, -0.129, -0.5, -0.415, -0.679, -2.542, -0.305, 0.013, 1.494, 0.532, -0.582, 0.189, 0.129, 0.037, 0.255, -0.962, -0.947, 0.202, -0.606, -0.879, -0.252, 0.261, 0.075, -0.873, -0.433, -0.834, -2.542, -0.327, -0.153, 0.532, 1.575, -1.007, 0.155, 0.597, -0.131, -0.045, -1.265, -1.194, 0.052, -0.915, -1.235, -0.131, 0.099, -0.121, -1.102, -0.8, -1.051, -2.542, 0.099, -0.741, -0.582, -1.007, 3.127, -0.777, -0.995, -0.499, -0.438, -0.174, -0.218, -0.884, -0.137, -0.189, -0.832, -0.091, -0.185, -0.446, -0.268, 0.016, -2.542, -0.126, 0.37, 0.189, 0.155, -0.777, 1.289, 0.494, -0.368, 0.251, -0.762, -0.636, 0.433, -0.259, -0.846, -0.241, 0.046, -0.01, -0.707, -0.444, -0.624, -2.542, -0.151, 0.117, 0.129, 0.597, -0.995, 0.494, 1.287, -0.38, -0.057, -0.967, -0.931, 0.342, -0.646, -1.116, -0.162, 0.007, -0.072, -0.97, -0.681, -0.758, -2.542, 0.041, -0.453, 0.037, -0.131, -0.499, -0.368, -0.38, 1.768, -0.42, -1.087, -1.022, -0.369, -0.758, -0.924, -0.333, 0.043, -0.342, -0.886, -0.863, -0.861, -2.542, -0.383, 0.178, 0.255, -0.045, -0.438, 0.251, -0.057, -0.42, 2.254, -0.781, -0.663, 0.033, -0.433, -0.29, -0.354, -0.099, -0.206, -0.217, 0.381, -0.664, -2.542, -0.26, -0.81, -0.962, -1.265, -0.174, -0.762, -0.967, -1.087, -0.781, 1.184, 0.585, -0.819, 0.423, 0.216, -0.799, -0.716, -0.29, -0.347, -0.258, 0.785, -2.542, -0.274, -0.701, -0.947, -1.194, -0.218, -0.636, -0.931, -1.022, -0.663, 0.585, 1.087, -0.799, 0.585, 0.39, -0.777, -0.73, -0.43, -0.123, -0.159, 0.343, -2.542, -0.236, 0.679, 0.202, 0.052, -0.884, 0.433, 0.342, -0.369, 0.033, -0.819, -0.799, 1.282, -0.5, -1.005, -0.175, 0.006, -0.035, -0.829, -0.564, -0.704, -2.542, -0.131, -0.46, -0.606, -0.915, -0.137, -0.259, -0.646, -0.758, -0.433, 0.423, 0.585, -0.5, 1.637, 0.285, -0.712, -0.406, -0.171, -0.135, -0.091, 0.219, -2.542, -0.451, -0.874, -0.879, -1.235, -0.189, -0.846, -1.116, -0.924, -0.29, 0.216, 0.39, -1.005, 0.285, 1.717, -0.837, -0.7, -0.528, 0.608, 0.885, 0.045, -2.542, -0.14, -0.311, -0.252, -0.131, -0.832, -0.241, -0.162, -0.333, -0.354, -0.799, -0.777, -0.175, -0.712, -0.837, 2.121, -0.015, -0.198, -0.777, -0.763, -0.588, -2.542, 0.191, -0.097, 0.261, 0.099, -0.091, 0.046, 0.007, 0.043, -0.099, -0.716, -0.73, 0.006, -0.406, -0.7, -0.015, 0.979, 0.427, -0.727, -0.505, -0.505, -2.542, 0.008, -0.129, 0.075, -0.121, -0.185, -0.01, -0.072, -0.342, -0.206, -0.29, -0.43, -0.035, -0.171, -0.528, -0.198, 0.427, 1.128, -0.659, -0.439, -0.062, -2.542, -0.586, -0.5, -0.873, -1.102, -0.446, -0.707, -0.97, -0.886, -0.217, -0.347, -0.123, -0.829, -0.135, 0.608, -0.777, -0.727, -0.659, 3.153, 0.763, -0.428, -2.542, -0.501, -0.415, -0.433, -0.8, -0.268, -0.444, -0.681, -0.863, 0.381, -0.258, -0.159, -0.564, -0.091, 0.885, -0.763, -0.505, -0.439, 0.763, 2.023, -0.287, -2.542, 0.015, -0.679, -0.834, -1.051, 0.016, -0.624, -0.758, -0.861, -0.664, 0.785, 0.343, -0.704, 0.219, 0.045, -0.588, -0.505, -0.062, -0.428, -0.287, 1.084, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, 3.235),
			nrow=21,
			ncol=21,
			dimnames=list(c(AA_STANDARD, "*"), c(AA_STANDARD, "*")))
		freqs <- setNames(c(0.0774, 0.0552, 0.0397, 0.0533, 0.0179, 0.0421, 0.0652, 0.0683, 0.0242, 0.0513, 0.0971, 0.0535, 0.0231, 0.0382, 0.0547, 0.0753, 0.0545, 0.0124, 0.029, 0.0655, 0.0019),
			c(AA_STANDARD, "*"))
		
		# optimized structure matrix
		structures <- PredictHEC(myXStringSet,
			type="probabilities")
		structureMatrix <- matrix(c(0.102, -0.119, -0.136, -0.119, 0.188, -0.16, -0.136, -0.16, 0.109),
			nrow=3) # order is H, E, C
	} else {
		# optimized substitution matrix
		SM <- matrix(c(0.763, -0.713, -0.513, -0.818, -0.713, 0.887, -0.618, -0.485, -0.513, -0.618, 0.83, -1.113, -0.818, -0.485, -1.113, 0.896),
			nrow=4,
			ncol=4,
			dimnames=list(DNA_BASES, DNA_BASES))
		freqs <- setNames(c(0.274, 0.229, 0.255, 0.241),
			DNA_BASES)
		
		structures <- NULL
		structureMatrix <- NULL
	}
	BG <- outer(freqs, freqs)
	
	if (type == 1L || type == 3L) { # tandem repeats
		if (verbose) {
			if (type > 1L) {
				cat("Detecting tandem repeats:\n")
				flush.console()
			}
			pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
			totW <- c(0, cumsum(as.numeric(width(myXStringSet))))
		}
		
		if (xtype == 3L) {
			size <- .Call("alphabetSizeReducedAA",
				myXStringSet,
				alphabet,
				PACKAGE="DECIPHER")
		} else {
			size <- .Call("alphabetSize",
				myXStringSet,
				PACKAGE="DECIPHER")
		}
		K <- as.integer(ceiling(log(maxPeriod/N, size)))
		if (xtype == 3L) {
			if (K > n)
				K <- n
			kmers <- .Call("enumerateSequenceReducedAA",
				myXStringSet,
				K,
				alphabet,
				FALSE, # mask repeats
				processors,
				PACKAGE="DECIPHER")
		} else {
			if (K > 16L)
				K <- 16L
			kmers <- .Call("enumerateSequence",
				myXStringSet,
				K,
				FALSE, # mask repeats
				processors,
				PACKAGE="DECIPHER")
		}
		
		o <- lapply(kmers, order, method="radix")
		d <- lapply(o, diff)
		d <- mapply(function(a, b) {
					a[a] <- c(b, NA)
					a
				},
			o,
			d,
			SIMPLIFY=FALSE)
		r <- mapply(function(x, y) {
				start <- NA # start of current repeat
				off <- mult*K # lookahead length
				values <- numeric(length(x)) # periodicity
				lengths <- numeric(length(x)) # total evidence
				i <- i2 <- 1L # current and next index
				while (i < length(x)) {
					if (!is.na(x[i]) &&
						!is.na(y[i]) &&
						x[i] > 0) {
						if (is.na(start)) {
							start <- i
							divisor <- x[i]
						} else if (x[i] < divisor) {
							if ((divisor/x[i]) %in% 1:3) {
								divisor <- x[i]
								if (i + 2*divisor <= length(y) &&
									!is.na(y[i + 2*divisor]) &&
									y[i + 2*divisor] == y[i]) {
									i2 <- i + divisor
								} else if (i + off <= length(x) &&
									!is.na(x[i + off]) &&
									x[i + off] == x[i]) {
									i2 <- i + off
								}
							} else {
								if (i2 > i && i2 < i + off) {
									i <- i2
									next
								}
								values[start] <- divisor
								lengths[start] <- i - start
								start <- i
								divisor <- x[i]
							}
						} else if (x[i] > divisor) {
							if ((x[i]/divisor) %in% 1:3) {
								if (i + 2*divisor <= length(y) &&
									!is.na(y[i + 2*divisor]) &&
									y[i + 2*divisor] == y[i]) {
									i2 <- i + divisor
								} else if (i + off <= length(x) &&
									!is.na(x[i + off]) &&
									x[i + off] == x[i]) {
									i2 <- i + off
								}
							} else {
								if (i2 > i && i2 < i + off) {
									i <- i2
									next
								}
								values[start] <- divisor
								lengths[start] <- i - start
								start <- i
								divisor <- x[i]
							}
						} else { # x[i] == divisor
							if (i + 2*divisor <= length(y) &&
								!is.na(y[i + 2*divisor]) &&
								y[i + 2*divisor] == y[i]) {
								i2 <- i + divisor
							} else if (i + off <= length(x) &&
								!is.na(x[i + off]) &&
								x[i + off] == x[i]) {
								i2 <- i + off
							}
						}
					} else if (!is.na(start)) {
						if (i2 > i && i2 < i + off) {
							i <- i2
							next
						}
						values[start] <- divisor
						lengths[start] <- i - start
						start <- NA
					}
					
					i <- i + 1L
				}
				list(values=values, lengths=lengths)
			},
			d,
			kmers,
			SIMPLIFY=FALSE)
		
		shift <- function(x, LnL, posL, posR, maxShifts) {
			# optimize the phase
			tempL <- posL
			tempR <- posR
			POSL <- posL - 1L
			POSR <- posR - 1L
			attempts <- 0L
			while (POSL[1L] >= 1L &&
				attempts <= maxShifts) {
				y <- .Call("replaceGaps",
					x,
					myXString,
					POSL,
					xtype,
					PACKAGE="DECIPHER")
				temp <- .score(y, gapCost, POSL, POSR)
				if (temp > LnL) {
					x <- y
					LnL <- temp
					posL <- POSL
					posR <- POSR
					attempts <- 0L
				} else {
					attempts <- attempts + 1L
				}
				POSL <- POSL - 1L
				POSR <- POSR - 1L
			}
			
			POSL <- tempL + 1L
			POSR <- tempR + 1L
			attempts <- 0L
			while (POSR[length(POSR)] <= l &&
				attempts <= maxShifts) {
				y <- .Call("replaceGaps",
					x,
					myXString,
					POSL,
					xtype,
					PACKAGE="DECIPHER")
				temp <- .score(y, gapCost, POSL, POSR)
				if (temp > LnL) {
					x <- y
					LnL <- temp
					posL <- POSL
					posR <- POSR
					attempts <- 0L
				} else {
					attempts <- attempts + 1L
				}
				POSL <- POSL + 1L
				POSR <- POSR + 1L
			}
			
			list(x, LnL, posL, posR)
		}
		
		result <- list()
		for (k in seq_along(r)) {
			myXString <- myXStringSet[[k]]
			l <- length(myXString)
			struct <- structures[[k]]
			values <- r[[k]]$values
			lengths <- r[[k]]$lengths
			
			if (xtype == 3L) {
				bg <- alphabetFrequency(myXString,
					as.prob=TRUE)[c(1:20, 27)]
				bg <- bg/sum(bg)
			} else {
				bg <- oligonucleotideFrequency(myXString,
					width=1,
					as.prob=TRUE)
			}
			sM <- SM + log(BG/outer(bg, bg))
			
			# only repeats occurring more frequently than expected
			w <- which(values/lengths <= maxPeriod/N &
				lengths <= maxPeriod)
			
			# reduce to the set of top ranked k-mer repeats
			visited <- integer(l)
			keep <- logical(length(w))
			o <- order(-lengths[w], values[w])
			for (i in seq_along(o)) {
				posL <- w[o[i]]
				posR <- posL + values[w[o[i]]] - 1L
				if (posR > l)
					posR <- l
				if (all(visited[posL:posR] > maxVisits))
					next
				visited[posL:posR] <- visited[posL:posR] + 1L
				keep[o[i]] <- TRUE
			}
			w <- w[keep]
			
			res <- vector("list", length(w))
			for (i in seq_along(w)) {
				posL <- seq(0,
					lengths[w[i]] + values[w[i]],
					values[w[i]])
				posL <- posL + w[i]
				posR <- posL + values[w[i]] - 1L
				keep <- posL <= l &
					posR <= l
				if (sum(keep) < length(keep)) {
					posL <- posL[keep]
					posR <- posR[keep]
				}
				
				delta <- as.integer(values[w[i]]/2)
				if (length(posR) > 1) { # extend unknown right bound
					posR[length(posR)] <- posR[length(posR)] + delta
					if (posR[length(posR)] > l)
						posR[length(posR)] <- l
				}
				
				if (length(posR) < 2L) {
					res[[i]] <- list(posL, posR, -Inf, k)
					if (verbose)
						setTxtProgressBar(pBar,
							(totW[k] + w[i])/totW[length(totW)])
					next
				}
				
				# align the repeats
				x <- extractAt(myXString, IRanges(posL, posR))
				ux <- unique(x)
				if (length(ux) > 1) {
					index <- match(x, ux)
					if (length(ux) == 2) {
						ux <- AlignProfiles(.subset(ux, 1),
							.subset(ux, 2),
							anchor=NA,
							processors=processors)
					} else {
						ux <- AlignSeqs(ux,
							iterations=0,
							refinements=0,
							anchor=NA,
							processors=processors,
							verbose=FALSE)
					}
					x <- .subset(ux, index)
					
					t <- TerminalChar(x)
					off <- min(t[-nrow(t), "trailingChar"])
					x <- subseq(x, end=width(x)[1L] - off)
					posR[length(posR)] <- posR[length(posR)] - off
				}
				
				# calculate the likelihood
				LnL <- .score(x, gapCost, posL, posR)
				if (posR[length(posR)] < posL[length(posL)] || LnL < 0) {
					res[[i]] <- list(posL, posR, LnL, k)
					if (verbose)
						setTxtProgressBar(pBar,
							(totW[k] + w[i])/totW[length(totW)])
					next
				}
				
				# try extending left
				start <- posL[1L] - values[w[i]]
				attempts <- 0L
				X <- x
				POSL <- posL
				POSR <- posR
				while (start >= 1L &&
					attempts <= maxFailures) {
					start <- start - delta
					if (start < 1)
						start <- 1L
					subseq <- .extract(myXStringSet,
						k,
						start,
						POSL[1L] - 1L)
					y <- AlignProfiles(subseq,
						X,
						anchor=NA,
						processors=processors)
					t <- TerminalChar(y)
					off <- min(t[-1L, "leadingChar"])
					y <- subseq(y, off + 1L)
					start <- start + off
					if (start >= POSL[1L])
						break # no overlap in alignment
					POSR <- c(POSL[1L] - 1L, POSR)
					POSL <- c(start, POSL)
					temp <- .score(y, gapCost, POSL, POSR)
					
					X <- y
					if (temp > LnL) {
						LnL <- temp
						posL <- POSL
						posR <- POSR
						x <- y
						attempts <- 0L
					} else {
						attempts <- attempts + 1L
					}
					start <- POSL[1L] - values[w[i]]
				}
				
				# try extending right
				end <- posR[length(posR)] + values[w[i]]
				attempts <- 0L
				X <- x
				POSL <- posL
				POSR <- posR
				while (end <= l &&
					attempts <= maxFailures) {
					end <- end + delta
					if (end > l)
						end <- l
					subseq <- .extract(myXStringSet,
						k,
						POSR[length(POSR)] + 1L,
						end)
					y <- AlignProfiles(X,
						subseq,
						anchor=NA,
						processors=processors)
					t <- TerminalChar(y)
					off <- min(t[-nrow(t), "trailingChar"])
					y <- subseq(y, end=width(y)[1L] - off)
					end <- end - off
					if (end <= POSR[length(POSR)])
						break # no overlap in alignment
					POSL <- c(POSL, POSR[length(POSR)] + 1L)
					POSR <- c(POSR, end)
					temp <- .score(y, gapCost, POSL, POSR)
					
					X <- y
					if (temp > LnL) {
						LnL <- temp
						posL <- POSL
						posR <- POSR
						x <- y
						attempts <- 0L
					} else {
						attempts <- attempts + 1L
					}
					end <- POSR[length(POSR)] + values[w[i]]
				}
				
				out <- shift(x, LnL, posL, posR, min(values[w[i]], maxShifts))
				x <- out[[1L]]
				LnL <- out[[2L]]
				posL <- out[[3L]]
				posR <- out[[4L]]
				
				# record the result
				res[[i]] <- list(posL, posR, LnL, k)
				if (verbose)
					setTxtProgressBar(pBar,
						(totW[k] + w[i])/totW[length(totW)])
			}
			
			result <- c(result, res)
		}
		
		result <- unique(result)
		if (length(result) > 0) {
			posL <- lapply(result, `[[`, 1L)
			posR <- lapply(result, `[[`, 2L)
			result <- data.frame(Index=sapply(result, `[[`, 4L),
				Begin=sapply(posL, `[`, 1L),
				End=sapply(posR, tail, n=1L),
				Left=I(posL),
				Right=I(posR),
				Score=sapply(result, `[[`, 3L))
			result <- result[result[, "Score"] >= minScore,]
			result <- result[order(result[, "Index"], result[, "Begin"]),]
		} else {
			result <- data.frame(Index=integer(),
				Begin=numeric(),
				End=numeric(),
				Left=I(numeric()),
				Right=I(numeric()),
				Score=numeric())
		}
		
		if (!allScores) {
			# pick highest score when overlapping
			if (nrow(result) > 1) {
				keep <- logical(nrow(result))
				keep [1L] <- TRUE
				i <- 1L
				j <- 2L
				while (j <= nrow(result)) {
					if (result[i, "Index"] == result[j, "Index"] &&
						result[i, "End"] >= result[j, "Begin"]) {
						if (result[i, "Score"] < result[j, "Score"]) {
							keep[i] <- FALSE
							keep[j] <- TRUE
							i <- j
						}
					} else {
						keep[j] <- TRUE
						i <- j
					}
					j <- j + 1L
				}
				result <- result[keep,]
			}
		}
		rownames(result) <- NULL
		
		if (verbose) {
			setTxtProgressBar(pBar, 1)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\n")
			if (type > 1L)
				cat("Detecting interspersed repeats:\n")
			flush.console()
		}
	}
	
	if (type > 1L) { # interspersed repeats
		dbConn <- dbConnect(SQLite(), ":memory:")
		on.exit(dbDisconnect(dbConn))
		Seqs2DB(myXStringSet,
			"XStringSet",
			dbConn,
			"1",
			verbose=FALSE)
		
		syn <- FindSynteny(dbConn,
			identifier=c("1", "1"),
			verbose=FALSE,
			...)
		
		# remove results within maxPeriod
		w <- which(syn[[2, 1]][, "index1"] == syn[[2, 1]][, "index2"] &
			abs(syn[[2, 1]][, "start1"] - syn[[2, 1]][, "start2"]) > maxPeriod &
			abs(syn[[2, 1]][, "end1"] - syn[[2, 1]][, "end2"]) > maxPeriod |
			syn[[2, 1]][, "index1"] != syn[[2, 1]][, "index2"])
		if (length(w) > 0) {
			hits <- mapply(seq,
				syn[[2, 1]][w, "first_hit"],
				syn[[2, 1]][w, "last_hit"],
				SIMPLIFY=FALSE)
			syn[[2, 1]] <- syn[[2, 1]][w,]
			syn[[1, 2]] <- syn[[1, 2]][unlist(hits),]
			hits <- cumsum(lengths(hits))
			syn[[2, 1]][, "first_hit"] <- c(1L, hits[-length(hits)] + 1L)
			syn[[2, 1]][, "last_hit"] <- hits
			
			ali <- AlignSynteny(syn,
				dbConn,
				processors=processors,
				verbose=verbose)
			
			ali <- ali[[1L]]
			res2 <- data.frame(syn[[2, 1]][, 1:8])
			if (length(ali) > 0) {
				bg <- oligonucleotideFrequency(myXStringSet,
					width=1,
					as.prob=TRUE,
					simplify.as="collapse")
				sM <- SM + log(BG/outer(bg, bg))
				
				res2[, "score"] <- sapply(ali,
					.score,
					gapCost=gapCost)
				res2 <- res2[res2[, "score"] >= minScore,]
				
				if (!allScores) {
					# eliminate overlapping repeats
					keep <- rep(TRUE, nrow(res2))
					o <- order(res2[, "score"], decreasing=TRUE)
					range1 <- IRanges(res2[, "start1"],
						res2[, "end1"])
					range2 <- IRanges(res2[, "start2"],
						res2[, "end2"])
					for (i in seq_along(o)[-1]) {
						w <- o[which(keep[o[seq_len(i - 1)]])]
						int <- intersect(range1[o[i],],
							range1[w[res2[o[i], "index1"] == res2[w, "index1"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
						int <- intersect(range1[o[i],],
							range2[w[res2[o[i], "index1"] == res2[w, "index2"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
						int <- intersect(range2[o[i],],
							range1[w[res2[o[i], "index2"] == res2[w, "index1"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
						int <- intersect(range2[o[i],],
							range2[w[res2[o[i], "index2"] == res2[w, "index2"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
					}
					res2 <- res2[keep,]
					rownames(res2) <- NULL
				}
			}
		} else {
			res2 <- data.frame(syn[[2, 1]][w, 1:8])
		}
		
		if (type == 2L) { # interspersed
			result <- res2
		} else { # both
			result <- list(result, res2)
		}
	}
	
	return(result)
}

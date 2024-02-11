DetectRepeats <- function(myXStringSet,
	type="tandem",
	minScore=10,
	allScores=FALSE,
	maxCopies=1000,
	maxPeriod=1000,
	maxFailures=3,
	maxShifts=5,
	alphabet=AA_REDUCED[[77]],
	useEmpirical=TRUE,
	correctBackground=TRUE,
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
	if (!is.numeric(maxCopies) || length(maxCopies) > 1)
		stop("maxCopies must be a numeric.")
	if (maxCopies < 2)
		stop("maxCopies must be at least two.")
	if (maxCopies != floor(maxCopies))
		stop("maxCopies must be a whole number.")
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
	if (!is.logical(correctBackground))
		stop("correctBackground must be a logical.")
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
		processors <- .detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	if (verbose)
		time.1 <- Sys.time()
	
	# default parameters
	N <- 10 # find k-mers on average N times per maxPeriod
	mult <- 2L # multiplier on K for lookahead length
	maxVisits <- 2L # maximum attempts to detect repeats in each position
	if (useEmpirical) {
		if (xtype == 3L) {
			if (correctBackground) {
				correction <- 0.32 # correction factor to standardize log-odds scores
			} else {
				correction <- 0.21 # correction factor to standardize log-odds scores
			}
			offset <- 0.06 # calibrate substitution matrix
			residues <- c(A=0.05, R=-0.02, N=0.05, D=0.06, C=0.29, Q=-0.07, E=-0.09, G=0.01, H=-0.79, I=0.27, L=0.1, K=-0.12, M=-0.31, F=0.21, P=-0.15, S=-0.19, T=0.03, W=0.28, Y=0.2, V=0.22) # log-odds in tandem repeat
			periods <- c(0.6281454, -4.807124, 2.698869, 1.994083, 0.9163926, 10.58584) # sigmoid fit for log-odds of periodicity relative to background
		} else {
			if (correctBackground) {
				correction <- 0.61 # correction factor to standardize log-odds scores
			} else {
				correction <- 0.33 # correction factor to standardize log-odds scores
			}
			periods <- c(0.6217107, -3.638397, 0.2537369, 0.1148134, 0.3390507, 11.63033) # sigmoid fit for log-odds of periodicity relative to background
		}
		lens <- c(NaN, 0.2, -2, 0, -0.5, 0.6, 0.9, 2.9, 1.1) # log-odds of repeat copy number relative to background
		gapCost <- c(-3, -0.2) # gap open and extension coefficients
		coef <- -0.65 # coefficient of the normalized quadratic column weight function (> -1)
	} else {
		if (xtype == 3L) {
			if (correctBackground) {
				correction <- 0.48 # correction factor to standardize log-odds scores
			} else {
				correction <- 0.19 # correction factor to standardize log-odds scores
			}
		} else {
			if (correctBackground) {
				correction <- 0.52 # correction factor to standardize log-odds scores
			} else {
				correction <- 0.23 # correction factor to standardize log-odds scores
			}
		}
		gapCost <- c(-2.6, -0.1) # gap open and extension coefficients
		coef <- 0 # coefficient of the normalized quadratic column weight function (> -1)
	}
	if (correctBackground)
		samples <- 5000 # total relative weight of background frequencies
	
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
				AAs <- alphabetFrequency(x, collapse=TRUE)[1:20]
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
		
		colScores <- ScoreAlignment(x,
			method="adjacent",
			type="scores",
			gapOpening=gapCost[1L],
			gapExtension=gapCost[2L],
			substitutionMatrix=sM,
			structures=structures,
			structureMatrix=structureMatrix,
			includeTerminalGaps=TRUE)
		
		# apply weight function to alignment columns
		w <- seq(-1, 1, length.out=length(colScores))
		w <- w*w
		w <- 1 + coef*w
		w <- w/mean(w) # normalize to mean of 1
		
		addScore + sum(w*colScores)
	}
	
	if (xtype == 3L) {
		# convert PFASUM40 to units of log-odds (PFASUM40*log(2)/3)
		SM <- matrix(c(0.8324524, -0.218532, -0.2657873, -0.2841384, 0.1137628, -0.1057223, -0.132703, 0.0656757, -0.3322948, -0.1996091, -0.2015672, -0.2129695, -0.07751118, -0.3568842, -0.1186668, 0.1685214, 0.01476403, -0.4825517, -0.4096673, 0.03146888, -2.54154, -0.218532, 1.324362, 0.05157015, -0.07413209, -0.6492536, 0.3704525, 0.1644838, -0.3646994, 0.1869765, -0.722242, -0.6231913, 0.6354253, -0.4040008, -0.7596027, -0.2274043, -0.05765252, -0.08515313, -0.4235649, -0.3423107, -0.5975795, -2.54154, -0.2657873, 0.05157015, 1.326718, 0.5145925, -0.5346764, 0.203716, 0.1608275, 0.06898547, 0.2330881, -0.8566086, -0.8444439, 0.227283, -0.5355428, -0.7694627, -0.1831468, 0.2529467, 0.08686867, -0.763259, -0.3633478, -0.7466928, -2.54154, -0.2841384, -0.07413209, 0.5145925, 1.451606, -0.8815099, 0.1895064, 0.5766811, -0.08309102, -0.01270192, -1.123522, -1.065003, 0.1085988, -0.8130616, -1.092175, -0.06352694, 0.1184415, -0.0733523, -0.9669403, -0.6929566, -0.9349516, -2.54154, 0.1137628, -0.6492536, -0.5346764, -0.8815099, 2.978956, -0.6824207, -0.8679589, -0.4107417, -0.3668135, -0.09125283, -0.1317846, -0.7753198, -0.05883087, -0.1113021, -0.709956, -0.07350826, -0.1492519, -0.3533838, -0.19876, 0.06574501, -2.54154, -0.1057223, 0.3704525, 0.203716, 0.1895064, -0.6824207, 1.120005, 0.4807669, -0.3028533, 0.2300036, -0.6820222, -0.5751562, 0.4190595, -0.2444557, -0.7449253, -0.1849837, 0.06706199, 0.01455609, -0.5997803, -0.368685, -0.5563199, -2.54154, -0.132703, 0.1644838, 0.1608275, 0.5766811, -0.8679589, 0.4807669, 1.167797, -0.3148448, -0.01814313, -0.8708528, -0.8331283, 0.362568, -0.5800775, -0.9887745, -0.1153397, 0.03599167, -0.03429346, -0.8595372, -0.5824169, -0.6889883, -2.54154, 0.0656757, -0.3646994, 0.06898547, -0.08309102, -0.4107417, -0.3028533, -0.3148448, 1.654577, -0.3493808, -0.928332, -0.8720658, -0.3017443, -0.6340564, -0.7735349, -0.2399156, 0.07392415, -0.2633959, -0.7408184, -0.7190016, -0.7276139, -2.54154, -0.3322948, 0.1869765, 0.2330881, -0.01270192, -0.3668135, 0.2300036, -0.01814313, -0.3493808, 2.122399, -0.692714, -0.589487, 0.05004523, -0.3718041, -0.2701368, -0.2695303, -0.07399346, -0.1642239, -0.1734601, 0.3390183, -0.5845483, -2.54154, -0.1996091, -0.722242, -0.8566086, -1.123522, -0.09125283, -0.6820222, -0.8708528, -0.928332, -0.692714, 1.05166, 0.5674103, -0.7364169, 0.4077265, 0.2621656, -0.692662, -0.6180274, -0.2506767, -0.2521496, -0.1753489, 0.7226753, -2.54154, -0.2015672, -0.6231913, -0.8444439, -1.065003, -0.1317846, -0.5751562, -0.8331283, -0.8720658, -0.589487, 0.5674103, 0.9706833, -0.7245641, 0.544987, 0.4017828, -0.6738084, -0.6305906, -0.3651673, -0.05576369, -0.09674602, 0.3510964, -2.54154, -0.2129695, 0.6354253, 0.227283, 0.1085988, -0.7753198, 0.4190595, 0.362568, -0.3017443, 0.05004523, -0.7364169, -0.7245641, 1.153068, -0.4585688, -0.8880428, -0.1207289, 0.03126094, -0.0091842, -0.7176846, -0.4753603, -0.6393243, -2.54154, -0.07751118, -0.4040008, -0.5355428, -0.8130616, -0.05883087, -0.2444557, -0.5800775, -0.6340564, -0.3718041, 0.4077265, 0.544987, -0.4585688, 1.370144, 0.2994396, -0.6052734, -0.3403353, -0.1373991, -0.06671542, -0.02880027, 0.2333827, -2.54154, -0.3568842, -0.7596027, -0.7694627, -1.092175, -0.1113021, -0.7449253, -0.9887745, -0.7735349, -0.2701368, 0.2621656, 0.4017828, -0.8880428, 0.2994396, 1.51541, -0.7086563, -0.5942004, -0.433321, 0.594391, 0.8069966, 0.1037295, -2.54154, -0.1186668, -0.2274043, -0.1831468, -0.06352694, -0.709956, -0.1849837, -0.1153397, -0.2399156, -0.2695303, -0.692662, -0.6738084, -0.1207289, -0.6052734, -0.7086563, 1.985312, 0.01431349, -0.1457342, -0.6387178, -0.6383886, -0.5052523, -2.54154, 0.1685214, -0.05765252, 0.2529467, 0.1184415, -0.07350826, 0.06706199, 0.03599167, 0.07392415, -0.07399346, -0.6180274, -0.6305906, 0.03126094, -0.3403353, -0.5942004, 0.01431349, 0.8363687, 0.3843674, -0.6161732, -0.4164255, -0.4361282, -2.54154, 0.01476403, -0.08515313, 0.08686867, -0.0733523, -0.1492519, 0.01455609, -0.03429346, -0.2633959, -0.1642239, -0.2506767, -0.3651673, -0.0091842, -0.1373991, -0.433321, -0.1457342, 0.3843674, 0.9569763, -0.5494924, -0.3483065, -0.05342432, -2.54154, -0.4825517, -0.4235649, -0.763259, -0.9669403, -0.3533838, -0.5997803, -0.8595372, -0.7408184, -0.1734601, -0.2521496, -0.05576369, -0.7176846, -0.06671542, 0.594391, -0.6387178, -0.6161732, -0.5494924, 2.96544, 0.7479405, -0.3300074, -2.54154, -0.4096673, -0.3423107, -0.3633478, -0.6929566, -0.19876, -0.368685, -0.5824169, -0.7190016, 0.3390183, -0.1753489, -0.09674602, -0.4753603, -0.02880027, 0.8069966, -0.6383886, -0.4164255, -0.3483065, 0.7479405, 1.820603, -0.2030748, -2.54154, 0.03146888, -0.5975795, -0.7466928, -0.9349516, 0.06574501, -0.5563199, -0.6889883, -0.7276139, -0.5845483, 0.7226753, 0.3510964, -0.6393243, 0.2333827, 0.1037295, -0.5052523, -0.4361282, -0.05342432, -0.3300074, -0.2030748, 0.956907, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, -2.54154, 3.003638),
			nrow=21,
			ncol=21,
			dimnames=list(c(AA_STANDARD, "*"), c(AA_STANDARD, "*")))
		if (useEmpirical)
			SM <- SM + offset
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
	if (correctBackground) {
		BG <- outer(freqs, freqs)
	} else {
		sM <- SM
	}
	
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
			size <- tapply(freqs[1:20], alphabet, sum)
			size <- exp(-sum(size*log(size)))
		} else {
			size <- exp(-sum(freqs*log(freqs)))
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
				FALSE, # mask low complexity regions
				processors,
				PACKAGE="DECIPHER")
		} else {
			if (K > 16L)
				K <- 16L
			kmers <- .Call("enumerateSequence",
				myXStringSet,
				K,
				FALSE, # mask repeats
				FALSE, # mask low complexity regions
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
							quotient <- divisor/x[i]
							if (quotient == 1 || quotient == 2 || quotient == 3) {
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
							quotient <- x[i]/divisor
							if (quotient == 1 || quotient == 2 || quotient == 3) {
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
		
		.shift <- function(x, LnL, posL, posR, maxShifts) {
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
			
			if (correctBackground) {
				if (xtype == 3L) {
					bg <- alphabetFrequency(myXString)[c(1:20, 27)]
				} else {
					bg <- oligonucleotideFrequency(myXString, width=1L)
				}
				bg <- bg + samples*freqs
				bg[bg == 0] <- 1 # add pseudocounts (needed if samples is zero)
				bg <- bg/sum(bg)
				sM <- SM + log(BG/outer(bg, bg))
			}
			
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
				if (.Call("all",
					visited[posL:posR] > maxVisits,
					PACKAGE="DECIPHER"))
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
				} else if (length(posR) > maxCopies) {
					length(posL) <- maxCopies
					length(posR) <- maxCopies
				}
				
				# align the repeats
				x <- .extractSet(myXString, posL, posR)
				ux <- unique(x)
				if (length(ux) > 1) {
					index <- match(x, ux)
					if (xtype == 3L) {
						rev_index <- match(ux, x)
						if (length(ux) == 2) {
							ux <- AlignProfiles(.subset(ux, 1),
								.subset(ux, 2),
								p.struct=struct[, posL[rev_index[1L]]:posR[rev_index[1L]], drop=FALSE],
								s.struct=struct[, posL[rev_index[2L]]:posR[rev_index[2L]], drop=FALSE],
								anchor=NA,
								processors=processors)
						} else {
							ux <- AlignSeqs(ux,
								iterations=0,
								refinements=0,
								structures=mapply(function(a, b)
										struct[, a:b, drop=FALSE],
									posL[rev_index],
									posR[rev_index],
									SIMPLIFY=FALSE),
								anchor=NA,
								processors=processors,
								verbose=FALSE)
						}
					} else {
						if (length(ux) == 2) {
							ux <- AlignProfiles(.subset(ux, 1),
								.subset(ux, 2),
								anchor=NA,
								processors=processors)
						} else {
							ux <- AlignSeqs(ux,
								iterations=0,
								refinements=0,
								useStructures=FALSE,
								anchor=NA,
								processors=processors,
								verbose=FALSE)
						}
					}
					x <- .subset(ux, index)
					
					t <- TerminalChar(x)
					off <- min(t[-nrow(t), "trailingChar"])
					x <- subseq(x, end=width(x)[1L] - off)
					posR[length(posR)] <- posR[length(posR)] - off
					off <- min(t[-1L, "leadingChar"])
					x <- subseq(x, off + 1L)
					posL[1L] <- posL[1L] + off
					
					if (width(x)[1L] == 0L) {
						res[[i]] <- list(posL, posR, -Inf, k)
						if (verbose)
							setTxtProgressBar(pBar,
								(totW[k] + w[i])/totW[length(totW)])
						next
					}
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
				
				prevLnL <- 0
				while (LnL > prevLnL) {
					# try extending left
					start <- posL[1L] - values[w[i]]
					attempts <- 0L
					X <- x
					POSL <- posL
					POSR <- posR
					while (start >= 1L &&
						attempts <= maxFailures &&
						length(POSL) < maxCopies) {
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
						attempts <= maxFailures &&
						length(POSL) < maxCopies) {
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
					
					if (LnL > prevLnL) {
						prevLnL <- LnL
						out <- .shift(x, LnL, posL, posR, min(diff(posL), maxShifts))
						x <- out[[1L]]
						LnL <- out[[2L]]
						posL <- out[[3L]]
						posR <- out[[4L]]
					} else { # no change since last attempt
						break # LnL == prevLnL
					}
				}
				
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
				Score=sapply(result, `[[`, 3L)*correction)
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
		if (!requireNamespace("RSQLite", quietly=TRUE))
			stop("Package 'RSQLite' must be installed.")
		dbConn <- dbConnect(dbDriver("SQLite"), ":memory:")
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
				if (correctBackground) {
					bg <- oligonucleotideFrequency(myXStringSet,
						width=1,
						simplify.as="collapse")
					bg <- bg + samples*freqs
					bg[bg == 0] <- 1 # add pseudocounts (needed if samples is zero)
					bg <- bg/sum(bg)
					sM <- SM + log(BG/outer(bg, bg))
				}
				
				res2[, "score"] <- sapply(ali,
					.score,
					gapCost=gapCost)*correction
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

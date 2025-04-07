InferSelection <- function(myDNAStringSet,
	readingFrame=1L,
	windowSize=NA,
	tolerance=5e-5,
	geneticCode=GENETIC_CODE,
	showPlot=FALSE,
	verbose=TRUE) {
	
	# error checking
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (length(myDNAStringSet) < 2L)
		stop("At least two sequences are required in myDNAStringSet.")
	uw <- unique(width(myDNAStringSet))
	if (length(uw) != 1L)
		stop("Sequences in myDNAStringSet must be the same width (aligned).")
	if (!is.numeric(readingFrame))
		stop("readingFrame must be a numeric.")
	if (length(readingFrame) != 1L)
		stop("readingFrame must be a single numeric.")
	if (!readingFrame %in% 1L:3L)
		stop("readingFrame must be 1, 2, or 3.")
	readingFrame <- as.integer(readingFrame)
	if (uw - readingFrame < 2)
		stop("myDNAStringSet must contain at least one codon.")
	if (length(windowSize) != 1L)
		stop("windowSize must be length one.")
	windows <- seq_len((uw - readingFrame + 1L) %/% 3L)
	if (is.na(windowSize)) {
		windows <- list(windows)
	} else if (!is.numeric(windowSize)) {
		stop("windowSize must be a numeric.")
	} else if (windowSize < 1) {
		stop("windowSize must be at least 1.")
	} else {
		windows <- split(windows,
			ceiling(seq_along(windows)/windowSize))
	}
	l <- length(windows)
	if (!is.character(geneticCode))
		stop("geneticCode must be a character vector.")
	if (length(geneticCode) != 64L)
		stop("geneticCode must include 64 codons.")
	if (any(is.na(geneticCode)))
		stop("geneticCode must not contain NA values.")
	if (is.null(names(geneticCode)))
		stop("geneticCode must have names.")
	if (any(nchar(geneticCode) != 1L) ||
		any(!geneticCode %in% c(AA_STANDARD, "*")))
		stop("geneticCode may only contain proteinogenic amino acids.")
	if (any(nchar(names(geneticCode)) != 3L) ||
		sum(grepl("[^ACTG]", names(geneticCode))) > 0L)
		stop("geneticCode must be named by codon.")
	if (any(duplicated(names(geneticCode))))
		stop("Duplicated codons in geneticCode.")
	if (!isTRUEorFALSE(showPlot))
		stop("showPlot must be a logical.")
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be a logical.")
	
	if (verbose)
		time.1 <- Sys.time()
	
	# initialize parameters
	iniShift <- 0.1 # initial shift for optimization (> 0)
	maxShift <- 1 # maximum shift per optimization iteration (> 0)
	minSlope <- 1e-2 # minimum slope to continue optimization (> 0)
	limits <- c(100, 10) # initial range of parameters (fold > 1)
	U <- 0.1 # initial theta (> 0)
	K <- 2 # initial kappa (> 0)
	omegas <- rep(-1, l) # initial log(omega) (> 0)
	minVal <- sqrt(.Machine$double.eps) # machine precision (> 0)
	.print <- function(val)
		formatC(val, digits=3, format="f")
	
	triplets <- strsplit(names(geneticCode), "", fixed=TRUE) # nucleotides per codon
	codons <- tapply(seq_along(geneticCode), geneticCode, c) # amino acid groups
	codons <- codons[-match("*", names(codons))]
	codons <- unlist(unname(codons))
	n <- length(codons)
	
	COUNTS <- sapply(unlist(windows),
		function(p)
			tabulate(match(subseq(myDNAStringSet,
						3L*(p - 1L) + readingFrame,
						3L*(p - 1L) + readingFrame + 2L),
					names(geneticCode)),
				64L)[codons],
		USE.NAMES=FALSE)
	
	# define components of codon substitution rate matrix
	rates <- kappa <- omega <- both <- matrix(0, n, n, dimnames=list(names(geneticCode)[codons], names(geneticCode)[codons]))
	for (i in 1L:(n - 1L)) {
		for (j in (i + 1L):n) {
			w <- which(triplets[[codons[i]]] != triplets[[codons[j]]])
			if (length(w) == 1L) {
				base1 <- triplets[[codons[i]]][w]
				base2 <- triplets[[codons[j]]][w]
				if ((base1 == "A" && base2 == "G") ||
					(base1 == "G" && base2 == "A") ||
					(base1 == "C" && base2 == "T") ||
					(base1 == "T" && base2 == "C")) { # transition
					if (geneticCode[codons[i]] == geneticCode[codons[j]]) { # synonymous
						kappa[i, j] <- kappa[j, i] <- 1
					} else { # non-synonymous
						both[i, j] <- both[j, i] <- 1
					}
				} else if (geneticCode[codons[i]] == geneticCode[codons[j]]) { # synonymous
					rates[i, j] <- rates[j, i] <- 1
				} else { # non-synonymous
					omega[i, j] <- omega[j, i] <- 1
				}
			}
		}
	}
	
	# equilibrium codon frequencies
	freqs <- oligonucleotideFrequency(myDNAStringSet,
		width=1,
		simplify.as="collapse")
	freqs[freqs == 0L] <- 1L
	freqs <- freqs/sum(freqs)
	freqs <- sapply(triplets[codons],
		function(x)
			prod(freqs[x]))
	freqs <- freqs/sum(freqs)
	
	.opt <- function(params, optOmegas=FALSE, getPvals=FALSE) {
		U <- params[1L] # theta
		K <- params[2L] # kappa
		
		.lik <- function(O, w) {
			# initialize rate matrix
			O <- exp(O) # ensure positive
			theta <- rates + K*kappa + O*omega + K*O*both
			
			# apply equilibrium frequencies
			sF <- sqrt(freqs)
			theta <- .Call("applyFreqs", theta, sF, PACKAGE="DECIPHER")
			
			# eigen decomposition of rate matrix
			D <- eigen(theta, TRUE)
			V <- D$vectors
			D <- D$values
			
			# compute conditional probabilities
			D <- 1/(1 + U/sum(freqs*diag(theta))*D)
			alpha <- .Call("conditionalProbs", V, D, sF, minVal, PACKAGE="DECIPHER")
			talpha <- t(alpha)
			
			tmu <- rowSums(alpha) - 1
			LnLini <- log(freqs) - rowSums(lgamma(alpha)) + lgamma(tmu) + log(tmu)
			
			# calculate likelihood
			likelihood <- 0
			LnLs <- rep(-Inf, length(positions))
			for (p in seq_along(positions)) {
				counts <- COUNTS[, positions[p]]
				
				N <- sum(counts)
				if (N == 0L)
					next # no protein coding codons
				LnLs <- LnLini + colSums(lgamma(counts + talpha)) - lgamma(N + tmu) - log(N + tmu)
				
				# add probabilities in log-space
				LnL <- -Inf
				for (i in seq_len(n)) { # each ancestral state
					diff <- LnL - LnLs[i]
					if (diff == 0) {
						LnL <- LnL + 0.6931471805599452862268 # log(2)
					} else if (diff < 0) {
						LnL <- LnLs[i] + log(1 + exp(diff))
					} else { # diff > 0
						LnL <- LnL + log(1 + exp(-diff))
					}
				}
				
				# add the multinomial coefficient
				likelihood <- likelihood + LnL + lgamma(N + 1L) - sum(lgamma(counts + 1L))
			}
			
			likelihood
		}
		
		total <- 0
		for (w in seq_len(l)) { # each window
			positions <- windows[[w]]
			
			O <- omegas[w] # initial value
			if (optOmegas) { # optimize omega
				delta <- iniShift
				shift <- Inf
				while (shift > 1e-4) {
					l0 <- .lik(O - delta, w)
					l1 <- .lik(O, w)
					l2 <- .lik(O + delta, w)
					
					df <- (l2 - l0)/(2*delta)
					if (abs(df) < minSlope)
						break
					ddf <- (l2 - 2*l1 + l0)/delta^2
					
					if (ddf == 0) {
						shift <- delta
					} else {
						shift <- abs(df/ddf)
						if (shift > maxShift)
							shift <- maxShift
					}
					if (df > 0) {
						O <- O + shift
					} else {
						O <- O - shift
					}
					delta <- delta/2
				}
				omegas[w] <<- O
			}
			
			likelihood <- .lik(O, w)
			if (getPvals) # compare to omega = 1
				pvals[w] <<- pchisq(2*(likelihood - .lik(1, w)),
					df=1L,
					lower.tail=FALSE)
			
			total <- total + likelihood
		}
		
		total
	}
	
	# optimize parameters
	lastLnL <- -Inf
	repeat {
		o <- optim(c(U, K),
			.opt,
			method="L-BFGS-B",
			lower=c(U, K)/limits,
			upper=c(U, K)*limits,
			control=list(fnscale=-1))
		
		if (verbose && interactive())
			cat("\rLnL = ",
				.print(o$value),
				sep="")
		
		U <- o$par[1L]
		K <- o$par[2L]
		bestLnL <- .opt(c(U, K), optOmegas=TRUE)
		
		if (lastLnL - bestLnL < bestLnL*tolerance) {
			if (verbose && interactive())
				cat("\rLnL = ",
					.print(bestLnL),
					sep="")
			lastLnL <- bestLnL
		} else {
			break
		}
		
		limits <- sqrt(limits) # narrow search limits
	}
	
	pvals <- numeric(l)
	bestLnL <- .opt(c(U, K), getPvals=TRUE)
	
	minPos <- sapply(windows, head, n=1L)
	maxPos <- sapply(windows, tail, n=1L)
	if (showPlot) {
		limits <- max(abs(omegas))
		plot(NA,
			xlim=c(minPos[1L], maxPos[length(maxPos)]),
			ylim=exp(c(-limits, limits)),
			xlab="Aligned codon position",
			ylab=expression(omega),
			log="y")
		abline(h=1, lty=2, col="gray")
		colors <- colorRampPalette(c("green3", "gray", "black"), bias=4)(101)
		rect(minPos - 1L, 1, maxPos, exp(omegas),
			col=colors[as.integer(pvals*100) + 1L],
			border=NA)
		legend("topright",
			c("Insignificant", "Significant"),
			fill=c("black", "green3"))
	}
	
	rngs <- ifelse(minPos == maxPos,
		minPos,
		paste(minPos, maxPos, sep="-"))
	results <- c(LogLikelihood=bestLnL,
		theta=U,
		kappa=K,
		setNames(exp(omegas), paste("omega", rngs)),
		setNames(pvals, paste("pvalue", rngs)))
	
	if (verbose) {
		cat("\rLnL = ",
			.print(bestLnL),
			sep="")
		time.2 <- Sys.time()
		cat("\n\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(results)
}

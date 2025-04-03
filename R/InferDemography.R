InferDemography <- function(x,
	readingFrame=NA,
	mu=1e-8,
	ploidy=1,
	informationCriterion="BIC",
	showPlot=FALSE,
	verbose=TRUE) {
	
	# error checking
	if (is.numeric(x)) {
		n <- 2L*(length(x) - 1L)
		if (is.na(x)[1L])
			stop("x[1] cannot be NA.")
		if (any(x < 0, na.rm=TRUE))
			stop("x cannot contain negative values.")
		if (n < 4L)
			stop("x must contain at least three values.")
	} else if (is(x, "DNAStringSet") || is(x, "RNAStringSet")) {
		n <- length(x)
		if (n < 4L)
			stop("x must contain at least four sequences.")
		uw <- unique(width(x))
		if (length(uw) != 1L)
			stop("Sequences in x must be the same width (aligned).")
		if (!(readingFrame %in% c(1L:3L, NA_integer_)))
			stop("readingFrame must be 1, 2, 3, or NA.")
	} else {
		stop("x must be a numeric vector, DNAStringSet, or RNAStringSet.")
	}
	if (!is.numeric(mu))
		stop("mu must be a numeric.")
	if (is.na(mu))
		stop("mu cannot be NA.")
	if (mu <= 0)
		stop("mu must be positive.")
	if (mu >= 1)
		stop("mu must be less than one.")
	if (!is.numeric(ploidy))
		stop("ploidy must be a numeric.")
	if (is.na(ploidy))
		stop("ploidy cannot be NA.")
	if (ploidy <= 0)
		stop("ploidy must be positive.")
	if (ploidy != floor(ploidy))
		stop("ploidy must be a whole number.")
	if (length(informationCriterion) != 1L)
		stop("Only one informationCriterion can be specified.")
	if (is.character(informationCriterion)) {
		ICs <- c("AIC", "BIC")
		informationCriterion <- pmatch(informationCriterion, ICs)
		if (is.na(informationCriterion))
			stop("Invalid informationCriterion.")
		if (informationCriterion == -1L)
			stop("Ambiguous informationCriterion.")
		informationCriterion <- -informationCriterion
	} else {
		if (!is.numeric(informationCriterion))
			stop("informationCriterion must be a numeric.")
		if (is.na(informationCriterion))
			stop("informationCriterion cannot be NA.")
		if (informationCriterion < 1L)
			stop("informationCriterion must be at least 1.")
	}
	if (!isTRUEorFALSE(showPlot))
		stop("showPlot must be a logical.")
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be a logical.")
	
	# initialize parameters
	MAX_ATTEMPTS <- 2L
	.print <- function(val)
		formatC(val, digits=3, format="f")
	
	if (verbose)
		time.1 <- Sys.time()
	
	if (is.numeric(x)) {
		L <- sum(x, na.rm=TRUE)
		SFS <- x[-1L] # first value is the number of monomorphic sites
	} else {
		# calculate the folded site frequency spectrum (SFS)
		L <- 0L # number of viable sites
		SFS <- numeric(floor(n/2)) # observed SFS
		if (is.na(readingFrame)) {
			w <- seq_len(uw)
		} else {
			w <- which((seq_len(uw) - readingFrame) %% 3L == 2L) # third codon position
		}
		c <- consensusMatrix(x, baseOnly=TRUE)
		for (i in w) {
			t <- c[, i]
			t <- t[t > 0L]
			if (length(t) <= 2L) {
				L <- L + 1L
				if (length(t) == 2L) { # dimorphic sites
					t <- min(t)
					SFS[t] <- SFS[t] + 1L
				}
			} # otherwise ignore under infinite sites assumption
		}
	}
	S <- sum(SFS, na.rm=TRUE)
	if (informationCriterion == -1L) { # AIC
		penalty <- 2
		informationCriterion <- Inf
	} else if (informationCriterion == -2L) { # BIC
		penalty <- log(length(SFS))
		informationCriterion <- Inf
	} else { # limited time intervals
		penalty <- -Inf
	}
	
	.Est <- function(N, k) {
		E <- numeric(length(SFS))
		for (r in seq_len(n - 1L)) {
			s <- seq_along(N)
			a <- n - k[s] + 1L
			b <- n - k[s + 1L] + 1L
			e <- sum(N*(choose(a, r) - choose(b, r)))
			e <- e*2*ploidy*mu*L/r/choose(n - 1L, r)
			if (r > length(SFS))
				r <- n - r # fold spectrum
			E[r] <- E[r] + e
		}
		E
	}
	
	.LnL <- function(Expect) {
		SFS <- c(L - S, SFS)
		Expect <- c(L - sum(Expect), Expect)
		sum(SFS*log(Expect) - Expect - lfactorial(SFS), na.rm=TRUE)
	}
	
	# Find best parameters
	.opt <- function(Ns, Levels) {
		# enforce min t >= 1 generation
		if (length(Levels) > 0L) {
			minN <- c(Levels*(Levels - 1L), n*(n - 1L))/(2L*ploidy)
		} else {
			minN <- n*(n - 1L)/(2L*ploidy)
		}
		w <- which(Ns < minN)
		if (length(w) > 0L)
			Ns[w] <- minN[w]
		
		c(optim(Ns,
				function(params) {
					if (sum(params < minN) > 0L) {
						-Inf
					} else {
						.LnL(.Est(params, c(2L, Levels, n + 1L)))
					}
				},
				control=list(fnscale=-1)),
			list(Levels=Levels))
	}
	
	Ns <- S/(2*mu*ploidy*L) # initial value
	while (sum(.Est(Ns, c(2L, n + 1L))) > L)
		Ns <- Ns/10
	Levels <- integer()
	
	best <- last <- suppressWarnings(.opt(Ns, Levels))
	if (verbose)
		cat("Intervals = 1: LnL =", .print(last$value))
	
	while (length(Levels) + 1L < informationCriterion) {
		temp <- c(2L, Levels, n + 1L)
		d <- diff(temp)
		w <- which.max(d)
		if (d[w] == 1L)
			break # no more levels possible
		
		# add new level
		Levels <- c((temp[w] + temp[w + 1L]) %/% 2L, Levels)
		Levels <- sort(Levels)
		Ns <- c(Ns[w], Ns)[order(c(w, seq_along(Ns)))]
		
		best <- suppressWarnings(.opt(Ns, Levels))
		if (verbose && interactive())
			cat("\nIntervals = ",
				length(Levels) + 1L,
				": LnL = ",
				.print(best$value),
				sep="")
		
		repeat { # optimize level bounds
			improved <- logical(length(Levels))
			for (j in seq_along(Levels)) { # try moving all levels
				attempts <- 1L
				while (Levels[j] - attempts > 2L &&
					(j == 1L || Levels[j] - attempts > Levels[j - 1L])) {
					temp <- Levels
					temp[j] <- temp[j] - attempts
					test <- suppressWarnings(.opt(Ns, temp))
					if (test$value > best$value) {
						best <- test
						Ns <- best$par
						Levels <- temp
						attempts <- 0L
						improved[j] <- TRUE
						if (verbose && interactive())
							cat("\rIntervals = ",
								length(Levels) + 1L,
								": LnL = ",
								.print(best$value),
								" ",
								sep="")
					} else if (attempts >= MAX_ATTEMPTS) {
						break
					}
					attempts <- attempts + 1L
				}
				
				if (!improved[j]) {
					attempts <- 1L
					while (Levels[j] + attempts <= n &&
						(j == length(Levels) || Levels[j] + attempts < Levels[j + 1L])) {
						temp <- Levels
						temp[j] <- temp[j] + attempts
						test <- suppressWarnings(.opt(Ns, temp))
						if (test$value > best$value) {
							best <- test
							Ns <- best$par
							Levels <- temp
							attempts <- 0L
							improved[j] <- TRUE
							if (verbose && interactive())
								cat("\rIntervals = ",
									length(Levels) + 1L,
									": LnL = ",
									.print(best$value),
									" ",
									sep="")
						} else if (attempts >= MAX_ATTEMPTS) {
							break
						}
						attempts <- attempts + 1L
					}
				}
			}
			if (sum(improved) == 0L)
				break
		}
		if (verbose && !interactive())
			cat("\nIntervals = ",
				length(Levels) + 1L,
				": LnL = ",
				.print(best$value),
				" ",
				sep="")
		
		if (best$value < last$value + penalty) {
			best <- last
			break
		} else {
			last <- best
		}
	}
	
	Expect <- .Est(best$par, c(2L, best$Levels, n + 1L))
	
	Ns <- c(best$par, best$par[length(best$par)])
	Levels <- c(2L, best$Levels, n)
	k <- n
	t <- 0
	Ts <- numeric(length(Ns))
	for (i in rev(seq_along(Ns))) {
		while (k >= Levels[i]) {
			t <- t + 2*ploidy/k/(k - 1L)*Ns[i]
			k <- k - 1L
		}
		Ts[i] <- t
	}
	
	if (showPlot) {
		layout(matrix(1L:2L))
		
		plot(SFS,
			xlab="Minor allele frequency",
			ylab="Count",
			main="Folded site frequency spectrum")
		lines(Expect)
		legend("topright",
			c("Observed", "Fitted"),
			pch=c(1, NA),
			lty=c(NA, 1))
		
		plot(rep(rev(Ts), each=2)[-c(1L, length(Ts)*2L)],
			rep(rev(Ns), each=2)[-1L:-2L],
			xlab="Time (generations)",
			ylab=expression(italic('N'['e'])),
			main="Inferred effective population size",
			type="l",
			log="xy")
	}
	
	results <- c(Intervals=length(Levels) - 1L,
		LogLikelihood=best$value,
		setNames(rev(Ts), paste("Time", rev(Levels))),
		setNames(rev(Ns)[-1L], paste("Ne", rev(Levels[-1L]))),
		setNames(c(L - S, SFS), paste("Observed", 0L:length(SFS))),
		setNames(c(L - sum(Expect), Expect), paste("Estimated", 0L:length(SFS))))
	
	if (verbose) {
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

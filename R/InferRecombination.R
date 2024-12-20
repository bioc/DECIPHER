InferRecombination <- function(x,
	readingFrame=NA,
	position=1:3,
	N=249,
	showPlot=FALSE,
	verbose=TRUE) {
	
	# error checking
	if (is(x, "XStringSet"))
		x <- list(x)
	if (is(x, "XStringSetList"))
		x <- as.list(x)
	if (!is.list(x))
		stop("x must be a list or XStringSet.")
	if (any(sapply(x, function(x) !(is(x, "DNAStringSet") || is(x, "RNAStringSet")))))
		stop("All elements of x must contain a DNAStringSet or RNAStringSet.")
	if (sum(position %in% c(1L:3L, NA_integer_)) != length(position))
		stop("position must be 1, 2, 3, or NA.")
	if (length(position) > 1L && sum(is.na(position)) != 0L)
		stop("position cannot be NA and a number.")
	if (any(duplicated(position)))
		stop("position cannot contain duplicated values.")
	if (!is.numeric(N))
		stop("N must be a numeric.")
	if (length(N) != 1L)
		stop("N must be a single numeric.")
	if (!is.finite(N))
		stop("N must be finite.")
	if (N < 9L)
		stop("N must be at least 9.")
	if (N != floor(N))
		stop("N must be a whole number.")
	if (N %% 3L != 0L)
		stop("N must be evenly divisible by 3.")
	if (!isTRUEorFALSE(showPlot))
		stop("showPlot must be a logical.")
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be a logical.")
	if (length(readingFrame) != 1L && length(readingFrame) != length(x))
		stop("readingFrame is the wrong length.")
	if (sum(readingFrame %in% c(NA_integer_, 1L:3L)) != length(readingFrame))
		stop("readingFrame must be 1, 2, 3, or NA.")
	readingFrame <- as.integer(readingFrame)
	if (length(readingFrame) == 1L)
		readingFrame <- rep(readingFrame, length(x))
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
	}
	
	ans <- .Call("correlationProfile",
		x,
		readingFrame,
		N,
		verbose,
		pBar,
		PACKAGE="DECIPHER")
	P <- ans[[1L]]
	C <- ans[[2L]]
	d <- ans[[3L]][1L]
	n <- ans[[3L]][2L]
	
	SSE <- function(params, pos, showPlot=FALSE, stats=FALSE) {
		Theta_s <- params[1L]
		Phi_s <- params[2L]
		f <- params[3L]
		
		if (is.na(pos)) {
			l <- seq(1L, length(P))
			d_s <- d/n
		} else {
			l <- seq(pos, length(P), 3L)
			d_s <- d/n*3*(sum(P[l])/sum(P)) # relative d_s for the position
		}
		w <- 2/3
		a <- 4/3
		
		if (Theta_s < 0 || Theta_s > d_s)
			return(Inf)
		if (f < 2 || f > 262144)
			return(Inf)
		if (Phi_s < 0 || Phi_s > 1)
			return(Inf)
		
		Theta_p <- (d_s*(1 + f*Phi_s*w + Theta_s*a) - Theta_s)/((1 - d_s*a)*(f*Phi_s*w + Theta_s*a) - d_s*a)
		Phi_p <- Theta_p*Phi_s/Theta_s
		
		c_s0 <- (1 + 2*Theta_s*a)/(1 + 2*Theta_s*a + Phi_s*w*(f + l))
		c_s1 <- (2*Phi_s*w*l)/(1 + 2*Theta_s*a + Phi_s*w*(f + l))
		c_s2 <- (Phi_s*w*(f - l))/(1 + 2*Theta_s*a + Phi_s*w*(f + l))
		
		d_p <- Theta_p/(1 + Theta_p*a)
		Q_p <- 2*((1 + Theta_p*a + Phi_p*w*l)/(1 + 2*Theta_p*a + 2*Phi_p*w*l))*d_p^2
		
		RHS <- c_s0*2*Theta_s/(1 + 2*Theta_s*a)*d_s + c_s1*d_s*d_p + c_s2*Q_p
		
		if (showPlot)
			lines(l, RHS/d_s, col=ifelse(is.na(pos), 1L, pos %% 3L + 5L))
		
		if (!stats && any(RHS <= 0))
			return(Inf)
		
		if (stats) {
			c(fragment=f,
				Theta_sample=Theta_s,
				Phi_sample=Phi_s,
				Theta_pool=Theta_p,
				Phi_pool=Phi_p,
				ratio=Phi_p/Theta_p,
				coverage=Phi_s*w*f/(1 + Theta_s*a + Phi_s*w*f),
				d_pool=d_p,
				d_clonal=Theta_s/(1 + a*Theta_s),
				d_sample=d_s,
				setNames(l, paste("Position", seq_along(l))),
				setNames(P[l]/C[l], paste("Profile", seq_along(l))),
				setNames(RHS, paste("Fitted", seq_along(l))))
		} else {
			sum(P[l]*(RHS/d_s - P[l]/C[l])^2) # weighted SSE
		}
	}
	
	if (showPlot) {
		if (length(position) == 1L && is.na(position)) {
			plot(P/C,
				ylab="P(l)",
				xlab="l (bp)")
		} else {
			l <- which(((seq_along(P) - 1L) %% 3L + 1L) %in% position)
			plot(l,
				P[l]/C[l],
				col=l %% 3L + 5L,
				ylab="P(l)",
				xlab="l (bp)")
			if (length(position) > 1L)
				legend("topright",
					paste("Position", position),
					pch=1,
					lty=1,
					col=position %% 3L + 5L,
					ncol=length(position))
		}
	}
	
	e <- expand.grid(Theta_s=2^(-10L:-30L),
		Phi_s=2^(-20L:0L),
		f=2^(0L:18L))
	
	results <- sapply(position,
		function(POS) {
			o <- numeric(nrow(e))
			for (i in seq_along(o))
				o[i] <- SSE(unname(unlist(e[i,])), pos=POS)
			params <- unname(unlist(e[which.min(o),]))
			o1 <- optim(params, SSE, control=list(maxit=1e4, reltol=1e-16), pos=POS)
			o2 <- optim(c(0.00001, 0.00005, 1000), SSE, control=list(maxit=1e4, reltol=1e-16), pos=POS)
			if (o1$value < o2$value) {
				o <- o1
			} else {
				o <- o2
			}
			SSE(o$par, showPlot=showPlot, stats=TRUE, pos=POS)
		})
	results <- as.matrix(results)
	colnames(results) <- paste("Position", position)
	
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
	
	return(results)
}

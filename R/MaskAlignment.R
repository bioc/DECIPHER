.centerPoint <- function(x, size) {
	
	# returns the center-point moving average of x
	# using a window size of size elements
	
	l <- length(x)
	index <- 1:l
	boundsL <- index - size
	boundsL[1:size] <- 1
	boundsL[(l - size + 1):l] <- boundsL[(l - size + 1):l] - 1:size
	boundsR <- index + size
	boundsR[(l - size + 1):l] <- l
	boundsR[1:size] <- boundsR[1:size] + size:1
	
	avgs <- numeric(length(x))
	for (i in index) {
		avgs[i] <- mean(x[boundsL[i]:boundsR[i]])
	}
	return(avgs)
}

MaskAlignment <- function(myXStringSet,
	type="sequences",
	windowSize=5,
	threshold=1.0,
	maxFractionGaps=0.2,
	includeTerminalGaps=FALSE,
	correction=FALSE,
	randomBackground=FALSE,
	showPlot=FALSE) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet") && !is(myXStringSet, "AAStringSet"))
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	TYPES <- c("sequences", "ranges", "values")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(windowSize))
		stop("windowSize must be a numeric.")
	if (floor(windowSize) != windowSize)
		stop("windowSize must be a whole number.")
	if (windowSize < 1)
		stop("windowSize must be greater than zero.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold < 0)
		stop("threshold must be at least zero.")
	u <- unique(width(myXStringSet))
	if (length(u) > 1)
		stop("myXStringSet must be aligned.")
	if (!is.numeric(maxFractionGaps))
		stop("maxFractionGaps must be a numeric.")
	if (maxFractionGaps > 1 || maxFractionGaps < 0)
		stop("maxFractionGaps must be between 0 and 1 inclusive.")
	if (!is.logical(includeTerminalGaps))
		stop("includeTerminalGaps must be a logical.")
	if (!is.logical(correction))
		stop("correction must be a logical.")
	if (!is.logical(randomBackground))
		stop("randomBackground must be a logical.")
	if (length(randomBackground) != 1L)
		stop("randomBackground must be a single logical.")
	if (!is.logical(showPlot))
		stop("showPlot must be a logical.")
	
	if (is(myXStringSet, "AAStringSet")) {
		pwm <- .Call("consensusProfileAA",
			myXStringSet,
			rep(1, length(myXStringSet)),
			NULL,
			PACKAGE="DECIPHER")
		cm <- pwm[24,]
		if (includeTerminalGaps) {
			w <- which(pwm[27,] > 0)
			if (length(w) > 0)
				cm[w] <- cm[w]*pwm[27, w]
			cm <- cm + 1 - pwm[27,]
		}
		a <- .Call("informationContentAA",
			pwm,
			length(myXStringSet),
			correction,
			randomBackground,
			PACKAGE="DECIPHER")
		if (showPlot)
			MAX <- max(4.321928094887362625798, # log2(20)
				a)
	} else {
		pwm <- .Call("consensusProfile",
			myXStringSet,
			rep(1, length(myXStringSet)),
			NULL,
			PACKAGE="DECIPHER")
		cm <- pwm[5,]
		if (includeTerminalGaps) {
			w <- which(pwm[8,] > 0)
			if (length(w) > 0)
				cm[w] <- cm[w]*pwm[8, w]
			cm <- cm + 1 - pwm[8,]
		}
		a <- .Call("informationContent",
			pwm,
			length(myXStringSet),
			correction,
			randomBackground,
			PACKAGE="DECIPHER")
		if (showPlot)
			MAX <- max(2, # log2(4)
				a)
	}
	
	gaps <- which(cm > maxFractionGaps)
	
	if (windowSize*2 + 1 > length(a) - length(gaps)) {
		a2 <- c <- numeric()
	} else {
		if (length(gaps) > 0) {
			a2 <- a[-gaps]
			c <- .centerPoint(a2, windowSize)
		} else {
			a2 <- a
			c <- .centerPoint(a2, windowSize)
		}
	}
	
	W <- which(c < threshold)
	if (type == 3L) { # "values"
		mask <- logical(length(a))
		if (length(W) > 0) {
			if (length(gaps) > 0) {
				mask[gaps] <- TRUE
				mask[-gaps][W] <- TRUE
			} else {
				mask[W] <- TRUE
			}
		} else {
			mask[gaps] <- TRUE
		}
		result <- data.frame(entropy=a,
			gaps=cm,
			mask=mask)
	} else {
		if (length(W) > 0) {
			if (length(W) > 1) {
				w <- which((W[2:length(W)] - 1) != W[1:(length(W) - 1)])
			} else {
				w <- integer()
			}
			index <- W[c(1, w + 1)]
			
			indicies <- list()
			below_threshold <- as.integer(a2 < threshold)
			rev_below_threshold <- rev(below_threshold)
			for (i in 1:length(index)) {
				rights <- .Call("multiMatch",
					below_threshold,
					1L,
					index[i])
				lefts <- length(below_threshold) - .Call("multiMatch",
					rev_below_threshold,
					1L,
					length(below_threshold) - index[i] + 1L) + 1
				indicies[[i]] <- c(lefts, rights)
			}
			
			W <- c(W, unlist(indicies))
			
			if (length(gaps) > 0) {
				W <- ((1:length(a))[-gaps])[W]
			}
			W <- c(W, gaps)
			W <- unique(sort(W))
			
			if (length(W) > 1) {
				w <- which((W[2:length(W)] - 1) != W[1:(length(W) - 1)])
			} else {
				w <- integer()
			}
			starts <- W[c(1, w + 1)]
			ends <- W[c(w, length(W))]
		} else if (length(gaps) > 0) {
			if (length(gaps) > 1) {
				w <- which((gaps[2:length(gaps)] - 1) != gaps[1:(length(gaps) - 1)])
			} else {
				w <- integer()
			}
			starts <- gaps[c(1, w + 1)]
			ends <- gaps[c(w, length(gaps))]
		} else {
			starts <- NULL
			ends <- NULL
		}
		
		if (type == 1L) {
			if (is(myXStringSet, "DNAStringSet")) {
				result <- DNAMultipleAlignment(myXStringSet,
					rowmask=as(IRanges(), "NormalIRanges"),
					colmask=as(IRanges(starts, ends), "NormalIRanges"))
			} else if (is(myXStringSet, "RNAStringSet")) {
				result <- RNAMultipleAlignment(myXStringSet,
					rowmask=as(IRanges(), "NormalIRanges"),
					colmask=as(IRanges(starts, ends), "NormalIRanges"))
			} else { # AAStringSet
				result <- AAMultipleAlignment(myXStringSet,
					rowmask=as(IRanges(), "NormalIRanges"),
					colmask=as(IRanges(starts, ends), "NormalIRanges"))
			}
		} else if (type == 2L) {
			result <- IRanges(starts, ends)
		}
	}
	
	if (showPlot) {
		if (length(gaps) > 0) {
			x <- (1:length(a))[-gaps]
		} else {
			x <- 1:length(a)
		}
		if (length(gaps) > 0) {
			org_mar <- par()$mar
			par(mar=org_mar + c(0, 0, 0, 2))
		}
		plot(x,
			c,
			type="l",
			xlim=c(1, u),
			ylim=c(0, MAX),
			ylab="Information Content (bits)",
			xlab="Column Position in Alignment")
		w <- which(!(W %in% gaps))
		if (length(w) > 0)
			points(W[w],
				a[W[w]],
				pch=20,
				col="red")
		if (length(W) < length(a)) {
			if (length(W) > 0) {
				x <- (1:length(a))[-W]
				y <- a[-W]
			} else {
				x <- 1:length(a)
				y <- a
			}
			points(x,
				y,
				pch=20,
				col="green")
		}
		if (length(gaps) > 0) {
			abline(h=MAX*maxFractionGaps,
				lty=2,
				col="orange")
			points(gaps,
				MAX*cm[gaps],
				pch=20,
				col="orange")
			axis(4,
				seq(0, MAX, length.out=5),
				labels=c("0", "25", "50", "75", "100"),
				col.axis="orange")
			mtext("Percent Gaps (%)", 4, line=3, col="orange")
			par(mar=org_mar)
		}
		abline(h=threshold,
			lty=3,
			col="red")
		text <- c("Kept Position",
			"Moving Avg.",
			"Masked Pos.",
			"Bit Threshold",
			ifelse(length(gaps) > 0, "Masked Gaps", ""),
			ifelse(length(gaps) > 0, "Gap threshold", ""))
		p <- par("usr")
		legend(0, p[4] + 0.12*(p[4] - p[3]),
			text,
			bg="white",
			lty=c(0, 1, 0, 3, 0, 2),
			pch=c(20, NA, 20, NA, 20, NA),
			col=c("green", "black", "red", "red", rep(ifelse(length(gaps) > 0, "orange", NA), 2)),
			ncol=3,
			bty='n',
			text.width=max(abs(strwidth(text)))*1.2,
			xpd=TRUE)
	}
	
	return(result)
}

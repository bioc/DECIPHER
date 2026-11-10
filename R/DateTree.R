DateTree <- function(x,
	time=1,
	penalty=0,
	epsilon=1e-6) {
	
	# error checking
	if (!is(x, "dendrogram"))
		stop("x must be a dendrogram object.")
	if (length(x) != 2L)
		stop("x must be rooted.")
	if (!is.numeric(time))
		stop("time must be numeric.")
	if (!is.numeric(penalty))
		stop("penalty must be numeric.")
	if (length(penalty) != 1L)
		stop("penalty must be a single value.")
	if (is.na(penalty))
		stop("penalty cannot be NA.")
	if (penalty < 0)
		stop("penalty cannot be negative.")
	if (!is.numeric(epsilon))
		stop("epsilon must be numeric.")
	if (length(epsilon) != 1L)
		stop("epsilon must be a single value.")
	if (is.na(epsilon))
		stop("epsilon cannot be NA.")
	if (epsilon <= 0)
		stop("epsilon must be positive.")
	dist <- attr(x, "height") - rapply(x, attr, which="height")
	if (any(dist < 0))
		stop("Leaves cannot be above the root.")
	index <- unlist(x)
	if (!is.numeric(index))
		stop("Leaf indices of x must be numeric.")
	if (length(index) != attr(x, "members"))
		stop("Number of leaves not equal to members of x.")
	if (max(index) > length(index) || # leaf indices exceed members of x
		min(index) < 1L || # leaf indices of x must be at least 1
		any(duplicated(index))) { # duplicate leaf indices in x
		INDEX <- index # original index
		REINDEX <- TRUE
		count <- 0L
		x <- dendrapply(x,
			function(y) {
				if (is.leaf(y)) {
					count <<- count + 1L
					y[1L] <- count
				}
				y
			})
		index <- seq_len(count)
	} else {
		REINDEX <- FALSE
		index <- as.integer(index)
		index <- order(index)
	}
	labs <- labels(x)
	if (length(time) == 0L) {
		stop("time cannot have zero length.")
	} else if (length(time) == 1L || # only the time of ancestor
		all(abs(time - time[1L]) < epsilon, na.rm=TRUE)) {
		if (is.na(time))
			stop("time cannot be NA.")
		slope <- mean(dist)
		intercept <- 0
		time <- setNames(rep(time, length(dist)), labs)
		rescale <- TRUE
	} else { # sampling time of tips
		if (is.null(names(time)))
			stop("time must be a named numeric vector.")
		time <- setNames(time[labs], labs)
		w <- which(!is.na(time))
		if (length(w) <= 2L)
			stop("Too few times matching labels of x.")
		
		# weight based on proportion of tree length
		# because measurements are dependent
		.getWeights <- function(y, w) {
			for (i in seq_along(y)) {
				l <- attr(y, "height") - attr(y[[i]], "height")
				if (l < epsilon)
					l <- epsilon
				if (is.leaf(y[[i]])) {
					weight[index[y[[i]][1L]]] <<- w + l
				} else {
					.getWeights(y[[i]], w + l/attr(y[[i]], "members"))
				}
			}
			NULL
		}
		weight <- setNames(numeric(length(labs)), labs)
		.getWeights(x, 0)
		# adjust weight by inverse distance from root
		# because mean = variance in Poisson process
		weight <- weight/pmax(dist, epsilon)
		# adjust weight by density of time points
		# because rarer times more informative
		kde <- density(time[w])
		weight <- weight/pmax(approx(kde$x, kde$y, time)$y, epsilon, na.rm=TRUE)
		# exclude missing tips from normalization
		weight <- weight/mean(weight[w])
		
		# regression to predict node root-to-node distance
		mean_times <- mean(weight*time, na.rm=TRUE)
		mean_dist <- mean(weight[w]*dist[w], na.rm=TRUE)
		slope <- sum(weight*(time - mean_times)*(dist - mean_dist), na.rm=TRUE)
		slope <- slope/sum(weight*(time - mean_times)^2, na.rm=TRUE)
		if (slope <= 0)
			stop("Unable to estimate root time.")
		intercept <- mean_dist - slope*mean_times
		if (length(w) < length(time)) # interpolate missing tip times
			time[-w] <- (dist[-w] - intercept)/slope
		rescale <- FALSE
	}
	
	.setRates <- function(y) {
		if (is.leaf(y)) {
			c(unname(time[index[y[1L]]]), 0)
		} else if (length(y) != 2L) {
			stop("x must be strictly bifurcating.")
		} else {
			# calculate balanced length below each branch
			count <<- count + 1L
			c <- count
			
			res1 <- .setRates(y[[1L]])
			res2 <- .setRates(y[[2L]])
			
			# split the results
			t1 <- res1[1L]
			t2 <- res2[1L]
			l1 <- res1[2L]
			l2 <- res2[2L]
			
			# accumulate branch lengths
			h <- attr(y, "height")
			l1 <- l1 + h - attr(y[[1L]], "height")
			l2 <- l2 + h - attr(y[[2L]], "height")
			if (l1 < epsilon)
				l1 <- epsilon
			if (l2 < epsilon)
				l2 <- epsilon
			
			# project root distance to each node
			d1 <- intercept + slope*t1
			d2 <- intercept + slope*t2
			
			# calculate scaled rate
			rate <- (d2 - d1 + sqrt((d1 - d2)^2 + 4*l1*l2))/(2*l2)
			rates[c] <<- rate
			
			if (t1 == t2) {
				c(t1, (l1/rate + l2*rate)/2)
			} else if (t1 < t2) {
				c(t1, l1/rate)
			} else {
				c(t2, l2*rate)
			}
		}
	}
	
	.setTimes <- function(y, rate) {
		if (is.leaf(y)) {
			h <- unname(time[index[y[1L]]])
			attr(y, "height") <- h
			list(h, 0, y)
		} else {
			count <<- count + 1L
			
			rate <- rate*c(rates[count], 1/rates[count])
			h1 <- attr(y[[1L]], "height")
			h2 <- attr(y[[2L]], "height")
			
			res1 <- .setTimes(y[[1L]], rate[1L])
			res2 <- .setTimes(y[[2L]], rate[2L])
			
			# split the results
			t1 <- res1[[1L]]
			t2 <- res2[[1L]]
			l1 <- res1[[2L]]
			l2 <- res2[[2L]]
			y[[1L]] <- res1[[3L]]
			y[[2L]] <- res2[[3L]]
			
			rate <- rate^power
			h <- attr(y, "height")
			if (t1 == t2) {
				l <- (l1 + (h - h1)/rate[1L] + l2 + (h - h2)/rate[2L])/2
				t <- t1
			} else if (t1 < t2) {
				l <- l1 + (h - h1)/rate[1L] # projected length to tip
				t <- t1
			} else {
				l <- l2 + (h - h2)/rate[2L] # projected length to tip
				t <- t2
			}
			if (l < epsilon)
				l <- epsilon
			
			# reverse regression formula to compute node time
			# intercept + slope*h = intercept + slope*t - l
			h <- t - l/slope
			attr(y, "height") <- h
			
			list(t, l, y)
		}
	}
	
	.resetTimes <- function(y, minH) {
		if (!is.leaf(y)) {
			h <- attr(y, "height")
			if (h < minH) {
				attr(y, "height") <- minH
			} else {
				minH <- h
			}
			y[[1L]] <- .resetTimes(y[[1L]], minH)
			y[[2L]] <- .resetTimes(y[[2L]], minH)
			if (rescale)
				attr(y, "height") <- b*attr(y, "height") + a
		}
		y
	}
	
	power <- 1/(1 + penalty) # regularization power
	rates <- numeric(attr(x, "members") - 1L)
	count <- 0L
	.setRates(x)
	
	count <- 0L
	z <- .setTimes(x, 1)[[3L]]
	
	if (rescale) { # scale heights from zero to time
		b <- time[1L]/(time[1L] - attr(z, "height"))
		a <- time[1L] - b*time[1L]
	}
	z <- .resetTimes(z, attr(z, "height"))
	
	if (REINDEX) {
		z <- dendrapply(z,
			function(y) {
				if (is.leaf(y)) {
					y[1L] <- INDEX[y[1L]]
				}
				y
			})
	}
	
	return(z)
}

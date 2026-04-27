# index for rows i of column j within dist object of size n
.index <- function (i, j, n) {
	I <- i
	J <- rep(j, length(i))
	flip <- which(i > j)
	if (length(flip) > 0L) {
		I[flip] <- J[flip]
		J[flip] <- i[flip]
	}
	I1 <- I - 1L
	n*I1 - (I*I1) %/% 2L + J - I
}

Zipline <- function(x,
	type="dendrogram",
	distance="length",
	weights=1,
	bootstraps=0,
	power=0.5,
	verbose=TRUE,
	...) {
	
	# error checking
	if (!is.list(x))
		stop("x must be a list.")
	l <- length(x) # number of trees
	if (l == 0L)
		stop("x cannot be an empty list.")
	if (sum(lengths(x) == 0L) > 0L)
		stop("x cannot contain empty list elements.")
	if (sum(sapply(x, function(x) is(x, "dendrogram"))) != l)
		stop("x can only contain objects of class 'dendrogram'.")
	TYPES <- c("dist", "dendrogram", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(weights))
		stop("weights must be a numeric.")
	if (length(weights) == 1L) {
		if (is.na(weights))
			stop("weights cannot be NA.")
		if (is.infinite(weights))
			stop("weights cannot be Inf.")
		if (weights <= 0)
			stop("weights must be positive.")
		weights <- rep(weights, l)
	} else if (length(weights) == l) {
		if (sum(is.na(weights)) > 0L)
			stop("weights cannot contain NA values.")
		if (sum(is.infinite(weights)) > 0L)
			stop("weights cannot contain Inf values.")
		if (sum(weights <= 0) > 0L)
			stop("weights can only contain positive values.")
	} else {
		stop("weights must be the same length as x.")
	}
	if (!is.numeric(bootstraps))
		stop("bootstraps must be a numeric.")
	if (length(bootstraps) != 1L)
		stop("bootstraps must be a single numeric.")
	if (is.na(bootstraps))
		stop("bootstraps cannot be NA.")
	if (bootstraps < 0L)
		stop("bootstraps must be at least zero.")
	if (floor(bootstraps) != bootstraps)
		stop("bootstraps must be a whole number.")
	if (bootstraps > 0L && type == 1L)
		stop("bootstraps must be zero when type is 'dist'.")
	bootstraps <- as.integer(bootstraps)
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(power))
		stop("power must be a numeric.")
	if (length(power) != 1L)
		stop("power must be a single numeric.")
	if (is.na(power))
		stop("power cannot be NA.")
	if (power <= 0)
		stop("power must be positive.")
	
	if (verbose) {
		if (type > 1L && bootstraps == 0L) {
			cat("Generating summary distance matrix:\n")
			flush.console()
		}
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	}
	
	# determine full set of species labels
	dupes <- logical(l)
	spp <- character()
	for (i in seq_len(l)) {
		labs <- as.character(labels(x[[i]]))
		if (sum(is.na(labs)) > 0L)
			stop("x[[", i, "]] contains NA labels.")
		if (sum(duplicated(labs)) > 0L)
			dupes[i] <- TRUE
		spp <- union(spp, labs)
	}
	spp <- sort(spp) # already unique
	N <- length(spp) # number of species
	L <- (N*(N - 1L)) %/% 2L # number of distances
	
	# compute mean patristic distances
	sums <- numeric(L)
	nums <- integer(L)
	for (i in seq_len(l)) {
		y <- Cophenetic(x[[i]], distance=distance)
		y[y < 0] <- 0 # nullify any negative distances
		m <- match(attr(y, "Labels"), spp)
		n <- length(m)
		if (n == N && sum(m != seq_along(m)) == 0L) {
			sums <- sums + y[]
			nums <- nums + 1L
		} else if (dupes[i]) {
			j <- 1L # column in y
			k <- 0L # starting index in y
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				w <- which(m[(j + 1L):n] != m[j])
				for (p in w) {
					sums[i2[p]] <- sums[i2[p]] + y[i1[p]]
					nums[i2[p]] <- nums[i2[p]] + 1L
				}
				j <- j + 1L
				k <- k + n - j + 1L
			}
		} else {
			j <- 1L # column in y
			k <- 0L # starting index in y
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				sums[i2] <- sums[i2] + y[i1]
				nums[i2] <- nums[i2] + 1L
				j <- j + 1L
				k <- k + n - j + 1L
			}
		}
		
		if (verbose && bootstraps == 0L)
			setTxtProgressBar(pBar, i/l/3)
	}
	rM <- sums/nums
	rm(sums, nums)
	
	# regress patristic distances onto average distance to determine relative rate per tree
	sM <- rM^power
	rM <- rM^(2*power)
	rS <- sum(rM)
	lS <- numeric(l)
	for (i in seq_len(l)) {
		y <- Cophenetic(x[[i]], distance=distance)
		y[y < 0] <- 0 # nullify any negative distances
		m <- match(attr(y, "Labels"), spp)
		n <- length(m)
		y <- y^power
		if (n == N && sum(m != seq_along(m)) == 0L) {
			lS[i] <- sum(sM*y)/rS
		} else if (dupes[i]) {
			j <- 1L # column in y
			k <- 0L # starting index in y
			s <- 0 # sum of distances
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				w <- which(m[(j + 1L):n] != m[j])
				lS[i] <- lS[i] + sum(sM[i2[w]]*y[i1[w]])
				s <- s + sum(rM[i2[w]])
				j <- j + 1L
				k <- k + n - j + 1L
			}
			lS[i] <- lS[i]/s
		} else {
			j <- 1L # column in y
			k <- 0L # starting index in y
			s <- 0 # sum of distances
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				lS[i] <- lS[i] + sum(sM[i2]*y[i1])
				s <- s + sum(rM[i2])
				j <- j + 1L
				k <- k + n - j + 1L
			}
			lS[i] <- lS[i]/s
		}
		
		if (verbose && bootstraps == 0L)
			setTxtProgressBar(pBar, (l + i)/l/3)
	}
	sS <- lS # lS^(power/power)
	lS <- lS^2 # lS^(2*power/power)
	rm(rM, sM)
	
	# regress patristic distances onto relative rates to summarize distance
	nominator <- numeric(L)
	denominator <- numeric(L)
	for (i in seq_len(l)) {
		y <- Cophenetic(x[[i]], distance=distance)
		y[y < 0] <- 0 # nullify any negative distances
		m <- match(attr(y, "Labels"), spp)
		n <- length(m)
		y <- y^power
		y <- weights[i]*sS[i]*y
		lW <- weights[i]*lS[i]
		if (n == N && sum(m != seq_along(m)) == 0L) {
			nominator <- nominator + y[]
			denominator <- denominator + lW
		} else if (dupes[i]) {
			j <- 1L # column in y
			k <- 0L # starting index in y
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				w <- which(m[(j + 1L):n] != m[j])
				for (p in w) {
					nominator[i2[p]] <- nominator[i2[p]] + y[i1[p]]
					denominator[i2[p]] <- denominator[i2[p]] + lW
				}
				j <- j + 1L
				k <- k + n - j + 1L
			}
		} else {
			j <- 1L # column in y
			k <- 0L # starting index in y
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				nominator[i2] <- nominator[i2] + y[i1]
				denominator[i2] <- denominator[i2] + lW
				j <- j + 1L
				k <- k + n - j + 1L
			}
		}
		
		if (verbose && bootstraps == 0L)
			setTxtProgressBar(pBar, (2*l + i)/l/3)
	}
	D <- (nominator/denominator)^(1/power)
	rm(nominator, denominator)
	
	class(D) <- "dist"
	attr(D, "Size") <- N
	attr(D, "Labels") <- spp
	attr(D, "Diag") <- TRUE
	attr(D, "Upper") <- TRUE
	
	if (verbose && bootstraps == 0L) {
		time.2 <- Sys.time()
		close(pBar)
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	if (type == 1L) { # dist
		return(D)
	} else {
		tree <- Treeline(myDistMatrix=D,
			type="dendrogram",
			verbose=verbose && bootstraps == 0L,
			...)
		
		if (bootstraps > 0L) {
			.partitions <- function(x) {
				if (is.leaf(x))
					return(NULL)
				x0 <- paste(sort(match(labels(x), spp)), collapse=" ")
				x1 <- .partitions(x[[1L]])
				x2 <- .partitions(x[[2L]])
				list(x0, x1, x2)
			}
			
			.labelEdges <- function(x) {
				if (!is.leaf(x)) {
					part <- paste(sort(match(labels(x), spp)), collapse=" ")
					attr(x, "edgetext") <- as.character(counts[part])
				}
				x
			}
			
			counts <- unlist(.partitions(tree))
			counts <- setNames(integer(length(counts)), counts)
			for (i in seq_len(bootstraps)) {
				if (verbose)
					setTxtProgressBar(pBar, i/(bootstraps + 1L))
				
				s <- sample(l, replace=TRUE)
				t <- tabulate(s)
				w <- which(t > 0L)
				temp <- Zipline(x[w],
					type="dendrogram",
					distance=distance,
					weights=weights[w]*t[w],
					bootstraps=0L,
					verbose=FALSE,
					...)
				temp <- unlist(.partitions(temp))
				m <- match(temp, names(counts))
				m <- m[!is.na(m)]
				counts[m] <- counts[m] + 1L
			}
			counts <- round(100*counts/bootstraps)
			
			tree <- dendrapply(tree, .labelEdges)
			attr(tree, "edgetext") <- NULL # remove text from (virtual) root branch
			
			if (verbose) {
				time.2 <- Sys.time()
				setTxtProgressBar(pBar, 1)
				close(pBar)
				cat("\n")
				print(round(difftime(time.2,
					time.1,
					units='secs'),
					digits=2))
				cat("\n")
			}
		}
		
		if (type == 2L) { # dendrogram
			return(tree)
		} else { # both
			return(list(D, tree))
		}
	}
}

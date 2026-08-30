.combineOrdered <- function(m, config=0L) {
	v <- vector("list", ncol(m))
	for (i in seq_along(v)) {
		part1 <- sort.int(m[1:2, i])
		part2 <- sort.int(m[3:4, i])
		if (config != 0L) {
			if (config == 1L) { # first NNI
				temp <- part1[2L]
				part1[2L] <- part2[1L]
				part2[1L] <- temp
			} else if (config == 2L) { # second NNI
				temp <- part1[1L]
				part1[1L] <- part2[1L]
				part2[1L] <- temp
			}
			part1 <- sort.int(part1)
			part2 <- sort.int(part2)
		} # else leave in current config
		if (part1[1L] > part2[1L]) {
			v[[i]] <- writeBin(c(part2, part1), raw())
		} else {
			v[[i]] <- writeBin(c(part1, part2), raw())
		}
	}
	v
}

.makeRaw <- function(shift, total) {
	bytes <- raw(total)
	byte <- (shift - 1L) %/% 8L
	bit <- shift - byte*8L - 1L
	byte <- byte + 1L
	bytes[byte] <- as.raw(2^bit)
	bytes
}

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

# collapse duplicate names based on minimum distance
.dereplicate <- function(y) {
	labs <- attr(y, "Labels")
	labs <- sort(unique(labs))
	m <- match(attr(y, "Labels"), labs)
	n <- length(m)
	N <- length(labs)
	alt <- rep(Inf, (N*(N - 1L)) %/% 2L)
	j <- 1L # column in y
	k <- 0L # starting index in y
	while (j < n) {
		i1 <- (k + 1L):(k + n - j) # index in y
		i2 <- .index(m[(j + 1L):n], m[j], N) # index in alt
		w <- which(m[(j + 1L):n] != m[j])
		for (p in w)
			if (y[i1[p]] < alt[i2[p]])
				alt[i2[p]] <- y[i1[p]]
		j <- j + 1L
		k <- k + n - j + 1L
	}
	attr(alt, "Labels") <- labs
	alt
}

Zipline <- function(x,
	type="dendrogram",
	distance="length",
	weights=1,
	bootstraps=0,
	power=0.5,
	concordanceVectors=FALSE,
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
	if (!is.numeric(power))
		stop("power must be a numeric.")
	if (length(power) != 1L)
		stop("power must be a single numeric.")
	if (is.na(power))
		stop("power cannot be NA.")
	if (power <= 0)
		stop("power must be positive.")
	if (length(concordanceVectors) == 0L)
		stop("concordanceVectors must have positive length.")
	if (any(is.na(concordanceVectors)))
		stop("concordanceVectors cannot be NA.")
	if (is.logical(concordanceVectors)) {
		if (concordanceVectors)
			concordanceWeights <- rep(1L, l)
	} else {
		if (is.numeric(concordanceVectors)) {
			if (length(concordanceVectors) != l)
				stop("concordanceVectors must be the same length as x.")
		} else if (is.character(concordanceVectors)) {
			if (length(concordanceVectors) != 1L)
				stop("concordanceVectors must be a single character string.")
		} else {
			stop("concordanceVectors must be a logical.")
		}
		concordanceWeights <- concordanceVectors
		concordanceVectors <- TRUE
	}
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be a logical.")
	
	if (verbose) {
		if (type > 1L && bootstraps == 0L) {
			cat("Generating summary distance matrix:\n")
			flush.console()
		} else if (bootstraps > 0L) {
			cat("Calculating bootstrap replicates:\n")
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
	nums <- numeric(L)
	for (i in seq_len(l)) {
		y <- Cophenetic(x[[i]], distance=distance)
		y[y < 0] <- 0 # nullify any negative distances
		if (dupes[i])
			y <- .dereplicate(y)
		labs <- attr(y, "Labels")
		n <- length(labs)
		m <- match(labs, spp)
		if (weights[i] != 1)
			y <- weights[i]*y
		if (n == N && sum(m != seq_along(m)) == 0L) {
			sums <- sums + y[]
			nums <- nums + weights[i]
		} else {
			j <- 1L # column in y
			k <- 0L # starting index in y
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				sums[i2] <- sums[i2] + y[i1]
				nums[i2] <- nums[i2] + weights[i]
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
		if (dupes[i])
			y <- .dereplicate(y)
		labs <- attr(y, "Labels")
		n <- length(labs)
		m <- match(labs, spp)
		y <- y^power
		if (weights[i] != 1)
			y <- y*weights[i]
		if (n == N && sum(m != seq_along(m)) == 0L) {
			lS[i] <- sum(sM*y)/(weights[i]*rS)
		} else {
			j <- 1L # column in y
			k <- 0L # starting index in y
			s <- 0 # sum of distances
			while (j < n) {
				i1 <- (k + 1L):(k + n - j) # index in y
				i2 <- .index(m[(j + 1L):n], m[j], N) # index in sums
				lS[i] <- lS[i] + sum(sM[i2]*y[i1])
				s <- s + weights[i]*sum(rM[i2])
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
		if (dupes[i])
			y <- .dereplicate(y)
		labs <- attr(y, "Labels")
		n <- length(labs)
		m <- match(labs, spp)
		y <- y^power
		if (weights[i] != 1) {
			y <- weights[i]*sS[i]*y
			lW <- weights[i]*lS[i]
		} else {
			y <- sS[i]*y
			lW <- lS[i]
		}
		if (n == N && sum(m != seq_along(m)) == 0L) {
			nominator <- nominator + y[]
			denominator <- denominator + lW
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
	
	if (type == 1L) # dist only
		return(D)
	tree <- Treeline(myDistMatrix=D,
		type="dendrogram",
		verbose=verbose && bootstraps == 0L,
		...)
	
	if (bootstraps > 0L) {
		.partitions <- function(x) {
			if (is.leaf(x))
				return(NULL)
			part <- spp %in% labels(x)
			s <- 2L*sum(part)
			if (s == N) { # record one partition
				if (part[1L])
					part <- !part
			} else if (s > N) { # record smaller partition
				part <- !part
			}
			part <- which(part)
			if (length(part) == 1L) {
				part <- NULL
			} else {
				part <- paste(part, collapse=" ")
			}
			x1 <- .partitions(x[[1L]])
			x2 <- .partitions(x[[2L]])
			list(part, x1, x2)
		}
		counts <- unique(unlist(.partitions(tree)))
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
			if (N != attr(temp, "members"))
				next # different label sets
			temp <- unlist(.partitions(temp))
			m <- match(temp, names(counts))
			m <- m[!is.na(m)]
			counts[m] <- counts[m] + 1L
		}
		counts <- round(counts/bootstraps,
			ceiling(log10(bootstraps)))
		
		.labelEdges <- function(x) {
			if (!is.leaf(x)) {
				part <- spp %in% labels(x)
				s <- 2L*sum(part)
				if (s == N) { # record one partition
					if (part[1L])
						part <- !part
				} else if (s > N) { # record smaller partition
					part <- !part
				}
				part <- which(part)
				if (length(part) > 1L) {
					part <- paste(part, collapse=" ")
					attr(x, "edgetext") <- as.character(counts[part])
				}
			}
			x
		}
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
	
	if (concordanceVectors) { # annotate internal branches with concordance vectors
		if (verbose) {
			cat("Computing concordance vectors:\n")
			flush.console()
			pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
			time.1 <- Sys.time()
		}
		
		# convert clade membership into a bit vector
		.recordParts <- function(y) {
			count <<- count + 1L
			c <- count
			if (is.leaf(y)) {
				m <- match(attr(y, "label"), spp)
				parts[[c]] <<- .makeRaw(m, total)
			} else {
				for (i in seq_along(y)) {
					m <- .recordParts(y[[i]])
					if (length(parts[[c]]) > 0L) {
						parts[[c]] <<- parts[[c]] |
							parts[[m]]
					} else {
						parts[[c]] <<- parts[[m]]
					}
				}
			}
			c
		}
		
		# record an attribute at each internal node
		.recordAttr <- function(y, root=FALSE) {
			count <<- count + 1L
			c <- count
			if (root) {
				if (length(y) == 1L) { # descend
					return(.recordAttr(y[[1L]], TRUE))
				} else {
					return(c(c, sapply(y, .recordAttr)))
				}
			}
			if (!is.leaf(y)) {
				branchWeight[count] <<- sum(as.numeric(attr(y, concordanceWeights)), na.rm=TRUE)
				sapply(y, .recordAttr)
			}
			c
		}
		
		# record index arrangment of quartet
		.arrangement <- function(y, root=FALSE) {
			count <<- count + 1L
			c <- count
			if (root) {
				if (length(y) == 1L) { # descend
					return(.arrangement(y[[1L]], TRUE))
				} else {
					return(c(c, sapply(y, .arrangement)))
				}
			}
			if (!is.leaf(y)) {
				if (length(y) == 2L) { # bifurcation
					counts <- integer(2L)
					for (i in 1:2) {
						counts[i] <- count + 1L
						m <- .arrangement(y[[i]])
						arrange[i, c] <<- m
					}
					arrange[3L, counts] <<- counts[2:1]
					arrange[4L, counts] <<- c
				} else {
					for (i in seq_along(y))
						.arrangement(y[[i]])
				}
			}
			c
		}
		
		# add attribute based on values in vectors
		.applyVectors <- function(y) {
			count <<- count + 1L
			if (!is.leaf(y)) {
				m <- match(count, internal, nomatch=0L)
				if (m != 0L)
					attr(y, "vector") <- vectors[, m]
				y[] <- lapply(y, .applyVectors)
			}
			y
		}
		
		# modify arrange to account for (virtual) root node
		.handleRoot <- function(counts) {
			if (length(counts) == 3L) { # bifurcation
				arrange[3L:4L, counts[2L:3L]] <<- NA_integer_
				arrange[, counts[1L]] <<- c(arrange[1L:2L, counts[2L]],
					arrange[1L:2L, counts[3L]])
			} else if (length(counts) == 4L) { # trifurcation
				arrange[3L:4L, counts[2L]] <<- counts[c(3L, 4L)]
				arrange[3L:4L, counts[3L]] <<- counts[c(2L, 4L)]
				arrange[3L:4L, counts[4L]] <<- counts[c(2L, 3L)]
			} # else unresolved
		}
		
		total <- N %/% 8L + 1L # number of bytes needed
		
		# record membership in partitions for summary tree
		parts <- vector("list", N*2L - 1L)
		count <- 0L
		.recordParts(tree) # modifies `parts`
		master <- parts
		length(master) <- count
		master <- c(master,
			lapply(master,
				function(x)
					xor(x, master[[1L]])))
		
		# record the arrangement of clades for summary tree
		arrange <- matrix(NA_integer_,
			4L,
			count)
		count <- 0L
		counts <- .arrangement(tree, TRUE) # modifies `arrange`
		arrange[4L,] <- arrange[4L,] + count
		.handleRoot(counts) # modifies `arrange`
		arrangement <- arrange
		internal <- which(colSums(arrangement) > 0L)
		arrangement <- arrangement[, internal, drop=FALSE]
		
		# initialize concordance vectors
		vectors <- matrix(0L,
			7L,
			length(internal),
			dimnames=list(c("CFsplit", "DF1split", "DF2split", "CFquartet", "DF1quartet", "DF2quartet", "decisive"),
				NULL))
		for (i in seq_len(l)) { # for each input tree
			# record membership in input tree partitions
			n <- attr(x[[i]], "members")
			parts <- vector("list", n*2L - 1L)
			count <- 0L
			.recordParts(x[[i]]) # modifies `parts`
			length(parts) <- count
			parts <- c(parts,
				lapply(parts,
					function(x)
						xor(x, parts[[1L]])))
			
			# determine partitions shared with master
			temp <- lapply(master,
				function(y) y & parts[[1L]])
			
			# record possible quartet arrangements in summary tree
			arrange <- arrangement
			arrange[] <- selfmatch(temp)[arrange]
			CFs <- .combineOrdered(arrange)
			DF1s <- .combineOrdered(arrange, config=1L)
			DF2s <- .combineOrdered(arrange, config=2L)
			
			# record nodes where all clades have representation
			present <- !temp %in% list(raw(total)) # decisive partition
			present <- matrix(present[arrangement],
				nrow(arrangement),
				ncol(arrangement))
			present <- colSums(present)
			w <- which(present == nrow(arrangement))
			if (is.character(concordanceWeights)) { # record weights at each branch
				branchWeight <- integer(count)
				count <- 0L
				counts <- .recordAttr(x[[i]], TRUE) # modifies `branchWeight`
				if (length(counts) == 3L) # root bifurcation
					branchWeight[counts[1L]] <- max(branchWeight[counts], na.rm=TRUE)
				branchWeight <- rep(branchWeight, 2L)
				# must use weight of 1 since no branch correspondence necessary
				vectors["decisive", w] <- vectors["decisive", w] + 1L
			} else {
				branchWeight <- rep(concordanceWeights[i], length(parts))
				vectors["decisive", w] <- vectors["decisive", w] + concordanceWeights[i]
			}
			
			# record possible partitions in summary tree
			part1 <- part2 <- part3 <- vector("list", ncol(arrange))
			for (j in w) {
				part1[[j]] <- temp[[arrange[1L, j]]] | temp[[arrange[2L, j]]] # original
				part2[[j]] <- temp[[arrange[1L, j]]] | temp[[arrange[3L, j]]] # first NNI
				part3[[j]] <- temp[[arrange[3L, j]]] | temp[[arrange[2L, j]]] # second NNI
			}
			m <- match(part1, parts, nomatch=0L)
			w <- which(m != 0L)
			vectors["CFsplit", w] <- vectors["CFsplit", w] + branchWeight[m[w]]
			m <- match(part2, parts, nomatch=0L)
			w <- which(m != 0L)
			vectors["DF1split", w] <- vectors["DF1split", w] + branchWeight[m[w]]
			m <- match(part3, parts, nomatch=0L)
			w <- which(m != 0L)
			vectors["DF2split", w] <- vectors["DF2split", w] + branchWeight[m[w]]
			
			# record clade arrangement for input tree quartets
			arrange <- matrix(NA_integer_,
				4L,
				length(parts))
			count <- 0L
			counts <- .arrangement(x[[i]], TRUE) # modifies `arrange`
			arrange[4L,] <- arrange[4L,] + count
			.handleRoot(counts) # modifies `arrange`
			
			# re-index to match master partitions
			present <- match(parts, temp, nomatch=0L)
			present <- matrix(present[arrange],
				nrow(arrange),
				ncol(arrange))
			
			# compute concordance vector
			w <- which(colSums(present != 0L) == nrow(present))
			temp <- .combineOrdered(present[, w, drop=FALSE])
			branchWeight <- branchWeight[w]
			m <- match(CFs, temp, nomatch=0L)
			w <- which(m != 0L)
			vectors["CFquartet", w] <- vectors["CFquartet", w] + branchWeight[m[w]]
			m <- match(DF1s, temp, nomatch=0L)
			w <- which(m != 0L)
			vectors["DF1quartet", w] <- vectors["DF1quartet", w] + branchWeight[m[w]]
			m <- match(DF2s, temp, nomatch=0L)
			w <- which(m != 0L)
			vectors["DF2quartet", w] <- vectors["DF2quartet", w] + branchWeight[m[w]]
			
			if (verbose)
				setTxtProgressBar(pBar, i/l)
		}
		
		# annotate tree with concordance vectors
		count <- 0L
		tree <- .applyVectors(tree)
		
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
	
	if (type == 2L) { # dendrogram only
		return(tree)
	} else { # both
		return(list(D, tree))
	}
}

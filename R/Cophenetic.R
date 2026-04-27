Cophenetic <- function(x,
	distance="length") {
	
	# error checking
	if (!is(x, "dendrogram"))
		stop("x must be an object of class 'dendrogram'.")
	l <- length(distance)
	if (l == 0L)
		stop("distance must have at least one element.")
	if (!is.character(distance))
		stop("distance must be a character string.")
	if (sum(is.na(distance)) > 0L)
		stop("distance cannot contain NA values.")
	
	DIST <- c("length", "edges")
	DIST <- match(distance, DIST)
	w <- which(is.na(DIST))
	if (l - length(w) > 1L) # more than one non-NA value
		stop("distance may contain one of either 'length' or 'edges'.")
	if (length(w) > 0L) {
		DIST[w] <- 3L # other distance
		if (l > 1L) {
			if (length(w) == l)
				stop("distance must contain either 'length' or 'edges'.")
			DIST <- sort(DIST) # place edge attributes last
		}
	}
	
	n <- attr(x, "members")
	if (n == 1L) {
		d <- dist(numeric())
		attr(d, "Diag") <- TRUE
		attr(d, "Upper") <- TRUE
		attr(d, "Labels") <- attr(x, "label")
		return(d)
	}
	
	u <- unlist(x)
	o <- order(u)
	u <- u[o]
	labs <- labels(x)
	labs <- labs[o]
	if (sum(seq_along(u) != u) != 0L) { # need to renumber
		x <- rapply(x,
			function(y) {
				y[] <- match(y[1L], u)
				y
			},
			how="replace")
	}
	
	index <- n - 1L
	C <- matrix(0L, index, 2L)
	H <- matrix(0, index, 2L)
	
	# convert dendrogram to matrix
	indices <- integer(index)
	stack <- vector("list", index)
	pos <- 1L
	stack[[pos]] <- x
	indices[pos] <- index
	while (pos > 0L) {
		y <- stack[[pos]]
		i <- indices[pos]
		pos <- pos - 1L # remove
		
		h <- attr(y, "height")
		while (length(y) == 1L)
			y <- y[[1L]] # descend
		
		if (DIST[1L] == 1L) { # length
			H[i, 1L] <- abs(h - attr(y[[1L]], "height"))
		} else {
			H[i, 1L] <- 1
		}
		for (j in seq_along(w)) {
			a <- attr(y[[1L]], distance[w[j]])
			if (length(a) == 0L) {
				H[i, 1L] <- 0
			} else {
				H[i, 1L] <- H[i, 1L]*sum(as.numeric(a), na.rm=TRUE)
			}
		}
		
		if (is.leaf(y[[1L]])) {
			C[i, 1L] <- -y[[1L]][1L]
		} else {
			index <- index - 1L
			C[i, 1L] <- index
			pos <- pos + 1L
			stack[pos] <- y[1L] # add
			indices[pos] <- index
		}
		if (length(y) == 2L) {
			if (DIST[1L] == 1L) { # length
				H[i, 2L] <- abs(h - attr(y[[2L]], "height"))
			} else {
				H[i, 2L] <- 1
			}
			for (j in seq_along(w)) {
				a <- attr(y[[2L]], distance[w[j]])
				if (length(a) == 0L) {
					H[i, 2L] <- 0
				} else {
					H[i, 2L] <- H[i, 2L]*sum(as.numeric(a), na.rm=TRUE)
				}
			}
			if (is.leaf(y[[2L]])) {
				C[i, 2L] <- -y[[2L]][1L]
			} else {
				index <- index - 1L
				C[i, 2L] <- index
				pos <- pos + 1L
				stack[pos] <- y[2L] # add
				indices[pos] <- index
			}
		} else { # length(y) > 2
			index <- index - 1L
			C[i, 2L] <- index
			pos <- pos + 1L
			y <- y[-1L]
			attr(y, "height") <- h
			stack[[pos]] <- y # add
			indices[pos] <- index
		}
	}
	if ((DIST[1L] == 2L) && # edges
		length(x) == 2L) # root is one edge
		H[nrow(H),] <- H[nrow(H),]/2
	
	d <- .Call("patristic", C, H, 1, PACKAGE="DECIPHER")
	
	class(d) <- "dist"
	attr(d, "Size") <- n
	attr(d, "Diag") <- TRUE
	attr(d, "Upper") <- TRUE
	attr(d, "Labels") <- labs
	
	return(d)
}

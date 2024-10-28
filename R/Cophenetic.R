Cophenetic <- function(x) {
	
	# error checking
	if (!is(x, "dendrogram"))
		stop("x must be an object of class 'dendrogram'.")
	
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
	labs <- rapply(x,
		function(x)
			attr(x, "label"))
	labs <- labs[o]
	x <- rapply(x,
		function(y) {
			y[] <- match(y[1L], u)
			y
		},
		how="replace")
	
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
			y <- stack[[1L]] # descend
		
		h1 <- attr(y[[1L]], "height")
		H[i, 1L] <- h - h1
		
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
			h2 <- attr(y[[2L]], "height")
			H[i, 2L] <- h - h2
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
	
	d <- .Call("patristic", C, H, 1, PACKAGE="DECIPHER")
	
	class(d) <- "dist"
	attr(d, "Size") <- n
	attr(d, "Diag") <- TRUE
	attr(d, "Upper") <- TRUE
	attr(d, "Labels") <- labs
	
	return(d)
}

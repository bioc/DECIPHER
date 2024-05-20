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
	d <- numeric(n*(n - 1)/2)
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
	
	.dist <- function(dend) {
		# initialize a stack of maximum length (n)
		stack <- vector("list", n)
		visit <- logical(n) # node already visited
		parent <- integer(n) # index of parent node
		index <- integer(n) # index in parent node
		pos <- 1L # current position in the stack
		stack[[pos]] <- dend
		while (pos > 0L) { # more nodes to visit
			if (visit[pos]) { # ascending tree
				visit[pos] <- FALSE # reset visit
				
				for (k in seq_along(stack[[pos]])) {
					h <- attr(stack[[pos]], "height") - attr(stack[[pos]][[k]], "height")
					I <- unlist(stack[[pos]][[k]])
					d <- .Call("cophenetic", # in-place change of d (requires previous temporary copy)
						I,
						n,
						d,
						h,
						PACKAGE="DECIPHER")
				}
				
				# replace self in parent
				if (parent[pos] > 0)
					stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
				pos <- pos - 1L # pop off of stack
			} else { # descending tree
				visit[pos] <- TRUE
				p <- pos
				for (i in seq_along(stack[[p]])) {
					if (!is.leaf(stack[[p]][[i]])) {
						# push subtree onto stack
						pos <- pos + 1L
						stack[[pos]] <- stack[[p]][[i]]
						parent[[pos]] <- p
						index[[pos]] <- i
					}
				}
			}
		}
		return(stack[[1L]])
	}
	.dist(x)
	class(d) <- "dist"
	attr(d, "Size") <- n
	attr(d, "Diag") <- TRUE
	attr(d, "Upper") <- TRUE
	attr(d, "Labels") <- labs
	
	return(d)
}

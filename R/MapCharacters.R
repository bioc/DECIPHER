MapCharacters <- function(x,
	refPositions=seq_len(nchar(attr(x, "state")[1])),
	labelEdges=FALSE,
	type="dendrogram",
	chars=LETTERS,
	ignoreIndels=TRUE) {
	
	# error checking
	if (!is(x, "dendrogram"))
		stop("x is not a dendrogram.")
	if (!is.logical(labelEdges))
		stop("labelEdges must be a logical.")
	if (!is.numeric(refPositions))
		stop("refPositions must be a numeric.")
	if (any(refPositions > nchar(attr(x, "state")[1])) ||
		any(refPositions < 1))
		stop("refPositions out of bounds.")
	if (is.unsorted(refPositions))
		stop("refPositions must be in ascending order.")
	if (length(chars) <= 1)
		stop("At least two chars must be specified.")
	if (any(is.na(chars)))
		stop("Chars cannot contain NA values.")
	if (any(nchar(chars) != 1))
		stop("All chars must be a single character.")
	gaps <- "-"
	if (any(chars %in% gaps))
		stop("Gap ('-') characters cannot be in chars.")
	if (!is.logical(ignoreIndels))
		stop("ignoreIndels must be a logical.")
	TYPES <- c("dendrogram", "table", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	
	.map <- function(x) {
		# initialize a stack of maximum length (l)
		stack <- vector("list", l)
		visit <- logical(l) # node already visited
		parent <- integer(l) # index of parent node
		index <- integer(l) # index in parent node
		pos <- 1L # current position in the stack
		stack[[pos]] <- x
		while (pos > 0L) { # more nodes to visit
			if (visit[pos]) { # ascending tree
				visit[pos] <- FALSE # reset visit
				
				for (i in seq_along(stack[[pos]])) {
					s1 <- attr(stack[[pos]], "state") # root ancestral state
					if (is.null(s1))
						stop("Missing 'state' attribute in x.")
					s1 <- s1[1]
					if (!is.character(s1))
						stop("Incorrect 'state' attribute in x.")
					s1 <- strsplit(s1, "", fixed=TRUE)[[1]]
					
					s2 <- attr(stack[[pos]][[i]], "state")
					if (is.null(s2))
						stop("Missing 'state' attribute in x.")
					s2 <- s2[1]
					if (!is.character(s2))
						stop("Incorrect 'state' attribute in x.")
					s2 <- strsplit(s2, "", fixed=TRUE)[[1]]
					if (length(s1) != length(s2))
						stop("Inconsistent 'state' attributes in x.")
					
					if (ignoreIndels) {
						nongaps <- which(!(s1 %in% gaps) &
							!(s2 %in% gaps))
					} else {
						nongaps <- which(!(s1 %in% gaps) |
							!(s2 %in% gaps))
					}
					w <- which(s1[nongaps] != s2[nongaps] &
						s1[nongaps] %in% chars &
						s2[nongaps] %in% chars)
					
					if (length(w) > 0) {
						ini <- s1[nongaps[w]]
						fin <- s2[nongaps[w]]
						
						ins <- which(s1[nongaps[w]] == "-")
						if (length(ins) > 0) {
							counts <- rep(1L, length(ins))
							remove <- logical(length(ins))
							start <- 1L
							for (j in seq_along(ins)[-1]) {
								if (w[ins[j]] == w[ins[j - 1]] + 1L) {
									counts[start] <- counts[start] + 1L
									remove[j] <- TRUE
								} else {
									start <- j
								}
							}
							ini[ins] <- "ins("
							fin[ins] <- paste(",",
								sapply(mapply(seq,
									ins,
									ins + counts - 1L,
									SIMPLIFY=FALSE),
									function(s) {
										paste(s2[nongaps[w[s]]],
										collapse="")
									}),
								")",
								sep="")
							remove <- which(remove)
							if (length(remove) > 0) {
								w <- w[-ins[remove]]
								ini <- ini[-ins[remove]]
								fin <- fin[-ins[remove]]
							}
						}
						
						del <- which(s2[nongaps[w]] == "-")
						if (length(del) > 0) {
							counts <- rep(1L, length(del))
							remove <- logical(length(del))
							start <- 1L
							for (j in seq_along(del)[-1]) {
								if (w[del[j]] == w[del[j - 1]] + 1L) {
									counts[start] <- counts[start] + 1L
									remove[j] <- TRUE
								} else {
									start <- j
								}
							}
							ini[del] <- "del("
							fin[del] <- paste(",",
								sapply(mapply(seq,
									del,
									del + counts - 1L,
									SIMPLIFY=FALSE),
									function(s) {
										paste(s1[nongaps[w[s]]],
										collapse="")
									}),
								")",
								sep="")
							remove <- which(remove)
							if (length(remove) > 0) {
								w <- w[-del[remove]]
								ini <- ini[-del[remove]]
								fin <- fin[-del[remove]]
							}
						}
						
						w <- match(nongaps[w], refPositions)
						keep <- which(!is.na(w))
						if (length(keep) > 0) {
							attr(stack[[pos]][[i]], "change") <- paste(ini[keep],
								w[keep],
								fin[keep],
								sep="")
							if (type > 1L) {
								res[[num]] <<- attr(stack[[pos]][[i]], "change")
								num <<- num + 1L
							}
						} else {
							attr(stack[[pos]][[i]], "change") <- character()
						}
						if (labelEdges)
							attr(stack[[pos]][[i]], "edgetext") <- length(keep)
					} else {
						attr(stack[[pos]][[i]], "change") <- character()
						if (labelEdges)
							attr(stack[[pos]][[i]], "edgetext") <- 0L
					}
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
	
	l <- attr(x, "members")
	if (type > 1L) {
		res <- vector("list",
			2*l - 2)
		num <- 1L
	}
	out <- .map(x)
	attr(out, "change") <- character() # empty set at root node
	if(type > 1L) {
		res <- table(unlist(res))
		res <- sort(res, decreasing=TRUE)
		if (type == 2L) { # table
			out <- res
		} else if (type == 3L) { # both
			out <- list(out, res)
		}
	}
	
	invisible(out)
}

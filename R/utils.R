.subset <- function(x, i) {
	ans <- new(class(x))
	g <- x@ranges@group[i]
	u <- unique(g)
	ans@pool <- x@pool[u]
	ans@ranges@group <- match(g, u)
	ans@ranges@start <- x@ranges@start[i]
	ans@ranges@width <- x@ranges@width[i]
	ans
}

.extract <- function(x, i, s, e) {
	ans <- new(class(x))
	g <- x@ranges@group[i]
	u <- unique(g)
	ans@pool <- x@pool[u]
	ans@ranges@group <- match(g, u)
	ans@ranges@start <- as.integer(x@ranges@start[i] + s - 1L)
	ans@ranges@width <- as.integer(e - s + 1L)
	ans
}

.extractSet <- function(x, s, e) {
	ans <- new(paste(class(x), "Set", sep=""))
	ans@pool@.link_to_cached_object_list <- list(x@shared@.link_to_cached_object)
	ans@pool@xp_list <- list(x@shared@xp)
	ans@ranges@group <- rep(1L, length(s))
	ans@ranges@start <- as.integer(s + x@offset)
	ans@ranges@width <- as.integer(e - s + 1L)
	ans
}

.replace <- function(x, y, index) {
	ans <- new(class(x))
	ans@pool <- c(x@pool,
		y@pool)
	g <- .Call("getPools",
		ans@pool@xp_list,
		PACKAGE="DECIPHER")
	w <- which(!duplicated(g))
	ans@pool <- ans@pool[w]
	m <- match(g, g[w])
	ans@ranges@group <- x@ranges@group
	ans@ranges@group[index] <- y@ranges@group + length(x@pool)
	g <- m[ans@ranges@group]
	u <- unique(g)
	ans@pool <- ans@pool[u]
	ans@ranges@group <- match(g, u)
	ans@ranges@start <- x@ranges@start
	ans@ranges@start[index] <- y@ranges@start
	ans@ranges@width <- x@ranges@width
	ans@ranges@width[index] <- y@ranges@width
	ans
}

.append <- function(x, y) {
	ans <- new(class(x))
	ans@pool <- c(x@pool,
		y@pool)
	ans@ranges@group <- c(x@ranges@group,
		y@ranges@group + length(x@pool))
	ans@ranges@start <- c(x@ranges@start,
		y@ranges@start)
	ans@ranges@width <- c(x@ranges@width,
		y@ranges@width)
	ans
}

.switch <- function(x) {
	if (is(x, "DNAStringSet")) {
		ans <- new("RNAStringSet")
	} else {
		ans <- new("DNAStringSet")
	}
	ans@pool <- x@pool
	ans@ranges@group <- x@ranges@group
	ans@ranges@start <- x@ranges@start
	ans@ranges@width <- x@ranges@width
	ans
}

.detectCores <- function()
	.Call("detectCores", PACKAGE="DECIPHER")

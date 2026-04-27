WriteDendrogram <- function(x,
	file="",
	quote="'",
	space=" ",
	internalLabels=TRUE,
	digits=10,
	append=FALSE,
	unroot=FALSE) {
	
	# error checking
	if (!isTRUEorFALSE(unroot))
		stop("unroot must be a logical.")
	if (!is(x, "dendrogram")) {
		if (is.list(x)) {
			if (length(x) == 0L)
				stop("x cannot be empty.")
			for (i in seq_along(x))
				if (!is(x[[i]], "dendrogram"))
					stop("x[[", i, "]] must be a dendrogram.")
			for (i in seq_along(x))
				WriteDendrogram(x[[i]],
					file=file,
					quote=quote,
					space=space,
					internalLabels=internalLabels,
					digits=digits,
					unroot=unroot,
					append=ifelse(i == 1L, append, TRUE))
			return(invisible(NULL))
		} else {
			stop("x must be a dendrogram.")
		}
	}
	if (!is.character(quote))
		stop("quote must be a character.")
	if (nchar(quote) > 1L)
		stop("quote must be a single character.")
	if (!is.character(space))
		stop("space must be a character.")
	if (nchar(space) != 1L)
		stop("space must be a single character.")
	if (!is.numeric(digits))
		stop("digits must be a numeric.")
	if (floor(digits) != digits)
		stop("digits must be a whole number.")
	if (digits < 1)
		stop("digits must be at least 1.")
	if (!isTRUEorFALSE(append))
		stop("append must be a logical.")
	if (!isTRUEorFALSE(internalLabels))
		stop("internalLabels must be a logical.")
	
	if (unroot) {
		if (attr(x, "members") < 3L) {
			warning("x contains too few members to unroot.")
		} else {
			if (length(x) == 2L && # rooted
				(length(x[[1L]]) != 1L || length(x[[2L]]) != 1L)) {
				if (length(x[[1L]]) == 1L) {
					poly <- 2L
				} else {
					poly <- 1L
				}
				deltaH <- attr(x, "height") - attr(x[[poly]], "height")
				attr(x, "height") <- attr(x[[poly]], "height")
				# add height to other side
				x[[3L - poly]] <- dendrapply(x[[3L - poly]],
					function(x) {
						attr(x, "height") <- attr(x, "height") - deltaH
						x
					})
				while (length(x[[poly]]) > 1L) {
					x[[length(x) + 1L]] <- x[[poly]][[length(x[[poly]])]]
					length(x[[poly]]) <- length(x[[poly]]) - 1L
				}
				x[[poly]] <- x[[poly]][[1L]]
			} else if (length(x) <= 2) {
				warning("x could not be unrooted.")
			}
		}
	}
	
	if (is.character(file)) {
		if (file == "") {
			file <- stdout()
		} else if (substring(file, 1L, 1L) == "|") {
			file <- pipe(substring(file, 2L), ifelse(append, "a", "w"))
			on.exit(close(file))
		} else {
			file <- file(file, ifelse(append, "a", "w"))
			on.exit(close(file))
		}
	}
	
	getLab <- function(LAB) {
		if (is.null(LAB))
			return("")
		if (space != " ")
			LAB <- gsub(" ", space, LAB, fixed=TRUE)
		if (quote != "") {
			LAB <- gsub(quote, '_', LAB, fixed=TRUE)
			LAB <- paste(quote, LAB, quote, sep="")
		}
		return(LAB)
	}
	
	.dendrogram2newick <- function(x, height=attr(x, "height"), root=TRUE) {
		if (is.leaf(x)) {
			cat(getLab(attr(x, "label")),
				":",
				formatC(height - attr(x, "height"),
					digits=digits,
					width=1,
					format="fg"),
				sep="",
				file=file,
				append=TRUE)
		} else {
			cat("(",
				file=file,
				append=TRUE)
			for (i in seq_along(x)) {
				.dendrogram2newick(x[[i]],
					attr(x, "height"),
					root=FALSE)
				if (i < length(x))
					cat(",",
						file=file,
						append=TRUE)
			}
			if (root) {
				cat(");\n",
					file=file,
					append=TRUE)
			} else {
				if (internalLabels) {
					edgetext <- attr(x, "edgetext")
					if (!is.null(edgetext))
						if (!is.numeric(edgetext))
							edgetext <- getLab(edgetext)
				} else {
					edgetext <- NULL
				}
				cat(")",
					edgetext,
					":",
					formatC(height - attr(x, "height"),
						digits=digits,
						width=1,
						format="fg"),
					sep="",
					file=file,
					append=TRUE)
			}
		}
	}
	
	.dendrogram2newick(x)
	invisible(NULL)
}

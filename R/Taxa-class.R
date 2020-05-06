.spoke <- function(x,
	labels=names(x),
	edges=200,
	radius=0.6,
	angle=45,
	init.angle=0,
	final.angle=270,
	col,
	cex=1,
	circle=TRUE,
	font=rep(1, length(x))) {
	labels <- as.graphicsAnnot(labels)
	x <- c(0, cumsum(x)/sum(x)*final.angle/360)
	dx <- diff(x)
	nx <- length(dx)
	p <- par(mar=c(0, 0, 0, 0))
	xlim <- ylim <- c(-1, 1)
	plot(NA,
		type="n",
		xlim=xlim,
		ylim=ylim,
		bty="n",
		xaxt="n",
		yaxt="n",
		xaxs="i",
		yaxs="i",
		main="",
		xlab="",
		ylab="",
		asp=1)
	dev.hold()
	on.exit(dev.flush())
	usr <- par("usr")
	xlim <- usr[1:2]
	ylim <- usr[3:4]
	deltax <- diff(xlim)/2
	deltay <- diff(ylim)/2
	angle <- rep(angle, nx)
	twopi <- -2*pi
	t2xy <- function(t) {
		t2p <- twopi*t + init.angle*pi/180
		list(x=radius*cos(t2p),
			y=radius*sin(t2p),
			an=t2p)
	}
	for (i in seq_len(nx)) {
		n <- max(2, floor(edges*dx[i]))
		P <- t2xy(seq.int(x[i], x[i + 1], length.out=n))
		if (circle) {
			polygon(c(P$x, 0),
				c(P$y, 0),
				angle=angle[i], 
				border=NULL,
				col=col[i])
		} else {
			polygon(c(P$x, 0),
				c(P$y, 0),
				angle=angle[i], 
				border=NA,
				col=col[i])
			segments(c(P$x[1],
					P$x[length(P$x)]),
				c(P$y[1],
					P$y[length(P$y)]),
				0,
				0)
		}
		P <- t2xy(mean(x[i + 0:1]))
		lab <- as.character(labels[i])
		dist <- min(abs(deltax/cos(P$an)),
			abs(deltay/sin(P$an)))
		dist <- dist - radius*1.11 # subtract offset from origin
		if (strwidth(labels[i], cex=cex) > dist) {
			for (n in (nchar(labels[i]) - 3):1) {
				labels[i] <- paste(substring(labels[i], 1, n),
					"...",
					sep="")
				if (strwidth(labels[i], cex=cex) <= dist)
					break
			}
		}
		if (!is.na(lab) && nzchar(lab)) {
			lines(c(1, 1.05)*P$x,
				c(1, 1.05)*P$y)
			text(1.1*P$x,
				1.1*P$y,
				labels[i],
				xpd=TRUE,
				srt=ifelse(P$x < 0,
					P$an/pi*180 + 180,
					P$an/pi*180),
				adj=ifelse(P$x < 0,
					1,
					0),
				cex=cex,
				font=font[i])
		}
	}
	par(p)
}

.plotIDs <- function(x, showRanks, n) {
	taxa <- lapply(x,
		function(x)
			x[[1]][-1])
	d <- which(!duplicated(taxa))
	u <- taxa[d]
	m <- match(taxa, u)
	tot <- tapply(n, m, sum)
	
	if (all(lengths(x[d]) > 2)) {
		hasRanks <- TRUE
		lastRank <- sapply(x[d],
			function(x) {
				if (startsWith(tail(x[[1]], n=1), "unclassified_")) {
					""
				} else {
					tail(x[[3]], n=1)
				}
			})
		rank <- lapply(x[d],
			function(x)
				x[[3]][-1])
		ranks <- unlist(lapply(rank,
				function(x)
					c("Root", x, "")),
			use.names=FALSE)
		all_ranks <- unlist(ranks)
		u_ranks <- unique(unlist(rank))
		u_all_ranks <- c("Root", u_ranks, "")
		l <- length(u_ranks)
		parent <- integer(l)
		for (i in seq_len(l)) {
			w <- which(all_ranks==u_ranks[i])
			m <- match(all_ranks[w - 1], u_all_ranks)
			parent[i] <- max(m)
		}
		temp <- c("Root", "")
		remaining <- seq_len(l)
		i <- 1L
		while (length(remaining) > 0) {
			w <- which(temp==u_all_ranks[parent[remaining[i]]])
			if (length(w) != 0) {
				temp <- c(temp[1:w],
					u_ranks[remaining[i]],
					temp[(w + 1):length(temp)])
				remaining <- remaining[-i]
			} else {
				i <- i + 1L
			}
			if (i > length(remaining))
				i <- 1L
		}
		u_ranks <- temp[-c(1, length(temp))]
	} else {
		hasRanks <- FALSE
		rank <- lapply(u,
			function(x) {
				paste("Level", seq_along(x))
			})
		u_ranks <- unique(unlist(rank))
		l <- length(u_ranks)
	}
	
	taxa <- mapply(function(x, y) {
			if (length(x) < l) {
				X <- character(l)
				m <- match(y, u_ranks)
				X[m] <- x
				last <- x[1]
				j <- 1L
				while (j <= l) {
					if (X[j]=="") {
						X[j] <- last
					} else {
						last <- X[j]
					}
					j <- j + 1L
				}
				return(X)
			} else {
				return(x)
			}
		},
		u,
		rank)
	idsTbl <- matrix(unlist(taxa),
		nrow=l)
	
	# reorder alphabetically by rank/group
	rev_order <- function(...)
		order(..., decreasing=TRUE)
	o <- do.call(rev_order,
		as.data.frame(t(idsTbl)))
	tot <- tot[o]
	idsTbl <- idsTbl[, o, drop=FALSE]
	if (hasRanks)
		lastRank <- lastRank[o]
	
	org_colors <- rainbow(length(o), s=0.6)
	colors <- col2rgb(org_colors)
	for (i in rev(seq_len(l))) {
		groups <- idsTbl[i,]
		counts <- tot
		t <- tapply(counts,
			groups,
			sum)
		t <- t[unique(groups)]
		t <- t/sum(t)
		m <- match(groups,
			names(t))
		cols <- character(length(t))
		for (j in seq_along(t)) {
			w <- which(m==j)
			cols[j] <- rgb(weighted.mean(colors[1, w],
					counts[w]),
				weighted.mean(colors[2, w],
					counts[w]),
				weighted.mean(colors[3, w],
					counts[w]),
				maxColorValue=255)
		}
		if (i==l) { # draw outer ring
			names(t)[t < 0.01] <- ""
			if (hasRanks) {
				font <- ifelse(lastRank %in% c("genus", "species"), 3, 1)
			} else {
				font <- 1
			}
			if (showRanks) {
				.spoke(t,
					col=cols,
					cex=0.8,
					font=font)
			} else {
				.spoke(t,
					col=cols,
					cex=0.8,
					final.angle=360,
					font=font)
				break
			}
		} else {
			p <- par(new=TRUE)
			names(t)[] <- ""
			.spoke(t,
				radius=0.6*i/l,
				col=cols,
				cex=0.8,
				circle=FALSE)
			p <- par(p)
		}
	}
	if (showRanks) {
		text(0,
			0.6*(seq_len(l) - 0.5)/l,
			u_ranks,
			pos=4)
		segments(0,
			0.6*seq_len(l)/l,
			0.02,
			0.6*seq_len(l)/l)
	}
	
	o <- order(idsTbl[nrow(idsTbl),])
	invisible(org_colors[o])
}

`[.Taxa` <- function(x, i) {
	ans <- NextMethod("[", x)
	if (class(x)[1] != "Taxa") {
		stop("x must be an object of class 'Taxa'.")
	} else if (class(x)[2]=="Train") {
		x <- unclass(x)[i]
	} else if (class(x)[2]=="Test") {
		x <- unclass(x)[i]
		class(x) <- c("Taxa", "Test")
	} else {
		stop("x has unrecognized class: ", class(x)[2])
	}
	return(x)
}

`c.Taxa` <- function(...) {
	ans <- NextMethod("c", unclass(...))
	if (class(as.list(substitute(...)))[2]=="Test")
		class(ans) <- c("Taxa", "Test")
	return(ans)
}

plot.Taxa <- function(x, y=NULL, showRanks=TRUE, n=NULL, ...) {
	if (class(x)[1] != "Taxa")
		stop("x must be an object of class 'Taxa'.")
	if (!is.logical(showRanks))
		stop("showRanks must be a logical.")
	if (!is.null(y)) {
		if (class(y)[1] != "Taxa")
			stop("y must be an object of class 'Taxa'.")
		if (class(x)[2]=="Train") {
			train <- x
			if (class(y)[2]!="Test")
				stop("y must be an object of class 'Test' when x is an object of class 'Train'.")
			test <- y
		} else if (class(x)[2]=="Test") {
			test <- x
			if (class(y)[2]!="Train")
				stop("y must be an object of class 'Train' when x is an object of class 'Test'.")
			train <- y
		} else {
			stop("x has unrecognized class: ", class(x)[2])
		}
		if (is.null(n)) {
			n <- rep(1, length(test))
		} else {
			if (length(n) != length(test))
				stop("n has the wrong length.")
			if (!is.numeric(n))
				stop("n must be a numeric.")
			if (any(n <= 0))
				stop("n must all be greater than zero.")
		}
		
		layout(matrix(1:2),
			heights=c(2, 1))
		
		if (length(test) > 0) {
			groups <- sapply(test,
				function(x)
					tail(x$taxon, n=1))
			upone <- sapply(test,
				function(x)
					tail(x$taxon, n=2)[1])
			t <- sort(unique(groups))
			
			# record the lowest taxon present in taxa
			m <- match(t, groups)
			groups <- t
			w <- which(!(groups %in% train$taxa))
			if (length(w) > 0)
				groups[w] <- upone[m[w]]
			
			# Create a pie chart
			colors <- .plotIDs(test, showRanks, n)
			
			# Create a taxonomic tree
			makeTree <- function(I) {
				if (length(train$children[[I]])==0) {
					z <- I
					attr(z, "members") <- 1L
					attr(z, "leaf") <- TRUE
				} else {
					l <- length(train$children[[I]])
					z <- vector("list", l)
					members <- 0L
					for (i in seq_len(l)) {
						z[[i]] <- makeTree(train$children[[I]][i])
						members <- members + attr(z[[i]], "members")
					}
					attr(z, "members") <- members
				}
				attr(z, "height") <- 1 - (train$levels[I] - 1)/m
				x <- match(train$taxa[I], groups)
				if (!is.na(x))
					attr(z, "edgePar") <- list(col=colors[x])
				
				return(z)
			}
			m <- max(train$levels)
			tree <- makeTree(1)
			attr(tree, "label") <- "Root"
			attr(tree, "height") <- 1
			attr(tree, "members") <- sum(sapply(tree, attr, "members"))
			tree <- .applyMidpoints(tree, length(train$levels))
			# edgePar for edge.root is disregarded, so add parent node
			tree <- list(tree)
			attr(tree, "height") <- 1 + 1/m
			attr(tree, "members") <- attr(tree[[1]], "members")
			attr(tree, "midpoint") <- attr(tree[[1]], "midpoint")
			class(tree) <- "dendrogram"
			
			p <- par(mar=c(0.1, 0.1, 2.1, 0.1))
			# plot taxonomic tree
			dev.hold()
			on.exit(dev.flush())
			plot(tree,
				main="Distribution on taxonomic tree",
				leaflab="none",
				yaxt="n")
			par(new=TRUE)
			# overlay colored branches on top
			tree <- dendrapply(tree,
				function(x) {
					if (is.null(attr(x, "edgePar")))
						attr(x, "edgePar") <- list(col=NA)
					x
				})
			plot(tree,
				main="",
				leaflab="none",
				yaxt="n")
			par(p)
		}
	} else if (class(x)[2]=="Train") {
		if (!is.null(n))
			stop("n must be NULL if x is of class 'Train'.")
		layout(matrix(c(1, 2, 3, 1, 2, 4), nrow=3))
		
		# Create a taxonomic tree
		makeTree <- function(I) {
			if (length(x$children[[I]])==0) {
				z <- I
				attr(z, "members") <- 1L
				attr(z, "leaf") <- TRUE
			} else {
				l <- length(x$children[[I]])
				z <- vector("list", l)
				members <- 0L
				for (i in seq_len(l)) {
					z[[i]] <- makeTree(x$children[[I]][i])
					members <- members + attr(z[[i]], "members")
				}
				attr(z, "members") <- members
			}
			attr(z, "height") <- 1 - (x$levels[I] - 1)/m
			if (x$taxonomy[I] %in% x$problemGroups)
				attr(z, "edgePar") <- list(col="magenta")
			
			return(z)
		}
		m <- max(x$levels)
		tree <- makeTree(1)
		attr(tree, "label") <- "Root"
		attr(tree, "height") <- 1
		attr(tree, "members") <- sum(sapply(tree, attr, "members"))
		if (x$taxonomy[1] %in% x$problemGroups)
			attr(tree, "edgePar") <- list(col="magenta")
		tree <- .applyMidpoints(tree, length(x$levels))
		# edgePar for edge.root is disregarded, so add parent node
		tree <- list(tree)
		attr(tree, "height") <- 1 + 1/m
		attr(tree, "members") <- attr(tree[[1]], "members")
		attr(tree, "midpoint") <- attr(tree[[1]], "midpoint")
		class(tree) <- "dendrogram"
		
		p <- par(mar=c(0.1, 0.1, 2.1, 0.1))
		# plot taxonomic tree
		dev.hold()
		on.exit(dev.flush())
		plot(tree,
			main="",
			leaflab="none",
			yaxt="n")
		par(new=TRUE)
		# overlay colored branches on top
		tree <- dendrapply(tree,
			function(x) {
				if (is.null(attr(x, "edgePar")))
					attr(x, "edgePar") <- list(col=NA)
				x
			})
		plot(tree,
			main="",
			leaflab="none",
			yaxt="n")
		mtext(expression(bold("Taxonomic tree (problem groups in "*phantom("magenta")*")")),
			col="black",
			cex=0.8,
			line=0.4)
		mtext(expression(bold(phantom("Taxonomic tree (problem groups in ")*"magenta"*phantom(")"))),
			col="magenta",
			cex=0.8,
			line=0.4)
		
		par(mar=c(3.1, 4.1, 2.1, 1.1))
		if (is.null(x$ranks)) {
			t <- table(x$levels)
			names(t) <- paste("Level", names(t))
			main <- "Frequency of taxonomic levels"
		} else {
			t <- table(x$ranks)
			w <- which(!duplicated(x$ranks))
			u <- x$ranks[w]
			p <- c(0,
				match(x$ranks[x$parents[w]],
					u))
			c <- match(x$ranks[sapply(x$children[w], `[`, 1)],
				u)
			c[is.na(c)] <- length(w) + 1L
			o <- order(p, c)
			t <- t[u[o]]
			main <- "Frequency of taxonomic ranks"
		}
		plot(t,
			main=main,
			ylab="Frequency")
		
		par(mar=c(4.1, 4.1, 2.1, 1.1))
		w <- which(lengths(x$children)==0)
		s <- lengths(x$sequences[w])
		if (diff(range(s)) >= 79) {
			s <- log10(s)
			useLog <- TRUE
		} else {
			useLog <- FALSE
		}
		hist(s,
			breaks=seq(0,
				ceiling(max(s)),
				length.out=10),
			main="Sequences per taxonomic label",
			xlab="Number of sequence representatives",
			xaxt="n")
		ticks <- axTicks(1)
		if (useLog) {
			axis(1, ticks, round(10^ticks))
		} else {
			axis(1, ticks, ticks)
		}
		
		idf <- x$IDFweights
		idf <- sort(idf)
		suppressWarnings(plot(idf,
			type="l",
			main="Inverse document frequency (IDF) weights",
			ylab="Weight",
			xlab=paste("Sorted ", x$K, "-mers", sep=""),
			log="y"))
		abline(h=1, lty=2)
		
		par(p)
	} else if (class(x)[2]=="Test") {
		if (is.null(n)) {
			n <- rep(1, length(x))
		} else {
			if (length(n) != length(x))
				stop("n must be the same length as x.")
			if (!is.numeric(n))
				stop("n must be a numeric.")
			if (any(n <= 0))
				stop("n must all be greater than zero.")
		}
		if (length(x) > 0)
			.plotIDs(x, showRanks, n)
	} else {
		stop("x has unrecognized class: ", class(x)[2])
	}
	
	invisible(x)
}

print.Taxa <- function(x, ...) {
	if (class(x)[1] != "Taxa") {
		stop("x must be an object of class 'Taxa'.")
	} else if (class(x)[2]=="Train") {
		cat(paste("  A training set of class 'Taxa'",
			"\n   * K-mer size: ",
			x$K,
			"\n   * Number of rank levels: ",
			max(x$levels) + 1L,
			"\n   * Total number of sequences: ",
			length(x$sequences[[1]]),
			"\n   * Number of taxonomic groups: ",
			sum(lengths(x$children)==0),
			"\n   * Number of problem groups: ",
			length(x$problemGroups),
			"\n   * Number of problem sequences: ",
			nrow(x$problemSequences),
			"\n",
			sep=""))
	} else if (class(x)[2]=="Test") {
		l <- length(x)
		cat(paste("  A test set of class 'Taxa' with length ",
			l,
			"\n",
			sep=""))
		if (l==0)
			return(x)
		if (l > 10) {
			w1 <- 1:5
			w2 <- (l - 4):l
		} else {
			w1 <- seq_len(l)
			w2 <- integer()
		}
		o <- options()$width
		width <- nchar(l) + 2
		header <- format("",
			width=width)
		header <- paste(header,
			" confidence ",
			sep="")
		p1 <- format(paste("[", w1, "]", sep=""),
			width=width,
			justify="right")
		p1 <- paste(p1,
			format(sapply(x[w1],
					function(x)
						tail(x$confidence,
							n=1)),
				width=10,
				nsmall=1,
				digits=3,
				justify="right"),
			"% ",
			sep="")
		pred <- sapply(x[w1],
			function(x)
				paste(x$taxon,
					collapse="; "))
		remaining <- o - width - 12
		if (is.null(names(x))) { # no names
			w <- which(nchar(pred) > remaining)
			if (length(w) > 0) {
				pred <- substring(pred,
					1,
					remaining - 3)
				pred <- paste(pred,
					"...",
					sep="")
			}
			header <- paste(header,
				"taxon",
				sep="")
			p1 <- paste(p1,
				format(pred,
					width=remaining,
					justify="left"),
				sep="")
		} else { # names
			remaining <- remaining - 21
			ns <- names(x)[w1]
			w <- which(nchar(ns) > 20)
			if (length(w) > 0) {
				ns <- substring(ns,
					1,
					17)
				ns <- paste(ns,
					"...",
					sep="")
			}
			p1 <- paste(p1,
				format(ns,
					width=20,
					justify="left"),
				sep="")
			w <- which(nchar(pred) > remaining)
			if (length(w) > 0) {
				pred <- substring(pred,
					1,
					remaining - 3)
				pred <- paste(pred,
					"...",
					sep="")
			}
			header <- paste(header,
				"name                 taxon",
				sep="")
			p1 <- paste(p1,
				" ",
				format(pred,
					width=remaining,
					justify="left"),
				sep="")
		}
		cat(header,
			p1,
			sep="\n")
		if (length(w2) > 0) {
			dots <- format("...",
				width=width,
				justify="right")
			dots <- paste(dots,
				"        ...",
				sep="")
			p2 <- format(paste("[", w2, "]", sep=""),
				width=width,
				justify="right")
			p2 <- paste(p2,
				format(sapply(x[w2],
						function(x)
							tail(x$confidence,
								n=1)),
					width=10,
					nsmall=1,
					digits=3,
					justify="right"),
				"% ",
				sep="")
			pred <- sapply(x[w2],
				function(x)
					paste(x$taxon,
						collapse="; "))
			if (is.null(names(x))) { # no names
				w <- which(nchar(pred) > remaining)
				dots <- paste(dots,
					" ...",
					sep="")
				if (length(w) > 0) {
					pred <- substring(pred,
						1,
						remaining - 3)
					pred <- paste(pred,
						"...",
						sep="")
				}
				p2 <- paste(p2,
					format(pred,
						width=remaining,
						justify="left"),
					sep="")
			} else { # names
				ns <- names(x)[w2]
				w <- which(nchar(ns) > 20)
				if (length(w) > 0) {
					ns <- substring(ns,
						1,
						17)
					ns <- paste(ns,
						"...",
						sep="")
				}
				p2 <- paste(p2,
					format(ns,
						width=20,
						justify="left"),
					sep="")
				w <- which(nchar(pred) > remaining)
				if (length(w) > 0) {
					pred <- substring(pred,
						1,
						remaining - 3)
					pred <- paste(pred,
						"...",
						sep="")
				}
				dots <- paste(dots,
					" ...                  ...",
					sep="")
				p2 <- paste(p2,
					" ",
					format(pred,
						width=remaining,
						justify="left"),
					sep="")
			}
			cat(dots, p2, sep="\n")
		}
	} else {
		stop("x has unrecognized class: ", class(x)[2])
	}
	
	invisible(x)
}

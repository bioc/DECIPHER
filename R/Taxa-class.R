.spoke <- function(x,
	labels=names(x),
	edges=200,
	radius=0.6,
	angle=45,
	init.angle=0,
	col,
	cex=1) {
	labels <- as.graphicsAnnot(labels)
	x <- c(0, cumsum(x)/sum(x))
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
		polygon(c(P$x, 0),
			c(P$y, 0),
			angle=angle[i], 
			border=NULL,
			col=col[i])
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
					P$an/pi*180+180,
					P$an/pi*180),
				adj=ifelse(P$x < 0,
					1,
					0),
				cex=cex)
		}
	}
	par(p)
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
	ans <- NextMethod("c", ...)
	if (class(as.list(substitute(...)))[2]=="Test")
		class(ans) <- c("Taxa", "Test")
	return(ans)
}

plot.Taxa <- function(x, y=NULL, ...) {
	if (class(x)[1] != "Taxa")
		stop("x must be an object of class 'Taxa'.")
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
		
		layout(matrix(1:2),
			heights=c(2, 1))
		
		if (length(test) > 0) {
			groups <- sapply(test,
				function(x)
					tail(x$taxon, n=1))
			upone <- sapply(test,
				function(x)
					tail(x$taxon, n=2)[1])
			t <- table(groups)
			t <- t/sum(t)
			
			# record the lowest taxon present in taxa
			m <- match(names(t), groups)
			groups <- names(t)
			w <- which(!(groups %in% train$taxa))
			if (length(w) > 0)
				groups[w] <- upone[m[w]]
			
			# Create a pie chart
			names(t)[t < 0.01] <- ""
			colors <- rainbow(length(t), s=0.6)
			.spoke(t,
				col=colors,
				cex=0.8)
			
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
				attr(z, "edgePar") <- list(col="red")
			
			return(z)
		}
		m <- max(x$levels)
		tree <- makeTree(1)
		attr(tree, "label") <- "Root"
		attr(tree, "height") <- 1
		attr(tree, "members") <- sum(sapply(tree, attr, "members"))
		if (x$taxonomy[1] %in% x$problemGroups)
			attr(tree, "edgePar") <- list(col="red")
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
		mtext(expression(bold("Taxonomic tree (problem groups in "*phantom("red")*")")),
			col="black",
			cex=0.8,
			line=0.4)
		mtext(expression(bold(phantom("Taxonomic tree (problem groups in ")*"red"*phantom(")"))),
			col="red",
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
		if (length(x) > 0) {
			groups <- sapply(x,
				function(x)
					tail(x$taxon, n=1))
			t <- table(groups)
			t <- t/sum(t)
			names(t)[t < 0.01] <- ""
			
			.spoke(t,
				col=rainbow(length(t), s=0.6),
				cex=0.8)
		}
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
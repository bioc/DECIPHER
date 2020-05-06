# below function modified from stats package
to.dendrogram <- function(object, states=NULL) {
	
	z <- list()
	oHgts <- object$lengths
	oHgt <- object$height
	nMerge <- length(oHgt)
	if (nMerge != nrow(object$merge))
		stop("'merge' and 'height' do not fit!")
	hMax <- oHgt[nMerge]
	
	count <- 1L
	one <- 1L
	two <- 2L
	for (k in seq_len(nMerge)) {
		x <- as.integer(object$merge[k, ])
		neg <- x < 0
		if (all(neg)) { # two leaves
			zk <- as.list(-x)
			attr(zk, "members") <- two
			objlabels <- object$labels[-x]
			attr(zk[[1L]], "label") <- objlabels[1L]
			attr(zk[[2L]], "label") <- objlabels[2L]
			attr(zk[[1L]], "members") <- attr(zk[[2L]], "members") <- one
			attr(zk[[1L]], "height") <- oHgt[k] - oHgts[k, 1]
			attr(zk[[2L]], "height") <- oHgt[k] - oHgts[k, 2]
			attr(zk[[1L]], "leaf") <- attr(zk[[2L]], "leaf") <- TRUE
		} else if (any(neg)) { # one leaf, one node
			X <- as.character(x)
			isL <- x[1L] < 0 # is leaf left?
			zk <-
			if (isL) {
				list(-x[1L], z[[X[2L]]])
			} else {
				list(z[[X[1L]]], -x[2L])
			}
			attr(zk, "members") <- attr(z[[X[1 + isL]]], "members") + one
			attr(zk[[2 - isL]], "members") <- one
			attr(zk[[2 - isL]], "height") <- oHgt[k] - oHgts[k, 2 - isL]
			attr(zk[[2 - isL]], "label") <- object$labels[-x[2 - isL]]
			attr(zk[[2 - isL]], "leaf") <- TRUE
			attr(zk[[1 + isL]], "state") <- states[count]
			count <- count + 1L
		} else { # two nodes
			x <- as.character(x)
			zk <- list(z[[x[1L]]], z[[x[2L]]])
			attr(zk, "members") <- attr(z[[x[1L]]], "members") + attr(z[[x[2L]]], "members")
			attr(zk[[1L]], "state") <- states[count]
			count <- count + 1L
			attr(zk[[2L]], "state") <- states[count]
			count <- count + 1L
		}
		attr(zk, "height") <- oHgt[k]
		k <- as.character(k)
		z[[k]] <- zk
	}
	z <- z[[k]]
	attr(z, "state") <- states[count]
	class(z) <- "dendrogram"
	z
}

.collapse <- function(dend, collapse, dim) {
	# initialize a stack of maximum length (dim)
	stack <- vector("list", dim)
	visit <- logical(dim) # node already visited
	parent <- integer(dim) # index of parent node
	index <- integer(dim) # index in parent node
	pos <- 1L # current position in the stack
	stack[[pos]] <- dend
	while (pos > 0L) { # more nodes to visit
		if (visit[pos]) { # ascending tree
			visit[pos] <- FALSE # reset visit
			
			if (!is.leaf(stack[[pos]][[1]])) {
				h1 <- attr(stack[[pos]][[1]], "height")
			} else {
				h1 <- -Inf
			}
			if (!is.leaf(stack[[pos]][[2]])) {
				h2 <- attr(stack[[pos]][[2]], "height")
			} else {
				h2 <- -Inf
			}
			
			h <- attr(stack[[pos]], "height")
			
			if ((h - h1) <= collapse || (h - h2) <= collapse) {
				# make multifurcating
				m1 <- attr(stack[[pos]][[1]], "members")
				m2 <- attr(stack[[pos]][[2]], "members")
				states <- c(attr(stack[[pos]][[1]], "state"),
					attr(stack[[pos]][[2]], "state"))
				m <- m1 + m2
				if ((h - h1) <= collapse && (h - h2) <= collapse) {
					l1 <- length(stack[[pos]][[1]])
					l2 <- length(stack[[pos]][[2]])
					x <- vector("list", l1 + l2)
					for (i in seq_len(l1))
						x[i] <- stack[[pos]][[1]][i]
					for (i in seq_len(l2))
						x[i + l1] <- stack[[pos]][[2]][i]
				} else if ((h - h1) <= collapse) {
					l <- length(stack[[pos]][[1]])
					x <- vector("list", l + 1)
					for (i in seq_len(l))
						x[i] <- stack[[pos]][[1]][i]
					x[l + 1] <- stack[[pos]][-1]
				} else if ((h - h2) <= collapse) {
					l <- length(stack[[pos]][[2]])
					x <- vector("list", l + 1)
					x[1] <- stack[[pos]][-2]
					for (i in seq_len(l))
						x[i + 1] <- stack[[pos]][[2]][i]
				}
				stack[[pos]] <- x
				
				attr(stack[[pos]], "height") <- h
				attr(stack[[pos]], "members") <- m
				attr(stack[[pos]], "state") <- unique(states)
				
				class(stack[[pos]]) <- "dendrogram"
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

.organizeClusters <- function(myClusters,
	dNames,
	o) {
	
	l <- length(dNames)
	clusters <- data.frame(cluster=integer(l),
		row.names=dNames)
	w <- which(myClusters[, 7] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 7]] <- as.integer(myClusters[w, 9])
	w <- which(myClusters[, 8] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 8]] <- as.integer(myClusters[w, 10])
	
	# order the cluster numbers to match
	# the order of the dendrogram
	temp <- 0
	l <- max(clusters$cluster)
	clustersTemp <- clusters
	v <- vector(mode="numeric",length=l)
	j <- 0
	for (i in 1:length(o)) {
		if (clusters$cluster[o[i]] != temp &
			length(which(v==clusters$cluster[o[i]]))==0) {
			temp <- clusters$cluster[o[i]]
			j <- j + 1
			v[j] <- temp
		}
	}
	for (k in 1:l) {
		w <- which(clusters$cluster == v[k])
		clustersTemp$cluster[w] <- k
	}
	clusters <- clustersTemp
	
	return(clusters)
}

.organizeClustersFast <- function(myClusters,
	dNames) {
	
	l <- length(dNames)
	clusters <- data.frame(cluster=integer(l),
		row.names=dNames)
	w <- which(myClusters[, 7] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 7]] <- as.integer(myClusters[w, 9])
	w <- which(myClusters[, 8] < 0)
	if (length(w) > 0)
		clusters$cluster[-1*myClusters[w, 8]] <- as.integer(myClusters[w, 10])
	
	return(clusters)
}

.rates <- function(alpha, nBins) {
	
	# Determine rates based on alpha and the number of bins
	# bins roots normalized to 1 of the General Laguerre Quadrature
	# first nBins elements are rates with mean 1
	# second nBins elements are probabilities with sum 1
	
	findRoots <- function(alpha, nBins) {
		
		# Determine rates based on Gamma's alpha and the number of bins
		# bins roots normalized to 1 of the General Laguerre Polynomial (GLP)
		
		coeff  <- integer(nBins + 1)
		for (i in 0:nBins) {
			n <- nBins + alpha
			k <- nBins - i
			coeff[i + 1] <- (-1)^i*choose(nBins + alpha, nBins - i)/factorial(i)
		}
		
		return(sort(Re(polyroot(coeff))))
	}
	
	roots <- findRoots(alpha - 1, nBins)
	
	Laguerre <- function(x, alpha, degree) {
		y <- 0
		for (i in 0:degree) {
			y <- y + (-1)^i*choose(degree + alpha, degree - i)*x^i/factorial(i)
		}
		return(y)
	}
	
	weights <- numeric(nBins)
	f <- prod(1 + (alpha - 1)/(1:nBins))
	
	for (i in 1:nBins) {
		weights[i] <- f*roots[i]/((nBins + 1)^2*Laguerre(roots[i],
			alpha - 1,
			nBins + 1)^2)
	}
	
	roots <- roots/alpha
	
	return(c(roots, weights))
}

.optimizeModel <- function(myClusters,
	model,
	myDNAStringSet,
	N,
	processors=1) {
	
	rates <- as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", model))
	model <- sub("([^+]*)(\\+G(\\d+))?", "\\1", model)
	
	if (model=="JC69" && is.na(rates)) {
		LnL <- .Call("clusterML",
			myClusters,
			myDNAStringSet,
			c(0.25, 0.25, 0.25, 0.25, 1, 1, 1, 1),
			integer(),
			numeric(),
			0,
			1L,
			processors,
			PACKAGE="DECIPHER")
		K <- 2*dim(myClusters)[1] - 1
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(NA, NA, NA, NA, NA, NA, NA, LnL, AICc, BIC))
	} else if (model=="JC69") { # rates is an integer
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(0.25, 0.25, 0.25, 0.25, 1, 1, .rates(params, rates)),
				integer(),
				numeric(),
				0,
				1L,
				processors,
				PACKAGE="DECIPHER")
		}
		o <- optimize(f, c(0.001, 500), tol=1e-4)
		K <- 2*dim(myClusters)[1]
		AICc <- 2*K + 2*o$objective + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*o$objective + K*log(N)
		return(c(NA, NA, NA, NA, NA, NA, o$minimum, o$objective, AICc, BIC))
	} else if (model=="K80") {
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(0.25, 0.25, 0.25, 0.25, params, params, 1, 1),
				integer(),
				numeric(),
				0,
				1L,
				processors,
				PACKAGE="DECIPHER")
		}
		o <- optimize(f, c(0, 10), tol=1e-4)
		
		if (!is.na(rates)) {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(0.25, 0.25, 0.25, 0.25, o$minimum, o$minimum, .rates(params, rates)),
					integer(),
					numeric(),
					0,
					1L,
					processors,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 1
		} else {
			a <- NA
			LnL <- o$objective
			K <- 2*dim(myClusters)[1]
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(NA, NA, NA, NA, rep(o$minimum, 2), a, LnL, AICc, BIC))
	} else if (model=="T92") {
		baseFreqs <- alphabetFrequency(myDNAStringSet,
			baseOnly=TRUE,
			as.prob=TRUE,
			collapse=TRUE)[1:4]
		baseFreqs <- ifelse(baseFreqs < 0.01, 0.01, baseFreqs)
		baseFreqs <- baseFreqs/sum(baseFreqs)
		baseFreqs <- c((baseFreqs[1] + baseFreqs[4])/(2*sum(baseFreqs)),
			(baseFreqs[2] + baseFreqs[3])/(2*sum(baseFreqs)))
		baseFreqs <- c(baseFreqs[1], baseFreqs[2], baseFreqs[2], baseFreqs[1])
#		f <- function(params) {
#			if (params[1] > 0.49)
#				return(Inf)
#			LnL <- .Call("clusterML",
#				myClusters,
#				myDNAStringSet,
#				c(params, rep((1 - 2*params)/2, 2), params, 1, 1, 1, 1),
#				integer(),
#				numeric(),
#				0,
#				1L,
#				processors,
#				PACKAGE="DECIPHER")
#		}
#		o <- nlminb(0.25,
#			f,
#			upper=0.49,
#			lower=0.01,
#			control=list(rel.tol=1e-4))
#		baseFreqs <- c(o$par, rep((1 - 2*o$par)/2, 2), o$par)
		
		f <- function(params) {
			LnL <- .Call("clusterML",
				myClusters,
				myDNAStringSet,
				c(baseFreqs, params, params, 1, 1),
				integer(),
				numeric(),
				0,
				1L,
				processors,
				PACKAGE="DECIPHER")
		}
		o <- nlminb(1,
			f,
			upper=10,
			lower=0.01,
			control=list(rel.tol=1e-4))
		
		if (!is.na(rates)) {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, o$par[1], o$par[1], .rates(params, rates)),
					integer(),
					numeric(),
					0,
					1L,
					processors,
					PACKAGE="DECIPHER")
			}
			a <- optimize(f, c(0.001, 500), tol=1e-4)
			LnL <- a$objective
			a <- a$minimum
			K <- 2*dim(myClusters)[1] + 2
		} else {
			a <- NA
			LnL <- o$objective
			K <- 2*dim(myClusters)[1] + 1
		}
		
		AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
		BIC <- 2*LnL + K*log(N)
		return(c(baseFreqs, o$par[1], o$par[1], a, LnL, AICc, BIC))
	} else {
		baseFreqs <- alphabetFrequency(myDNAStringSet,
			baseOnly=TRUE,
			as.prob=TRUE,
			collapse=TRUE)[1:4]
		baseFreqs <- ifelse(baseFreqs < 0.01, 0.01, baseFreqs)
		baseFreqs <- baseFreqs/sum(baseFreqs)
#	f <- function(params) {
#			if (any(params > 0.8) ||
#				any(params < 0.05) ||
#				sum(params) > 0.95)
#				return(Inf)
#			LnL <- .Call("clusterML",
#				myClusters,
#				myDNAStringSet,
#				c(params[1], params[2], params[3], 1 - sum(params), 1, 1, 1, 1),
#				integer(),
#				numeric(),
#				0,
#				1L,
#				processors,
#				PACKAGE="DECIPHER")
#		}
#		o <- optim(baseFreqs[1:3],
#			f,
#			control=list(reltol=1e-3))
#		baseFreqs[1:3] <- o$par
#		baseFreqs[4] <- 1 - sum(baseFreqs[1:3])
		
		if (model=="F81") {
			if (!is.na(rates)) {
				f <- function(params) {
					LnL <- .Call("clusterML",
						myClusters,
						myDNAStringSet,
						c(baseFreqs, 1, 1, .rates(params, rates)),
						integer(),
						numeric(),
						0,
						1L,
						processors,
						PACKAGE="DECIPHER")
				}
				a <- optimize(f, c(0.001, 500), tol=1e-4)
				LnL <- a$objective
				a <- a$minimum
				K <- 2*dim(myClusters)[1] + 3
			} else {
				K <- 2*dim(myClusters)[1] + 2
				a <- NA
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, 1, 1, 1, 1),
					integer(),
					numeric(),
					0,
					1L,
					processors,
					PACKAGE="DECIPHER")
			}
			
			AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
			BIC <- 2*LnL + K*log(N)
			return(c(baseFreqs, NA, NA, a, LnL, AICc, BIC))
		} else if (model=="HKY85") {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, params, params, 1, 1),
					integer(),
					numeric(),
					0,
					1L,
					processors,
					PACKAGE="DECIPHER")
			}
			o <- nlminb(1,
				f,
				upper=10,
				lower=0.01,
				control=list(rel.tol=1e-4))
			
			if (!is.na(rates)) {
				f <- function(params) {
					LnL <- .Call("clusterML",
						myClusters,
						myDNAStringSet,
						c(baseFreqs, o$par[1], o$par[1], .rates(params, rates)),
						integer(),
						numeric(),
						0,
						1L,
						processors,
						PACKAGE="DECIPHER")
				}
				a <- optimize(f, c(0.001, 500), tol=1e-4)
				LnL <- a$objective
				a <- a$minimum
				K <- 2*dim(myClusters)[1] + 4
			} else {
				a <- NA
				LnL <- o$objective
				K <- 2*dim(myClusters)[1] + 3
			}
			
			AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
			BIC <- 2*LnL + K*log(N)
			return(c(baseFreqs, o$par[1], o$par[1], a, LnL, AICc, BIC))
		} else if (model=="TN93") {
			f <- function(params) {
				LnL <- .Call("clusterML",
					myClusters,
					myDNAStringSet,
					c(baseFreqs, params[1], params[2], 1, 1),
					integer(),
					numeric(),
					0,
					1L,
					processors,
					PACKAGE="DECIPHER")
			}
			o <- nlminb(c(1, 1),
				f,
				upper=c(10, 10),
				lower=c(0.01, 0.01),
				control=list(rel.tol=1e-4))
			
			if (!is.na(rates)) {
				f <- function(params) {
					LnL <- .Call("clusterML",
						myClusters,
						myDNAStringSet,
						c(baseFreqs, o$par[1], o$par[2], .rates(params, rates)),
						integer(),
						numeric(),
						0,
						1L,
						processors,
						PACKAGE="DECIPHER")
				}
				a <- optimize(f, c(0.001, 500), tol=1e-4)
				LnL <- a$objective
				a <- a$minimum
				K <- 2*dim(myClusters)[1] + 5
			} else {
				a <- NA
				LnL <- o$objective
				K <- 2*dim(myClusters)[1] + 4
			}
			
			AICc <- 2*K + 2*LnL + 2*K*(K + 1)/(N - K - 1)
			BIC <- 2*LnL + K*log(N)
			return(c(baseFreqs, o$par[1], o$par[2], a, LnL, AICc, BIC))
		}
	}
}

.simultaneousBrent <- function(f, # scalar function
	a, # lower bounds
	b, # best guesses
	c, # upper bounds
	tol=1e-5, # accuracy of params change in x
	overallTol=1e0, # accuracy of overall scalar
	relTol=1e-2) { # accuracy of params change in y
	
	# performs optimization using Brent's method
	# simultaneously optimizes nearly-independent parameters
	
	phi <- (3 - sqrt(5))/2
	
	x <- w <- v <- b
	fx <- numeric(length(x))
	baseline <- f(x)
	fx[] <- baseline
	fw <- fv <- fx
	
	b <- ifelse(a > c, a, c)
	
	W <- (1:length(x))[(c - a) > tol]
	if (length(W)==0)
		return(x)
	
	e <- xm <- tol1 <- tol2 <- numeric(length(W))
	q <- r <- p <- numeric(length(W))
	
	iteration <- 0L
	maxIterations <- 1000L
	delta <- 0
	f_delta <- rep(Inf, length(W))
	while(iteration < maxIterations &&
		(iteration < 3 || delta > overallTol)) {
		iteration <- iteration + 1L
		
		xm[W] <- (a[W] + b[W])/2
		tol1[W] <- tol*abs(x[W]) + 1e-10
		tol2[W] <- 2*(tol1[W])
		
		w1 <- which(abs(x[W] - xm[W]) <= (tol2[W] - (b[W] - a[W])/2) |
			f_delta < relTol)
		if (length(w1) > 0) {
			W <- W[-w1]
			if (length(W)==0)
				break
		}
		
		w1 <- which(abs(e[W]) > tol1[W])
		w2 <- which(!(1:length(W) %in% w1))
		d <- numeric(length(W))
		
		if (length(w1) > 0) {
			r[W[w1]] <- (x[W[w1]] - w[W[w1]])*(fx[W[w1]] - fv[W[w1]])
			q[W[w1]] <- (x[W[w1]] - v[W[w1]])*(fx[W[w1]] - fw[W[w1]])
			p[W[w1]] <- (x[W[w1]] - v[W[w1]])*q[W[w1]] - (x[W[w1]] - w[W[w1]])*r[W[w1]]
			q[W[w1]] <- 2*(q[W[w1]] - r[W[w1]])
			p[W[w1]] <- ifelse(q[W[w1]] > 0,  p[W[w1]], -p[W[w1]])
			q[W[w1]] <- abs(q[W[w1]])
			
			etemp <- e[W[w1]]
			e[W[w1]] <- d[W[w1]]
			
			w3 <- which(abs(p[W[w1]]) >= abs(q[W[w1]]*etemp/2) |
				p[W[w1]] <= q[W[w1]]*(a[W[w1]] - x[W[w1]]) |
				p[W[w1]] >= q[W[w1]]*(b[W[w1]] - x[W[w1]]))
			w4 <- which(!(1:length(w1) %in% w3))
			if (length(w3) > 0) {
				e[W[w1[w3]]] <- ifelse(x[W[w1[w3]]] >= xm[W[w1[w3]]],
					a[W[w1[w3]]] - x[W[w1[w3]]],
					b[W[w1[w3]]] - x[W[w1[w3]]])
				d[W[w1[w3]]] <- phi*e[W[w1[w3]]]
			}
			if (length(w4) > 0) {
				d[W[w1[w4]]] <- p[W[w1[w4]]]/q[W[w1[w4]]]
				u <- x[W[w1[w4]]] + d[W[w1[w4]]]
				d[W[w1[w4]]] <- ifelse(u - a[W[w1[w4]]] < tol2[W[w1[w4]]] |
						b[W[w1[w4]]] - u < tol2[W[w1[w4]]],
					ifelse(xm[W[w1[w4]]] - x[W[w1[w4]]] >= 0,
						abs(tol1[W[w1[w4]]]),
						-abs(tol1[W[w1[w4]]])),
					d[W[w1[w4]]])
			}
		}
		if (length(w2) > 0) {
			e[W[w2]] <- ifelse(x[W[w2]] >= xm[W[w2]],
				a[W[w2]] - x[W[w2]],
				b[W[w2]] - x[W[w2]])
			d[W[w2]] <- phi*e[W[w2]]
		}
		
		# perform one function call per iteration
		u <- ifelse(abs(d[W]) >= tol1[W],
			x[W] + d[W],
			x[W] + ifelse(d[W] > 0,
				abs(tol1[W]),
				-abs(tol1[W])))
		fu <- f(x, W, u) # provide all alternative parameters
		
		newBaseline <- fu[1]
		fu <- fu[2:length(fu)]
		f_delta <- abs(fu - fx[W])
		offset <- fx[W] - newBaseline
		w1 <- which(offset > 0)
		if (length(w1) > 0) {
			fx[W[w1]] <- fx[W[w1]] - offset[w1]
			fw[W[w1]] <- fw[W[w1]] - offset[w1]
			fv[W[w1]] <- fv[W[w1]] - offset[w1]
		}
		delta <- abs(baseline - newBaseline)
		baseline <- newBaseline
		
		w1 <- which(fu <= fx[W])
		w2 <- which(!(1:length(W) %in% w1))
		if (length(w1) > 0){
			w3 <- which(u[w1] >= x[W[w1]])
			if (length(w3) > 0)
				a[W[w1[w3]]] <- x[W[w1[w3]]]
			w4 <- which(u[w1] < x[W[w1]])
			if (length(w4) > 0)
				b[W[w1[w4]]] <- x[W[w1[w4]]]
			v[W[w1]] <- w[W[w1]]
			w[W[w1]] <- x[W[w1]]
			x[W[w1]] <- u[w1]
			fv[W[w1]] <- fw[W[w1]]
			fw[W[w1]] <- fx[W[w1]]
			fx[W[w1]] <- fu[w1]
		}
		if (length(w2) > 0) {
			w3 <- which(u[w2] < x[W[w2]])
			if (length(w3) > 0)
				a[W[w2[w3]]] <- u[w2[w3]]
			w4 <- which(u[w2] >= x[W[w2]])
			if (length(w4) > 0)
				b[W[w2[w4]]] <- u[w2[w4]]
			
			w3 <- which(fu[w2] <= fw[W[w2]] | w[W[w2]]==x[W[w2]])
			if (length(w3) > 0) {
				v[W[w2[w3]]] <- w[W[w2[w3]]]
				w[W[w2[w3]]] <- u[w2[w3]]
				fv[W[w2[w3]]] <- fw[W[w2[w3]]]
				fw[W[w2[w3]]] <- fu[w2[w3]]
			}
			
			w4 <- which(fu[w2] <= fv[W[w2]] | v[W[w2]]==x[W[w2]] | v[W[w2]]==w[W[w2]])
			if (length(w4) > 0) {
				v[W[w2[w4]]] <- u[w2[w4]]
				fv[W[w2[w4]]] <- fu[w2[w4]]
			}
		}
	}
	
	return(x)
}

.reorderClusters <- function(myClusters, all=FALSE) {
	
	# order clusters by branching pattern
	repeat {
		a <- apply(as.matrix(myClusters[, 7:8], ncol=2), 1, max)
		w <- which(a > 1:dim(myClusters)[1])[1]
		if (is.na(w))
			break
		temp <- myClusters[w, c(4, 5, 7, 8)]
		myClusters[w:(a[w] - 1), c(4, 5, 7, 8)] <- myClusters[(w + 1):a[w], c(4, 5, 7, 8)]
		myClusters[a[w], c(4, 5, 7, 8)] <- temp
		w1 <- which(myClusters[w:dim(myClusters)[1], 7:8] %in% (w + 1):a[w])
		w2 <- which(myClusters[(a[w] + 1):dim(myClusters)[1], 7:8] %in% w)
		if (length(w1) > 0)
			myClusters[w:dim(myClusters)[1], 7:8][w1] <- myClusters[w:dim(myClusters)[1], 7:8][w1] - 1
		if (length(w2) > 0)
			myClusters[(a[w] + 1):dim(myClusters)[1], 7:8][w2] <- a[w]
	}
	
	if (all) { # also renumber columns 1 to 3
		count <- 0L
		myClusters[, 1:2] <- myClusters[, 7:8]
		for (i in 1:dim(myClusters)[1]) {
			if (myClusters[i, 7] > 0 && myClusters[i, 8] > 0) {
				myClusters[i, 1] <- myClusters[myClusters[i, 7], 3]
				myClusters[i, 2] <- myClusters[myClusters[i, 8], 3]
				count <- count + 1L
				myClusters[i, 3] <- count
			} else if (myClusters[i, 7] > 0) {
				myClusters[i, 1] <- myClusters[myClusters[i, 7], 3]
				myClusters[i, 3] <- myClusters[i, 1]
			} else if (myClusters[i, 8] > 0) {
				myClusters[i, 2] <- myClusters[myClusters[i, 8], 3]
				myClusters[i, 3] <- myClusters[i, 2]
			} else {
				count <- count + 1L
				myClusters[i, 3] <- count
			}
		}
	}
	
	return(myClusters)
}

.swapBranches <- function(myClusters, r1, c1, r2, c2) {
	
	# swap branch [r1, c1] with [r2, c2]
	temp <- myClusters[r1, c1]
	myClusters[r1, c1] <- myClusters[r2, c2]
	myClusters[r2, c2] <- temp
	temp <- myClusters[r1, c1 - 3]
	myClusters[r1, c1 - 3] <- myClusters[r2, c2 - 3]
	myClusters[r2, c2 - 3] <- temp
	
	myClusters <- .reorderClusters(myClusters)
	
	return(myClusters)
}

.NNI <- function(myClusters, bestLnL, NNIs, maximizeLikelihood, tol=1e-1) {
	
	if (dim(myClusters)[1]==1)
		return(myClusters)
	
	# perform rounds of nearest-neighbor interchanges
	W <- dim(myClusters)[1]:2
	while (length(W) > 0) {
		i <- W[1]
		W <- W[-1]
		
		if (i==dim(myClusters)[1]) { # swap nodes
			if (all(myClusters[i, 7:8] > 0)) { # two nodes
				myClustersTemp1 <- .swapBranches(myClusters,
					myClusters[i, 8], 7,
					myClusters[i, 7], 7)
				tempLnL1 <- maximizeLikelihood(myClustersTemp1, NNIs + 1, tol)
				myClustersTemp2 <- .swapBranches(myClusters,
					myClusters[i, 8], 7,
					myClusters[i, 7], 8)
				tempLnL2 <- maximizeLikelihood(myClustersTemp2, NNIs + 1, tol)
				if (tempLnL1 < bestLnL - tol &&
					tempLnL1 < tempLnL2) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp1
					bestLnL <- tempLnL1
				} else if (tempLnL2 < bestLnL - tol) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp2
					bestLnL <- tempLnL2
				}
			}
		} else {
			w <- which(myClusters[i, 7:8] > 0)
			if (length(w)==2) {
				myClustersTemp1 <- .swapBranches(myClusters,
					i, 8,
					myClusters[i, 7], 7)
				tempLnL1 <- maximizeLikelihood(myClustersTemp1, NNIs + 1, tol)
				myClustersTemp2 <- .swapBranches(myClusters,
					i, 8,
					myClusters[i, 7], 8)
				tempLnL2 <- maximizeLikelihood(myClustersTemp2, NNIs + 1, tol)
				myClustersTemp3 <- .swapBranches(myClusters,
					i, 7,
					myClusters[i, 8], 7)
				tempLnL3 <- maximizeLikelihood(myClustersTemp3, NNIs + 1, tol)
				myClustersTemp4 <- .swapBranches(myClusters,
					i, 7,
					myClusters[i, 8], 8)
				tempLnL4 <- maximizeLikelihood(myClustersTemp4, NNIs + 1, tol)
				if (tempLnL1 < bestLnL - tol &&
					tempLnL1 < tempLnL2 &&
					tempLnL1 < tempLnL3 &&
					tempLnL1 < tempLnL4) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp1
					bestLnL <- tempLnL1
				} else if (tempLnL2 < bestLnL - tol &&
					tempLnL2 < tempLnL3 &&
					tempLnL2 < tempLnL4) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp2
					bestLnL <- tempLnL2
				} else if (tempLnL3 < bestLnL - tol &&
					tempLnL3 < tempLnL4) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp3
					bestLnL <- tempLnL3
				} else if (tempLnL4 < bestLnL - tol) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp4
					bestLnL <- tempLnL4
				}
			} else if (length(w)==1) {
				myClustersTemp1 <- .swapBranches(myClusters,
					i, ifelse(w==1, 8, 7),
					myClusters[i, w + 6], 7)
				tempLnL1 <- maximizeLikelihood(myClustersTemp1, NNIs + 1, tol)
				myClustersTemp2 <- .swapBranches(myClusters,
					i, ifelse(w==1, 8, 7),
					myClusters[i, w + 6], 8)
				tempLnL2 <- maximizeLikelihood(myClustersTemp2, NNIs + 1, tol)
				if (tempLnL1 < bestLnL - tol &&
					tempLnL1 < tempLnL2) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp1
					bestLnL <- tempLnL1
				} else if (tempLnL2 < bestLnL - tol) {
					W <- c(W, i)
					NNIs <- NNIs + 1
					myClusters <- myClustersTemp2
					bestLnL <- tempLnL2
				}
			}
		}
	}
	
	# return the highest likelihood tree
	return(myClusters)
}

MODELS <- c("JC69",
	"JC69+G4",
	"K80",
	"K80+G4",
	"F81",
	"F81+G4",
	"HKY85",
	"HKY85+G4",
	"T92",
	"T92+G4",
	"TN93",
	"TN93+G4")

.splitClusters <- function(x, y) {
	
	clusterNum <- 0L
	X <- integer(length(x))
	u.y <- unique(y)
	for (i in u.y) {
		w.y <- which(y==i)
		u.x <- unique(x[w.y])
		for (j in u.x) {
			clusterNum <- clusterNum + 1L
			w.x <- which(x[w.y]==j)
			X[w.y[w.x]] <- clusterNum
		}
	}
	return(X)
}

.root <- function(x1, root) {
	
	# if root is zero then midpoint root the tree
	# otherwise outgroup root based on root index
	# (note: output columns 1:3 are uncorrected)
	
	n <- nrow(x1)
	if (root==0) { # midpoint root
		# find the leaf at minimum height
		r1 <- which(x1[, 7] < 0)
		h1 <- x1[r1, 6] - x1[r1, 4]
		r2 <- which(x1[, 8] < 0)
		h2 <- x1[r2, 6] - x1[r2, 5]
		r <- c(r1, r2) # row number of leaf
		z <- rep(c(7L, 8L), c(length(r1), length(r2)))
		h <- c(h1, h2) # height of leaf
		
		# reorder by sequence number
		o <- order(x1[cbind(r, z)],
			decreasing=TRUE)
		h <- h[o]
		minH <- which.min(h) # index of lowest leaf
	} else { # outgroup root
		w <- which(x1[n, 7:8]==-root)
		if (length(w) > 0) { # already outgroup rooted
			# extend the root node
			x1[n, 6] <- x1[n, 6] + x1[n, 3 + w]
			x1[n, 3 + w] <- 0
			return(x1)
		}
		minH <- root # index of root
	}
	
	# find most distant leaf from minH
	longest <- numeric(n) # length of longest path
	merged <- logical(n) # whether merged yet
	index <- numeric(n) # column back to minH
	for (i in seq_len(n)) {
		b1 <- x1[i, 7]
		if (b1 < 0) { # merged with leaf
			if (b1==-minH) {
				merged[i] <- TRUE
				b1 <- NA_real_
			} else {
				l1 <- x1[i, 4]
			}
		} else { # merged with node
			if (merged[b1]) {
				merged[i] <- TRUE
				b1 <- NA_real_
			} else {
				l1 <- longest[b1] + x1[i, 4]
			}
		}
		
		b2 <- x1[i, 8]
		if (b2 < 0) { # merged with leaf
			if (b2==-minH) {
				merged[i] <- TRUE
				b2 <- NA_real_
			} else {
				l2 <- x1[i, 5]
			}
		} else { # merged with node
			if (merged[b2]) {
				merged[i] <- TRUE
				b2 <- NA_real_
			} else {
				l2 <- longest[b2] + x1[i, 5]
			}
		}
		
		if (is.na(b1)) { # b1 contains minH
			longest[i] <- l2
			# leave index[i] at zero
		} else if (is.na(b2)) { # b2 contains minH
			longest[i] <- l1
			index[i] <- 1
		} else if (l1 >= l2) {
			longest[i] <- l1
			# index[i] not needed
		} else { # l2 > l1
			longest[i] <- l2
			# index[i] not needed
		}
	}
	
	if (root==0) { # determine height of the midpoint
		w <- which(merged)
		longest <- longest + x1[, 6] - h[minH]
		m <- w[which.max(longest[w])]
		midH <- longest[m]/2
		if (isTRUE(all.equal(x1[n, 6], midH)))
			return(x1) # already midpoint rooted
		
		# find the edge containing the midpoint
		lowH <- x1[m, 6] - x1[m, 4 + index[m]]
		while (lowH > midH) { # descend the tree
			m <- x1[m, 7 + index[m]]
			lowH <- x1[m, 6] - x1[m, 4 + index[m]]
		}
	} else { # root at tip of outgroup
		w <- which(x1[, 7:8]==-root, arr.ind=TRUE)
		midH <- x1[w[1], 6] - x1[w[1], 3 + w[2]]
		m <- w[1]
	}
	
	# invert and lower nodes above rotation point
	.dropH <- function(i, delta) {
		stack <- integer(n)
		pos <- 1L
		stack[pos] <- i
		while (pos > 0) {
			i <- stack[pos]
			x1[i, 6] <<- x1[i, 6] - delta
			pos <- pos - 1L
			if (x1[i, 7] > 0) {
				pos <- pos + 1L
				stack[pos] <- x1[i, 7]
			}
			if (x1[i, 8] > 0) {
				pos <- pos + 1L
				stack[pos] <- x1[i, 8]
			}
		}
	}
	up <- integer(n) # pointers up tree
	w <- which(x1[, 7:8] > 0, arr.ind=TRUE)
	up[x1[, 7:8][w]] <- w[, "row"]
	remove <- logical(n) # replaced nodes
	x2 <- x1 # new rooted tree
	count <- n # row in x2
	# make new root node
	delta <- x1[m, 6] - midH
	x2[count, 4:10] <- c(x1[m, 4 + index[m]] - delta,
		delta,
		midH,
		x1[m, 7 + index[m]],
		count - 1,
		x1[m, 9 + index[m]],
		-1)
	if (up[m]) {
		while (up[m]) {
			count <- count - 1
			delta <- x1[m, 6] - midH
			remove[m] <- TRUE
			x2[count, 4:10] <- c(x1[m, 5 - index[m]],
				x1[up[m], 4 + index[up[m]]],
				midH - delta,
				x1[m, 8 - index[m]],
				count - 1,
				x1[m, 10 - index[m]],
				-1)
			if (x1[m, 8 - index[m]] > 0)
				.dropH(x1[m, 8 - index[m]], 2*delta)
			m <- up[m]
		}
		delta <- x1[m, 6] - midH
		x2[count, 5] <- sum(x1[m, 4:5])
	}
	remove[m] <- TRUE
	keep <- which(!remove)
	x2[count, 8] <- x1[m, 8 - index[m]]
	if (x2[count, 8] > 0)
		x2[count, 8] <- match(x2[count, 8], keep)
	x2[count, 10] <- x1[m, 10 - index[m]]
	if (x1[m, 8 - index[m]] > 0)
		.dropH(x1[m, 8 - index[m]], 2*delta)
	if (length(keep) > 0) {
		x2[1:(count - 1),] <- x1[keep,]
		w <- which(x2[1:(count - 1), 7:8] > 0)
		if (length(w) > 0)
			x2[1:(count - 1), 7:8][w] <- match(x2[1:(count - 1), 7:8][w], keep)
		w <- which(x2[n:count, 7] %in% keep)
		x2[n:count, 7][w] <- match(x2[n:count, 7][w], keep)
	}
	
	if (root > 0) {
		x2[, 6] <- x2[, 6] - min(x2[, 6] - x2[, 4], x2[, 6] - x2[, 5])
	}
	
	return(x2)
}

.applyMidpoints <- function(dend, dim) {
	# initialize a stack of maximum length (dim)
	stack <- vector("list", dim)
	visit <- logical(dim) # node already visited
	parent <- integer(dim) # index of parent node
	index <- integer(dim) # index in parent node
	pos <- 1L # current position in the stack
	stack[[pos]] <- dend
	while (pos > 0L) { # more nodes to visit
		if (visit[pos]) { # ascending tree
			visit[pos] <- FALSE # reset visit
			
			members <- sapply(stack[[pos]],
				function(x) {
					m <- attr(x, "members")
					if (is.null(m)) {
						return(1L)
					} else {
						return(m)
					}
				})
			
			l <- length(stack[[pos]])
			if (is.leaf(stack[[pos]][[1]]) && is.leaf(stack[[pos]][[l]])) {
				attr(stack[[pos]], "midpoint") <- (sum(members) - 1)/2
			} else if (is.leaf(stack[[pos]][[1]])) {
				attr(stack[[pos]], "midpoint") <- (sum(members[-l]) + attr(stack[[pos]][[l]], "midpoint"))/2
			} else if (is.leaf(stack[[pos]][[l]])) {
				attr(stack[[pos]], "midpoint") <- (attr(stack[[pos]][[1]], "midpoint") + sum(members[-l]))/2
			} else {
				attr(stack[[pos]], "midpoint") <- (sum(members[-l]) + attr(stack[[pos]][[1]], "midpoint") + attr(stack[[pos]][[l]], "midpoint"))/2
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

IdClusters <- function(myDistMatrix=NULL,
	method="UPGMA",
	cutoff=-Inf,
	showPlot=FALSE,
	type="clusters",
	myXStringSet=NULL,
	model=MODELS,
	collapse=0,
	reconstruct=FALSE,
	root=0,
	processors=1,
	verbose=TRUE) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	METHODS <- c("NJ","UPGMA", "ML", "complete", "single", "WPGMA", "inexact")
	method <- pmatch(method, METHODS)
	if (is.na(method) && !is.null(myDistMatrix)) {
		stop("Invalid method.  Choose either ML, NJ, complete, single, WPGMA, or UPGMA.")
	} else if (is.na(method)) {
		stop("Invalid method.")
	}
	if (method==-1 && !is.null(myDistMatrix)) {
		stop("Ambiguous method.  Choose either ML, NJ, complete, single, WPGMA, or UPGMA.")
	} else if (method==-1) {
		stop("Ambiguous method.")
	}
	if (!is.logical(reconstruct) &&
		!is.numeric(reconstruct))
		stop("reconstruct must be a logical or numeric.")
	if (is.numeric(reconstruct)) {
		if (reconstruct <= 0)
			stop("reconstruct must be be greater than zero.")
		if (reconstruct > 1)
			stop("reconstruct can be at most one.")
	}
	if (method==3 ||
		(type > 1 && reconstruct)) {
		if (length(model) < 1)
			stop("No model(s) specified.")
		if (!is.character(model))
			stop("model must be a character vector.")
		w <- which(!(model %in% MODELS))
		if (length(w) > 0) {
			submodels <- sub("([^+]*)(\\+G(\\d+))?", "\\1", model[w])
			if (!all(submodels %in% sub("\\+G(\\d+)$", "", MODELS)))
				stop(paste("Available models are:",
					paste(MODELS, collapse=", "),
					collapse=" "))
			rates <- as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", model[w]))
			if (any(is.na(rates)) || any(floor(rates)!=rates))
				stop("The number rates in the discrete Gamma distribution (i.e., +G4) should be an integer value.")
			if (any(rates > 10))
				stop("Up to 10 rates are allowed in the discrete Gamma distribution (i.e., +G10).")
			if (any(rates < 2))
				stop("A minimum of two rates are required for the discrete Gamma distribution (i.e., +G2).")
		}
		model <- unique(model)
	}
	if (!is.numeric(cutoff))
		stop("cutoff must be a numeric.")
	if (is.integer(cutoff))
		cutoff <- as.numeric(cutoff)
	if (!is.logical(showPlot))
		stop("showPlot must be a logical.")
	TYPES <- c("clusters", "dendrogram", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (method==7 && showPlot)
		stop("showPlot must be FALSE if method is 'inexact'")
	if (method==7 && type > 1)
		stop("type must be 'clusters' when method is 'inexact'")
	if (length(cutoff) > 1 && type > 1)
		warning("More than one cutoff specified when type is ", TYPES[type], ".")
	if (method==7 && any(cutoff < 0))
		stop("cutoff must be at least zero when method is 'inexact'.")
	if (method==7 && any(cutoff >= 1))
		stop("cutoff must be less than one when method is 'inexact'.")
	ASC <- TRUE
	if (is.unsorted(cutoff)) {
		if (is.unsorted(rev(cutoff))) {
			stop("cutoff must be sorted.")
		} else {
			ASC <- FALSE
		}
	}
	if (!is.numeric(collapse))
		stop("collapse must be a numeric.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors)!=processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	if (method != 7) {
		if (is(myDistMatrix, "matrix")) {
			dim <- dim(myDistMatrix)
			if (dim[2]!=dim[1])
				stop("myDistMatrix is not square.")
			dim <- dim[1]
		} else if (is(myDistMatrix, "dist")) {
			dim <- attr(myDistMatrix, "Size")
		} else {
			stop(paste("myDistMatrix must be a matrix for method '", METHODS[method], "'.", sep=""))
		}
		if (dim < 2)
			stop("myDistMatrix is too small.")
		if (typeof(myDistMatrix)=="integer")
			myDistMatrix[] <- as.numeric(myDistMatrix)
		
		if (!is.numeric(root))
			stop("root must be a numeric.")
		if (length(root) != 1)
			stop("root must be a single numeric.")
		if (floor(root) != root)
			stop("root must be an integer.")
		if (root < 0)
			stop("root must be at least 0.")
		if (root > dim)
			stop(paste("root cannot be greater than ", dim, ".", sep=""))
	}
	if (method == 3 ||
		(reconstruct && type > 1)) {
		if (is(myXStringSet, "DNAStringSet")) {
			typeX <- 1L
		} else if (is(myXStringSet, "RNAStringSet")) {
			typeX <- 2L
		} else {
			stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
		}
		if (length(myXStringSet)!=dim)
			stop("myDistMatrix must have as many rows as the number of sequences.")
		if (length(unique(width(myXStringSet))) != 1)
			stop("All sequences in myXStringSet must be the same width (aligned).")
		if (method==3 &&
			!is.null(attr(myDistMatrix, "correction")) &&
			attr(myDistMatrix, "correction") == "none")
			warning('myDistMatrix should probably have a correction for method="ML".')
	}
	
	if (method == 7) {
		if (is(myXStringSet, "DNAStringSet")) {
			typeX <- 1L
		} else if (is(myXStringSet, "RNAStringSet")) {
			typeX <- 2L
		} else if (is(myXStringSet, "AAStringSet")) {
			typeX <- 3L
		} else {
			stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
		}
		a <- vcountPattern("-", myXStringSet)
		if (any(a > 0))
			stop("Gap characters ('-') must be removed before inexact clustering.")
		a <- vcountPattern("+", myXStringSet)
		if (any(a > 0))
			stop("Mask characters ('+') must be removed before inexact clustering.")
		a <- vcountPattern(".", myXStringSet)
		if (any(a > 0))
			stop("Unknown characters ('.') must be removed before inexact clustering.")
		if (all(width(myXStringSet)==0L))
			stop("All sequences in myXStringSet are zero width.")
		
		if (verbose) {
			lastValue <- 0
			pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
		}
		if (typeX==3L) { # AAStringSet
			wordSize <- ceiling(log(100*quantile(width(myXStringSet), 0.99),
				.Call("alphabetSizeReducedAA",
					myXStringSet,
					0:19,
					PACKAGE="DECIPHER")))
			if (wordSize > 7)
				wordSize <- 7
			if (wordSize < 1)
				wordSize <- 1
			words <- 20^wordSize
		} else { # DNAStringSet or RNAStringSet
			wordSize <- ceiling(log(100*quantile(width(myXStringSet), 0.99),
				.Call("alphabetSize",
					myXStringSet,
					PACKAGE="DECIPHER")))
			if (wordSize > 15)
				wordSize <- 15
			if (wordSize < 1)
				wordSize <- 1
			words <- 4^wordSize
		}
		
		l <- length(myXStringSet)
		if (l==0)
			stop("myXStringSet contains no sequences.")
		lc <- length(cutoff)
		if (is.null(names(myXStringSet))) {
			dNames <- 1:l
		} else {
			dNames <- names(myXStringSet)
			w <- which(duplicated(dNames))
			if (length(w) > 0) {
				warning("Duplicated names of myXStringSet appended with index.")
				dNames[w] <- paste(dNames[w],
					w,
					sep="_")
			}
		}
		if (lc > 1) {
			cNames <- paste("cluster",
				gsub("\\.", "_", cutoff),
				METHODS[method],
				sep="")
		} else {
			cNames <- "cluster"
		}
		c <- matrix(0L,
			nrow=l,
			ncol=lc,
			dimnames=list(dNames,
				cNames))
		
		# identify duplicated sequences
		x <- selfmatch(myXStringSet)
		u <- unique(x)
		l <- length(u)
		t <- tabulate(x, length(x))[u]
		
		cutoff <- cutoff/2 # cluster radius is half the diameter
		# probability of chance match is (1 - cutoff)^N
		# N is wordSize - 1 to correct for chance matches
		threshold <- (1 - cutoff)^(wordSize - 1)
		for (i in seq_len(lc)) {
			if (!ASC && i > 1) {
				o <- u[order(c[u, i - 1],
					width(myXStringSet)[u],
					t, # frequency
					decreasing=TRUE)]
			} else {
				o <- u[order(width(myXStringSet)[u],
					t, # frequency
					decreasing=TRUE)]
			}
			
			if (typeX==3L) { # AAStringSet
				v <- .Call("enumerateSequenceAA",
					.subset(myXStringSet, o[1:ifelse(l > 999, 999, l)]),
					wordSize,
					PACKAGE="DECIPHER")
			} else { # DNAStringSet or RNAStringSet
				v <- .Call("enumerateSequence",
					.subset(myXStringSet, o[1:ifelse(l > 999, 999, l)]),
					wordSize,
					PACKAGE="DECIPHER")
			}
			v <- lapply(v,
				sort.int,
				method="radix")
			
			cluster_num <- 1L
			offset <- 0L
			c[o[cluster_num], i] <- cluster_num
			
			nGroups <- 1L
			seeds.reps <- v[1L]
			seeds.nums <- list(v[1L])
			seeds.seqs <- list(.subset(myXStringSet, o[1]))
			seeds.clust <- list(1L)
			seeds.index <- list(o[1])
			
			index <- 1L
			for (j in seq_along(o)[-1]) {
				if (verbose) {
					value <- round(((i - 1)*l + j)/lc/l, 2)
					if (value > lastValue) {
						lastValue <- value
						setTxtProgressBar(pBar, value)
					}
				}
				
				if (j %% 1000 == 0) {
					index <- 1L
					if (typeX==3L) { # AAStringSet
						v <- .Call("enumerateSequenceAA",
							.subset(myXStringSet, o[j:ifelse(j + 999 > l, l, j + 999)]),
							wordSize,
							PACKAGE="DECIPHER")
					} else { # DNAStringSet or RNAStringSet
						v <- .Call("enumerateSequence",
							.subset(myXStringSet, o[j:ifelse(j + 999 > l, l, j + 999)]),
							wordSize,
							PACKAGE="DECIPHER")
					}
					v <- lapply(v,
						sort.int,
						method="radix")
				} else {
					index <- index + 1L
				}
				
				if (!ASC && i > 1) {
					if (c[o[j], i - 1] != c[o[j - 1], i - 1]) {
						# different clusters in last cutoff
						offset <- offset + cluster_num
						cluster_num <- 1L
						c[o[j], i] <- cluster_num + offset
						nGroups <- 1L
						seeds.reps <- v[index]
						seeds.nums <- list(v[index])
						seeds.seqs <- list(.subset(myXStringSet, o[j]))
						seeds.clust <- list(cluster_num)
						seeds.index <- list(o[j])
						next
					} else if (c[o[j], i] > 0) { # cluster pre-assigned
						c[o[j], i] <- c[c[o[j], i], i]
						next
					}
				}
				
				# Split the sequences into "alignable" groups
				m <- .Call("matchListsDual",
					v[index],
					seeds.reps,
					FALSE, # verbose
					NULL, # pBar
					processors,
					PACKAGE="DECIPHER")
				
				# Determine to which group does the sequence belong
				# expected number of occurrences of a random word
				probs <- lengths(seeds.reps)/words
				probs <- ifelse(probs > 1, 1, probs)
				# probability of >= `q` matches in `size` trials
				prob <- 1 - pbinom(q=m*length(v[[index]]),
					size=length(v[[index]]),
					prob=probs)
				group <- which(prob < 0.01) # < 1% likelihood by chance
				
				if (length(group)==0) {
					# form a new group
					cluster_num <- cluster_num + 1L
					c[o[j], i] <- cluster_num
					nGroups <- nGroups + 1L
					seeds.reps[nGroups] <- v[index]
					seeds.nums[[nGroups]] <- v[index]
					seeds.seqs[[nGroups]] <- .subset(myXStringSet, o[j])
					seeds.clust[[nGroups]] <- cluster_num
					seeds.index[[nGroups]] <- o[j]
				} else { # part of an existing group
					if (length(group) > 1)
						group <- group[which.max(m[group])]
					
					# Try to place into a cluster based on k-mers
					m <- .Call("matchListsDual",
						v[index],
						seeds.nums[[group]],
						FALSE, # verbose
						NULL, # pBar
						processors,
						PACKAGE="DECIPHER")
					
					# percent identity >= threshold
					w <- which(m >= threshold[i])
					if (length(w) > 0) {
						w <- w[which.max(m[w])]
						c[o[j], i] <- seeds.clust[[group]][w] + offset
						if (!ASC && i < lc) {
							cols <- (i + 1):lc
							cols <- cols[which(m[w] >= threshold[cols])]
							if (length(cols) > 0) # assign forward
								c[o[j], cols] <- seeds.index[[group]][w]
						}
					} else { # need to align to determine cluster
						pattern <- .subset(myXStringSet, o[j])
						
						weights <- m^3 # emphasize closest sequences
						weights <- weights/mean(weights)
						if (any(!is.finite(weights)))
							weights <- 1
						
						# pairwise alignment
						temp <- AlignProfiles(pattern,
							seeds.seqs[[group]],
							s.weight=weights,
							processors=processors)
						d <- .Call("distMatrix",
							temp,
							typeX,
							FALSE, # includeTerminalGaps
							FALSE, # penalizeGapGapMatches
							TRUE, # penalizeGapLetterMatches
							FALSE, # full matrix
							1L, # matrix
							FALSE, # verbose
							NULL, # pBar
							processors,
							PACKAGE="DECIPHER")[-1L]
						w <- which(d <= cutoff[i])
						if (length(w) > 0) {
							w <- w[which.min(d[w])]
							c[o[j], i] <- seeds.clust[[group]][w] + offset
							if (!ASC && i < lc) {
								cols <- (i + 1):lc
								cols <- cols[which(d[w] <= cutoff[cols])]
								if (length(cols) > 0) # assign forward
									c[o[j], cols] <- seeds.index[[group]][w]
							}
						} else { # form a new cluster
							cluster_num <- cluster_num + 1L
							c[o[j], i] <- cluster_num + offset
							seeds.nums[[group]] <- c(v[index],
								seeds.nums[[group]])
							seeds.seqs[[group]] <- temp
							seeds.clust[[group]] <- c(cluster_num,
								seeds.clust[[group]])
							seeds.index[[group]] <- c(o[j],
								seeds.index[[group]])
						}
					}
				}
			}
			
			c[, i] <- c[x, i]
		}
		myClusters <- as.data.frame(c)
	} else {
		# initialize a progress bar
		if (verbose) {
			if (method==3) {
				cat("Constructing initial neighbor-joining tree:\n")
				flush.console()
			}
			pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		} else {
			pBar <- NULL
		}
		
		myClusters <- .Call("cluster",
			myDistMatrix,
			ifelse(method==3, -Inf, cutoff[1]),
			ifelse(method==3, 1L, method),
			ifelse(is(myDistMatrix, "matrix"),
				dim,
				-dim),
			verbose,
			pBar,
			processors,
			PACKAGE="DECIPHER")
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
		}
		
		if (method==3 ||
			(reconstruct && type > 1)) {
			myClusters <- .reorderClusters(myClusters,
				all=method != 3)
			
			m <- matrix(NA,
				nrow=length(model),
				ncol=10,
				dimnames=list(model,
					c("FreqA", "FreqC", "FreqG", "FreqT",
						"A2G", "C2T", "alpha",
						"-LnL", "AICc", "BIC")))
			
			if (!(length(model)==1 && model=="JC69")) {
				a <- apply(consensusMatrix(myXStringSet),
					2,
					function(x)
						sum(x[1:15] > 0)) > 1
				N <- sum(a)
				
				for (i in seq_along(model)) {
					m[model[i],] <- .optimizeModel(myClusters,
						model[i],
						myXStringSet,
						N,
						processors=processors)
					if (verbose)
						cat("\n", model[i],
							":", paste(rep(" ",
									max(nchar(rownames(m))) - nchar(model[i]) + 1),
								collapse=""),
							"-ln(L)=", round(m[model[i], "-LnL"], 0),
							", AICc=", round(m[model[i], "AICc"], 0),
							", BIC=", round(m[model[i], "BIC"], 0),
							sep="")
				}
			}
			
			if (length(model) > 1) { # choose the best model
				w <- which.min(m[,"BIC"])
				m <- m[w,, drop=FALSE]
				model <- model[w]
				if (verbose)
					cat("\n\nThe selected model was:  ",
						model,
						"\n",
						sep="")
			}
			
			.giveParams <- function(model_params) {
				w <- which(is.na(model_params))
				if (length(w) > 0) {
					model_params[w] <- c(0.25, 0.25, 0.25, 0.25, 1, 1, NA)[w]
				}
				if (is.na(model_params[7])) {
					model_params <- c(model_params[1:6], 1, 1)
				} else {
					model_params <- c(model_params[1:6], .rates(model_params[7], as.integer(sub("([^+]*)(\\+G(\\d+))?", "\\3", rownames(m)))))
				}
			}
			model_params <- .giveParams(as.numeric(m[1:7]))
			
			if (method==3) {
				# given myClusters return adjusted heights
				adjustTreeHeights <- function(myClusters) {
					cumHeight <- numeric(max(myClusters[, 3]))
					for (i in 1:dim(myClusters)[1]) {
						if (myClusters[i, 1] < 0 && myClusters[i, 2] < 0) {
							cumHeight[myClusters[i, 3]] <- max(myClusters[i, 4], myClusters[i, 5])
							myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
						} else if (myClusters[i, 1] > 0 && myClusters[i, 2] > 0) {
							cumHeight[myClusters[i, 3]] <- max(myClusters[i, 4] + cumHeight[myClusters[i, 1]],
								myClusters[i, 5] + cumHeight[myClusters[i, 2]])
							myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
						} else if (myClusters[i, 1] > 0) {
							cumHeight[myClusters[i, 3]] <- cumHeight[myClusters[i, 1]] + myClusters[i, 4]
							if (myClusters[i, 5] > cumHeight[myClusters[i, 3]])
								cumHeight[myClusters[i, 3]] <- myClusters[i, 5]
							myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
						} else {
							cumHeight[myClusters[i, 3]] <- cumHeight[myClusters[i, 2]] + myClusters[i, 5]
							if (myClusters[i, 4] > cumHeight[myClusters[i, 3]])
								cumHeight[myClusters[i, 3]] <- myClusters[i, 4]
							myClusters[i, 6] <- cumHeight[myClusters[i, 3]]
						}
					}
					
					myClusters <- .Call("adjustHeights",
						myClusters,
						PACKAGE="DECIPHER")
					return(myClusters)
				}
				
				# print progress of likelihood maximization
				.startLnL <- Inf
				.bestLnL <- Inf
				.NNIs <- 0L
				if (verbose) {
					cat("\nMaximizing Likelihood of Tree:\n")
					flush.console()
					printLine <- function(value, percentComplete) {
						cat("\r-ln(Likelihood) = ",
							formatC(round(value, 0),
								digits=0,
								format="f"),
							" (",
							formatC(abs(round(100*(value - .startLnL)/.startLnL,
									2)),
								digits=2,
								format="f"),
							"% improvement), ",
							.NNIs,
							ifelse(.NNIs==1, " NNI  ", " NNIs "),
							sep="")
						flush.console()
						invisible(value)
					}
				}
				
				# given branch lengths return -LnL
				maximizeLikelihood <- function(x, branches=integer(), lengths=numeric()) {
					myClusters[, 4:5] <- x
					LnL <- .Call("clusterML",
						myClusters,
						myXStringSet,
						model_params,
						branches,
						lengths,
						0,
						1L,
						processors,
						PACKAGE="DECIPHER")
					
					w <- which.min(LnL)
					if (LnL[w] < .bestLnL) {
						if (w > 1) {
							x[branches[w - 1]] <- lengths[w - 1]
							params <<- x
						} else {
							params <<- x
						}
						.bestLnL <<- LnL[w]
						if (verbose) {
							if (is.infinite(.startLnL))
								.startLnL <<- LnL[1]
							printLine(LnL[w])
						}
					}
					
					return(LnL)
				}
				
				maximizeLikelihood2 <- function(myClusters, NNIs, tol) {
					LnL <- .Call("clusterML",
						myClusters,
						myXStringSet,
						model_params,
						integer(),
						numeric(),
						0,
						1L,
						processors,
						PACKAGE="DECIPHER")
					
					if (LnL < .bestLnL - tol) {
						.bestLnL <<- LnL
						if (verbose) {
							.NNIs <<- NNIs
							if (is.infinite(.startLnL))
								.startLnL <<- LnL
							printLine(LnL)
						}
					}
					
					return(LnL)
				}
				
				# maximize likelihood of tree
				index <- rep(TRUE, 2*dim(myClusters)[1])
				repeat {
					currentLnL <- .bestLnL
					currentNNIs <- .NNIs
					params <- as.numeric(myClusters[, 4:5])
					tempClusters <- paste(myClusters[, 7], myClusters[, 8])
					.simultaneousBrent(maximizeLikelihood,
						ifelse(index, 0, params),
						params,
						ifelse(index, 10*params, params))
					myClusters[, 4:5] <- params
					
					myClusters <- .NNI(myClusters,
						.bestLnL,
						.NNIs,
						maximizeLikelihood2)
					
					if (abs(.bestLnL - currentLnL) < 1e0 &&
						currentNNIs==.NNIs) {
						temp_params <- model_params
					} else if (model != "JC69") {
						m[1,] <- .optimizeModel(myClusters,
							rownames(m),
							myXStringSet,
							N,
							processors=processors)
						temp_params <- .giveParams(as.numeric(m[1:7]))
					} else {
						temp_params <- .giveParams(as.numeric(m[1:7]))
					}
					
					if ((abs(.bestLnL - currentLnL) < 1e0 || # negligible improvement in likelihood
						currentNNIs==.NNIs) && # no new nearest neighbor interchanges (NNIs)
						all(abs(temp_params - model_params) < 0.01)) # negligible change in model_params
						break
					model_params <- temp_params
					index <- rep(!(paste(myClusters[, 7], myClusters[, 8]) %in% tempClusters), 2)
				}
				
				myClusters <- .reorderClusters(myClusters, all=TRUE)
				myClusters <- adjustTreeHeights(myClusters)
				myClusters <- .Call("reclusterNJ",
					myClusters,
					cutoff[1],
					PACKAGE="DECIPHER")
				
				if (verbose) {
					printLine(.bestLnL, 100)
					cat("\n")
					flush.console()
				}
			}
			
			if (verbose) {
				params <- formatC(round(m, 3),
					digits=3,
					format="f")
				cat(ifelse(grepl("NA", params[1], fixed=TRUE),
						"",
						"\nModel parameters:"),
					ifelse(grepl("NA", params[1], fixed=TRUE),
						"",
						ifelse(grepl("T92", rownames(m), fixed=TRUE),
							paste("\nFrequency(A) = Frequency(T) = ", params[1],
								"\nFrequency(C) = Frequency(G) = ", params[2],
								sep=""),
							paste("\nFrequency(A) = ", params[1],
								"\nFrequency(C) = ", params[2],
								"\nFrequency(G) = ", params[3],
								"\nFrequency(T) = ", params[4],
								sep=""))),
					ifelse(grepl("NA", params[5], fixed=TRUE),
						"",
						ifelse(grepl("TN93", rownames(m), fixed=TRUE),
							paste("\nRate A <-> G = ", params[5],
								"\nRate C <-> T = ", params[6],
								"\nTransversion rates = 1",
								sep=""),
							paste("\nTransition rates = ", params[5],
								"\nTransversion rates = 1",
								sep=""))),
					ifelse(grepl("NA", params[7], fixed=TRUE),
						"",
						paste("\nAlpha = ", params[7], sep="")),
					ifelse(grepl("NA", params[1], fixed=TRUE),
						"",
						"\n"),
					sep="")
			}
			params <- m[, 1:7]
		} else {
			model <- NULL
			params <- NULL
		}
		
		if (showPlot || type > 1) {
			if (method==1 ||
				method==3 ||
				root > 0)
				myClusters <- .root(myClusters, root)
			
			# create a dendrogram object
			myClustersList <- list()
			dNames <- labels(myDistMatrix)
			if (is.null(dNames)) {
				myClustersList$labels <- 1:(dim(myClusters)[1] + 1)
			} else {
				if (is(dNames, "list")) {
					myClustersList$labels <- dNames[[1]]
				} else {
					myClustersList$labels <- dNames
				}
				
				w <- which(duplicated(myClustersList$labels))
				if (length(w) > 0) {
					warning("Duplicated labels in dendrogram appended with index.")
					myClustersList$labels[w] <- paste(myClustersList$labels[w],
						w,
						sep="_")
				}
			}
			
			myClustersList$merge <- matrix(myClusters[,7:8], ncol=2)
			myClustersList$height <- matrix(myClusters[,6], ncol=1)
			myClustersList$lengths <- matrix(myClusters[,4:5], ncol=2)
			if (dim > 100) {
				fontSize <- 0.6
			} else if (dim > 70) {
				fontSize <- 0.7
			} else if (dim > 40) {
				fontSize <- 0.8
			} else {
				fontSize <- 0.9
			}
			if (dim > 300) {
				leaves <- "none"
			} else {
				leaves <- "perpendicular"
			}
			
			if (reconstruct && type > 1) {
				# perform ancestral state reconstruction
				states <- .Call("clusterML",
					myClusters,
					myXStringSet,
					model_params,
					integer(),
					numeric(),
					reconstruct,
					typeX,
					processors,
					PACKAGE="DECIPHER")
				d <- to.dendrogram(myClustersList,
					states)
				myXStringSet <- unname(as.character(myXStringSet))
				d <- rapply(d,
					function(x) {
						attr(x, "state") <- myXStringSet[x]
						x
					},
					how="replace")
			} else {
				d <- to.dendrogram(myClustersList)
			}
			
			# convert bifurcating tree to multifurcating
			if (collapse >= 0)
				d <- .collapse(d, collapse, dim)
			
			attr(d, "method") <- METHODS[method]
			attr(d, "model") <- model
			attr(d, "parameters") <- params
			
			# specify the order of clusters that
			# will match the plotted dendrogram
			orderDendrogram <- order.dendrogram(d)
			
			c <- .organizeClusters(myClusters, myClustersList$labels, orderDendrogram)
			
			if (is.finite(cutoff[1])) {
				# create a visibily different vector of colors
				cl <- colors()
				v1 <- c(117,254,73,69,152,51,26,450,503,596,652,610,563,552,97)
				r <- cl[v1]
				
				# color edges by cluster
				.colEdge <- function(dend) {
					# initialize a stack of maximum length (dim)
					stack <- vector("list", dim)
					visit <- logical(dim) # node already visited
					parent <- integer(dim) # index of parent node
					index <- integer(dim) # index in parent node
					pos <- 1L # current position in the stack
					stack[[pos]] <- dend
					while (pos > 0L) { # more nodes to visit
						if (visit[pos]) { # ascending tree
							visit[pos] <- FALSE # reset visit
							
							for (i in seq_along(stack[[pos]])) {
								if (is.leaf(stack[[pos]][[i]])) {
									a <- attributes(stack[[pos]][[i]])
									num <- c$cluster[which(myClustersList$labels==as.character(a$label))]
									attr(stack[[pos]][[i]], "edgePar") <- list(col=r[num %% 15 + 1])
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
				d <- .colEdge(d)
				
				.reorder <- function(dend) {
					# initialize a stack of maximum length (dim)
					stack <- vector("list", dim)
					visit <- logical(dim) # node already visited
					parent <- integer(dim) # index of parent node
					index <- integer(dim) # index in parent node
					pos <- 1L # current position in the stack
					stack[[pos]] <- dend
					while (pos > 0L) { # more nodes to visit
						if (visit[pos]) { # ascending tree
							visit[pos] <- FALSE # reset visit
							
							members <- lapply(stack[[pos]], unlist)
							# sort tree by ascending cluster number
							o <- sort.list(sapply(members,
									function(x)
										min(c[x, 1])))
							stack[[pos]][] <- stack[[pos]][o]
							
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
				d <- .reorder(d)
			}
			
			# add midpoints to the tree
			d <- .applyMidpoints(d, dim)
		}
		if (type==1 || type==3) {
			dNames <- labels(myDistMatrix)
			if (is.null(dNames)) {
				dNames <- 1:(dim(myClusters)[1] + 1)
			} else {
				if (is(dNames, "list"))
					dNames <- dNames[[1]]
				
				w <- which(duplicated(dNames))
				if (length(w) > 0) {
					warning("Duplicated labels in myDistMatrix appended with index.")
					dNames[w] <- paste(dNames[w],
						w,
						sep="_")
				}
			}
			
			if (type==1) # do not number clusters by order of appearance
				c <- .organizeClustersFast(myClusters, dNames)
			
			if (length(cutoff) > 1) {
				names(c) <- paste("cluster",
					gsub("\\.", "_", cutoff[1]),
					METHODS[method],
					sep="")
				for (i in 2:length(cutoff)) {
					if (method==1 ||
						method==3 ||
						root > 0) {
						myClusters <- .Call("reclusterNJ",
							myClusters,
							cutoff[i],
							PACKAGE="DECIPHER")
					} else { # ultrametric
						myClusters <- .Call("reclusterUPGMA",
							myClusters,
							cutoff[i],
							PACKAGE="DECIPHER")
					}
					x <- .organizeClustersFast(myClusters, dNames)
					if ((method==1 ||
						method==3 ||
						root > 0) &&
						!ASC) # ensure clusters are subsets
						x[, 1] <- .splitClusters(x[, 1], c[, dim(c)[2]])
					names(x) <- paste("cluster",
						gsub("\\.", "_", cutoff[i]),
						METHODS[method],
						sep="")
					c <- cbind(c, x)
				}
			}
			myClusters <- c
		}
		
		if (showPlot)
			plot(d,
				horiz=FALSE,
				leaflab=leaves,	
				nodePar=list(lab.cex=fontSize, pch = NA))
	}
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	if (type==1) {
		return(myClusters)
	} else if (type==2) {
		return(d)
	} else {
		return(list(myClusters, d))
	}
}

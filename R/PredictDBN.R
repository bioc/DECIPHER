.getWeights <- function(weight,
	myXStringSet,
	verbose,
	processors) {
	
	l <- length(myXStringSet)
	
	if (length(weight) == 1) {
		if (l == 1L || is.numeric(weight)) {
			weight <- rep(1, l)
		} else if (is.na(weight)) {
			if (verbose) {
				cat("Determining distance matrix based on alignment:\n")
				flush.console()
			}
			d <- DistanceMatrix(myXStringSet,
				type="dist",
				penalizeGapLetterMatches=FALSE,
				verbose=verbose,
				processors=processors)
			
			if (verbose) {
				cat("Determining sequence weights:\n")
				flush.console()
			}
			suppressWarnings(guideTree <- TreeLine(myDistMatrix=d,
				method="UPGMA",
				verbose=verbose,
				processors=processors))
			
			.weights <- function(guideTree) {
				# initialize a stack of maximum length (l)
				stack <- vector("list", l)
				visit <- logical(l) # node already visited
				parent <- integer(l) # index of parent node
				index <- integer(l) # index in parent node
				pos <- 1L # current position in the stack
				stack[[pos]] <- guideTree
				while (pos > 0L) { # more nodes to visit
					if (visit[pos]) { # ascending tree
						visit[pos] <- FALSE # reset visit
						dend <- stack[[pos]]
						
						# align subtrees
						treeLengths <- numeric(length(dend))
						inherit <- attr(dend, "inherit")
						if (length(dend) > 1) {
							for (i in seq_len(length(dend))) {
								h <- attr(dend[[i]], "treeLength")
								if (!is.null(h))
									treeLengths[i] <- h
							}
							
							h <- attr(dend, "height")
							members <- vector("list", length(dend))
							for (i in seq_len(length(dend))) {
								m <- unlist(dend[i])
								treeLengths[i] <- treeLengths[i] + h - height[m[1]]
								weight[m] <<- weight[m] + (h - height[m])/length(m)
								height[m] <<- h
								members[[i]] <- m
							}
						} else if (is.leaf(dend)) {
							height[unlist(dend)] <- attr(dend, "height")
						} else { # inherit from subtree
							treeLengths[1] <- attr(dend[[1]], "treeLength")
						}
						
						attr(stack[[pos]], "treeLength") <- sum(treeLengths)
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
				
				return(attr(stack[[1]], "treeLength"))
			}
			
			weight <- height <- numeric(l)
			.weights(guideTree)
			weight <- weight/mean(weight)
		} else {
			stop("weight must be a numeric or NA.")
		}
	} else {
		if (length(weight) != l)
			stop("Length of weight must equal one or the length of the myXStringSet.")
		if (any(is.na(weight)))
			stop("weight must be a numeric.")
		if (!isTRUE(all.equal(1, mean(weight))))
			stop("The mean of weight must be 1.")
	}
	
	return(weight)
}

.parseDBN <- function(dbn) {
	s <- strsplit(dbn, "")[[1]]
	n <- sum(s != "." & s != "-")
	c1 <- c2 <- c3 <- integer(n/2)
	i <- 0L
	parens <- square <- curly <- straight <- integer()
	for (j in seq_along(s)) {
		if (s[j] != "." && s[j] != "-") {
			if (s[j] == "(") {
				parens <- c(parens, j)
			} else if (s[j] == ")") {
				i <- i + 1L
				c1[i] <- parens[length(parens)]
				c2[i] <- j
				length(parens) <- length(parens) - 1L
			} else if (s[j] == "[") {
				square <- c(square, j)
			} else if (s[j] == "]") {
				i <- i + 1L
				c1[i] <- square[length(square)]
				c2[i] <- j
				c3[i] <- 1L
				length(square) <- length(square) - 1L
			} else if (s[j] == "{") {
				curly <- c(curly, j)
			} else if (s[j] == "}") {
				i <- i + 1L
				c1[i] <- curly[length(curly)]
				c2[i] <- j
				c3[i] <- 2L
				length(curly) <- length(curly) - 1L
			} else if (s[j] == "<") {
				straight <- c(straight, j)
			} else if (s[j] == ">") {
				i <- i + 1L
				c1[i] <- straight[length(straight)]
				c2[i] <- j
				c3[i] <- 3L
				length(straight) <- length(straight) - 1L
			} else {
				stop("Unrecognized character in structure: ", s[j])
			}
		}
	}
	
	o <- order(c3, c1, c2)
	matrix(c(c1[o], c2[o], c3[o]),
		ncol=3,
		dimnames=list(NULL, c("(", ")", "order")))
}

PredictDBN <- function(myXStringSet,
	type="states",
	minOccupancy=0.4,
	impact=c(1, 1.2, 0.4, -1),
	avgProdCorr=1,
	slope=2,
	shift=1.3,
	threshold=0.4,
	pseudoknots=1,
	weight=NA,
	useFreeEnergy=TRUE,
	deltaGrules=NULL,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!(is(myXStringSet, "DNAStringSet") ||
		is(myXStringSet, "RNAStringSet")))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	l <- length(myXStringSet)
	if (l == 0L)
		stop("myXStringSet must contain at least one sequence.")
	if (!is.logical(useFreeEnergy))
		stop("useFreeEnergy must be a logical.")
	if (l == 1L && !useFreeEnergy)
		stop("useFreeEnergy must be TRUE when there is only one sequence in myXStringSet.")
	u <- unique(width(myXStringSet))
	if (length(u) != 1)
		stop("Sequences in myXStringSet must be the same width (aligned).")
	TYPES <- c("states", "pairs", "scores", "structures", "search", "evidence")
	if (length(type) == 0)
		stop("No type specified.")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (u == 0L) { # zero width sequences
		if (type == 1L) {
			return("")
		} else if (type == 2L) {
			return(matrix(nrow=0, ncol=3, dimnames=list(NULL, c("(", ")", "order"))))
		} else if (type == 3L) {
			return(matrix(nrow=3, ncol=0, dimnames=list(c(".", "(", ")"), NULL)))
		} else if (type == 6L) {
			return(matrix(nrow=0, ncol=3, dimnames=list(NULL, c("(", ")", "score"))))
		} else {
			x <- matrix(0, nrow=3, ncol=0, dimnames=list(c(".", "(", ")"), NULL))
			x <- lapply(seq_len(l),
				function(...)
					return(x))
			return(x)
		}
	}
	if (!is.numeric(minOccupancy))
		stop("minOccupancy must be a numeric.")
	if (minOccupancy < 0 || minOccupancy > 1)
		stop("minOccupancy must be between zero and one.")
	if (!is.numeric(impact) || length(impact) != 4)
		stop("impact must be a numeric vector of length four.")
	if (impact[4] > 0)
		stop("The impact of inconsistent pairs (element four) can be at most zero.")
	if (any(impact[1:2] < 0))
		stop("The impact of A/T and G/C base pairs must be at least zero.")
	if (!is.numeric(avgProdCorr))
		stop("avgProdCorr must be a numeric.")
	if (avgProdCorr < 0)
		stop("avgProdCorr must be at least zero.")
	if (!is.numeric(slope))
		stop("slope must be a numeric.")
	if (slope <= 0)
		stop("slope must be greater than zero.")
	if (!is.numeric(shift))
		stop("shift must be a numeric.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold < 0 || threshold > 1)
		stop("threshold must be between zero and one.")
	if (!is.numeric(pseudoknots))
		stop("pseudoknots must be a numeric.")
	if (pseudoknots < 0 || pseudoknots > 3)
		stop("pseudoknots must be between zero and three.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors) != processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	weight <- .getWeights(weight,
		myXStringSet,
		verbose,
		processors)
	
	if (useFreeEnergy) {
		# initialize a progress bar
		if (verbose) {
			time.1 <- Sys.time()
			cat("Computing Free Energies:\n")
			pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
		}
		
		if (is.null(deltaGrules)) {
			data("deltaGrulesRNA", envir=environment(), package="DECIPHER")
			deltaGrules <- deltaGrulesRNA
		} else {
			if (!is.numeric(deltaGrules))
				stop("deltaGrules must be numeric.")
			if (length(deltaGrules) != 390625L)
				stop("deltaGrules must be of dimensions 5 x 5 x 5 x 5 x 5 x 5 x 5 x 5.")
		}
		
		pos <- as.matrix(myXStringSet)
		pos <- lapply(seq_len(nrow(pos)),
			function(x)
				which(!(pos[x,] %in% c("-", "."))))
		
		noGaps <- RemoveGaps(myXStringSet)
		J <- W <- rep(list(integer()), u)
		for (i in seq_along(myXStringSet)) {
			pals <- findPalindromes(noGaps[[i]],
				min.armlength=4,
				max.looplength=500,
				min.looplength=3,
				max.mismatch=1,
				allow.wobble=TRUE)
			
			dna <- DNAStringSet(pals)
			s1 <- as.character(subseq(dna, 1L, 1L))
			e1 <- as.character(subseq(dna, -1L))
			w <- which((s1 == "A" & e1 == "T") |
				(s1 == "T" & e1 == "A") |
				(s1 == "C" & e1 == "G") |
				(s1 == "G" & e1 == "C") |
				(s1 == "G" & e1 == "T") |
				(s1 == "T" & e1 == "G"))
			pals <- pals[w]
			dna <- dna[w]
			
			arms <- palindromeArmLength(dna,
				max.mismatch=0,
				allow.wobble=TRUE)
			arms1 <- palindromeArmLength(dna,
				max.mismatch=1,
				allow.wobble=TRUE)
			w <- which(arms + 1L != arms1)
			if (length(w) > 0)
				arms[w] <- arms1[w] # mismatch not at end
			max_arms <- as.integer((width(pals) - 3)/2)
			w <- which(arms > max_arms)
			if (length(w) > 0)
				arms[w] <- max_arms[w]
			
			s1 <- start(pals)
			s2 <- start(pals) + arms - 1L
			e1 <- end(pals) - arms + 1L
			e2 <- end(pals)
			p <- pos[[i]]
			
			w <- which(width(dna) - 3L > 2L*arms)
			arms[w] <- arms[w] + 1L # include end bases
			
			dG <- .Call("calculateHairpinDeltaG",
				dna,
				arms,
				deltaGrulesRNA,
				PACKAGE="DECIPHER")
			o <- order(dG)
			o <- o[dG[o] < 0]
			dG <- dG/mean(dG[o][seq_len(sum(cumsum(arms[o]) < length(p)))])
			dG <- 1 - exp(-dG) # convert to probability
			
			for (j in seq_along(o)) {
				p1 <- s1[o[j]]:s2[o[j]]
				p2 <- e2[o[j]]:e1[o[j]]
				p1 <- p[p1]
				p2 <- p[p2]
				for (k in seq_along(p1)) {
					m <- .Call("firstHit",
						p2[k],
						J[[p1[k]]],
						PACKAGE="DECIPHER")
					if (is.na(m)) {
						w <- .Call("multiMatchLower",
							J[[p1[k]]],
							p2[k],
							length(J[[p1[k]]]),
							PACKAGE="DECIPHER")
						if (length(w) == 0L) {
							J[[p1[k]]] <- c(p2[k], J[[p1[k]]])
							W[[p1[k]]] <- c(dG[o[j]]*weight[i], W[[p1[k]]])
						} else {
							w <- seq_len(w)
							J1 <- J[[p1[k]]]
							x <- logical(length(J1))
							x[w] <- TRUE
							y <- !x
							J2 <- J1[y]
							J1 <- J1[x]
							J[[p1[k]]] <- c(J1, p2[k], J2)
							W1 <- W[[p1[k]]]
							W2 <- W1[y]
							W1 <- W1[x]
							W[[p1[k]]] <- c(W1, dG[o[j]]*weight[i], W2)
						}
					} else {
						W[[p1[k]]][m] <- W[[p1[k]]][m] + dG[o[j]]*weight[i]
					}
				}
			}
			
			if (verbose)
				setTxtProgressBar(pBar, i/l)
		}
		
		if (l > 1L) {
			t <- TerminalChar(myXStringSet)
			w <- numeric(u)
			for (i in seq_along(myXStringSet))
				if (t[i, 3L] > 0)
					w[(t[i, 1L] + 1):(u - t[i, 2L])] <- w[(t[i, 1L] + 1):(u - t[i, 2L])] + weight[i]
			for (i in seq_along(J))
				W[[i]] <- W[[i]]/pmin(w[J[[i]]], w[i])
		}
		patterns <- list(J, W)
		
		if (verbose) {
			close(pBar)
			time.2 <- Sys.time()
			cat("\n")
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\n")
		}
	} else {
		patterns <- NULL
	}
	
	# initialize a progress bar
	if (verbose) {
		time.1 <- Sys.time()
		cat("Predicting RNA Secondary Structures:\n")
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
	} else {
		pBar <- NULL
	}
	
	ans <- .Call("predictDBN",
		myXStringSet,
		type,
		minOccupancy,
		impact,
		avgProdCorr,
		slope,
		shift,
		weight,
		pseudoknots,
		threshold,
		patterns,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (type == 2L) {
		ans <- .parseDBN(ans)
	} else if (type == 3L) {
		rownames(ans) <- c(".", "(", ")")
	} else if (type == 4L) {
		ans <- lapply(ans,
			function(x) {
				rownames(x) <- c(".", "(", ")")
				return(x)
			})
	} else if (type == 6L) {
		colnames(ans) <- c("(", ")", "score")
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(ans)
}

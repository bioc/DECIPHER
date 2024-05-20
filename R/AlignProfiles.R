AlignProfiles <- function(pattern,
	subject,
	p.weight=1,
	s.weight=1,
	p.struct=NULL,
	s.struct=NULL,
	perfectMatch=2,
	misMatch=-1,
	gapOpening=-12,
	gapExtension=-3,
	gapPower=-1,
	terminalGap=-4,
	restrict=c(-1000, 2, 10),
	anchor=0.7,
	normPower=c(1, 0),
	standardize=TRUE,
	substitutionMatrix=NULL,
	structureMatrix=NULL,
	processors=1) {
	
	# error checking
	if (is(pattern, "DNAStringSet")) {
		type <- 1L
		if (!is(subject, "DNAStringSet"))
			stop("pattern and subject must be of the same class.")
	} else if (is(pattern, "RNAStringSet")) {
		type <- 2L
		if (!is(subject, "RNAStringSet"))
			stop("pattern and subject must be of the same class.")
	} else if (is(pattern, "AAStringSet")) {
		type <- 3L
		if (!is(subject, "AAStringSet"))
			stop("pattern and subject must be of the same class.")
	} else {
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	if (length(subject) < 1)
		stop("At least one sequence is required in the subject.")
	w.p <- unique(width(pattern))
	if (length(w.p) != 1)
		stop("Sequences in pattern must be the same width (aligned).")
	w.s <- unique(width(subject))
	if (length(w.s) != 1)
		stop("Sequences in subject must be the same width (aligned).")
	if (!is.numeric(p.weight))
		stop("p.weight must be a numeric.")
	if (length(p.weight) != 1 && length(p.weight) != length(pattern))
		stop("Length of p.weight must equal one or the length of the pattern.")
	if (length(p.weight) == 1) {
		p.weight <- rep(1, length(pattern))
	} else {
		if (!isTRUE(all.equal(1, mean(p.weight))))
			stop("The mean of p.weight must be 1.")
	}
	if (!is.numeric(s.weight))
		stop("s.weight must be a numeric.")
	if (length(s.weight) != 1 && length(s.weight) != length(subject))
		stop("Length of s.weight must equal one or the length of the subject.")
	if (length(s.weight) == 1) {
		s.weight <- rep(1, length(subject))
	} else {
		if (!isTRUE(all.equal(1, mean(s.weight))))
			stop("The mean of s.weight must be 1.")
	}
	if (is.null(p.struct) != is.null(s.struct))
		stop("Both p.struct and s.struct must be specified.")
	if (!is.null(p.struct)) {
		if (is.matrix(p.struct)) {
			if (dim(p.struct)[2] != w.p)
				stop("The number of columns in p.struct does not match the width of the pattern.")
		} else if (is.list(p.struct)) {
			if (length(p.struct) != length(pattern))
				stop("p.struct is not the same length as the pattern.")
		} else {
			stop("p.struct must be a matrix or list.")
		}
	}
	if (!is.null(s.struct)) {
		if (is.matrix(s.struct)) {
			if (dim(s.struct)[2] != w.s)
				stop("The number of columns in s.struct does not match the width of the subject.")
		} else if (is.list(s.struct)) {
			if (length(s.struct) != length(subject))
				stop("s.struct is not the same length as the subject.")
		} else {
			stop("s.struct must be a matrix or list.")
		}
	}
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	gapOpening <- gapOpening/2 # split into gap opening and closing
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!is.numeric(gapPower))
		stop("gapPower must be a numeric.")
	if (!is.numeric(terminalGap))
		stop("terminalGap must be a numeric.")
	if (length(terminalGap) > 2 || length(terminalGap) < 1)
		stop("Length of terminalGap must be 1 or 2.")
	if (any(is.infinite(terminalGap)))
		stop("terminalGap must be finite.")
	if (length(terminalGap) == 1)
		terminalGap[2] <- terminalGap[1]
	if (!is.numeric(restrict))
		stop("restrict must be a numeric.")
	if (length(restrict) != 3)
		stop("Length of restrict must be 3.")
	if (restrict[1] >= 0)
		stop("restrict[1] must be less than zero.")
	if (restrict[2] < 0)
		stop("restrict[2] must be at least zero.")
	if (restrict[3] <= 0)
		stop("restrict[3] must be greater than zero.")
	if (floor(restrict[3]) != restrict[3])
		stop("restrict[3] must be a whole number.")
	restrict <- as.double(restrict)
	if (!is.numeric(anchor) && !is.na(anchor))
		stop("anchor must be numeric.")
	if (is.matrix(anchor)) {
		if (any(is.na(anchor)))
			stop("The matrix anchor must not contain NA values.")
		if (dim(anchor)[1] != 4)
			stop("The matrix anchor must have four rows.")
		if (any(anchor[2,] < anchor[1,] |
			anchor[4,] < anchor[3,]) ||
			any(anchor[2, -ncol(anchor)] >= anchor[1, -1] |
			anchor[4, -ncol(anchor)] >= anchor[3, -1]))
			stop("Anchors must be ascending order.")
		if (dim(anchor)[2] > 0) {
			if (anchor[1, 1] < 1)
				stop("The first anchor is outside the range of the pattern.")
			if (anchor[3, 1] < 1)
				stop("The first anchor is outside the range of the subject.")
			if (anchor[2, ncol(anchor)] > w.p)
				stop("The last anchor is outside the range of the pattern.")
			if (anchor[4, ncol(anchor)] > w.s)
				stop("The last anchor is outside the range of the subject.")
		}
	} else {
		if (is.numeric(anchor) && anchor <= 0)
			stop("anchor must be greater than zero.")
		if (is.numeric(anchor) && anchor > 1)
			stop("anchor must be less than or equal to one.")
	}
	if (!is.numeric(normPower))
		stop("normPower must be a numeric.")
	if (any(normPower < 0))
		stop("normPower must be at least zero.")
	if (length(normPower) < 2) {
		normPower <- rep(normPower, 2)
	} else if (length(normPower) > 2) {
		stop("Length of normPower must be 1 or 2.")
	}
	normPower <- as.double(normPower)
	if (!isTRUEorFALSE(standardize))
		stop("standardize must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors) != processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- .detectCores()
	} else {
		processors <- as.integer(processors)
	}
	if (length(pattern) > 2147483647)
		stop(paste("Length of pattern (",
			length(pattern),
			") longer than the maximum allowable length (2,147,483,647).",
			sep=""))
	if (length(subject) > 2147483647)
		stop(paste("Length of subject (",
			length(subject),
			") longer than the maximum allowable length (2,147,483,647).",
			sep=""))
	
	if (type == 3L) { # AAStringSet
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.null(substitutionMatrix)) {
			# use PFASUM50
			subMatrix <- matrix(c(4.1181,-1.1516,-1.3187,-1.4135,0.4271,-0.5467,-0.6527,0.1777,-1.6582,-1.1243,-1.1843,-1.0235,-0.5685,-1.9515,-0.6072,0.8284,0.0361,-2.5368,-2.1701,0.0661,-11,-1.1516,6.341,0.0543,-0.6628,-3.2085,1.6006,0.5067,-1.961,0.7706,-3.5053,-3.0357,2.938,-1.9894,-3.7846,-1.3455,-0.4194,-0.5594,-2.1629,-1.7957,-2.9403,-11,-1.3187,0.0543,6.4672,2.3024,-2.5179,0.8192,0.5566,0.1585,1.104,-4.1629,-4.0977,0.8743,-2.6216,-3.805,-1.0904,1.1291,0.3253,-3.7763,-1.874,-3.6076,-11,-1.4135,-0.6628,2.3024,6.8156,-4.358,0.6705,2.582,-0.5667,-0.196,-5.475,-5.1661,0.226,-3.9595,-5.3456,-0.5662,0.4273,-0.5218,-4.7691,-3.4644,-4.5477,-11,0.4271,-3.2085,-2.5179,-4.358,13.5349,-3.3641,-4.3086,-2.1614,-1.8945,-0.7546,-0.9453,-3.8239,-0.5923,-0.8182,-3.6019,-0.3927,-0.801,-1.9317,-1.1607,0.0673,-11,-0.5467,1.6006,0.8192,0.6705,-3.3641,5.5795,2.1372,-1.5923,1.0862,-3.3001,-2.7545,1.872,-1.1216,-3.6631,-1.0426,0.1982,-0.0434,-3.061,-1.9214,-2.6993,-11,-0.6527,0.5067,0.5566,2.582,-4.3086,2.1372,5.5684,-1.6462,-0.2488,-4.1849,-4.0275,1.4821,-2.7964,-4.8311,-0.7028,0.0283,-0.312,-4.1969,-2.9489,-3.281,-11,0.1777,-1.961,0.1585,-0.5667,-2.1614,-1.5923,-1.6462,7.6508,-1.8185,-4.7058,-4.4215,-1.5991,-3.2786,-3.9992,-1.4409,0.184,-1.4823,-3.8328,-3.7343,-3.7264,-11,-1.6582,0.7706,1.104,-0.196,-1.8945,1.0862,-0.2488,-1.8185,9.7543,-3.3812,-2.8685,0.1425,-1.8724,-1.2545,-1.5333,-0.4285,-0.8896,-0.9385,1.6476,-2.8729,-11,-1.1243,-3.5053,-4.1629,-5.475,-0.7546,-3.3001,-4.1849,-4.7058,-3.3812,5.1229,2.5319,-3.5454,1.8309,0.9346,-3.4603,-3.0985,-1.2543,-1.5006,-1.117,3.3961,-11,-1.1843,-3.0357,-4.0977,-5.1661,-0.9453,-2.7545,-4.0275,-4.4215,-2.8685,2.5319,4.7049,-3.4581,2.5303,1.687,-3.365,-3.1578,-1.8626,-0.5308,-0.6881,1.4829,-11,-1.0235,2.938,0.8743,0.226,-3.8239,1.872,1.4821,-1.5991,0.1425,-3.5454,-3.4581,5.5476,-2.164,-4.3516,-0.7583,0.0275,-0.1516,-3.5889,-2.4422,-3.0453,-11,-0.5685,-1.9894,-2.6216,-3.9595,-0.5923,-1.1216,-2.7964,-3.2786,-1.8724,1.8309,2.5303,-2.164,7.0856,1.2339,-3.0823,-1.7587,-0.7402,-0.5841,-0.3946,0.9477,-11,-1.9515,-3.7846,-3.805,-5.3456,-0.8182,-3.6631,-4.8311,-3.9992,-1.2545,0.9346,1.687,-4.3516,1.2339,7.4322,-3.6222,-3.0316,-2.2851,2.6305,3.8302,0.1942,-11,-0.6072,-1.3455,-1.0904,-0.5662,-3.6019,-1.0426,-0.7028,-1.4409,-1.5333,-3.4603,-3.365,-0.7583,-3.0823,-3.6222,9.1796,-0.0652,-0.8587,-3.3634,-3.3006,-2.5443,-11,0.8284,-0.4194,1.1291,0.4273,-0.3927,0.1982,0.0283,0.184,-0.4285,-3.0985,-3.1578,0.0275,-1.7587,-3.0316,-0.0652,4.2366,1.8491,-3.1454,-2.1838,-2.1839,-11,0.0361,-0.5594,0.3253,-0.5218,-0.801,-0.0434,-0.312,-1.4823,-0.8896,-1.2543,-1.8626,-0.1516,-0.7402,-2.2851,-0.8587,1.8491,4.8833,-2.8511,-1.8993,-0.2699,-11,-2.5368,-2.1629,-3.7763,-4.7691,-1.9317,-3.061,-4.1969,-3.8328,-0.9385,-1.5006,-0.5308,-3.5889,-0.5841,2.6305,-3.3634,-3.1454,-2.8511,13.6485,3.3017,-1.851,-11,-2.1701,-1.7957,-1.874,-3.4644,-1.1607,-1.9214,-2.9489,-3.7343,1.6476,-1.117,-0.6881,-2.4422,-0.3946,3.8302,-3.3006,-2.1838,-1.8993,3.3017,8.7568,-1.2438,-11,0.0661,-2.9403,-3.6076,-4.5477,0.0673,-2.6993,-3.281,-3.7264,-2.8729,3.3961,1.4829,-3.0453,0.9477,0.1942,-2.5443,-2.1839,-0.2699,-1.851,-1.2438,4.6928,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,14),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (is.character(substitutionMatrix)) {
			subMatrix <- .getSubMatrix(substitutionMatrix)
		} else if (is.matrix(substitutionMatrix)) {
			if (nrow(substitutionMatrix) != ncol(substitutionMatrix))
				stop("substitutionMatrix is not square.")
			if (sum(!(AAs %in% rownames(substitutionMatrix))) > 0L ||
				sum(!(AAs %in% colnames(substitutionMatrix))) > 0L)
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			stop("Invalid substitutionMatrix must be NULL, a character string, or a matrix.")
		}
		if (nrow(subMatrix) != length(AAs) || sum(rownames(subMatrix) != AAs) > 0L)
			subMatrix <- subMatrix[AAs, AAs]
		if (!is.double(subMatrix))
			mode(subMatrix) <- "numeric"
	} else {
		bases <- c("A", "C", "G",
			ifelse(type == 2L, "U", "T"))
		if (!is.null(substitutionMatrix)) {
			if (is.matrix(substitutionMatrix)) {
				if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
					any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
					stop("substitutionMatrix is incomplete.")
				substitutionMatrix <- substitutionMatrix[bases, bases]
				substitutionMatrix <- as.numeric(substitutionMatrix)
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		} else if (type == 2L && missing(perfectMatch) && missing(misMatch)) {
			sM <- matrix(c(11, 4, 5, 4, 4, 12, 4, 5, 5, 4, 12, 4, 4, 5, 4, 11),
				nrow=4,
				ncol=4,
				dimnames=list(bases, bases))
		}
	}
	
	if (type == 3L) {
		consensusProfile <- "consensusProfileAA"
	} else {
		consensusProfile <- "consensusProfile"
	}
	if (is.null(p.struct)) {
		p.profile <- .Call(consensusProfile,
			pattern,
			p.weight,
			NULL,
			PACKAGE="DECIPHER")
		s.profile <- .Call(consensusProfile,
			subject,
			s.weight,
			NULL,
			PACKAGE="DECIPHER")
	} else {
		if (is.null(structureMatrix)) {
			if (type == 3L) {
				# assume structures from PredictHEC
				structureMatrix <- matrix(c(3, 2, -1, 2, 12, -4, -1, -4, 1),
					nrow=3) # order is H, E, C
			} else {
				structureMatrix <- matrix(c(7, -3, -3, -3, 11, -8, -3, -8, 11),
					nrow=3) # order is ., (, )
			}
		} else {
			# assume structures and matrix are ordered the same
			if (!is.double(structureMatrix))
				stop("structureMatrix must be contain numerics.")
			if (!is.matrix(structureMatrix))
				stop("structureMatrix must be a matrix.")
			if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
				stop("structureMatrix is not square.")
		}
		
		if (is.list(p.struct)) {
			if (dim(structureMatrix)[1] != dim(p.struct[[1]])[1])
				stop("Dimensions of structureMatrix are incompatible with p.struct.")
			
			p.profile <- .Call(consensusProfile,
				pattern,
				p.weight,
				p.struct,
				PACKAGE="DECIPHER")
		} else { # p.struct is a matrix
			if (dim(structureMatrix)[1] != dim(p.struct)[1])
				stop("Dimensions of structureMatrix are incompatible with p.struct.")
			
			p.profile <- .Call(consensusProfile,
				pattern,
				p.weight,
				NULL,
				PACKAGE="DECIPHER")
			
			p.profile <- rbind(p.profile, p.struct)
		}
		if (is.list(s.struct)) {
			if (dim(structureMatrix)[1] != dim(s.struct[[1]])[1])
				stop("Dimensions of structureMatrix are incompatible with s.struct.")
			
			s.profile <- .Call(consensusProfile,
				subject,
				s.weight,
				s.struct,
				PACKAGE="DECIPHER")
		} else { # s.struct is a matrix
			if (dim(structureMatrix)[1] != dim(s.struct)[1])
				stop("Dimensions of structureMatrix are incompatible with s.struct.")
			
			s.profile <- .Call(consensusProfile,
				subject,
				s.weight,
				NULL,
				PACKAGE="DECIPHER")
			
			s.profile <- rbind(s.profile, s.struct)
		}
	}
	
	f <- function(p.profile, s.profile, tGaps=terminalGap) {
		p.d <- dim(p.profile)[2]
		s.d <- dim(s.profile)[2]
		size <- as.numeric(p.d)*as.numeric(s.d)
		if (size > 2147483647) # maximum when indexing by signed integer
			stop(paste("Alignment larger (",
				format(size, big.mark=","),
				") than the maximum allowable size (2,147,483,647).",
				sep=""))
		
		if (type == 3) { # AAStringSet
			if (is.null(p.struct)) {
				t <- .Call("alignProfilesAA",
					p.profile,
					s.profile,
					subMatrix,
					numeric(),
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					restrict,
					standardize,
					processors,
					PACKAGE="DECIPHER")
			} else {
				t <- .Call("alignProfilesAA",
					p.profile,
					s.profile,
					subMatrix,
					structureMatrix,
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					restrict,
					standardize,
					processors,
					PACKAGE="DECIPHER")
			}
		} else { # DNAStringSet or RNAStringSet
			if (is.null(p.struct)) {
				t <- .Call("alignProfiles",
					p.profile,
					s.profile,
					type,
					substitutionMatrix,
					numeric(),
					perfectMatch,
					misMatch,
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					restrict,
					standardize,
					processors,
					PACKAGE="DECIPHER")
			} else {
				t <- .Call("alignProfiles",
					p.profile,
					s.profile,
					type,
					substitutionMatrix,
					structureMatrix,
					perfectMatch,
					misMatch,
					gapOpening,
					gapExtension,
					gapPower,
					normPower,
					tGaps[1],
					tGaps[2],
					restrict,
					standardize,
					processors,
					PACKAGE="DECIPHER")
			}
		}
	}
	
	if (is.na(anchor[1])) { # don't use anchors
		inserts <- f(p.profile, s.profile)
	} else { # use anchors
		if (is.matrix(anchor)) {
			anchors <- anchor
		} else { # find anchors
			if (type == 3L) { # AAStringSet
				wordSize <- 7
			} else {
				wordSize <- 15
			}
			l <- min(length(pattern), length(subject))
			o.p <- order(p.weight, decreasing=TRUE)
			o.s <- order(s.weight, decreasing=TRUE)
			if (type == 3L) { # AAStringSet
				num.p <- .Call("enumerateGappedSequenceAA",
					pattern,
					wordSize,
					o.p[1:l],
					PACKAGE="DECIPHER")
				num.s <- .Call("enumerateGappedSequenceAA",
					subject,
					wordSize,
					o.s[1:l],
					PACKAGE="DECIPHER")
			} else {
				num.p <- .Call("enumerateGappedSequence",
					pattern,
					wordSize,
					o.p[1:l],
					PACKAGE="DECIPHER")
				num.s <- .Call("enumerateGappedSequence",
					subject,
					wordSize,
					o.s[1:l],
					PACKAGE="DECIPHER")
			}
			
			anchors <- .Call("matchRanges",
				num.p,
				num.s,
				wordSize,
				w.p,
				anchor,
				PACKAGE="DECIPHER")
			if (ncol(anchors) > 1L) {
				# reject unevenly spaced adjacent anchors
				d <- diff(anchors[1L,])/diff(anchors[3L,])
				d <- d >= 0.8 & d <= 1.25 # within rejection ratio
				d <- c(FALSE, d) | c(d, FALSE)
				anchors <- anchors[, d, drop=FALSE]
			}
		}
		
		numAnchors <- dim(anchors)[2]
		if (numAnchors == 0) {
			inserts <- f(p.profile, s.profile)
		} else {
			max.p <- which.max(p.weight)
			max.s <- which.max(s.weight)
			
			if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				1L, anchors[1, 1],
				1L, anchors[3, 1],
				max.p,
				max.s,
				PACKAGE="DECIPHER")) {
				inserts <- f(p.profile[, 1L:anchors[1, 1], drop=FALSE],
					s.profile[, 1L:anchors[3, 1], drop=FALSE],
					tGaps=c(terminalGap[1], -1e9))
			} else {
				inserts <- list(integer(), integer(),
					integer(), integer())
			}
			
			n <- 2L
			while (n <= numAnchors) { # align regions between anchors
				if (!.Call("firstSeqsEqual",
					pattern,
					subject,
					anchors[2, n - 1], anchors[1, n],
					anchors[4, n - 1], anchors[3, n],
					max.p,
					max.s,
					PACKAGE="DECIPHER")) {
					temp <- f(p.profile[, anchors[2, n - 1]:anchors[1, n], drop=FALSE],
						s.profile[, anchors[4, n - 1]:anchors[3, n], drop=FALSE],
						tGaps=c(-1e9, -1e9))
					inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, n - 1] - 1L)
					inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, n - 1] - 1L)
					inserts[[2]] <- c(inserts[[2]], temp[[2]])
					inserts[[4]] <- c(inserts[[4]], temp[[4]])
				}
				n <- n + 1L
			}
			
			n <- 1L
			while (n <= numAnchors) { # align anchor regions
				if (!.Call("firstSeqsGapsEqual",
					pattern,
					subject,
					anchors[1, n], anchors[2, n],
					anchors[3, n], anchors[4, n],
					type,
					max.p,
					max.s,
					PACKAGE="DECIPHER")) {
					temp <- .Call("firstSeqsPosEqual",
						pattern,
						subject,
						anchors[1, n], anchors[2, n],
						anchors[3, n], anchors[4, n],
						type,
						max.p,
						max.s,
						PACKAGE="DECIPHER")
					if (length(temp) == 4) {
						inserts[[1]] <- c(inserts[[1]], temp[[1]])
						inserts[[3]] <- c(inserts[[3]], temp[[3]])
					} else { # number of sites differs
						temp <- f(p.profile[, anchors[1, n]:anchors[2, n], drop=FALSE],
							s.profile[, anchors[3, n]:anchors[4, n], drop=FALSE],
							tGaps=c(-1e9, -1e9))
						inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[1, n] - 1L)
						inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[3, n] - 1L)
					}
					inserts[[2]] <- c(inserts[[2]], temp[[2]])
					inserts[[4]] <- c(inserts[[4]], temp[[4]])
				}
				n <- n + 1L
			}
			
			end.p <- anchors[2, numAnchors] == w.p
			end.s <- anchors[4, numAnchors] == w.s
			if (end.p && !end.s) { # need to add gaps after pattern
				inserts[[1]] <- c(inserts[[1]], w.p + 1L)
				inserts[[2]] <- c(inserts[[2]],
					w.s - anchors[4, numAnchors])
			} else if (end.s && !end.p) { # need to add gaps after subject
				inserts[[3]] <- c(inserts[[3]], w.s + 1L)
				inserts[[4]] <- c(inserts[[4]],
					w.p - anchors[2, numAnchors])
			} else if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				anchors[2, numAnchors], w.p,
				anchors[4, numAnchors], w.s,
				max.p,
				max.s,
				PACKAGE="DECIPHER")) { # need to align
				temp <- f(p.profile[, anchors[2, numAnchors]:w.p, drop=FALSE],
					s.profile[, anchors[4, numAnchors]:w.s, drop=FALSE],
					tGaps=c(-1e9, terminalGap[2]))
				inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, numAnchors] - 1L)
				inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, numAnchors] - 1L)
				inserts[[2]] <- c(inserts[[2]], temp[[2]])
				inserts[[4]] <- c(inserts[[4]], temp[[4]])
			} # else don't do anything
		}
	}
	
	ns.p <- names(pattern)
	ns.s <- names(subject)
	
	if (length(inserts[[1]]) > 0) {
		o <- order(inserts[[1]])
		pattern <- .Call("insertGaps",
			pattern,
			as.integer(inserts[[1]][o]),
			as.integer(inserts[[2]][o]),
			type,
			processors,
			PACKAGE="DECIPHER")
	}
	if (length(inserts[[3]]) > 0) {
		o <- order(inserts[[3]])
		subject <- .Call("insertGaps",
			subject,
			as.integer(inserts[[3]][o]),
			as.integer(inserts[[4]][o]),
			type,
			processors,
			PACKAGE="DECIPHER")
	}
	
	result <- .append(pattern, subject)
	
	if (!(is.null(ns.p) && is.null(ns.s))) {
		if (is.null(ns.p))
			ns.p <- rep(NA_character_, length(pattern))
		if (is.null(ns.s))
			ns.s <- rep(NA_character_, length(subject))
		names(result) <- c(ns.p, ns.s)
	}
	
	return(result)
}

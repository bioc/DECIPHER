.chainGenes <- function(allstarts,
	totScore,
	codScore,
	sameScores,
	oppoScores,
	maxOverlapSame,
	maxOverlapOpposite,
	minScore,
	includeLength,
	alternateScores,
	startScores,
	maxScore) {
	
	maxFracOverlap <- 0.25 # max overlap with adjacent genes
	maxFracOverlapInclude <- 0.01 # max overlap for forcing inclusion
	ORFanCutoff <- 300L # length cutoff for ORFans
	ORFanDensity <- 0.2 # allowed density of ORFans
	ORFanScore <- seq(minScore[2],
		minScore[1],
		length.out=ORFanDensity*100) # score cutoffs for ORFans
	ORFanDist <- 10L # minimum distance from nearest non-ORFan
	minScore <- minScore[1]
	
	alllengths <- allstarts[, 4L] - allstarts[, 3L] + 1L
	w <- which(totScore >= minScore)
	topstarts <- allstarts[w,, drop=FALSE]
	o <- order(topstarts[, 1L], topstarts[, 3L])
	topstarts <- topstarts[o,]
	toplengths <- alllengths[w[o]]
	topScore <- totScore[w[o]]
	# prevent choosing short ORFs over extending a longer ORF
	topScore <- topScore - log(toplengths) # analogous to multiple testing correction
	
	maxL <- (length(sameScores) - 1L)/2L
	if (maxL == maxOverlapSame &&
		maxL == maxOverlapOpposite) {
		scoreIntergenic <- TRUE
	} else {
		scoreIntergenic <- FALSE
	}
	
	indices <- .Call("chainGenes",
		topstarts,
		topScore,
		toplengths,
		scoreIntergenic,
		maxOverlapSame,
		maxOverlapOpposite,
		maxFracOverlap,
		sameScores,
		oppoScores,
		PACKAGE="DECIPHER")
	
	# eliminate ORFans (small genes in low density regions)
	remove <- logical(length(indices))
	currIndex <- 0L
	i <- 1L
	while (i < length(indices)) {
		if (currIndex != topstarts[indices[i], 1L]) {
			currIndex <- topstarts[indices[i], 1L]
			last <- 0L
		}
		if (toplengths[indices[i]] < ORFanCutoff &&
			topstarts[indices[i], 3L] - last >= ORFanDist) {
			j <- i + 1L
			while (j < length(indices) &&
				topstarts[indices[j], 1L] == topstarts[indices[i], 1L]) {
				if (toplengths[indices[j]] < ORFanCutoff) {
					j <- j + 1L
				} else {
					break
				}
			}
			if (topstarts[indices[j], 1L] == topstarts[indices[i], 1L]) {
				if (topstarts[indices[j], 3L] - topstarts[indices[j - 1L], 4L] < ORFanDist)
					j <- j - 1L
				if (j > i) {
					r <- i:(j - 1L)
					density <- sum(toplengths[indices[r]])/(topstarts[indices[j], 3L] - last)
					if (density < ORFanDensity)
						remove[r] <- topScore[indices[r]] < ORFanScore[ceiling(density*100)]
				}
			}
			i <- j
		}
		last <- topstarts[indices[i], 4L]
		i <- i + 1L
	}
	remove <- which(remove)
	if (length(remove) > 0)
		indices <- indices[-remove]
	
	indices <- w[o[indices]]
	
	missing <- which(alllengths >= includeLength &
		maxScore < minScore &
		alternateScores >= minScore &
		startScores > 0)
	currIndex <- 0
	while (length(missing) > 0) {
		if (currIndex != allstarts[missing[1L], 1L]) {
			currIndex <- allstarts[missing[1L], 1L]
			W <- which(allstarts[, 1L] == allstarts[missing[1L], 1L])
			begins <- allstarts[W, 3L]
			ends <- allstarts[W, 4L]
		}
		if (allstarts[missing[1], 2L] == 1L) { # reverse strand
			w <- W[begins == allstarts[missing[1L], 3L]]
		} else { # forward strand
			w <- W[ends == allstarts[missing[1L], 4L]]
		}
		if (all(totScore[w] < minScore)) { # could not be included already
			# pick max scoring (long) ORF to include
			i <- w[alllengths[w] >= includeLength]
			i <- i[which.max(totScore[i] - codScore[i])]
			
			# check for complete overlap
			overlap <- which(allstarts[indices, 3L] <= allstarts[i, 3L] &
				allstarts[indices, 4L] >= allstarts[i, 4L] &
				allstarts[indices, 1L] == allstarts[i, 1L])
			
			if (length(overlap) == 0) { # check beginning overlap
				overlap <- which(allstarts[indices, 3L] > allstarts[i, 3L] &
					allstarts[indices, 3L] < allstarts[i, 4L] &
					allstarts[indices, 1L] == allstarts[i, 1L])
				if (length(overlap) > 0) {
					delta <- allstarts[i, 4L] - allstarts[indices[overlap], 3L] + 1
					overlap <- which(delta > alllengths[i]*maxFracOverlapInclude |
						delta > alllengths[indices[overlap]]*maxFracOverlapInclude)
				}
			}
			
			if (length(overlap) == 0) { # check ending overlap
				overlap <- which(allstarts[indices, 4L] > allstarts[i, 3L] &
					allstarts[indices, 4L] < allstarts[i, 4L] &
					allstarts[indices, 1L] == allstarts[i, 1L])
				if (length(overlap) > 0) {
					delta <- allstarts[indices[overlap], 4L] - allstarts[i, 3L] + 1
					overlap <- which(delta > alllengths[i]*maxFracOverlapInclude |
						delta > alllengths[indices[overlap]]*maxFracOverlapInclude)
				}
			}
			
			if (length(overlap) == 0)
				indices <- c(indices, i)
		}
		missing <- missing[!(missing %in% w)]
	}
	
	return(sort(indices))
}

.extractSD <- function(dna) {
	
	up_region <- DNAString("WWGGDVDTGATCCWKCYGCASVTTCMSSTACRGMW")
	maxL <- 60 # length of upstream region to align
	SD <- DNAString("AAGGAGGT")
	
	dna <- xscat(dna, paste(rep("-", maxL), collapse=""))
	dna <- unlist(dna)
	dna <- xscat(dna,
		paste(rep("-", maxL), sep="", collapse=""),
		reverseComplement(dna))
	
	views1 <- matchPattern("ACCTTGTTACGACTT",
		dna,
		max.mismatch=2,
		with.indels=TRUE,
		fixed="subject")
	
	views2 <- matchPattern("GACGGGCGGTGTGTRCAA",
		dna,
		max.mismatch=3,
		with.indels=TRUE,
		fixed="subject")
	
	keep <- logical(length(views1))
	for (j in seq_along(views1)) {
		L <- start(views2) - end(views1)[j]
		keep[j] <- any(L >= 50 & L <= 150)
	}
	
	views <- views1[keep]
	
	starts <- start(views) - maxL
	starts <- starts[starts > 0]
	if (length(starts) > 0) {
		dna <- extractAt(dna,
			IRanges(starts,
				width=maxL))
		
		p <- pairwiseAlignment(RemoveGaps(dna),
			up_region,
			substitutionMatrix=nucleotideSubstitutionMatrix(),
			gapOpening=5,
			gapExtension=0,
			type="local-global")
		
		p <- DNAStringSet(pattern(p))
		p <- subseq(p, 1, length(SD))
		p <- RemoveGaps(p)
		p <- unique(p)
		
		if (length(p) == 1 &&
			width(p) == length(SD))
			SD <- p[[1]]
	}
	
	return(SD)
}

.extractRBS <- function(dna,
	deltaGrulesRNA,
	upstreamWidth,
	pattern,
	precision=2,
	spacing=3) {
	
	mapping <- diag(4) # don't penalize mismatches
	dimnames(mapping) <- list(DNA_BASES,
			DNA_BASES)
	mapping["G", "G"] <- 4
	mapping["C", "C"] <- 4
	
	w <- which(width(dna) == upstreamWidth)
	p <- pairwiseAlignment(dna[w],
		pattern,
		type="local-global",
		substitutionMatrix=mapping,
		gapOpening=10,
		gapExtension=4)
	pattern <- as.character(pattern(p))
	subject <- as.character(subject(p))
	
	dG <- .Call("calculateDeltaG",
		pattern,
		subject,
		deltaGrulesRNA,
		PACKAGE="DECIPHER")
	dG <- round(dG/precision)*precision
	
	pos <- (start(pattern(p)) + end(pattern(p)))/2
	pos <- round(pos/spacing)*spacing
	
	motif <- rep(NA_character_,
		length(dna))
	motif[w] <- paste(dG, pos)
	
	return(motif)
}

.performFolding <- function(foldLeft,
	foldRight,
	deltaGrulesRNA) {
	
	mapping <- diag(4)
	mapping[mapping == 0] <- -4
	dimnames(mapping) <- list(DNA_BASES,
			DNA_BASES)
	mapping["G", "G"] <- 6
	mapping["C", "C"] <- 6
	
	x <- xscat(foldLeft, foldRight)
	d <- duplicated(x)
	w <- which(!d)
	m <- match(x, x[w])
	
	p <- pairwiseAlignment(foldLeft[w],
		foldRight[w],
		type="overlap",
		substitutionMatrix=mapping,
		gapOpening=5,
		gapExtension=1)
	pattern <- as.character(pattern(p))
	subject <- as.character(subject(p))
	
	dG <- .Call("calculateDeltaG",
		pattern,
		subject,
		deltaGrulesRNA,
		PACKAGE="DECIPHER")
	
	return(dG[m])
}

.logisticRegression <- function(response,
	predictor,
	indices,
	remove,
	penalty=7.879, # qchisq(0.005, 1, lower.tail=FALSE)
	interactions=FALSE,
	bins=25,
	samples) {
	
	other <- seq_len(length(response))[-c(indices, remove)]
	w <- c(indices,
		sample(other,
			min(length(other),
				samples)))
	
	# remove any remaining NA values from indices
	sub_data <- cbind(response=response[w],
		predictor[w,, drop=FALSE])
	sub_data <- na.omit(sub_data)
	
	# fit a logistic regression model
	if (interactions) { # include pairwise interactions
		suppressWarnings(model <- glm(response ~ .^2,
			binomial,
			data=sub_data,
			control=glm.control(1e-6)))
	} else {
		suppressWarnings(model <- glm(response ~ .,
			binomial,
			data=sub_data,
			control=glm.control(1e-6)))
	}
	
#	# perform stepwise feature selection
#	suppressWarnings(model <- step(model,
#		k=penalty,
#		trace=0,
#		direction="both"))
	
	# convert probabilies into scores
	suppressWarnings(p <- predict(model,
		predictor,
		"response"))
	p <- .bincode(p,
		c(0,
			quantile(p[indices],
				seq(1/bins,
					1 - 1/bins,
					1/bins),
				na.rm=TRUE),
			1))
	fg <- tabulate(p[indices], bins)
	fg <- ifelse(fg == 0, 1, fg)
	fg <- fg/sum(fg)
	bg <- tabulate(p[-indices], bins)
	bg <- ifelse(bg == 0, 1, bg)
	bg <- bg/sum(bg)
	score <- log(fg/bg)
	score[p]
}

FindGenes <- function(myDNAStringSet,
	geneticCode=getGeneticCode("11"),
	minGeneLength=60,
	includeGenes=NULL,
	allowEdges=TRUE,
	allScores=FALSE,
	showPlot=FALSE,
	processors=1,
	verbose=TRUE) {
	
	# set default parameters
	threshold <- c(0.995, 0.999) # lower/upper cutoff quantiles for background scoring
	noiseCutoff <- 150L # cutoff length for background scoring
	backgroundLength <- 200L # background for lenScr
	initialCodons <- 8L # initial codon frequencies
	terminationCodons <- 6L # termination codon frequencies
	maxL <- 1000L # maximum gene overlap
	maxOverlapSame <- 60L # region for spline fitting overlap
	maxOverlapOpposite <- 200L # region for spline fitting overlap
	upstreamWidth <- 20L # upstream nucleotides for motif search
	upstreamMotifs <- 30L # upstream of upstreamWidth
	upstreamKmerSize <- 4L # must be <= 8
	RBSk <- 6L # length of possible RBS motifs
	foldingWidth <- c(15L, 25L) # nucleotides either side of start_offset and stop_offset
	start_offset <- seq(-10, 20, 10) # folding sites relative to start
	stop_offset <- seq(0, 30, 15) # folding sites relative to stop
	temp <- 37 # degrees Celsius
	ions <- 1 # Na+ equivalents
	length_params <- c(-4, 1400, 0.25, -10, 200) # initial parameters for sigmoid fitting the distribution of gene lengths
	includeMultiplier <- 2 # multiplier on x_intercept to seek inclusion despite negative score
	startCodonDist <- c("ATG"=85, "GTG"=10, "TTG"=4.9996, "CTG"=0.0001, "ATA"=0.0001, "ATT"=0.0001, "ATC"=0.0001) # prior distribution of initial codon frequencies
	signals <- matrix(c(3, 8, 1, 9, 12, 10, 4, 0, 2, 6, 7, 11),
		nrow=3)
	codonFreqCutoff <- c(0.10, 0.15, 0.25) # maximum average distance of codon frequencies within models
	improvement <- 20 # score increase required to choose an alternate codon model
	lenScrMultiplier <- 2 # multiplier on positive length scores
	staScrMultiplier <- 2 # multiplier on start scores
	intergenicCount <- 30 # number of intergenic observations to use a precise score
	maxIterations <- 19L # maximum number of iterations (must be > 1)
	includeFrac <- 1 # threshold of replicates required for inclusion in indices
	intercept <- 1.5 # log10(number) of predicted background genes at x_intercept
	samples <- c(1e3, 1e5) # minimum and maximum number of background samples for logistic regression
	maxD <- 4L # maximum coupling order (must be > 1 if couplingModel = TRUE)
	maxSD <- 0.05 # maximum standard deviation in coupling score when presumably uncoupled
	maxNormSD <- 1 # maximum normalized standard deviation in coupling score when presumably uncoupled
	codons <- c('AAA', 'AAC', 'AAG', 'AAT',
		'ACA', 'ACC', 'ACG', 'ACT',
		'AGA', 'AGC', 'AGG', 'AGT',
		'ATA', 'ATC', 'ATG', 'ATT',
		'CAA', 'CAC', 'CAG', 'CAT',
		'CCA', 'CCC', 'CCG', 'CCT',
		'CGA', 'CGC', 'CGG', 'CGT',
		'CTA', 'CTC', 'CTG', 'CTT',
		'GAA', 'GAC', 'GAG', 'GAT',
		'GCA', 'GCC', 'GCG', 'GCT',
		'GGA', 'GGC', 'GGG', 'GGT',
		'GTA', 'GTC', 'GTG', 'GTT',
		'TAA', 'TAC', 'TAG', 'TAT',
		'TCA', 'TCC', 'TCG', 'TCT',
		'TGA', 'TGC', 'TGG', 'TGT',
		'TTA', 'TTC', 'TTG', 'TTT')
	coefs <- matrix(c(1.085, -1.072, -0.701, 0.313, 1.058, -0.901, -0.749, 1.063, 1.18, -0.251, 0.184, 0.782, 0.553, -1.911, -0.353, 0.44, 0.898, -1.152, -0.523, -0.268, 1.933, 0, -0.065, 1.95, 0.166, -0.386, -0.564, 0.665, 0.201, -1.359, -1.796, 0.398, 1.587, -0.55, -0.119, 1.408, 2.355, 0.073, 0.051, 2.504, 2.111, 0.415, 1.19, 1.961, 1.015, -1.328, -0.857, 1.182, 0, -1.595, 0, -0.407, 1.09, -0.724, -1.462, 1.102, 0, -0.978, 0.091, -0.081, 1.142, -2.179, -0.179, 0.099, -1.574, 2.453, 1.888, -0.48, -2.906, 2.087, 1.137, -3.124, -3.688, 0.247, -1.39, -2.683, -1.961, 4.657, 1.569, -0.527, -2.012, 1.5, 1.381, 0.023, -4.899, -0.711, -0.11, -4.922, -2.039, 0.707, 0.198, -1.851, -1.732, 2.633, 4.284, -0.95, -1.714, 1.732, 1.038, -1.715, -4.357, 0.676, 0.26, -4.841, -4.426, -0.118, -2.677, -3.609, -2.691, 2.779, 2.279, -2.311, 0, 3.014, 0, 0.78, -3.426, 0.771, 2.107, -3.403, 0, 0.186, -0.529, -2.207, -2.956, 4.664, 0.223, -0.015),
		ncol=2,
		dimnames=list(codons,
			c("Intercept", "Slope")))
	
	# error checking
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	if (sum(width(myDNAStringSet)) == 0)
		stop("myDNAStringSet must contain sequences.")
	a <- vcountPattern("-", myDNAStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') must be removed from myDNAStringSet.")
	a <- vcountPattern(".", myDNAStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') must be removed from myDNAStringSet.")
	if (!is.character(geneticCode))
		stop("geneticCode must be a character vector.")
	if (length(geneticCode) != 64)
		stop("geneticCode must contain 64 codons.")
	if (is.null(names(geneticCode)))
		stop("geneticCode must be a named character vector.")
	startCodons <- names(geneticCode[which(geneticCode == "M")])
	startCodons <- c(startCodons,
		attr(geneticCode, "alt_init_codons"))
	if (length(startCodons) == 0)
		stop("geneticCode must contain at least one start codon.")
	if (length(startCodons) >= 10)
		stop("geneticCode must contain fewer than 10 start codons.")
	startIndex <- match(startCodons, codons) - 1L
	if (any(is.na(startIndex)))
		stop("geneticCode may only contain standard 3-letter start codons.")
	stopIndex <- names(geneticCode[which(geneticCode == "*")])
	stopIndex <- match(stopIndex, codons) - 1L
	if (length(stopIndex) == 0)
		stop("No stop codons (*) found in the geneticCode.")
	if (length(stopIndex) >= 10)
		stop("geneticCode must contain fewer than 10 stop codons.")
	if (any(startIndex %in% stopIndex))
		stop("Start codons cannot also be stop codons in the geneticCode.")
	gC <- geneticCode
	geneticCode <- geneticCode[codons]
	if (any(is.na(geneticCode)))
		stop("geneticCode is missing one or more codons.")
	geneticCode <- match(geneticCode, AA_STANDARD) - 1L
	if (!is.numeric(minGeneLength))
		stop("minGeneLength must be a numeric.")
	if (length(minGeneLength) != 1)
		stop("minGeneLength must contain a single integer.")
	if (minGeneLength < (2 + initialCodons)*3)
		stop("minGeneLength must be at least ",
			(2 + max(initialCodons, terminationCodons))*3,
			".")
	if (minGeneLength > 0.9*noiseCutoff)
		stop("minGeneLength can be at most ",
			floor(0.9*noiseCutoff),
			".")
	if (minGeneLength != floor(minGeneLength))
		stop("minGeneLength must be a whole number.")
	minGeneLength <- as.integer(minGeneLength)
	if (!is.null(includeGenes) &&
		!is(includeGenes, "Genes"))
		stop("includeGenes must be a Genes object.")
	if (!is.logical(allowEdges))
		stop("allowEdges must be a logical.")
	if (!is.logical(allScores))
		stop("allScores must be a logical.")
	if (!is.logical(showPlot))
		stop("showPlot must be a logical.")
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
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Iter  Models Start Motif  Init  Fold UpsNt  Term   RBS  Auto  Stop Genes\n1          1")
		flush.console()
		iter <- 1L
	}
	
	# find all possible genes
	ORFs <- .Call("getORFs",
		myDNAStringSet,
		startIndex,
		stopIndex,
		minGeneLength,
		allowEdges,
		PACKAGE="DECIPHER")
	colnames(ORFs) <- c("Index", "Strand", "Begin", "End")
	includeScores <- includeGenes[, "TotalScore"]
	annotations <- attr(includeGenes, "annotations")
	includeGene <- includeGenes[, "Gene"]
	includeGenes <- includeGenes[, colnames(ORFs)]
	mode(includeGenes) <- "integer" # coerce to integer matrix
	
	if (nrow(ORFs) == 0)
		stop("myDNAStringSet does not contain any open reading frames.")
	
	longest <- .Call("longestORFs",
		ORFs,
		PACKAGE="DECIPHER")
	l <- ORFs[, 4] - ORFs[, 3] + 1L
	
	# determine the background intergenic length distribution
	o <- longest[l[longest] >= noiseCutoff]
	o <- o[order(ORFs[o, 1], ORFs[o, 3])]
	intergenic <- ORFs[o[-1], 3] - ORFs[o[-length(o)], 4] + 1L
	sameStrand <- ORFs[o[-1], 2] == ORFs[o[-length(o)], 2]
	sameIndex <- ORFs[o[-1], 1] == ORFs[o[-length(o)], 1]
	r <- -maxL:maxL # symmetric length range of scores
	bg_same <- tabulate(intergenic[sameStrand & sameIndex] + maxL + 1L, length(r))
	first <- maxL - maxOverlapSame
	bg_same[seq_len(first)] <- 0L
	bg_same <- bg_same/sum(bg_same)
	bg_oppo <- tabulate(intergenic[!sameStrand & sameIndex] + maxL + 1L, length(r))
	bg_oppo[seq_len(first)] <- 0L
	bg_oppo <- bg_oppo/sum(bg_oppo)
	same_scores <- oppo_scores <- rep(0, length(r))
	
	# fit the background ORF length distribution
	t <- tabulate(l[longest])
	x <- which(t[seq_len(backgroundLength)] > 0)
	weight <- sqrt(t[x])
	logt <- log10(t)
	# fit logt = slope*x + y_intercept
	mx <- mean(x)
	mt <- mean(logt[x])
	slope <- sum(weight*(x - mx)*(logt[x] - mt))/sum(weight*(x - mx)^2)
	y_intercept <- mt - slope*mx
	x_intercept <- (intercept - y_intercept)/slope
	x_intercept <- as.integer(x_intercept)
	
	.PDF2 <- function(x, p)
		-((exp(p[1]*log(x/p[2]))+1)^(-p[3]-1)*(exp(p[4]*log(x/p[5]))+1)^(-p[3]-1)*(((p[4]*p[3]+p[1]*p[3])*exp(p[1]*log(x/p[2]))+p[4]*p[3])*exp(p[4]*log(x/p[5]))+p[1]*p[3]*exp(p[1]*log(x/p[2]))))/x
	.CDF2 <- function(x, p) {
		sig1 <- (1 + exp(p[1]*log(x/p[2])))^-p[3]
		sig2 <- (1 + exp(p[4]*log(x/p[5])))^-p[3]
		sig1*sig2
	}
	.fitSigmoid <- function(x, y, ini) {
		weights <- diff(c(0, log(x)))
		.SSE <- function(p) {
			# avoid conditions with numeric errors
			pdf <- .PDF2(x, p)
			if (any(is.nan(pdf) |
				is.infinite(pdf) |
				pdf < 0))
				return(sum(y^2))
			
			expected <- .CDF2(x, p)
			
			sum(weights*(y - expected)^2)
		}
		o <- suppressWarnings(optim(ini,
			.SSE,
			control=list(maxit=1e4)))
		
		return(o$par)
	}
	
	# fit the distribution of approximate gene lengths
	s <- seq(3, length(t), 3)
	bg <- 10^(slope*seq_along(t) + y_intercept)
	observed <- t[s] - bg[s] # expected number of genes
	observed[seq_len(backgroundLength/3)] <- 0
	observed[observed < 0] <- 0
	observed <- cumsum(observed)
	observed <- observed/observed[length(observed)]
	o <- .fitSigmoid(s,
		observed,
		length_params)
	
	# calculate log-odds model of gene length
	fg <- .PDF2(s, o)
	fg <- fg/sum(fg)
	bg <- bg[s]
	bg <- bg/sum(bg)
	fg[fg == 0] <- min(fg[fg > 0])
	lenScr <- log(fg/bg)
	w <- which(is.infinite(lenScr))
	if (length(w) > 0)
		lenScr[w] <- (lenScr[w[1L] - 1] - lenScr[w[1L] - 2])*seq_along(w) + lenScr[w[1L] - 1]
	lenScr[lenScr > 0] <- lenScrMultiplier*lenScr[lenScr > 0]
	lenScr <- lenScr[l/3]
	
	# compute preliminary codon scores from GC-content
	GC <- oligonucleotideFrequency(myDNAStringSet,
		width=1,
		step=1,
		as.prob=TRUE,
		simplify.as="collapse")
	GC <- GC["C"] + GC["G"]
	codon_scores <- coefs[, 1] + coefs[, 2]*GC
	preScr <- .Call("scoreCodonModel",
		myDNAStringSet,
		ORFs,
		codon_scores,
		PACKAGE="DECIPHER")
#	preScr <- preScr + lenScr
	
	# compute single codon log-odds model
	codon_scores <- .Call("codonModel",
		myDNAStringSet,
		ORFs,
		stopIndex,
		x_intercept,
		preScr, # >= 0 is considered
		PACKAGE="DECIPHER")
	codScr <- .Call("scoreCodonModel",
		myDNAStringSet,
		ORFs,
		codon_scores,
		PACKAGE="DECIPHER")
	orf_scores <- lenScr + codScr
	
	# add preliminary start codon scoring
	starts <- .Call("getRegion",
		myDNAStringSet,
		ORFs,
		3L,
		3L,
		TRUE,
		PACKAGE="DECIPHER")
	bgStarts <- table(starts)[startCodons]
	bgStarts[is.na(bgStarts)] <- 0L
	bgStarts <- bgStarts/sum(bgStarts)
	fgStarts <- startCodonDist[startCodons]
	fgStarts[is.na(fgStarts)] <- 0.01 # add missing start codons
	fgStarts <- fgStarts/sum(fgStarts)
	start_scores <- log(fgStarts/bgStarts)
	staScr <- unname(start_scores[starts])
	if (allowEdges) # potential non-canonical start
		staScr[is.na(staScr)] <- 0
	staScr <- staScrMultiplier*staScr
	orf_scores <- orf_scores + staScr
	
	# chain high-scoring ORFs into putative genes
	cutoff <- quantile(orf_scores[longest[which(l[longest] <= noiseCutoff)]],
		threshold)
	if (length(includeScores) > 0) {
		indices <- .chainGenes(rbind(ORFs,
				includeGenes),
			c(orf_scores,
				includeScores),
			c(codScr,
				rep(0, length(includeScores))),
			sameScores=same_scores,
			oppoScores=oppo_scores,
			maxOverlapSame=maxOverlapSame,
			maxOverlapOpposite=maxOverlapOpposite,
			minScore=cutoff,
			includeLength=includeMultiplier*x_intercept,
			alternateScores=c(preScr,
				rep(0, length(includeScores))),
			startScores=c(staScr,
				rep(0, length(includeScores))),
			maxScore=c(.Call("maxPerORF",
					ORFs,
					orf_scores,
					PACKAGE="DECIPHER"),
				includeScores))
	} else {
		indices <- .chainGenes(ORFs,
			orf_scores,
			codScr,
			sameScores=same_scores,
			oppoScores=oppo_scores,
			maxOverlapSame=maxOverlapSame,
			maxOverlapOpposite=maxOverlapOpposite,
			minScore=cutoff,
			includeLength=includeMultiplier*x_intercept,
			alternateScores=preScr,
			startScores=staScr,
			maxScore=.Call("maxPerORF",
				ORFs,
				orf_scores,
				PACKAGE="DECIPHER"))
	}
	
	if (verbose) {
		deltaSta <- mean(staScr[indices], na.rm=TRUE) - mean(staScr[-indices], na.rm=TRUE)
		show1 <- formatC(deltaSta,
			width=6,
			digits=2,
			format="f")
		cat(show1)
		cat("                 _                              ")
		show2 <- formatC(length(indices),
			width=6,
			digits=0,
			format="f")
		cat(show2)
		flush.console()
	}
	
	if (length(includeScores) > 0)
		indices <- indices[indices <= nrow(ORFs)]
	
	# split the upstream region into k-mers
	upstream <- .Call("getRegion",
		myDNAStringSet,
		ORFs,
		upstreamWidth,
		0L,
		TRUE,
		PACKAGE="DECIPHER")
	
	data("deltaHrulesRNA", envir=environment(), package="DECIPHER")
	data("deltaSrulesRNA", envir=environment(), package="DECIPHER")
	deltaSrulesRNA <- deltaSrulesRNA + 0.368*log(ions)/1000
	deltaGrulesRNA <- deltaHrulesRNA - (273.15 + temp)*deltaSrulesRNA
	
	# Shine-Delgarno free energy
	SD <- .extractSD(myDNAStringSet)
	upstream <- DNAStringSet(upstream)
	dG_RBS <- .extractRBS(upstream,
		deltaGrulesRNA,
		upstreamWidth,
		SD)
	
	# compute mRNA folding free energy
	dG_Fold <- matrix(NA_real_,
		nrow=nrow(ORFs),
		ncol=length(start_offset) + length(stop_offset))
	if (verbose) {
		count <- 2L
		symbols <- c("_", "\\", "|", "/")
		cat("\r1          1")
		cat(show1)
		cat("                 \\                              ")
		cat(show2)
		flush.console()
	}
	for(i in seq_len(ncol(dG_Fold))) {
		# tile subsequences 3 bases apart
		if (i <= length(start_offset)) {
			foldLeft <- .Call("getRegion",
				myDNAStringSet,
				ORFs,
				foldingWidth[1],
				start_offset[i],
				TRUE, # relative to start
				PACKAGE="DECIPHER")
			
			foldRight <- .Call("getRegion",
				myDNAStringSet,
				ORFs,
				-foldingWidth[1],
				start_offset[i],
				TRUE, # relative to start
				PACKAGE="DECIPHER")
			
			w <- which(nchar(foldLeft) == foldingWidth[1] &
				nchar(foldRight) == foldingWidth[1])
		} else {
			foldLeft <- .Call("getRegion",
				myDNAStringSet,
				ORFs,
				foldingWidth[2],
				stop_offset[i - length(start_offset)],
				FALSE, # relative to stop
				PACKAGE="DECIPHER")
			
			foldRight <- .Call("getRegion",
				myDNAStringSet,
				ORFs,
				-foldingWidth[2],
				stop_offset[i - length(start_offset)],
				FALSE, # relative to stop
				PACKAGE="DECIPHER")
			
			w <- which(nchar(foldLeft) == foldingWidth[2] &
				nchar(foldRight) == foldingWidth[2])
		}
		foldLeft <- DNAStringSet(unlist(foldLeft))
		foldRight <- DNAStringSet(unlist(foldRight))
		foldRight <- reverseComplement(foldRight)
		
		dG_Fold[w, i] <- .performFolding(foldLeft[w],
			foldRight[w],
			deltaGrulesRNA)
		
		if (verbose) {
			count <- count + 1L
			if (count > length(symbols))
				count <- 1L
			cat("\r1          1")
			cat(show1)
			cat("                 ")
			cat(symbols[count])
			cat("                              ")
			cat(show2)
			flush.console()
		}
	}
	reject <- which(rowSums(is.na(dG_Fold)) > 0)
	dG_Fold <- as.data.frame(dG_Fold)
	colnames(dG_Fold) <- c(paste("start",
			formatC(start_offset, flag="+"),
			sep=""),
		paste("stop",
			formatC(stop_offset, flag="+"),
			sep=""))
	
	.getGenes <- function(indices,
		multiModel=FALSE,
		startCodonModel=FALSE,
		initialCodons=0L,
		terminationCodons=0L,
		scoreRBS=FALSE,
		scoreAutocorr=FALSE,
		upstreamMotifs=0L,
		scoreFolding=FALSE,
		scoreIntergenicLength=FALSE,
		runLengthModel=FALSE,
		stopCodonModel=FALSE,
		scoreUpstream=FALSE,
		couplingModel=TRUE,
		allScores=FALSE,
		showPlot=FALSE) {
		
		genes <- ORFs[indices,, drop=FALSE]
		
		# compute gene length model
		observed <- tabulate(l[indices], s[length(s)])
		observed <- observed[s]
		observed <- cumsum(observed)
		observed <- observed/observed[length(observed)]
		# reset initial condition for quicker convergence
		length_params <<- .fitSigmoid(s,
			observed,
			length_params)
		fg <- .PDF2(s, length_params)
		fg <- fg/sum(fg)
		fg[fg == 0] <- min(fg[fg > 0])
		lenScr <- log(fg/bg)
		w <- which(is.infinite(lenScr))
		if (length(w) > 0)
			lenScr[w] <- (lenScr[w[1L] - 1] - lenScr[w[1L] - 2])*seq_along(w) + lenScr[w[1L] - 1]
		lenScr[lenScr > 0] <- lenScrMultiplier*lenScr[lenScr > 0]
		lenScr <- lenScr[l/3]
		
		if (showPlot) {
			layout(matrix(c(1, 3, 2, 4),
				nrow=2))
			plot(logt,
				xlab='ORF length (nucleotides)',
				ylab=expression('log'[10]*'(Frequency)'),
				xlim=c(minGeneLength, length(logt)),
				pch=46,
				log="x")
			abline(a=y_intercept,
				b=slope,
				h=intercept,
				v=x_intercept,
				col=c(6, 5, 5),
				lty=c(1, 2, 2),
				untf=TRUE)
			legend("topright",
				c("Background", "Intercepts", "ORF"),
				col=c(6, 5, 1),
				lty=c(1, 2, NA),
				pch=c(NA, NA, 15),
				bg="white")
			
			plot(s,
				100*diff(c(0, observed)),
				xlab="Length (nucleotides)",
				ylab="Frequency (%)",
				xlim=c(minGeneLength,
					max(1e4,
						s[length(s)])),
				pch=46,
				col="darkgreen",
				log="x")
			lines(s,
				100*fg,
				lwd=2,
				col="darkgreen")
			lines(s,
				100*bg,
				col="magenta")
			legend("topright",
				c("Background", "Estimated", "Gene"),
				col=c("magenta", "darkgreen", "darkgreen"),
				lty=c(1, 1, NA),
				lwd=c(1, 2, NA),
				pch=c(NA, NA, 15),
				bg="white")
		}
		
		# compute codon log-odds models
		codFreqs <- .Call("codonFrequencies",
			myDNAStringSet,
			ORFs,
			indices,
			PACKAGE="DECIPHER")
		freqs <- colMeans(codFreqs)
		minSize <- quantile(freqs[freqs > 0], 0.1)
		minSize <- c(90*minSize^-1, # rare codons ~30 times
			0.3*minSize^-2) # rare dicodons ~0.1 times
		if (multiModel) {
			if (nrow(codFreqs) > 10000) # use a subset
				codFreqs <- codFreqs[seq_len(10000),]
			dists <- .Call("dist",
				codFreqs,
				processors,
				PACKAGE="DECIPHER")
			attr(dists, "Size") <- nrow(codFreqs)
			class(dists) <- "dist"
			clust <- TreeLine(myDistMatrix=dists,
				type="clusters",
				method="UPGMA",
				cutoff=codonFreqCutoff,
				processors=processors,
				verbose=FALSE)
			clusts <- list()
			for (i in seq_along(codonFreqCutoff))
				clusts <- c(clusts,
					lapply(seq_len(max(clust[, i])),
						function(x)
							which(clust[, i] == x)))
			size <- sapply(clusts,
				function(x)
					sum(l[indices[x]]))
			keep <- size >= minSize[1]
			clusts <- clusts[keep]
			size <- size[keep]
			if (length(clusts) > 1 ||
				length(clusts) == 0 ||
				(length(clusts) == 1 &&
				length(clusts[[1]]) < length(indices))) {
				clusts <- c(list(seq_along(indices)),
					clusts)
				size <- c(sum(l[indices]),
					size)
			}
		} else {
			clusts <- list(seq_along(indices))
			size <- sum(l[indices])
		}
		for (i in seq_along(clusts)) {
			if (size[i] >= minSize[2]) { # use dicodon model
				codon_scores <- .Call("dicodonModel",
					myDNAStringSet,
					genes[clusts[[i]],, drop=FALSE],
					stopIndex,
					PACKAGE="DECIPHER")
				#dimnames(codon_scores) <- list(codons, codons)
				# codon_scores matrix is [next, prev]
			} else { # use unicodon model
				codon_scores <- .Call("unicodonModel",
					myDNAStringSet,
					genes[clusts[[i]],, drop=FALSE],
					stopIndex,
					PACKAGE="DECIPHER")
			}
			tempScr <- .Call("scoreCodonModel",
				myDNAStringSet,
				ORFs,
				codon_scores,
				PACKAGE="DECIPHER")
			
			if (i == 1L) {
				codScr <- tempScr
				
				# add run length into model
				if (runLengthModel) {
					run_scores <- .Call("runLengthModel",
						myDNAStringSet,
						genes,
						codon_scores,
						PACKAGE="DECIPHER")
					runScr <- .Call("scoreRunLengthModel",
						myDNAStringSet,
						ORFs,
						codon_scores,
						run_scores,
						PACKAGE="DECIPHER")
					orf_scores <- orf_scores + runScr
				}
			} else if (i == 2L) {
				altScr <- tempScr
				if (allScores)
					models <- rep(i, length(altScr))
			} else { # record the best alternative codScr
				w <- which(.Call("maxPerORF", ORFs, tempScr, PACKAGE="DECIPHER") > .Call("maxPerORF", ORFs, altScr, PACKAGE="DECIPHER"))
				if (length(w) > 0) {
					altScr[w] <- tempScr[w]
					if (allScores)
						models[w] <- i
				}
			}
		}
		if (length(clusts) > 1) {
#			w <- which(altScr > 0 &
#				altScr - improvement >= codScr)
			w <- which(.Call("maxPerORF", ORFs, altScr, PACKAGE="DECIPHER") - improvement > .Call("maxPerORF", ORFs, codScr, PACKAGE="DECIPHER"))
			if (length(w) > 0) {
				codScr[w] <- altScr[w]
				if (allScores)
					models[-w] <- 1L
			} else {
				models <- rep(1L, length(codScr))
			}
		} else if (allScores) {
			models <- rep(1L, length(codScr))
		}
		orf_scores <- lenScr + codScr
		
		# add amino acid coupling into model
		if (couplingModel) {
			coupling_scores <- .Call("couplingModel",
				myDNAStringSet,
				ORFs,
				indices,
				geneticCode,
				maxD + 20L,
				PACKAGE="DECIPHER")
			# nullify score for rare couplings (assume uncoupled past maxD)
			sdevs <- apply(coupling_scores[, -seq_len(maxD), drop=FALSE], 1, sd)
			norm_sdevs <- abs(sdevs/rowMeans(coupling_scores[, -seq_len(maxD), drop=FALSE]))
			coupling_scores[norm_sdevs > maxNormSD & sdevs > maxSD,] <- 0
			coupling_scores <- coupling_scores[, seq_len(maxD), drop=FALSE]
			di <- size[1] >= minSize[2] # dicodon model used
			if (di) {
				coupling_scores[, 1] <- 0 # coupling already included in dicodon scoring
				coupling_scores <- coupling_scores/(maxD - 1)
			} else {
				coupling_scores <- coupling_scores/maxD
			}
			couScr <- .Call("scoreCouplingModel",
				myDNAStringSet,
				ORFs,
				coupling_scores,
				geneticCode,
				PACKAGE="DECIPHER")
			orf_scores <- orf_scores + couScr
		}
		
		# add start codon preferences into model
		if (startCodonModel) {
			start_scores <- .Call("startCodonModel",
				myDNAStringSet,
				ORFs,
				indices,
				startIndex,
				PACKAGE="DECIPHER")
			staScr <- .Call("scoreStartCodonModel",
				myDNAStringSet,
				ORFs,
				start_scores,
				PACKAGE="DECIPHER")
			staScr <- staScrMultiplier*staScr
		}
		deltaSta <- mean(staScr[indices]) - mean(staScr[-indices])
		if (deltaSta > 0) {
			orf_scores <- orf_scores + staScr
		} else if (verbose) {
			deltaSta <- "-"
		}
		
		# add initial codon bias into model
		if (initialCodons) {
			ini_codon_scores <- .Call("initialCodonModel",
				myDNAStringSet,
				ORFs,
				indices,
				initialCodons,
				PACKAGE="DECIPHER")
			iniScr <- .Call("scoreInitialCodonModel",
				myDNAStringSet,
				ORFs,
				ini_codon_scores,
				PACKAGE="DECIPHER")
			deltaIni <- mean(iniScr[indices]) - mean(iniScr[-indices])
			if (deltaIni < 0) {
				deltaIni <- "-"
			} else {
				orf_scores <- orf_scores + iniScr
			}
		}
		
		# add termination codon bias into model
		if (terminationCodons) {
			ter_codon_scores <- .Call("terminationCodonModel",
				myDNAStringSet,
				ORFs,
				indices,
				terminationCodons,
				PACKAGE="DECIPHER")
			terScr <- .Call("scoreTerminationCodonModel",
				myDNAStringSet,
				ORFs,
				ter_codon_scores,
				PACKAGE="DECIPHER")
			deltaTer <- mean(terScr[indices]) - mean(terScr[-indices])
			if (deltaTer < 0) {
				deltaTer <- "-"
			} else {
				orf_scores <- orf_scores + terScr
			}
		}
		
		# add stop codon preferences into model
		if (stopCodonModel) {
			stop_scores <- .Call("stopCodonModel",
				myDNAStringSet,
				ORFs,
				indices,
				stopIndex,
				PACKAGE="DECIPHER")
			stoScr <- .Call("scoreStopCodonModel",
				myDNAStringSet,
				ORFs,
				stop_scores,
				PACKAGE="DECIPHER")
			deltaSto <- mean(stoScr[indices]) - mean(stoScr[-indices])
			if (deltaSto < 0) {
				deltaSto <- "-"
			} else {
				orf_scores <- orf_scores + stoScr
			}
		}
		
		# add autocorrelation into model
		if (scoreAutocorr) {
			autocorr_scores <- .Call("autocorrelationModel",
				myDNAStringSet,
				ORFs,
				indices,
				geneticCode,
				PACKAGE="DECIPHER")
			autScr <- .Call("scoreAutocorrelationModel",
				myDNAStringSet,
				ORFs,
				autocorr_scores,
				geneticCode,
				PACKAGE="DECIPHER")
			deltaAut <- mean(autScr[indices]) - mean(autScr[-indices])
			if (deltaAut < 0) {
				deltaAut <- "-"
			} else {
				orf_scores <- orf_scores + autScr
			}
		}
		
		# add upstream nucleotide bias model
		if (scoreUpstream) {
			up_nuc_scores <- .Call("nucleotideBiasModel",
				myDNAStringSet,
				ORFs,
				indices,
				upstreamWidth,
				PACKAGE="DECIPHER")
			upsScr <- .Call("scoreNucleotideBiasModel",
				myDNAStringSet,
				ORFs,
				up_nuc_scores,
				PACKAGE="DECIPHER")
			deltaUps <- mean(upsScr[indices]) - mean(upsScr[-indices])
			if (deltaUps < 0) {
				deltaUps <- "-"
			} else {
				orf_scores <- orf_scores + upsScr
			}
		}
		
		# add upstream motifs into model
		if (upstreamMotifs) {
			upstream_motif_scores <- .Call("upstreamMotifModel",
				myDNAStringSet,
				ORFs,
				indices,
				upstreamWidth + 1L,
				upstreamMotifs,
				upstreamKmerSize,
				PACKAGE="DECIPHER")
			motScr <- .Call("scoreUpstreamMotifModel",
				myDNAStringSet,
				ORFs,
				upstream_motif_scores,
				upstreamWidth + 1L,
				upstreamMotifs,
				upstreamKmerSize,
				PACKAGE="DECIPHER")
			deltaMot <- mean(motScr[indices]) - mean(motScr[-indices])
			if (deltaMot < 0) {
				deltaMot <- "-"
			} else {
				orf_scores <- orf_scores + motScr
			}
		}
		
		if (scoreFolding ||
			scoreRBS) {
			response <- integer(nrow(ORFs))
			response[indices] <- 1L
		}
		
		# add folding around start into model
		if (scoreFolding) {
			folScr <- .logisticRegression(response,
				dG_Fold,
				indices,
				reject,
#				interactions=TRUE,
				samples=ifelse(allScores,
					samples[2], # last iteration
					samples[1]))
			
			if (length(reject) > 0)
				folScr[reject] <- 0
			
			deltaFol <- mean(folScr[indices]) - mean(folScr[-indices])
			if (deltaFol < 0) {
				deltaFol <- "-"
			} else {
				orf_scores <- orf_scores + folScr
			}
		}
		
		# add ribosome binding sites into model
		if (scoreRBS) {
			t <- unique(dG_RBS)
			t <- setNames(rep(1,
					length(t)),
				t)
			
			obs <- t
			temp <- table(dG_RBS[indices])
			obs[names(temp)] <- temp
			obs[is.na(names(obs))] <- 0
			obs <- obs/sum(obs)
			
			bg <- t
			temp <- table(dG_RBS[-indices])
			bg[names(temp)] <- temp
			bg[is.na(names(bg))] <- 0
			bg <- bg/sum(bg)
			
			rbsScr <- log(obs/bg)
			rbsScr <- rbsScr[dG_RBS]
			rbsScr[is.na(dG_RBS)] <- 0
			
			# learn alternative RBS sites
			fg <- oligonucleotideFrequency(upstream[indices],
				RBSk,
				simplify.as="collapse")
			fg <- fg/sum(fg)
			bg <- oligonucleotideFrequency(upstream[-indices],
				RBSk,
				simplify.as="collapse")
			bg[bg == 0] <- 1L # add pseudocounts
			bg <- bg/sum(bg)
			score <- log(fg/bg)
			score <- sort(score,
				decreasing=TRUE)
			score <- score[score >= 1]
			if (length(score) > 0) {
				if (length(score) > 1000)
					score <- score[1:1000]
				
				motif <- names(score)
				motif <- DNAStringSet(motif)
				
				dists <- DistanceMatrix(motif,
					processors=processors,
					verbose=FALSE)
				clusts <- TreeLine(myDistMatrix=dists,
					type="clusters",
					method="complete",
					cutoff=1/RBSk, # 1 mismatch
					processors=processors,
					verbose=FALSE)
				clusts <- clusts$cluster
				o <- tapply(score,
					clusts,
					sum)
				o <- order(o,
					decreasing=TRUE)
				if (length(o) > 20)
					o <- o[1:20]
				
				a <- oligonucleotideFrequency(upstream,
					1, # base frequencies
					simplify.as="collapse",
					as.prob=TRUE)
				PWMs <- lapply(o,
					function(x) {
						w <- which(clusts == x)
						
						# weight motifs approximately by score
						m <- rep(motif[w],
							score[w])
						
						pwm <- PWM(m,
							prior.params=a)
						
						ns <- ConsensusSequence(m)
						ns <- as.character(ns)
						
						list(pwm, ns)
					})
				
				hits <- vector("list",
					length(PWMs) + 1L)
				names(hits) <- c(sapply(PWMs,
						`[`,
						2L),
					"dG_SD")
				PWMs <- lapply(PWMs,
					`[[`,
					1L)
				
				begin <- cumsum(width(upstream)) - width(upstream)
				w <- which(width(upstream) != upstreamWidth)
				if (length(w) > 0)
					begin <- begin[-w]
				merged <- unlist(upstream)
				pos <- seq(0L,
					upstreamWidth - RBSk,
					1L)
				for (i in seq_along(PWMs)) {
					temp <- .Call("scoreTopPWM",
						PWMs[[i]],
						merged,
						begin,
						pos,
						processors,
						PACKAGE="DECIPHER")
					if (length(w) > 0) {
						scores <- numeric(length(upstream))
						scores[-w] <- temp
					} else {
						scores <- temp
					}
					hits[[i]] <- scores
				}
				hits[[length(hits)]] <- rbsScr
				hits <- data.frame(hits)
				
				rbsScr <- .logisticRegression(response,
					hits,
					indices,
					w,
					samples=ifelse(allScores,
						samples[2], # last iteration
						samples[1]))
				if (length(w) > 0)
					rbsScr[w] <- 0
			}
			
			deltaRbs <- mean(rbsScr[indices]) - mean(rbsScr[-indices])
			if (deltaRbs < 0) {
				deltaRbs <- "-"
			} else {
				orf_scores <- orf_scores + rbsScr
			}
		}
		
		# add intergenic length into model
		if (scoreIntergenicLength) {
			o <- order(ORFs[indices, 1], ORFs[indices, 3])
			intergenic <- ORFs[indices[o][-1], 3] - ORFs[indices[o][-length(o)], 4] + 1L
			sameStrand <- ORFs[indices[o][-1], 2] == ORFs[indices[o][-length(o)], 2]
			sameIndex <- ORFs[indices[o][-1], 1] == ORFs[indices[o][-length(o)], 1]
			obs_same <- tabulate(intergenic[sameStrand & sameIndex] + maxL + 1L, length(r))
			obs_same[seq_len(first)] <- 0L
			# accumulate background
			sum <- 0
			for (i in (maxL - 1):(first + 1)) {
				if (obs_same[i] == 0) {
					sum <- sum + bg_same[i]
				} else if (sum > 0) {
					bg_same[i] <- bg_same[i] + sum
					sum <- 0
				}
			}
			obs_oppo <- tabulate(intergenic[!sameStrand & sameIndex] + maxL + 1L, length(r))
			obs_oppo[seq_len(first)] <- 0L
			# accumulate background
			sum <- 0
			for (i in (maxL - 1):(first + 1)) {
				if (obs_oppo[i] == 0) {
					sum <- sum + bg_oppo[i]
				} else if (sum > 0) {
					bg_oppo[i] <- bg_oppo[i] + sum
					sum <- 0
				}
			}
			
			obs_sum <- sum(obs_same, obs_oppo)
			odds_same <- log(obs_same/obs_sum/bg_same)
			odds_same[is.nan(odds_same) | is.infinite(odds_same)] <- 0
			t <- try(same_scores <- predict(smooth.spline(r,
					odds_same,
					obs_same,
					df=20))[[2]],
				silent=TRUE)
			if (is(t, "try-error")) {
				same_scores <- numeric(length(odds_same))
			} else {
				sub <- which(obs_same >= intergenicCount)
				if (length(sub) > 0)
					same_scores[sub] <- odds_same[sub]
			}
			
			odds_oppo <- log(obs_oppo/obs_sum/bg_oppo)
			odds_oppo[is.nan(odds_oppo) | is.infinite(odds_oppo)] <- 0
			t <- try(oppo_scores <- predict(smooth.spline(r,
					odds_oppo,
					obs_oppo,
					df=20))[[2]],
				silent=TRUE)
			if (is(t, "try-error")) {
				oppo_scores <- numeric(length(odds_oppo))
			} else {
				sub <- which(obs_oppo >= intergenicCount)
				if (length(sub) > 0)
					oppo_scores[sub] <- odds_oppo[sub]
			}
			
			maxOverlapSame <- maxL
			maxOverlapOpposite <- maxL
			if (showPlot) {
				halfMaxObs <- max(obs_same, obs_oppo)/2
				plot(r,
					odds_same,
					xlab="Intergenic length (nucleotides)",
					ylab="Log-odds",
					col="darkblue",
					cex=obs_same/halfMaxObs,
					ylim=c(-5, 5),
					xlim=c(-200, 200))
				points(r,
					odds_oppo,
					cex=obs_oppo/halfMaxObs,
					pch=0,
					col="darkred")
				lines(r,
					same_scores,
					col="darkblue")
				lines(r,
					oppo_scores,
					lty=3,
					lwd=2,
					col="darkred")
				legend("topleft",
					c("Same", "Oppo."),
					col=c("darkblue", "darkred"),
					lty=c(1, 3),
					pch=c(1, 0),
					lwd=c(1, 2),
					title=expression(underline("Strand")),
					bg="white")
			}
		}
		
		cutoff <- quantile(orf_scores[l <= noiseCutoff],
			threshold)
		cutoff <- cutoff + same_scores[maxL + 1L] # add typical intergenic score
		
		# chain high-scoring ORFs into genes
		if (length(includeScores) > 0) {
			indices <- .chainGenes(rbind(ORFs,
					includeGenes),
				c(orf_scores,
					includeScores),
				c(codScr,
					rep(0, length(includeScores))),
				sameScores=same_scores,
				oppoScores=oppo_scores,
				maxOverlapSame=maxOverlapSame,
				maxOverlapOpposite=maxOverlapOpposite,
				minScore=cutoff,
				includeLength=includeMultiplier*x_intercept,
				alternateScores=c(preScr,
					rep(0, length(includeScores))),
				startScores=c(staScr,
					rep(0, length(includeScores))),
				maxScore=c(.Call("maxPerORF",
						ORFs,
						orf_scores,
						PACKAGE="DECIPHER"),
					includeScores))
		} else {
			indices <- .chainGenes(ORFs,
				orf_scores,
				codScr,
				sameScores=same_scores,
				oppoScores=oppo_scores,
				maxOverlapSame=maxOverlapSame,
				maxOverlapOpposite=maxOverlapOpposite,
				minScore=cutoff,
				includeLength=includeMultiplier*x_intercept,
				alternateScores=preScr,
				startScores=staScr,
				maxScore=.Call("maxPerORF",
					ORFs,
					orf_scores,
					PACKAGE="DECIPHER"))
		}
		
		if (showPlot) {
			plot(l[longest],
				orf_scores[longest],
				xlab="ORF length (nucleotides)",
				ylab="Total score",
				log="x",
				col="#00000033",
				pch=46)
			points(l[indices],
				orf_scores[indices],
				col="darkgreen",
				pch=46)
			abline(h=cutoff[1],	
				lty=2,
				col="orange")
			legend("topleft",
				c("ORF", "Gene", "Threshold"),
				pch=c(15, 15, NA),
				lty=c(NA, NA, 2),
				col=c("black", "darkgreen", "orange"),
				bg="white")
		}
		
		if (verbose) {
			cat("\r")
			iter <<- iter + 1L
			show <- formatC(iter,
				width=-6,
				digits=0,
				format="f")
			cat(show)
			show <- formatC(length(size),
				width=6,
				digits=0,
				format="f")
			cat(show)
			show <- formatC(deltaSta,
				width=6,
				digits=2,
				format="f")
			cat(show)
			if (upstreamMotifs) {
				show <- formatC(deltaMot,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			if (initialCodons) {
				show <- formatC(deltaIni,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			if (scoreFolding) {
				show <- formatC(deltaFol,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			if (scoreUpstream) {
				show <- formatC(deltaUps,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			if (terminationCodons) {
				show <- formatC(deltaTer,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			if (scoreRBS) {
				show <- formatC(deltaRbs,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			if (scoreAutocorr) {
				show <- formatC(deltaAut,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			if (stopCodonModel) {
				show <- formatC(deltaSto,
					width=6,
					digits=2,
					format="f")
				cat(show)
			} else {
				cat("      ")
			}
			show <- formatC(length(indices),
				width=6,
				digits=0,
				format="f")
			cat(show)
			flush.console()
		}
		
		if (allScores) {
			ans <- cbind(ORFs,
				TotalScore=orf_scores,
				LengthScore=lenScr,
				CodingScore=codScr,
				CodonModel=models)
			if (couplingModel)
				ans <- cbind(ans,
					CouplingScore=couScr)
			if (startCodonModel)
				ans <- cbind(ans,
					StartScore=staScr)
			if (stopCodonModel)
				ans <- cbind(ans,
					StopScore=stoScr)
			if (initialCodons > 0)
				ans <- cbind(ans,
					InitialCodonScore=iniScr)
			if (terminationCodons > 0)
				ans <- cbind(ans,
					TerminationCodonScore=terScr)
			if (scoreRBS)
				ans <- cbind(ans,
					RibosomeBindingSiteScore=rbsScr)
			if (scoreAutocorr)
				ans <- cbind(ans,
					AutocorrelationScore=autScr)
			if (scoreUpstream)
				ans <-cbind(ans,
					UpstreamNucleotideScore=upsScr)
			if (upstreamMotifs)
				ans <- cbind(ans,
					UpstreamMotifScore=motScr)
			if (scoreFolding)
				ans <- cbind(ans,
					FoldingScore=folScr)
			if (runLengthModel)
				ans <- cbind(ans,
					RunLengthScore=runScr)
			return(list(ans, indices))
		} else {
			return(indices)
		}
	}
	
	params <- list()
	for (j in seq_len(ncol(signals))) {
		for (i in seq_len(nrow(signals))) {
			if (signals[i, j] == 1) {
				params <- c(params,
					list(multiModel=TRUE))
			} else if (signals[i, j] == 2) {
				params <- c(params,
					list(startCodonModel=TRUE))
			} else if (signals[i, j] == 3) {
				params <- c(params,
					list(initialCodons=initialCodons))
			} else if (signals[i, j] == 4) {
				params <- c(params,
					list(terminationCodons=terminationCodons))
			} else if (signals[i, j] == 5) {
				params <- c(params,
					list(runLengthModel=TRUE))
			} else if (signals[i, j] == 6) {
				params <- c(params,
					list(scoreRBS=TRUE))
			} else if (signals[i, j] == 7) {
				params <- c(params,
					list(scoreAutocorr=TRUE))
			} else if (signals[i, j] == 8) {
				params <- c(params,
					list(upstreamMotifs=upstreamMotifs))
			} else if (signals[i, j] == 9) {
				params <- c(params,
					list(scoreFolding=TRUE))
			} else if (signals[i, j] == 10) {
				params <- c(params,
					list(scoreIntergenicLength=TRUE))
			} else if (signals[i, j] == 11) {
				params <- c(params,
					list(stopCodonModel=TRUE))
			} else if (signals[i, j] == 12) {
				params <- c(params,
					list(scoreUpstream=TRUE))
			}
		}
		indices <- do.call(.getGenes,
			c(params, list(indices)))
		if (length(includeScores) > 0)
			indices <- indices[indices <= nrow(ORFs)]
	}
	
	bootstraps <- integer(nrow(ORFs) + length(includeScores))
	bootstraps[indices] <- 1L
	prev_indices <- indices
	for (i in 2L:maxIterations) {
		if (length(includeScores) > 0) {
			indices <- do.call(.getGenes,
				c(params, list(indices[indices <= nrow(ORFs)])))
		} else {
			indices <- do.call(.getGenes,
				c(params, list(indices)))
		}
		
		bootstraps[indices] <- bootstraps[indices] + 1L
		indices <- which(bootstraps >= i*includeFrac)
		if (length(prev_indices) == length(indices) &&
			all(prev_indices == indices))
			break
		prev_indices <- indices
	}
	
	if (length(includeScores) > 0)
		indices <- indices[indices <= nrow(ORFs)]
	params <- c(params,
		list(indices=indices, # indices no longer changes
			allScores=TRUE,
			showPlot=showPlot))
	ans <- do.call(.getGenes,
		params)
	indices <- ans[[2]]
	bootstraps[indices] <- bootstraps[indices] + 1L
	ans <- ans[[1]]
	bootstraps[ans[, "TotalScore"] <= 0] <- 0L
	if (length(includeScores) > 0) {
		if (allScores) {
			gene <- numeric(nrow(ans))
			gene[indices[indices <= nrow(ORFs)]] <- 1
			
			indices <- indices[indices > nrow(ORFs)] - nrow(ORFs)
			if (length(indices) > 0) {
				includeGene[-indices] <- 0
				includeGene[-indices] <- -includeGene[-indices]
			} else {
				includeGene[] <- 0
				includeGene[] <- -includeGene[]
			}
			gene <- c(gene, includeGene)
			
			includeGenes <- matrix(c(includeGenes,
					includeScores,
					rep(NA_real_,
						length(includeScores)*(ncol(ans) - ncol(includeGenes) - 1L))),
				nrow=length(includeScores))
			
			ans <- cbind(rbind(ans, includeGenes),
				FractionReps=bootstraps/(i + 1),
				Gene=gene)
		} else {
			subindices <- indices[indices > nrow(ORFs)] - nrow(ORFs)
			includeGenes <- matrix(c(includeGenes[subindices,],
					includeScores[subindices],
					rep(NA_real_,
						length(subindices)*(ncol(ans) - ncol(includeGenes) - 1L))),
				nrow=length(subindices))
			
			ans <- ans[indices[indices <= nrow(ORFs)],]
			ans <- cbind(rbind(ans,
					includeGenes),
				FractionReps=bootstraps[indices]/(i + 1),
				Gene=c(rep(1, nrow(ans)),
					includeGene[subindices]))
		}
	} else {
		if (allScores) {
			gene <- numeric(nrow(ans))
			gene[indices] <- 1
			ans <- cbind(ans,
				FractionReps=bootstraps/(i + 1),
				Gene=gene)
		} else {
			ans <- ans[indices,]
			ans <- cbind(ans,
				FractionReps=bootstraps[indices]/(i + 1),
				Gene=rep(1, nrow(ans)))
		}
	}
	
	o <- order(ans[, "Index"], ans[, "Begin"])
	ans <- ans[o,]
	rownames(ans) <- seq_len(nrow(ans))
	
	class(ans) <- "Genes"
	attr(ans, "widths") <- setNames(width(myDNAStringSet),
		names(myDNAStringSet))
	attr(ans, "geneticCode") <- gC
	attr(ans, "minGeneLength") <- minGeneLength
	attr(ans, "allowEdges") <- allowEdges
	attr(ans, "annotations") <- c(annotations,
		setNames(0:1,
			c(NA_character_, "Protein coding gene")))
	
	if (verbose) {
		cat("\n\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(ans)
}

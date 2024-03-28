AlignPairs <- function(pattern,
	subject,
	pairs=NULL,
	type="values",
	perfectMatch=2,
	misMatch=-2,
	gapOpening=-16,
	gapExtension=-1.2,
	substitutionMatrix=NULL,
	bandWidth=50,
	dropScore=-100,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (is(pattern, "DNAStringSet")) {
		if (!is(subject, "DNAStringSet"))
			stop("pattern and subject must be of the same class.")
		xtype <- 1L
	} else if (is(pattern, "RNAStringSet")) {
		if (!is(subject, "RNAStringSet"))
			stop("pattern and subject must be of the same class.")
		xtype <- 2L
	} else if (is(pattern, "AAStringSet")) {
		xtype <- 3L
		if (!is(subject, "AAStringSet"))
			stop("pattern and subject must be of the same class.")
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "J", "X", "*")
	} else {
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	if (length(pattern) < 1)
		stop("At least one sequence is required in the pattern.")
	if (length(subject) < 1)
		stop("At least one sequence is required in the subject.")
	if (is.null(pairs)) {
		if (length(pattern) != length(subject))
			stop("pattern and subject must be the same length when pairs is NULL.")
		pairs <- data.frame(Pattern=seq_along(pattern), Subject=seq_along(subject))
		pairs$Position <- rep(list(matrix(integer(), nrow=4L)), length(pattern))
	} else {
		if (!is.data.frame(pairs))
			stop("pairs must be a data frame.")
		if (is.na(match("Pattern", colnames(pairs))))
			stop("pairs must contain a column named 'Pattern'.")
		if (is.na(match("Subject", colnames(pairs))))
			stop("pairs must contain a column named 'Subject'.")
		if (!is.numeric(pairs$Pattern))
			stop("The 'Pattern' column of pairs must be a numeric.")
		if (!is.numeric(pairs$Subject))
			stop("The 'Subject' column of pairs must be a numeric.")
		rng <- range(pairs$Pattern)
		if (rng[1L] < 1L)
			stop("pairs contains a 'Pattern' index less than 1.")
		if (rng[1L] > length(pattern))
			stop("pairs contains a 'Pattern' index greater than the length of pattern.")
		rng <- range(pairs$Subject)
		if (rng[1L] < 1L)
			stop("pairs contains a 'Subject' index less than 1.")
		if (rng[1L] > length(subject))
			stop("pairs contains a 'Subject' index greater than the length of subject.")
		if (is.null(pairs$Position)) {
			pairs$Position <- rep(list(matrix(integer(), nrow=4L)), length(pairs$Pattern))
		} else {
			if (!is.list(pairs$Position))
				stop("The 'Position' column of pairs must be a list.")
		}
	}
	TYPES <- c("values", "sequences", "both")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (length(gapOpening) > 1L)
		stop("gapOpening must be a single number.")
	if (is.na(gapOpening))
		stop("gapOpening cannot be NA.")
	if (gapOpening > 0)
		stop("gapOpening can be at most zero.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (length(gapExtension) > 1L)
		stop("gapExtension must be a single number.")
	if (is.na(gapExtension))
		stop("gapExtension cannot be NA.")
	if (gapExtension > 0)
		stop("gapExtension can be at most zero.")
	if (is.null(substitutionMatrix)) {
		if (xtype == 3L) { # AAStringSet
			# use PFASUM50
			substitutionMatrix <- matrix(c(4.1181,-1.1516,-1.3187,-1.4135,0.4271,-0.5467,-0.6527,0.1777,-1.6582,-1.1243,-1.1843,-1.0235,-0.5685,-1.9515,-0.6072,0.8284,0.0361,-2.5368,-2.1701,0.0661,-4.3724,-3.6123,-4.1617,-13.4661,-11,-1.1516,6.341,0.0543,-0.6628,-3.2085,1.6006,0.5067,-1.961,0.7706,-3.5053,-3.0357,2.938,-1.9894,-3.7846,-1.3455,-0.4194,-0.5594,-2.1629,-1.7957,-2.9403,-3.3391,-2.0471,-6.2059,-13.5008,-11,-1.3187,0.0543,6.4672,2.3024,-2.5179,0.8192,0.5566,0.1585,1.104,-4.1629,-4.0977,0.8743,-2.6216,-3.805,-1.0904,1.1291,0.3253,-3.7763,-1.874,-3.6076,-0.1568,-2.3424,-7.1221,-13.3773,-11,-1.4135,-0.6628,2.3024,6.8156,-4.358,0.6705,2.582,-0.5667,-0.196,-5.475,-5.1661,0.226,-3.9595,-5.3456,-0.5662,0.4273,-0.5218,-4.7691,-3.4644,-4.5477,0.2269,-1.0448,-8.2793,-13.568,-11,0.4271,-3.2085,-2.5179,-4.358,13.5349,-3.3641,-4.3086,-2.1614,-1.8945,-0.7546,-0.9453,-3.8239,-0.5923,-0.8182,-3.6019,-0.3927,-0.801,-1.9317,-1.1607,0.0673,-6.4678,-6.9273,-3.8728,-13.6649,-11,-0.5467,1.6006,0.8192,0.6705,-3.3641,5.5795,2.1372,-1.5923,1.0862,-3.3001,-2.7545,1.872,-1.1216,-3.6631,-1.0426,0.1982,-0.0434,-3.061,-1.9214,-2.6993,-2.2648,-0.6904,-5.9511,-13.2661,-11,-0.6527,0.5067,0.5566,2.582,-4.3086,2.1372,5.5684,-1.6462,-0.2488,-4.1849,-4.0275,1.4821,-2.7964,-4.8311,-0.7028,0.0283,-0.312,-4.1969,-2.9489,-3.281,-1.1782,-0.5894,-7.0859,-13.4713,-11,0.1777,-1.961,0.1585,-0.5667,-2.1614,-1.5923,-1.6462,7.6508,-1.8185,-4.7058,-4.4215,-1.5991,-3.2786,-3.9992,-1.4409,0.184,-1.4823,-3.8328,-3.7343,-3.7264,-3.2391,-4.6258,-7.5259,-13.9053,-11,-1.6582,0.7706,1.104,-0.196,-1.8945,1.0862,-0.2488,-1.8185,9.7543,-3.3812,-2.8685,0.1425,-1.8724,-1.2545,-1.5333,-0.4285,-0.8896,-0.9385,1.6476,-2.8729,-2.5873,-2.6952,-6.0537,-13.4251,-11,-1.1243,-3.5053,-4.1629,-5.475,-0.7546,-3.3001,-4.1849,-4.7058,-3.3812,5.1229,2.5319,-3.5454,1.8309,0.9346,-3.4603,-3.0985,-1.2543,-1.5006,-1.117,3.3961,-7.8602,-6.8291,-0.617,-13.4445,-11,-1.1843,-3.0357,-4.0977,-5.1661,-0.9453,-2.7545,-4.0275,-4.4215,-2.8685,2.5319,4.7049,-3.4581,2.5303,1.687,-3.365,-3.1578,-1.8626,-0.5308,-0.6881,1.4829,-7.673,-6.5019,-0.9662,-13.7126,-11,-1.0235,2.938,0.8743,0.226,-3.8239,1.872,1.4821,-1.5991,0.1425,-3.5454,-3.4581,5.5476,-2.164,-4.3516,-0.7583,0.0275,-0.1516,-3.5889,-2.4422,-3.0453,-2.4826,-1.3665,-6.4906,-13.3896,-11,-0.5685,-1.9894,-2.6216,-3.9595,-0.5923,-1.1216,-2.7964,-3.2786,-1.8724,1.8309,2.5303,-2.164,7.0856,1.2339,-3.0823,-1.7587,-0.7402,-0.5841,-0.3946,0.9477,-6.3316,-5.0858,-0.7189,-13.2071,-11,-1.9515,-3.7846,-3.805,-5.3456,-0.8182,-3.6631,-4.8311,-3.9992,-1.2545,0.9346,1.687,-4.3516,1.2339,7.4322,-3.6222,-3.0316,-2.2851,2.6305,3.8302,0.1942,-7.6136,-7.3522,-1.58,-13.5461,-11,-0.6072,-1.3455,-1.0904,-0.5662,-3.6019,-1.0426,-0.7028,-1.4409,-1.5333,-3.4603,-3.365,-0.7583,-3.0823,-3.6222,9.1796,-0.0652,-0.8587,-3.3634,-3.3006,-2.5443,-3.7844,-3.828,-6.4005,-13.7659,-11,0.8284,-0.4194,1.1291,0.4273,-0.3927,0.1982,0.0283,0.184,-0.4285,-3.0985,-3.1578,0.0275,-1.7587,-3.0316,-0.0652,4.2366,1.8491,-3.1454,-2.1838,-2.1839,-2.2562,-2.9067,-6.1355,-13.3514,-11,0.0361,-0.5594,0.3253,-0.5218,-0.801,-0.0434,-0.312,-1.4823,-0.8896,-1.2543,-1.8626,-0.1516,-0.7402,-2.2851,-0.8587,1.8491,4.8833,-2.8511,-1.8993,-0.2699,-3.1362,-3.2086,-4.6244,-13.3413,-11,-2.5368,-2.1629,-3.7763,-4.7691,-1.9317,-3.061,-4.1969,-3.8328,-0.9385,-1.5006,-0.5308,-3.5889,-0.5841,2.6305,-3.3634,-3.1454,-2.8511,13.6485,3.3017,-1.851,-7.313,-6.7322,-3.8695,-13.7409,-11,-2.1701,-1.7957,-1.874,-3.4644,-1.1607,-1.9214,-2.9489,-3.7343,1.6476,-1.117,-0.6881,-2.4422,-0.3946,3.8302,-3.3006,-2.1838,-1.8993,3.3017,8.7568,-1.2438,-5.7064,-5.5317,-3.844,-13.5781,-11,0.0661,-2.9403,-3.6076,-4.5477,0.0673,-2.6993,-3.281,-3.7264,-2.8729,3.3961,1.4829,-3.0453,0.9477,0.1942,-2.5443,-2.1839,-0.2699,-1.851,-1.2438,4.6928,-7.1173,-6.052,-0.6976,-13.4507,-11,-4.3724,-3.3391,-0.1568,0.2269,-6.4678,-2.2648,-1.1782,-3.2391,-2.5873,-7.8602,-7.673,-2.4826,-6.3316,-7.6136,-3.7844,-2.2562,-3.1362,-7.313,-5.7064,-7.1173,0.7734,-4.5572,-10.7422,-16.4853,-11,-3.6123,-2.0471,-2.3424,-1.0448,-6.9273,-0.6904,-0.5894,-4.6258,-2.6952,-6.8291,-6.5019,-1.3665,-5.0858,-7.3522,-3.828,-2.9067,-3.2086,-6.7322,-5.5317,-6.052,-4.5572,-0.1769,-9.6217,-16.3932,-11,-4.1617,-6.2059,-7.1221,-8.2793,-3.8728,-5.9511,-7.0859,-7.5259,-6.0537,-0.617,-0.9662,-6.4906,-0.7189,-1.58,-6.4005,-6.1355,-4.6244,-3.8695,-3.844,-0.6976,-10.7422,-9.6217,-0.918,-16.6105,-11,-13.4661,-13.5008,-13.3773,-13.568,-13.6649,-13.2661,-13.4713,-13.9053,-13.4251,-13.4445,-13.7126,-13.3896,-13.2071,-13.5461,-13.7659,-13.3514,-13.3413,-13.7409,-13.5781,-13.4507,-16.4853,-16.3932,-16.6105,-20.4804,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,14),
				nrow=length(AAs),
				ncol=length(AAs),
				dimnames=list(AAs, AAs))
		} else {
			if (!is.numeric(perfectMatch))
				stop("perfectMatch must be a numeric.")
			if (!is.numeric(misMatch))
				stop("misMatch must be a numeric.")
			substitutionMatrix <- nucleotideSubstitutionMatrix(match=perfectMatch,
				mismatch=misMatch,
				type=ifelse(xtype == 1L, "DNA", "RNA"))
		}
	} else if (is.character(substitutionMatrix)) {
		if (xtype == 3L) { # AAStringSet
			substitutionMatrix <- .getSubMatrix(substitutionMatrix)
		} else {
			stop("substitutionMatrix must be NULL or a matrix.")
		}
	} else if (is(substitutionMatrix, "matrix")) {
		if (nrow(substitutionMatrix) != ncol(substitutionMatrix))
			stop("substitutionMatrix must be a square matrix.")
		if (is.null(rownames(substitutionMatrix)))
			stop("substitutionMatrix must have named rows.")
		if (is.null(colnames(substitutionMatrix)))
			stop("substitutionMatrix must have named columns.")
		if (any(rownames(substitutionMatrix) != colnames(substitutionMatrix)))
			stop("Row and column names of substitutionMatrix must be the same.")
	} else {
		stop("Invalid substitutionMatrix.")
	}
	if (!is.numeric(bandWidth))
		stop("bandWidth must be a number.")
	if (length(bandWidth) != 1L)
		stop("bandWidth must be a single number.")
	if (is.na(bandWidth))
		stop("bandWidth cannot be NA.")
	if (bandWidth < 4)
		stop("bandWidth must be at least 4.")
	if (bandWidth > 10000)
		stop("bandWidth can be at most 10000.")
	if (bandWidth != floor(bandWidth))
		stop("bandWidth must be a whole number.")
	if (!is.numeric(dropScore))
		stop("dropScore must be a number.")
	if (length(dropScore) != 1L)
		stop("dropScore must be a single number.")
	if (is.na(dropScore))
		stop("dropScore cannot be NA.")
	if (dropScore >= gapOpening)
		stop("dropScore must be less than ", gapOpening, ".")
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be TRUE or FALSE.")
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
	
	# initialize variables
	lkup <- paste(rownames(substitutionMatrix), collapse="")
	if (xtype == 1L) {
		lkup <- DNAStringSet(lkup)
		matchMatrix <- nucleotideSubstitutionMatrix(type="DNA") > 0
	} else if (xtype == 2L) {
		lkup <- RNAStringSet(lkup)
		matchMatrix <- nucleotideSubstitutionMatrix(type="RNA") > 0
	} else { # xtype == 3L
		lkup <- AAStringSet(lkup)
		matchMatrix <- matrix(FALSE, nrow=length(AAs), ncol=length(AAs), dimnames=list(AAs, AAs))
		diag(matchMatrix) <- TRUE
		matchMatrix["B", c("N", "D")] <- matchMatrix[c("N", "D"), "B"] <- TRUE
		matchMatrix["Z", c("Q", "E")] <- matchMatrix[c("Q", "E"), "Z"] <- TRUE
		matchMatrix["J", c("I", "L")] <- matchMatrix[c("I", "L"), "J"] <- TRUE
		matchMatrix["X", AAs[-length(AAs)]] <- matchMatrix[AAs[-length(AAs)], "X"] <- TRUE
	}
	m <- match(rownames(substitutionMatrix), rownames(matchMatrix))
	if (any(is.na(m)))
		stop("Unexpected letters in substitutionMatrix: '", paste(rownames(substitutionMatrix)[which(is.na(m))], collapse="', '"), "'.")
	matchMatrix <- matchMatrix[rownames(substitutionMatrix), colnames(substitutionMatrix)]
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
	} else {
		pBar <- NULL
	}
	
	ans <- .Call("alignPairs",
		pattern,
		subject,
		pairs$Pattern,
		pairs$Subject,
		pairs$Position,
		bandWidth,
		gapOpening,
		gapExtension,
		dropScore,
		substitutionMatrix,
		matchMatrix,
		lkup,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (type != 2L) {
		results <- data.frame(Pattern=pairs$Pattern,
			PatternStart=ans[[1L]],
			PatternEnd=ans[[2L]],
			Subject=pairs$Subject,
			SubjectStart=ans[[3L]],
			SubjectEnd=ans[[4L]],
			Matches=ans[[5L]],
			Mismatches=ans[[6L]],
			AlignmentLength=ans[[7L]],
			Score=ans[[8L]])
		results$PatternGapPosition <- ans[[9L]]
		results$PatternGapLength <- ans[[10L]]
		results$SubjectGapPosition <- ans[[11L]]
		results$SubjectGapLength <- ans[[12L]]
	}
	
	if (type > 1L) {
		gaps <- tabulate(unlist(c(ans[[10L]], ans[[12L]])))
		gaps <- lapply(seq_along(gaps)*(gaps > 0),
			function(l)
				if (l > 0) {
					paste(rep("-", l), collapse="")
				} else {
					""
				})
		patterns <- replaceAt(subseq(pattern[pairs$Pattern],
				ans[[1L]],
				ans[[2L]]),
			ans[[9L]],
			lapply(ans[[10L]],
				function(x)
					unlist(gaps[x])))
		subjects <- replaceAt(subseq(subject[pairs$Subject],
				ans[[3L]],
				ans[[4L]]),
			ans[[11L]],
			lapply(ans[[12L]],
				function(x)
					unlist(gaps[x])))
		if (type == 2L) {
			results <- list(patterns, subjects)
		} else if (type == 3L) {
			results <- list(results, patterns, subjects)
		}
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(results)
}

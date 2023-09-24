AdjustAlignment <- function(myXStringSet,
	perfectMatch=5,
	misMatch=0,
	gapLetter=-3,
	gapOpening=-0.1,
	gapExtension=0,
	substitutionMatrix=NULL,
	shiftPenalty=-0.2,
	threshold=0.1,
	weight=1,
	processors=1) {
	
	# error checking
	if (is(myXStringSet, "DNAStringSet")) {
		type <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		type <- 2L
	} else if (is(myXStringSet, "AAStringSet")) {
		type <- 3L
	} else {
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	if (length(myXStringSet) < 2)
		return(myXStringSet)
	u <- unique(width(myXStringSet))
	if (length(u) != 1)
		stop("Sequences in myXStringSet must be the same width (aligned).")
	if (u < 4) # no changes can be made
		return(myXStringSet)
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (perfectMatch < 0)
		stop("perfectMatch must be at least zero.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (misMatch > 0)
		stop("misMatch must be less than or equal to zero.")
	if (!is.numeric(gapLetter))
		stop("gapLetter must be a numeric.")
	if (gapLetter > 0)
		stop("gapLetter must be less than or equal to zero.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (gapOpening > 0)
		stop("gapOpening must be less than or equal to zero.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (gapExtension > 0)
		stop("gapExtension must be less than or equal to zero.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold < 0)
		stop("threshold must be at least zero.")
	if (!is.numeric(shiftPenalty))
		stop("shiftPenalty must be a numeric.")
	if (shiftPenalty > 0)
		stop("shiftPenalty must be less than or equal to zero.")
	if (!is.numeric(weight))
		stop("weight must be a numeric.")
	if (length(weight) != 1 && length(weight) != length(myXStringSet))
		stop("Length of weight must equal one or the length of the myXStringSet.")
	if (length(weight) == 1) {
		weight <- rep(1, length(myXStringSet))
	} else {
		if (!isTRUE(all.equal(1, mean(weight))))
			stop("The mean of weight must be 1.")
	}
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
	
	if (type == 3L) { # AAStringSet
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (is.null(substitutionMatrix)) {
			# use PFASUM50
			substitutionMatrix <- matrix(c(4.1181,-1.1516,-1.3187,-1.4135,0.4271,-0.5467,-0.6527,0.1777,-1.6582,-1.1243,-1.1843,-1.0235,-0.5685,-1.9515,-0.6072,0.8284,0.0361,-2.5368,-2.1701,0.0661,-11,-1.1516,6.341,0.0543,-0.6628,-3.2085,1.6006,0.5067,-1.961,0.7706,-3.5053,-3.0357,2.938,-1.9894,-3.7846,-1.3455,-0.4194,-0.5594,-2.1629,-1.7957,-2.9403,-11,-1.3187,0.0543,6.4672,2.3024,-2.5179,0.8192,0.5566,0.1585,1.104,-4.1629,-4.0977,0.8743,-2.6216,-3.805,-1.0904,1.1291,0.3253,-3.7763,-1.874,-3.6076,-11,-1.4135,-0.6628,2.3024,6.8156,-4.358,0.6705,2.582,-0.5667,-0.196,-5.475,-5.1661,0.226,-3.9595,-5.3456,-0.5662,0.4273,-0.5218,-4.7691,-3.4644,-4.5477,-11,0.4271,-3.2085,-2.5179,-4.358,13.5349,-3.3641,-4.3086,-2.1614,-1.8945,-0.7546,-0.9453,-3.8239,-0.5923,-0.8182,-3.6019,-0.3927,-0.801,-1.9317,-1.1607,0.0673,-11,-0.5467,1.6006,0.8192,0.6705,-3.3641,5.5795,2.1372,-1.5923,1.0862,-3.3001,-2.7545,1.872,-1.1216,-3.6631,-1.0426,0.1982,-0.0434,-3.061,-1.9214,-2.6993,-11,-0.6527,0.5067,0.5566,2.582,-4.3086,2.1372,5.5684,-1.6462,-0.2488,-4.1849,-4.0275,1.4821,-2.7964,-4.8311,-0.7028,0.0283,-0.312,-4.1969,-2.9489,-3.281,-11,0.1777,-1.961,0.1585,-0.5667,-2.1614,-1.5923,-1.6462,7.6508,-1.8185,-4.7058,-4.4215,-1.5991,-3.2786,-3.9992,-1.4409,0.184,-1.4823,-3.8328,-3.7343,-3.7264,-11,-1.6582,0.7706,1.104,-0.196,-1.8945,1.0862,-0.2488,-1.8185,9.7543,-3.3812,-2.8685,0.1425,-1.8724,-1.2545,-1.5333,-0.4285,-0.8896,-0.9385,1.6476,-2.8729,-11,-1.1243,-3.5053,-4.1629,-5.475,-0.7546,-3.3001,-4.1849,-4.7058,-3.3812,5.1229,2.5319,-3.5454,1.8309,0.9346,-3.4603,-3.0985,-1.2543,-1.5006,-1.117,3.3961,-11,-1.1843,-3.0357,-4.0977,-5.1661,-0.9453,-2.7545,-4.0275,-4.4215,-2.8685,2.5319,4.7049,-3.4581,2.5303,1.687,-3.365,-3.1578,-1.8626,-0.5308,-0.6881,1.4829,-11,-1.0235,2.938,0.8743,0.226,-3.8239,1.872,1.4821,-1.5991,0.1425,-3.5454,-3.4581,5.5476,-2.164,-4.3516,-0.7583,0.0275,-0.1516,-3.5889,-2.4422,-3.0453,-11,-0.5685,-1.9894,-2.6216,-3.9595,-0.5923,-1.1216,-2.7964,-3.2786,-1.8724,1.8309,2.5303,-2.164,7.0856,1.2339,-3.0823,-1.7587,-0.7402,-0.5841,-0.3946,0.9477,-11,-1.9515,-3.7846,-3.805,-5.3456,-0.8182,-3.6631,-4.8311,-3.9992,-1.2545,0.9346,1.687,-4.3516,1.2339,7.4322,-3.6222,-3.0316,-2.2851,2.6305,3.8302,0.1942,-11,-0.6072,-1.3455,-1.0904,-0.5662,-3.6019,-1.0426,-0.7028,-1.4409,-1.5333,-3.4603,-3.365,-0.7583,-3.0823,-3.6222,9.1796,-0.0652,-0.8587,-3.3634,-3.3006,-2.5443,-11,0.8284,-0.4194,1.1291,0.4273,-0.3927,0.1982,0.0283,0.184,-0.4285,-3.0985,-3.1578,0.0275,-1.7587,-3.0316,-0.0652,4.2366,1.8491,-3.1454,-2.1838,-2.1839,-11,0.0361,-0.5594,0.3253,-0.5218,-0.801,-0.0434,-0.312,-1.4823,-0.8896,-1.2543,-1.8626,-0.1516,-0.7402,-2.2851,-0.8587,1.8491,4.8833,-2.8511,-1.8993,-0.2699,-11,-2.5368,-2.1629,-3.7763,-4.7691,-1.9317,-3.061,-4.1969,-3.8328,-0.9385,-1.5006,-0.5308,-3.5889,-0.5841,2.6305,-3.3634,-3.1454,-2.8511,13.6485,3.3017,-1.851,-11,-2.1701,-1.7957,-1.874,-3.4644,-1.1607,-1.9214,-2.9489,-3.7343,1.6476,-1.117,-0.6881,-2.4422,-0.3946,3.8302,-3.3006,-2.1838,-1.8993,3.3017,8.7568,-1.2438,-11,0.0661,-2.9403,-3.6076,-4.5477,0.0673,-2.6993,-3.281,-3.7264,-2.8729,3.3961,1.4829,-3.0453,0.9477,0.1942,-2.5443,-2.1839,-0.2699,-1.851,-1.2438,4.6928,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,14),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		} else if (is.character(substitutionMatrix)) {
			if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")))
				stop("Invalid substitutionMatrix.")
		}
		if (is.matrix(substitutionMatrix)) {
			if (any(!(AAs %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(AAs %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix
		} else {
			subMatrix <- eval(parse(text=data(list=substitutionMatrix, envir=environment(), package=ifelse(substitutionMatrix == "MIQS", "DECIPHER", "Biostrings"))))
		}
		subMatrix <- subMatrix[AAs, AAs]
		subMatrix <- as.numeric(subMatrix)
		
		functionCall <- "shiftGapsAA"
	} else { # DNAStringSet or RNAStringSet
		if (is.matrix(substitutionMatrix)) {
			bases <- c("A", "C", "G",
				ifelse(type == 2L, "U", "T"))
			if (any(!(bases %in% dimnames(substitutionMatrix)[[1]])) ||
				any(!(bases %in% dimnames(substitutionMatrix)[[2]])))
				stop("substitutionMatrix is incomplete.")
			subMatrix <- substitutionMatrix[bases, bases]
		} else if (is.null(substitutionMatrix)) {
			subMatrix <- matrix(misMatch,
				nrow=4, ncol=4)
			diag(subMatrix) <- perfectMatch
		} else {
			stop("substitutionMatrix must be NULL or a matrix.")
		}
		
		functionCall <- "shiftGaps"
	}
	
	# add gaps to both ends of the sequences
	ns <- names(myXStringSet)
	myXStringSet <- .Call("insertGaps",
		myXStringSet,
		c(1L, u + 1L),
		c(1L, 1L),
		type,
		processors,
		PACKAGE="DECIPHER")
	
	# adjust the alignment
	myXStringSet <- .Call(functionCall,
		myXStringSet, # in-place change of myXStringSet (requires previous temporary copy)
		as.numeric(subMatrix),
		gapOpening,
		gapExtension,
		gapLetter,
		shiftPenalty,
		threshold,
		weight,
		PACKAGE="DECIPHER")
	
	# remove all 100% gap columns
	myXStringSet <- .Call("removeCommonGaps",
		myXStringSet,
		type,
		FALSE, # includeMask
		processors,
		PACKAGE="DECIPHER")
	
	names(myXStringSet) <- ns
	
	return(myXStringSet)
}

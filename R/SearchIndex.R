SearchIndex <- function(pattern,
	invertedIndex,
	subject=NULL,
	type="one",
	minScore=NA,
	scoreOnly=FALSE,
	sepCost=-0.4,
	gapCost=-2.5,
	maskRepeats=TRUE,
	maskLCRs=TRUE,
	dropScore=-10,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is(invertedIndex, "InvertedIndex"))
		stop("invertedIndex must be an InvertedIndex object.")
	TYPES <- c("all", "one", "top")
	type <- pmatch(type[1], TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (length(minScore) != 1L)
		stop("minScore must be a single numeric.")
	if (!is.na(minScore) && !is.numeric(minScore))
		stop("minScore must be a numeric.")
	if (!isTRUEorFALSE(scoreOnly))
		stop("scoreOnly must be TRUE or FALSE.")
	if (length(sepCost) != 1L)
		stop("sepCost must be a single numeric.")
	if (is.na(sepCost) || !is.numeric(sepCost))
		stop("sepCost must be a numeric.")
	if (length(gapCost) != 1L)
		stop("gapCost must be a single numeric.")
	if (is.na(gapCost) || !is.numeric(gapCost))
		stop("gapCost must be a numeric.")
	if (!isTRUEorFALSE(maskRepeats))
		stop("maskRepeats must be TRUE or FALSE.")
	if (!isTRUEorFALSE(maskLCRs))
		stop("maskLCRs must be TRUE or FALSE.")
	if (!isTRUEorFALSE(verbose))
		stop("verbose must be TRUE or FALSE.")
	if (is.na(dropScore) || !is.numeric(dropScore))
		stop("dropScore must be a numeric.")
	if (length(dropScore) != 1L)
		stop("dropScore must be a single numeric.")
	if (dropScore > 0)
		stop("dropScore must be less than or equal to zero.")
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
	
	K <- invertedIndex$k
	step <- invertedIndex$step
	alphabet <- invertedIndex$alphabet
	freqs <- invertedIndex$frequency
	num <- invertedIndex$count
	len <- invertedIndex$length
	loc <- invertedIndex$location
	ind <- invertedIndex$index
	
	if (length(alphabet) == 20L) {
		if (!is(pattern, "AAStringSet"))
			stop("pattern must be an AAStringSet.")
		xtype <- 3L
	} else if (length(alphabet) == 4L) {
		if (is(pattern, "DNAStringSet")) {
			xtype <- 1L
		} else if (is(pattern, "RNAStringSet")) {
			xtype <- 2L
		} else {
			stop("pattern must be a DNAStringSet or RNAStringSet.")
		}
	} else {
		stop("Corrupted invertedIndex.")
	}
	if (!is.null(subject)) {
		if (length(subject) != length(len))
			stop("subject is incompatible with the invertedIndex.")
		if (xtype == 1L) {
			if (!is(subject, "DNAStringSet"))
				stop("subject must be a DNAStringSet.")
			subMatrix <- nucleotideSubstitutionMatrix(match=1.386294, # log(4)
				mismatch=-5.545177, # -4*log(4) allows 20% distance
				type="DNA")
			subMatrix["A", "A"] <- -log(freqs["A"])
			subMatrix["C", "C"] <- -log(freqs["C"])
			subMatrix["G", "G"] <- -log(freqs["G"])
			subMatrix["T", "T"] <- -log(freqs[!(names(freqs) %in% c("A", "C", "G"))])
			chars <- DNAStringSet(paste(rownames(subMatrix), collapse=""))
		} else if (xtype == 2L) {
			if (!is(subject, "RNAStringSet"))
				stop("subject must be a RNAStringSet.")
			subMatrix <- nucleotideSubstitutionMatrix(match=1.386294, # log(4)
				mismatch=-5.545177, # -4*log(4) allows 20% distance
				type="RNA")
			subMatrix["A", "A"] <- -log(freqs["A"])
			subMatrix["C", "C"] <- -log(freqs["C"])
			subMatrix["G", "G"] <- -log(freqs["G"])
			subMatrix["U", "U"] <- -log(freqs[!(names(freqs) %in% c("A", "C", "G"))])
			chars <- RNAStringSet(paste(rownames(subMatrix), collapse=""))
		} else { # xtype == 3L
			if (!is(subject, "AAStringSet"))
				stop("subject must be an AAStringSet.")
			# use the full PFASUM90 matrix calibrated to expectation (matches ~= -log(freqs))
			subMatrix <- matrix(c(2.1531,-0.8109,-0.7551,-0.8281,0.0688,-0.4222,-0.4138,-0.083,-0.979,-0.7883,-0.8611,-0.6089,-0.5443,-1.2266,-0.4539,0.3974,-0.0709,-1.5202,-1.3301,-0.1153,-1.9435,-1.5628,-1.9793,-5.2371,-4.9659,-0.8109,3.0611,-0.2391,-0.7324,-1.6635,0.5087,-0.1707,-1.2636,0.1364,-1.8352,-1.5685,1.2628,-1.14,-2.0327,-1.0208,-0.4922,-0.5567,-1.2355,-1.1356,-1.5974,-1.6558,-1.0345,-2.8099,-5.2764,-4.9659,-0.7551,-0.2391,3.2811,0.8499,-1.2776,0.1577,-0.0353,-0.1961,0.4775,-2.0728,-2.0216,0.139,-1.3669,-1.9526,-0.8666,0.3896,-0.0149,-1.9669,-1.0592,-1.8046,0.3604,-1.1079,-3.1867,-5.1994,-4.9659,-0.8281,-0.7324,0.8499,3.1735,-2.2112,0.0287,1.0668,-0.5565,-0.3268,-2.7662,-2.5827,-0.2546,-2.0885,-2.7204,-0.6335,-0.0641,-0.5175,-2.4706,-1.8456,-2.2981,0.4893,-0.3887,-3.7953,-5.2878,-4.9659,0.0688,-1.6635,-1.2776,-2.2112,5.6006,-1.7383,-2.3082,-1.24,-1.1362,-0.6809,-0.7673,-2.0293,-0.5984,-0.7222,-1.9964,-0.1872,-0.4273,-1.1734,-0.8035,-0.1612,-2.9049,-3.2221,-1.8803,-5.3412,-4.9659,-0.4222,0.5087,0.1577,0.0287,-1.7383,2.9961,0.8029,-1.0276,0.4654,-1.6674,-1.3027,0.6752,-0.5506,-1.8842,-0.7148,-0.142,-0.255,-1.619,-1.1273,-1.3761,-1.0625,0.1206,-2.5766,-5.149,-4.9659,-0.4138,-0.1707,-0.0353,1.0668,-2.3082,0.8029,2.7201,-1.025,-0.3998,-2.1262,-2.0174,0.354,-1.5047,-2.5066,-0.6275,-0.2714,-0.414,-2.2105,-1.6607,-1.6411,-0.4516,0.1855,-3.2034,-5.2489,-4.9659,-0.083,-1.2636,-0.1961,-0.5565,-1.24,-1.0276,-1.025,3.3406,-1.1352,-2.5059,-2.3459,-1.0166,-1.8755,-2.1928,-1.1142,-0.1185,-0.9926,-2.0906,-2.1353,-2.0053,-1.5434,-2.172,-3.5503,-5.4503,-4.9659,-0.979,0.1364,0.4775,-0.3268,-1.1362,0.4654,-0.3998,-1.1352,4.3656,-1.7771,-1.4664,-0.1425,-1.1296,-0.6338,-0.9922,-0.4026,-0.6345,-0.6464,0.7216,-1.5571,-1.091,-1.1744,-2.7225,-5.2275,-4.9659,-0.7883,-1.8352,-2.0728,-2.7662,-0.6809,-1.6674,-2.1262,-2.5059,-1.7771,2.5416,0.9454,-1.8231,0.6983,0.0484,-1.873,-1.6653,-0.6662,-1.0454,-0.8996,1.4591,-3.5888,-3.0891,-0.0165,-5.2127,-4.9659,-0.8611,-1.5685,-2.0216,-2.5827,-0.7673,-1.3027,-2.0174,-2.3459,-1.4664,0.9454,2.2779,-1.6988,1.0466,0.501,-1.782,-1.653,-1.0118,-0.576,-0.6119,0.4308,-3.4726,-2.8646,-0.0818,-5.3335,-4.9659,-0.6089,1.2628,0.139,-0.2546,-2.0293,0.6752,0.354,-1.0166,-0.1425,-1.8231,-1.6988,2.7764,-1.1255,-2.2758,-0.6652,-0.2555,-0.2942,-1.8695,-1.4346,-1.5386,-1.2259,-0.667,-2.8904,-5.2074,-4.9659,-0.5443,-1.14,-1.3669,-2.0885,-0.5984,-0.5506,-1.5047,-1.8755,-1.1296,0.6983,1.0466,-1.1255,3.6813,0.2825,-1.7701,-1.0421,-0.4541,-0.6095,-0.5404,0.2329,-2.8964,-2.2349,-0.2219,-5.1362,-4.9659,-1.2266,-2.0327,-1.9526,-2.7204,-0.7222,-1.8842,-2.5066,-2.1928,-0.6338,0.0484,0.501,-2.2758,0.2825,3.5507,-2.0426,-1.6262,-1.3124,0.9849,1.6772,-0.2784,-3.504,-3.3967,-0.8011,-5.2912,-4.9659,-0.4539,-1.0208,-0.8666,-0.6335,-1.9964,-0.7148,-0.6275,-1.1142,-0.9922,-1.873,-1.782,-0.6652,-1.7701,-2.0426,4.0183,-0.2166,-0.6586,-1.9473,-1.8794,-1.4024,-1.8721,-1.8049,-2.9616,-5.4065,-4.9659,0.3974,-0.4922,0.3896,-0.0641,-0.1872,-0.142,-0.2714,-0.1185,-0.4026,-1.6653,-1.653,-0.2555,-1.0421,-1.6262,-0.2166,2.3526,0.8047,-1.7097,-1.2892,-1.1841,-1.0067,-1.3689,-2.8036,-5.1774,-4.9659,-0.0709,-0.5567,-0.0149,-0.5175,-0.4273,-0.255,-0.414,-0.9926,-0.6345,-0.6662,-1.0118,-0.2942,-0.4541,-1.3124,-0.6586,0.8047,2.6474,-1.6094,-1.2019,-0.1755,-1.4363,-1.5,-2.0191,-5.1852,-4.9659,-1.5202,-1.2355,-1.9669,-2.4706,-1.1734,-1.619,-2.2105,-2.0906,-0.6464,-1.0454,-0.576,-1.8695,-0.6095,0.9849,-1.9473,-1.7097,-1.6094,5.8259,1.2574,-1.2247,-3.3889,-3.1146,-1.8833,-5.3907,-4.9659,-1.3301,-1.1356,-1.0592,-1.8456,-0.8035,-1.1273,-1.6607,-2.1353,0.7216,-0.8996,-0.6119,-1.4346,-0.5404,1.6772,-1.8794,-1.2892,-1.2019,1.2574,4.0647,-0.9588,-2.6193,-2.5909,-1.8604,-5.3107,-4.9659,-0.1153,-1.5974,-1.8046,-2.2981,-0.1612,-1.3761,-1.6411,-2.0053,-1.5571,1.4591,0.4308,-1.5386,0.2329,-0.2784,-1.4024,-1.1841,-0.1755,-1.2247,-0.9588,2.3577,-3.2213,-2.685,-0.2505,-5.222,-4.9659,-1.9435,-1.6558,0.3604,0.4893,-2.9049,-1.0625,-0.4516,-1.5434,-1.091,-3.5888,-3.4726,-1.2259,-2.8964,-3.504,-1.8721,-1.0067,-1.4363,-3.3889,-2.6193,-3.2213,0.9658,-1.7961,-4.6614,-6.3968,-4.9659,-1.5628,-1.0345,-1.1079,-0.3887,-3.2221,0.1206,0.1855,-2.172,-1.1744,-3.0891,-2.8646,-0.667,-2.2349,-3.3967,-1.8049,-1.3689,-1.5,-3.1146,-2.5909,-2.685,-1.7961,0.6161,-4.0915,-6.3578,-4.9659,-1.9793,-2.8099,-3.1867,-3.7953,-1.8803,-2.5766,-3.2034,-3.5503,-2.7225,-0.0165,-0.0818,-2.8904,-0.2219,-0.8011,-2.9616,-2.8036,-2.0191,-1.8833,-1.8604,-0.2505,-4.6614,-4.0915,0.1589,-6.4332,-4.9659,-5.2371,-5.2764,-5.1994,-5.2878,-5.3412,-5.149,-5.2489,-5.4503,-5.2275,-5.2127,-5.3335,-5.2074,-5.1362,-5.2912,-5.4065,-5.1774,-5.1852,-5.3907,-5.3107,-5.222,-6.3968,-6.3578,-6.4332,-7.9243,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,-4.9659,6.1118),
				nrow=25L)
			chars <- AAStringSet("ARNDCQEGHILKMFPSTWYVBZJX*")
		}
	} else {
		subMatrix <- numeric()
		chars <- NULL
	}
	
	# initialize a progress bar
	if (verbose) {
		pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	} else {
		pBar <- NULL
	}
	
	L <- length(num)
	maxSep <- as.integer(sqrt(L))
	
	if (xtype == 3L) {
		kmers <- .Call("enumerateSequenceReducedAA",
			pattern,
			K,
			alphabet,
			maskRepeats,
			maskLCRs,
			1L, # left is fast moving side
			processors,
			PACKAGE="DECIPHER")
	} else {
		kmers <- .Call("enumerateSequence",
			pattern,
			K,
			maskRepeats,
			maskLCRs,
			1L, # left is fast moving side
			processors,
			PACKAGE="DECIPHER")
	}
	
	ans <- .Call("searchIndex",
		kmers, # query k-mers [0 to length(num)]
		K, # wordSize
		step, # separation between k-mers (>= 1 and <= K)
		-log(freqs), # -log of normalized letter frequencies
		num, # count
		loc, # location
		ind, # index
		len, # positions
		sepCost, # sepC
		gapCost, # gapC
		type, # output
		sum(len) + 1, # size of target database
		minScore, # minimum score or NA to calculate
		scoreOnly, # FALSE to output anchor positions
		pattern, # optional query sequences
		subject, # optional target sequences
		subMatrix, # substitution matrix
		chars, # concatenated row/column names of subMatrix
		dropScore, # to terminate extension
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (!scoreOnly) {
		pos <- ans[[4L]]
		ans <- ans[1:3]
	}
	names(ans) <- c("Pattern", "Subject", "Score")
	ans <- data.frame(ans)
	if (!scoreOnly)
		ans$Position <- pos
	
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

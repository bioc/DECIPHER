FindSynteny <- function(dbFile,
	tblName="Seqs",
	identifier="",
	useFrames=TRUE,
	alphabet=AA_REDUCED[[172]],
	geneticCode=GENETIC_CODE,
	sepCost=-3,
	sepPower=0.5,
	gapCost=-12,
	gapPower=0.5,
	shiftCost=0,
	codingCost=0,
	maxSep=150,
	maxGap=15,
	minScore=100,
	N=10,
	dropScore=-6,
	maskRepeats=TRUE,
	maskLCRs=TRUE,
	allowOverlap=FALSE,
	storage=0.5,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.logical(useFrames))
		stop("useFrames must be a logical.")
	if (!is.numeric(sepCost))
		stop("sepCost must be a single numeric.")
	if (sepCost > 0)
		stop("sepCost must be less than or equal to zero.")
	if (!is.numeric(sepPower))
		stop("sepPower must be a single numeric.")
	if (sepPower <= 0)
		stop("sepPower must be greater than zero.")
	if (!is.numeric(gapCost))
		stop("gapCost must be a single numeric.")
	if (gapCost > 0)
		stop("gapCost must be less than or equal to zero.")
	if (!is.numeric(gapPower))
		stop("gapPower must be a single numeric.")
	if (gapPower <= 0)
		stop("gapPower must be greater than zero.")
	if (!is.numeric(shiftCost))
		stop("shiftCost must be a single numeric.")
	if (shiftCost > 0)
		stop("shiftCost must be less than or equal to zero.")
	if (!is.numeric(codingCost))
		stop("codingCost must be a single numeric.")
	if (codingCost > 0)
		stop("codingCost must be less than or equal to zero.")
	if (!is.numeric(maxSep))
		stop("maxSep must be a single numeric.")
	if (maxSep <= 0)
		stop("maxSep must be greater than zero.")
	if (!is.numeric(maxGap))
		stop("maxGap must be a single numeric.")
	if (maxGap < 0)
		stop("maxGap must be at least zero.")
	if (!is.numeric(dropScore))
		stop("dropScore must be a single numeric.")
	if (dropScore > 0)
		stop("dropScore must be less than or equal to zero.")
	if (is.na(N) || !is.numeric(N))
		stop("N must be a numeric.")
	if (length(N) > 1L)
		stop("N must be a single numeric.")
	if (N < 1)
		stop("N must be at least one.")
	if (!is.numeric(minScore))
		stop("minScore must be a single numeric.")
	if (minScore <= 0)
		stop("minScore must be greater than zero.")
	if (!is.logical(maskRepeats))
		stop("maskRepeats must be a logical.")
	if (!is.logical(maskLCRs))
		stop("maskLCRs must be a logical.")
	if (!is.logical(allowOverlap))
		stop("allowOverlap must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.numeric(storage))
		stop("storage must be a numeric.")
	if (storage < 0)
		stop("storage must be at least zero.")
	storage <- storage*1e9 # convert to bytes
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
	
	if (verbose)
		time.1 <- Sys.time()
	
	# initialize database
	if (is.character(dbFile)) {
		if (!requireNamespace("RSQLite", quietly=TRUE))
			stop("Package 'RSQLite' must be installed.")
		dbConn <- dbConnect(dbDriver("SQLite"), dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn <- dbFile
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	
	ids <- dbGetQuery(dbConn,
		paste("select distinct",
			dbQuoteIdentifier(dbConn, "identifier"),
			"from",
			dbQuoteIdentifier(dbConn, tblName)))$identifier
	if (identifier[1] == "") {
		identifier <- ids
	} else {
		w <- which(!(identifier %in% ids))
		if (length(w) > 0)
			stop("identifier(s) not found in ",
				tblName,
				": ",
				paste(identifier[w],
					collapse=", "))
	}
	l <- length(identifier)
	if (l < 2)
		stop("At least two identifers are required in ",
			tblName,
			".")
	
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
	f <- function(gC) {
		m <- match(codons, names(gC))
		if (any(is.na(m)))
			stop("geneticCode must contain all 64 codons.")
		AAStringSet(paste(gC[m],
			collapse=""))
	}
	if (is.list(geneticCode)) {
		if (length(geneticCode) != l)
			stop("The list geneticCode must have as many items as the number of identifiers.")
		if (!is.null(names(geneticCode))) {
			m <- match(identifier, names(geneticCode))
			if (any(is.na(m)))
				stop("All identifiers must be present in the names of the list geneticCode.")
			geneticCode <- geneticCode[m]
			geneticCode <- lapply(geneticCode, f)
		}
	} else { # all identifiers use the same geneticCode
		geneticCode <- rep(list(f(geneticCode)), l)
	}
	
	if (!is.character(alphabet))
		stop("alphabet must be a character vector.")
	if (any(alphabet == ""))
		stop("No elements of alphabet can be empty.")
	r <- strsplit(alphabet, "", fixed=TRUE)
	alphabet <- setNames(rep(0L, 20),
		AA_STANDARD)
	for (i in seq_along(r)) {
		w <- which(!(r[[i]] %in% AA_STANDARD))
		if (length(w) > 0)
			stop("Unrecognized letter(s) found in alphabet:  ",
				paste(r[[i]][w], collapse=", "),
				".")
		w <- which(alphabet[r[[i]]] != 0L)
		if (length(w) > 0)
			stop("Repeated amino acids found in alphabet:  ",
				paste(r[[i]][w], collapse=", "),
				".")
		alphabet[r[[i]]] <- i
	}
	w <- which(alphabet == 0L)
	if (length(w) > 0)
		stop("Standard amino acids missing from alphabet:  ",
			paste(names(w), collapse=", "),
			".")
	sizeAA <- max(alphabet)
	if (sizeAA == 1)
		stop("More than one grouping of amino acids is required in the alphabet.")
	n <- as.integer(floor(log(4294967295, sizeAA)))
	alphabet <- alphabet - 1L
	
	# subMatrixDNA gets multiplied by log(size) to calibrate for DNA distribution
	# -3 allows 25% transition rate and -4 allows 20% transversion rate
	subMatrixDNA <- matrix(c(1, -4, -3, -4, -4, 1, -4, -3, -3, -4, 1, -4, -4, -3, -4, 1),
		nrow=4L)
	lettersDNA <- DNAStringSet("ACGT") # concatenated row/column names
	# subMatrixAA uses PFASUM90 in units of log-odds (PFASUM90*log(2)/3)
	# replace letter/* with dropScore to prevent extension beyond matched reading frames
	subMatrixAA <- matrix(c(2.1531,-0.8109,-0.7551,-0.8281,0.0688,-0.4222,-0.4138,-0.083,-0.979,-0.7883,-0.8611,-0.6089,-0.5443,-1.2266,-0.4539,0.3974,-0.0709,-1.5202,-1.3301,-0.1153,dropScore,-0.8109,3.0611,-0.2391,-0.7324,-1.6635,0.5087,-0.1707,-1.2636,0.1364,-1.8352,-1.5685,1.2628,-1.14,-2.0327,-1.0208,-0.4922,-0.5567,-1.2355,-1.1356,-1.5974,dropScore,-0.7551,-0.2391,3.2811,0.8499,-1.2776,0.1577,-0.0353,-0.1961,0.4775,-2.0728,-2.0216,0.139,-1.3669,-1.9526,-0.8666,0.3896,-0.0149,-1.9669,-1.0592,-1.8046,dropScore,-0.8281,-0.7324,0.8499,3.1735,-2.2112,0.0287,1.0668,-0.5565,-0.3268,-2.7662,-2.5827,-0.2546,-2.0885,-2.7204,-0.6335,-0.0641,-0.5175,-2.4706,-1.8456,-2.2981,dropScore,0.0688,-1.6635,-1.2776,-2.2112,5.6006,-1.7383,-2.3082,-1.24,-1.1362,-0.6809,-0.7673,-2.0293,-0.5984,-0.7222,-1.9964,-0.1872,-0.4273,-1.1734,-0.8035,-0.1612,dropScore,-0.4222,0.5087,0.1577,0.0287,-1.7383,2.9961,0.8029,-1.0276,0.4654,-1.6674,-1.3027,0.6752,-0.5506,-1.8842,-0.7148,-0.142,-0.255,-1.619,-1.1273,-1.3761,dropScore,-0.4138,-0.1707,-0.0353,1.0668,-2.3082,0.8029,2.7201,-1.025,-0.3998,-2.1262,-2.0174,0.354,-1.5047,-2.5066,-0.6275,-0.2714,-0.414,-2.2105,-1.6607,-1.6411,dropScore,-0.083,-1.2636,-0.1961,-0.5565,-1.24,-1.0276,-1.025,3.3406,-1.1352,-2.5059,-2.3459,-1.0166,-1.8755,-2.1928,-1.1142,-0.1185,-0.9926,-2.0906,-2.1353,-2.0053,dropScore,-0.979,0.1364,0.4775,-0.3268,-1.1362,0.4654,-0.3998,-1.1352,4.3656,-1.7771,-1.4664,-0.1425,-1.1296,-0.6338,-0.9922,-0.4026,-0.6345,-0.6464,0.7216,-1.5571,dropScore,-0.7883,-1.8352,-2.0728,-2.7662,-0.6809,-1.6674,-2.1262,-2.5059,-1.7771,2.5416,0.9454,-1.8231,0.6983,0.0484,-1.873,-1.6653,-0.6662,-1.0454,-0.8996,1.4591,dropScore,-0.8611,-1.5685,-2.0216,-2.5827,-0.7673,-1.3027,-2.0174,-2.3459,-1.4664,0.9454,2.2779,-1.6988,1.0466,0.501,-1.782,-1.653,-1.0118,-0.576,-0.6119,0.4308,dropScore,-0.6089,1.2628,0.139,-0.2546,-2.0293,0.6752,0.354,-1.0166,-0.1425,-1.8231,-1.6988,2.7764,-1.1255,-2.2758,-0.6652,-0.2555,-0.2942,-1.8695,-1.4346,-1.5386,dropScore,-0.5443,-1.14,-1.3669,-2.0885,-0.5984,-0.5506,-1.5047,-1.8755,-1.1296,0.6983,1.0466,-1.1255,3.6813,0.2825,-1.7701,-1.0421,-0.4541,-0.6095,-0.5404,0.2329,dropScore,-1.2266,-2.0327,-1.9526,-2.7204,-0.7222,-1.8842,-2.5066,-2.1928,-0.6338,0.0484,0.501,-2.2758,0.2825,3.5507,-2.0426,-1.6262,-1.3124,0.9849,1.6772,-0.2784,dropScore,-0.4539,-1.0208,-0.8666,-0.6335,-1.9964,-0.7148,-0.6275,-1.1142,-0.9922,-1.873,-1.782,-0.6652,-1.7701,-2.0426,4.0183,-0.2166,-0.6586,-1.9473,-1.8794,-1.4024,dropScore,0.3974,-0.4922,0.3896,-0.0641,-0.1872,-0.142,-0.2714,-0.1185,-0.4026,-1.6653,-1.653,-0.2555,-1.0421,-1.6262,-0.2166,2.3526,0.8047,-1.7097,-1.2892,-1.1841,dropScore,-0.0709,-0.5567,-0.0149,-0.5175,-0.4273,-0.255,-0.414,-0.9926,-0.6345,-0.6662,-1.0118,-0.2942,-0.4541,-1.3124,-0.6586,0.8047,2.6474,-1.6094,-1.2019,-0.1755,dropScore,-1.5202,-1.2355,-1.9669,-2.4706,-1.1734,-1.619,-2.2105,-2.0906,-0.6464,-1.0454,-0.576,-1.8695,-0.6095,0.9849,-1.9473,-1.7097,-1.6094,5.8259,1.2574,-1.2247,dropScore,-1.3301,-1.1356,-1.0592,-1.8456,-0.8035,-1.1273,-1.6607,-2.1353,0.7216,-0.8996,-0.6119,-1.4346,-0.5404,1.6772,-1.8794,-1.2892,-1.2019,1.2574,4.0647,-0.9588,dropScore,-0.1153,-1.5974,-1.8046,-2.2981,-0.1612,-1.3761,-1.6411,-2.0053,-1.5571,1.4591,0.4308,-1.5386,0.2329,-0.2784,-1.4024,-1.1841,-0.1755,-1.2247,-0.9588,2.3577,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,dropScore,6.1118),
		nrow=21L)
	lettersAA <- AAStringSet("ARNDCQEGHILKMFPSTWYV*") # concatenated row/column names
	
	results <- matrix(data.frame(),
		nrow=l,
		ncol=l,
		dimnames=list(identifier, identifier))
	empty_upper <- matrix(integer(),
		nrow=0,
		ncol=8,
		dimnames=list(NULL,
			c("index1",
				"index2",
				"strand",
				"width",
				"start1",
				"start2",
				"frame1",
				"frame2")))
	empty_lower <- matrix(integer(),
		nrow=0,
		ncol=10,
		dimnames=list(NULL,
			c("index1",
				"index2",
				"strand",
				"score",
				"start1",
				"start2",
				"end1",
				"end2",
				"first_hit",
				"last_hit")))
	
	store <- rep(list(list(S=list(),
			E=list(`nt`=list(),
				`nt_rc`=list(),
				`aa`=list(list(), list(), list()),
				`aa_rc`=list(list(), list(), list())),
			O=list(`nt`=list(),
				`nt_rc`=list(),
				`aa`=list(list(), list(), list()),
				`aa_rc`=list(list(), list(), list())))),
		l)
	
	if (verbose) {
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
		tot <- l*(l - 1)/2
		its <- 0
	}
	for (g1 in 1:(l - 1)) {
		# remove unnecessary items from store
		store[g1][[1L]][["E"]][["nt_rc"]] <- list()
		store[g1][[1L]][["E"]][["aa_rc"]] <- list()
		store[g1][[1L]][["O"]][["nt_rc"]] <- list()
		store[g1][[1L]][["O"]][["aa_rc"]] <- list()
		
		s1 <- store[g1][[1L]][["S"]]
		if (length(s1) == 0) {
			s1 <- SearchDB(dbConn,
				tblName=tblName,
				identifier=identifier[g1],
				type="DNAStringSet",
				removeGaps="all",
				processors=processors,
				verbose=FALSE)
		} else {
			s1 <- s1[[1L]]
			store[g1][[1L]][["S"]] <- list() # clear s1 from store
		}
		w1 <- width(s1)
		INDEX1 <- seq_along(s1)
		results[g1, g1] <- list(w1)
		
		w <- which(w1 <= 32)
		if (length(w) > 0) {
			s1 <- s1[-w]
			INDEX1 <- INDEX1[-w]
		}
		if (length(s1) == 0) {
			for (g2 in (g1 + 1):l) {
				results[g1, g2][[1]] <- empty_upper
				results[g2, g1][[1]] <- empty_lower
			}
			
			next
		}
		
		WIDTH1 <- cumsum(width(s1))	
		if (length(s1) > 1) {
			seq1 <- DNAStringSet(unlist(s1))
		} else {
			seq1 <- s1
		}
		
		# determine the entropy equivalent alphabet size
		size <- .Call("alphabetSize",
			seq1,
			PACKAGE="DECIPHER")
		
		for (g2 in (g1 + 1):l) {
			s2 <- store[g2][[1L]][["S"]]
			if (length(s2) == 0) {
				s2 <- SearchDB(dbConn,
					tblName=tblName,
					identifier=identifier[g2],
					type="DNAStringSet",
					removeGaps="all",
					processors=processors,
					verbose=FALSE)
				if ((object.size(s2) + object.size(store) + object.size(results)) < storage)
					store[g2][[1L]][["S"]][[1L]] <- s2
			} else {
				s2 <- s2[[1L]]
			}
			w2 <- width(s2)
			INDEX2 <- seq_along(s2)
			
			w <- which(w2 <= 32)
			if (length(w) > 0) {
				s2 <- s2[-w]
				INDEX2 <- INDEX2[-w]
			}
			if (length(s2) == 0) {
				results[g1, g2][[1]] <- empty_upper
				results[g2, g1][[1]] <- empty_lower
				
				next
			}
			
			WIDTH2 <- cumsum(width(s2))
			if (length(s2) > 1) {
				seq2 <- DNAStringSet(unlist(s2))
			} else {
				seq2 <- s2
			}
			
			# calculate the k-mer size likely to occur 1/N times by chance
			M <- max(WIDTH1[length(WIDTH1)], WIDTH2[length(WIDTH2)])*N
			K <- as.integer(log(M, size))
			if (K < 1L) {
				K <- 1L
			} else if (K > 16L) {
				K <- 16L
			}
			
			E1 <- store[g1][[1L]][["E"]][["nt"]][K][[1L]]
			if (is.null(E1)) {
				E1 <- .Call("enumerateSequence",
					seq1,
					K,
					maskRepeats,
					maskLCRs,
					integer(), # mask numerous k-mers
					1L, # left is fast moving side
					processors,
					PACKAGE="DECIPHER")[[1]]
				for (i in which(WIDTH1 > (K - 2) & WIDTH1 < length(E1)))
					E1[(WIDTH1[i] - (K - 2)):WIDTH1[i]] <- NA
				if ((object.size(E1) + object.size(store) + object.size(results)) < storage)
					store[g1][[1L]][["E"]][["nt"]][[K]] <- E1
			}
			
			e2 <- store[g2][[1L]][["E"]][["nt"]][K][[1L]]
			if (is.null(e2)) {
				e2 <- .Call("enumerateSequence",
					seq2,
					K,
					maskRepeats,
					maskLCRs,
					integer(), # mask numerous k-mers
					1L, # left is fast moving side
					processors,
					PACKAGE="DECIPHER")[[1]]
				for (i in which(WIDTH2 > (K - 2) & WIDTH2 < length(e2)))
					e2[(WIDTH2[i] - (K - 2)):WIDTH2[i]] <- NA
				if ((object.size(e2) + object.size(store) + object.size(results)) < storage)
					store[g2][[1L]][["E"]][["nt"]][[K]] <- e2
			}
			
			O1 <- store[g1][[1L]][["O"]][["nt"]][K][[1L]]
			if (is.null(O1)) {
				O1 <- order(E1, method="radix", na.last=FALSE) - 1L
				if ((object.size(O1) + object.size(store) + object.size(results)) < storage)
					store[g1][[1L]][["O"]][["nt"]][[K]] <- O1
			}
			
			o2 <- store[g2][[1L]][["O"]][["nt"]][K][[1L]]
			if (is.null(o2)) {
				o2 <- order(e2, method="radix", na.last=FALSE) - 1L
				if ((object.size(o2) + object.size(store) + object.size(results)) < storage)
					store[g2][[1L]][["O"]][["nt"]][[K]] <- o2
			}
			
			if (identifier[g1] != identifier[g2]) {
				# match E1 to e2
				m <- .Call("intMatchOnce",
					E1,
					e2,
					O1,
					o2,
					PACKAGE="DECIPHER")
			} else {
				# match E1 to itself
				m <- .Call("intMatchSelfOnce",
					E1,
					O1,
					PACKAGE="DECIPHER")
			}
			m <- .Call("fillOverlaps", m, K, PACKAGE="DECIPHER") # in-place change of m
			r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
			w <- which(r@values == 1)
			widths <- r@lengths[w] + K
			ends <- cumsum(r@lengths)[w] + K
			
			ends1 <- ends
			starts1 <- ends1 - widths + 1L
			ends2 <- m[ends - K] + K
			starts2 <- ends2 - widths + 1L
			
			# re-index the sequences by contig
			index1 <- .Call("indexByContig",
				starts1,
				ends1,
				seq_along(starts1),
				INDEX1,
				WIDTH1)
			starts1 <- index1[[2]]
			ends1 <- index1[[3]]
			index1 <- index1[[1]]
			o <- order(starts2, method="radix", na.last=FALSE)
			index2 <- .Call("indexByContig",
				starts2,
				ends2,
				o,
				INDEX2,
				WIDTH2)
			starts2 <- index2[[2]]
			ends2 <- index2[[3]]
			index2 <- index2[[1]]
			
			# score and extend hits
			subScore <- .Call("extendMatches",
				seq1,
				seq2,
				starts1,
				ends1,
				index1,
				starts2,
				ends2,
				index2,
				WIDTH1,
				WIDTH2,
				subMatrixDNA*log(size), # scale subMatrix
				lettersDNA,
				dropScore,
				processors,
				PACKAGE="DECIPHER")
			starts1 <- subScore[[2L]]
			ends1 <- subScore[[3L]]
			starts2 <- subScore[[4L]]
			ends2 <- subScore[[5L]]
			subScore <- subScore[[1L]]
			
			if (useFrames) {
				x.s <- list(starts1)
				x.e <- list(ends1)
				x.i <- list(index1)
				y.s <- list(starts2)
				y.e <- list(ends2)
				y.i <- list(index2)
				x.f <- y.f <- list(rep(0L, length(starts1)))
				weights <- list(subScore)
				
				for (rF1 in 1:3) {
					t1 <- .Call("basicTranslate",
						s1,
						geneticCode[[g1]],
						rep(rF1, length(s1)),
						PACKAGE="DECIPHER")
					width1 <- cumsum(width(t1))
					if (length(t1) > 1) {
						t1 <- AAStringSet(unlist(t1))
					}
					
					# determine the optimal alphabet size
					if (rF1 == 1) {
						sizeAA <- .Call("alphabetSizeReducedAA",
							t1,
							alphabet,
							PACKAGE="DECIPHER")
						N_AA <- as.integer(log(M/3, sizeAA))
						if (N_AA < 1L) {
							N_AA <- 1L
						} else if (N_AA > n) {
							N_AA <- n
						}
					}
					
					e1 <- store[g1][[1L]][["E"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(e1)) {
						e1 <- .Call("enumerateSequenceReducedAA",
							t1,
							N_AA,
							alphabet,
							maskRepeats,
							maskLCRs,
							integer(), # mask numerous k-mers
							1L, # left is fast moving side
							processors,
							PACKAGE="DECIPHER")[[1]]
						for (i in which(width1 > (N_AA - 2) & width1 < length(e1)))
							e1[(width1[i] - (N_AA - 2)):width1[i]] <- NA
						if ((object.size(e1) + object.size(store) + object.size(results)) < storage)
							store[g1][[1L]][["E"]][["aa"]][[rF1]][[N_AA]] <- e1
					}
					
					o1 <- store[g1][[1L]][["O"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(o1)) {
						o1 <- order(e1, method="radix", na.last=FALSE) - 1L
						if ((object.size(o1) + object.size(store) + object.size(results)) < storage)
							store[g1][[1L]][["O"]][["aa"]][[rF1]][[N_AA]] <- o1
					}
					
					for (rF2 in 1:3) {
						t2 <- .Call("basicTranslate",
							s2,
							geneticCode[[g2]],
							rep(rF2, length(s2)),
							PACKAGE="DECIPHER")
						width2 <- cumsum(width(t2))
						if (length(t2) > 1) {
							t2 <- AAStringSet(unlist(t2))
						}
						
						e2 <- store[g2][[1L]][["E"]][["aa"]][[rF2]][N_AA][[1L]]
						if (is.null(e2)) {
							e2 <- .Call("enumerateSequenceReducedAA",
								t2,
								N_AA,
								alphabet,
								maskRepeats,
								maskLCRs,
								integer(), # mask numerous k-mers
								1L, # left is fast moving side
								processors,
								PACKAGE="DECIPHER")[[1]]
							for (i in which(width2 > (N_AA - 2) & width2 < length(e2)))
								e2[(width2[i] - (N_AA - 2)):width2[i]] <- NA
							if ((object.size(e2) + object.size(store) + object.size(results)) < storage)
								store[g2][[1L]][["E"]][["aa"]][[rF2]][[N_AA]] <- e2
						}
						
						o2 <- store[g2][[1L]][["O"]][["aa"]][[rF2]][N_AA][[1L]]
						if (is.null(o2)) {
							o2 <- order(e2, method="radix", na.last=FALSE) - 1L
							if ((object.size(o2) + object.size(store) + object.size(results)) < storage)
								store[g2][[1L]][["O"]][["aa"]][[rF2]][[N_AA]] <- o2
						}
						
						if (identifier[g1] != identifier[g2]) {
							# match e1 to e2
							m <- .Call("intMatchOnce",
								e1,
								e2,
								o1,
								o2,
								PACKAGE="DECIPHER")
						} else {
							# match E1 to itself
							m <- .Call("intMatchSelfOnce",
								e1,
								o1,
								PACKAGE="DECIPHER")
						}
						m <- .Call("fillOverlaps", m, N_AA, PACKAGE="DECIPHER") # in-place change of m
						r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
						w <- which(r@values == 1)
						widths <- r@lengths[w] + N_AA
						ends <- cumsum(r@lengths)[w] + N_AA
						
						ends1 <- ends
						starts1 <- ends1 - widths + 1L
						ends2 <- m[ends - N_AA] + N_AA
						starts2 <- ends2 - widths + 1L
						
						# convert from amino acid to nucleotide positions
						starts1 <- (starts1 - 1L)*3L + rF1
						ends1 <- ends1*3L + rF1 - 1L
						starts2 <- (starts2 - 1L)*3L + rF2
						ends2 <- ends2*3L + rF2 - 1L
						
						# re-index the sequences by contig
						index1 <- .Call("indexByContig",
							starts1,
							ends1,
							seq_along(starts1),
							INDEX1,
							width1*3L)
						starts1 <- index1[[2]]
						ends1 <- index1[[3]]
						index1 <- index1[[1]]
						o <- order(starts2, method="radix", na.last=FALSE)
						index2 <- .Call("indexByContig",
							starts2,
							ends2,
							o,
							INDEX2,
							width2*3L)
						starts2 <- index2[[2]]
						ends2 <- index2[[3]]
						index2 <- index2[[1]]
						
						# convert from nucleotide to amino acid positions
						starts1 <- (starts1 - rF1) %/% 3L + 1L
						ends1 <- (ends1 + 1L - rF1) %/% 3L
						starts2 <- (starts2 - rF2) %/% 3L + 1L
						ends2 <- (ends2 + 1L - rF2) %/% 3L
						
						# score and extend hits
						subScore <- .Call("extendMatches",
							t1,
							t2,
							starts1,
							ends1,
							index1,
							starts2,
							ends2,
							index2,
							width1,
							width2,
							subMatrixAA,
							lettersAA,
							dropScore,
							processors,
							PACKAGE="DECIPHER")
						starts1 <- subScore[[2L]]
						ends1 <- subScore[[3L]]
						starts2 <- subScore[[4L]]
						ends2 <- subScore[[5L]]
						subScore <- subScore[[1L]]
						
						# convert from amino acid to nucleotide positions
						starts1 <- (starts1 - 1L)*3L + rF1
						ends1 <- ends1*3L + rF1 - 1L
						starts2 <- (starts2 - 1L)*3L + rF2
						ends2 <- ends2*3L + rF2 - 1L
						
						x.s <- c(x.s, list(starts1))
						x.e <- c(x.e, list(ends1))
						x.i <- c(x.i, list(index1))
						y.s <- c(y.s, list(starts2))
						y.e <- c(y.e, list(ends2))
						y.i <- c(y.i, list(index2))
						x.f <- c(x.f, list(rep(rF1, length(starts1))))
						y.f <- c(y.f, list(rep(rF2, length(starts1))))
						weights <- c(weights, list(subScore))
					}
				}
				
				x.s <- unlist(x.s)
				x.e <- unlist(x.e)
				x.i <- unlist(x.i)
				y.s <- unlist(y.s)
				y.e <- unlist(y.e)
				y.i <- unlist(y.i)
				x.f <- unlist(x.f)
				y.f <- unlist(y.f)
				weights <- unlist(weights)
				maxW <- max(WIDTH1[length(WIDTH1)],
					WIDTH2[length(WIDTH2)])*10 # equivalent to ~10 searches
			} else {
				x.s <- starts1
				x.e <- ends1
				x.i <- index1
				y.s <- starts2
				y.e <- ends2
				y.i <- index2
				x.f <- y.f <- rep(0L, length(starts1))
				weights <- subScore
				maxW <- max(WIDTH1[length(WIDTH1)],
					WIDTH2[length(WIDTH2)])
			}
			
			# order by increasing sequence index in g1
			o <- order(x.i, x.s)
			x.s <- x.s[o]
			x.e <- x.e[o]
			x.i <- x.i[o]
			y.s <- y.s[o]
			y.e <- y.e[o]
			y.i <- y.i[o]
			x.f <- x.f[o]
			y.f <- y.f[o]
			weights <- weights[o]
			
			chains <- .Call("chainSegments",
				x.s,
				x.e,
				x.i,
				x.f,
				y.s,
				y.e,
				y.i,
				y.f,
				weights,
				sepCost,
				sepPower,
				gapCost,
				gapPower,
				shiftCost,
				codingCost,
				maxSep,
				maxGap,
				order(x.i, x.e) - 1L,
				minScore,
				maxW,
				allowOverlap,
				PACKAGE="DECIPHER")
			
			scores <- chains[[2]]
			chains <- chains[[1]]
			used <- unlist(chains,
				use.names=FALSE)
			count <- 0L
			for (i in seq_along(chains)) {
				chains[[i]] <- (count + 1L):(count + length(chains[[i]]))
				count <- count + length(chains[[i]])
			}
			results[g2, g1][[1]] <- list(chain=chains,
				score=scores)
			
			result1 <- matrix(c(x.i[used],
					y.i[used],
					rep(0L, length(used)),
					x.e[used] - x.s[used] + 1L,
					x.s[used],
					y.s[used],
					x.f[used],
					y.f[used]),
				nrow=length(used),
				ncol=8,
				dimnames=list(NULL,
					c("index1",
						"index2",
						"strand",
						"width",
						"start1",
						"start2",
						"frame1",
						"frame2")))
			
			if (verbose) {
				its <- its + 0.5
				setTxtProgressBar(pBar, its/tot)
			}
			
			s3 <- reverseComplement(s2)
			if (length(s3) > 1) {
				seq2 <- DNAStringSet(unlist(s3))
			} else {
				seq2 <- s3
			}
			
			e2 <- store[g2][[1L]][["E"]][["nt_rc"]][K][[1L]]
			if (is.null(e2)) {
				e2 <- .Call("enumerateSequence",
					seq2,
					K,
					maskRepeats,
					maskLCRs,
					integer(), # mask numerous k-mers
					1L, # left is fast moving side
					processors,
					PACKAGE="DECIPHER")[[1]]
				for (i in which(WIDTH2 > (K - 2) & WIDTH2 < length(e2)))
					e2[(WIDTH2[i] - (K - 2)):WIDTH2[i]] <- NA
				if ((object.size(e2) + object.size(store) + object.size(results)) < storage)
					store[g2][[1L]][["E"]][["nt_rc"]][[K]] <- e2
			}
			
			o2 <- store[g2][[1L]][["O"]][["nt_rc"]][K][[1L]]
			if (is.null(o2)) {
				o2 <- order(e2, method="radix", na.last=FALSE) - 1L
				if ((object.size(o2) + object.size(store) + object.size(results)) < storage)
					store[g2][[1L]][["O"]][["nt_rc"]][[K]] <- o2
			}
			
			# match E1 to e2
			m <- .Call("intMatchOnce",
				E1,
				e2,
				O1,
				o2,
				PACKAGE="DECIPHER")
			m <- .Call("fillOverlaps", m, K, PACKAGE="DECIPHER") # in-place change of m
			r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
			w <- which(r@values == 1)
			widths <- r@lengths[w] + K
			ends <- cumsum(r@lengths)[w] + K
			
			ends1 <- ends
			starts1 <- ends1 - widths + 1L
			ends2 <- m[ends - K] + K
			starts2 <- ends2 - widths + 1L
			
			# re-index the sequences by contig
			index1 <- .Call("indexByContig",
				starts1,
				ends1,
				seq_along(starts1),
				INDEX1,
				WIDTH1)
			starts1 <- index1[[2]]
			ends1 <- index1[[3]]
			index1 <- index1[[1]]
			o <- order(starts2, method="radix", na.last=FALSE)
			index2 <- .Call("indexByContig",
				starts2,
				ends2,
				o,
				INDEX2,
				WIDTH2)
			starts2 <- index2[[2]]
			ends2 <- index2[[3]]
			index2 <- index2[[1]]
			
			# score and extend hits
			subScore <- .Call("extendMatches",
				seq1,
				seq2,
				starts1,
				ends1,
				index1,
				starts2,
				ends2,
				index2,
				WIDTH1,
				WIDTH2,
				subMatrixDNA*log(size), # scale subMatrix
				lettersDNA,
				dropScore,
				processors,
				PACKAGE="DECIPHER")
			starts1 <- subScore[[2L]]
			ends1 <- subScore[[3L]]
			starts2 <- subScore[[4L]]
			ends2 <- subScore[[5L]]
			subScore <- subScore[[1L]]
			
			# correct for reverse complement positioning
			w <- w2[index2]
			temp <- w - starts2 + 1L
			starts2 <- w - ends2 + 1L
			ends2 <- temp
			
			if (useFrames) {
				x.s <- list(starts1)
				x.e <- list(ends1)
				x.i <- list(index1)
				y.s <- list(starts2)
				y.e <- list(ends2)
				y.i <- list(index2)
				x.f <- y.f <- list(rep(0L, length(starts1)))
				weights <- list(subScore)
				
				for (rF1 in 1:3) {
					t1 <- .Call("basicTranslate",
						s1,
						geneticCode[[g1]],
						rep(rF1, length(s1)),
						PACKAGE="DECIPHER")
					width1 <- cumsum(width(t1))
					if (length(t1) > 1) {
						t1 <- AAStringSet(unlist(t1))
					}
					
					e1 <- store[g1][[1L]][["E"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(e1)) {
						e1 <- .Call("enumerateSequenceReducedAA",
							t1,
							N_AA,
							alphabet,
							maskRepeats,
							maskLCRs,
							integer(), # mask numerous k-mers
							1L, # left is fast moving side
							processors,
							PACKAGE="DECIPHER")[[1]]
						for (i in which(width1 > (N_AA - 2) & width1 < length(e1)))
							e1[(width1[i] - (N_AA - 2)):width1[i]] <- NA
						if ((object.size(e1) + object.size(store) + object.size(results)) < storage)
							store[g1][[1L]][["E"]][["aa"]][[rF1]][[N_AA]] <- e1
					}
					
					o1 <- store[g1][[1L]][["O"]][["aa"]][[rF1]][N_AA][[1L]]
					if (is.null(o1)) {
						o1 <- order(e1, method="radix", na.last=FALSE) - 1L
						if ((object.size(o1) + object.size(store) + object.size(results)) < storage)
							store[g1][[1L]][["O"]][["aa"]][[rF1]][[N_AA]] <- o1
					}
					
					for (rF2 in 1:3) {
						t2 <- .Call("basicTranslate",
							s3,
							geneticCode[[g2]],
							rep(rF2, length(s3)),
							PACKAGE="DECIPHER")
						width2 <- cumsum(width(t2))
						if (length(t2) > 1) {
							t2 <- AAStringSet(unlist(t2))
						}
						
						e2 <- store[g2][[1L]][["E"]][["aa_rc"]][[rF2]][N_AA][[1L]]
						if (is.null(e2)) {
							e2 <- .Call("enumerateSequenceReducedAA",
								t2,
								N_AA,
								alphabet,
								maskRepeats,
								maskLCRs,
								integer(), # mask numerous k-mers
								1L, # left is fast moving side
								processors,
								PACKAGE="DECIPHER")[[1]]
							for (i in which(width2 > (N_AA - 2) & width2 < length(e2)))
								e2[(width2[i] - (N_AA - 2)):width2[i]] <- NA
							if ((object.size(e2) + object.size(store) + object.size(results)) < storage)
								store[g2][[1L]][["E"]][["aa_rc"]][[rF2]][[N_AA]] <- e2
						}
						
						o2 <- store[g2][[1L]][["O"]][["aa_rc"]][[rF2]][N_AA][[1L]]
						if (is.null(o2)) {
							o2 <- order(e2, method="radix", na.last=FALSE) - 1L
							if ((object.size(o2) + object.size(store) + object.size(results)) < storage)
								store[g2][[1L]][["O"]][["aa_rc"]][[rF2]][[N_AA]] <- o2
						}
						
						# match e1 to e2
						m <- .Call("intMatchOnce",
							e1,
							e2,
							o1,
							o2,
							PACKAGE="DECIPHER")
						m <- .Call("fillOverlaps", m, N_AA, PACKAGE="DECIPHER") # in-place change of m
						r <- Rle(.Call("intDiff", m, PACKAGE="DECIPHER"))
						w <- which(r@values == 1)
						widths <- r@lengths[w] + N_AA
						ends <- cumsum(r@lengths)[w] + N_AA
						
						ends1 <- ends
						starts1 <- ends1 - widths + 1L
						ends2 <- m[ends - N_AA] + N_AA
						starts2 <- ends2 - widths + 1L
						
						# convert from amino acid to nucleotide positions
						starts1 <- (starts1 - 1L)*3L + rF1
						ends1 <- ends1*3L + rF1 - 1L
						starts2 <- (starts2 - 1L)*3L + rF2
						ends2 <- ends2*3L + rF2 - 1L
						
						# re-index the sequences by contig
						index1 <- .Call("indexByContig",
							starts1,
							ends1,
							seq_along(starts1),
							INDEX1,
							width1*3L)
						starts1 <- index1[[2]]
						ends1 <- index1[[3]]
						index1 <- index1[[1]]
						o <- order(starts2, method="radix", na.last=FALSE)
						index2 <- .Call("indexByContig",
							starts2,
							ends2,
							o,
							INDEX2,
							width2*3L)
						starts2 <- index2[[2]]
						ends2 <- index2[[3]]
						index2 <- index2[[1]]
						
						# convert from nucleotide to amino acid positions
						starts1 <- (starts1 - rF1) %/% 3L + 1L
						ends1 <- (ends1 + 1L - rF1) %/% 3L
						starts2 <- (starts2 - rF2) %/% 3L + 1L
						ends2 <- (ends2 + 1L - rF2) %/% 3L
						
						# score and extend hits
						subScore <- .Call("extendMatches",
							t1,
							t2,
							starts1,
							ends1,
							index1,
							starts2,
							ends2,
							index2,
							width1,
							width2,
							subMatrixAA,
							lettersAA,
							dropScore,
							processors,
							PACKAGE="DECIPHER")
						starts1 <- subScore[[2L]]
						ends1 <- subScore[[3L]]
						starts2 <- subScore[[4L]]
						ends2 <- subScore[[5L]]
						subScore <- subScore[[1L]]
						
						# convert from amino acid to nucleotide positions
						starts1 <- (starts1 - 1L)*3L + rF1
						ends1 <- ends1*3L + rF1 - 1L
						starts2 <- (starts2 - 1L)*3L + rF2
						ends2 <- ends2*3L + rF2 - 1L
						
						# correct for reverse complement positioning
						w <- w2[index2]
						temp <- w - starts2 + 1L
						starts2 <- w - ends2 + 1L
						ends2 <- temp
						
						x.s <- c(x.s, list(starts1))
						x.e <- c(x.e, list(ends1))
						x.i <- c(x.i, list(index1))
						y.s <- c(y.s, list(starts2))
						y.e <- c(y.e, list(ends2))
						y.i <- c(y.i, list(index2))
						x.f <- c(x.f, list(rep(rF1, length(starts1))))
						y.f <- c(y.f, list(rep(rF2, length(starts1))))
						weights <- c(weights, list(subScore))
					}
				}
				
				x.s <- unlist(x.s)
				x.e <- unlist(x.e)
				x.i <- unlist(x.i)
				y.s <- unlist(y.s)
				y.e <- unlist(y.e)
				y.i <- unlist(y.i)
				x.f <- unlist(x.f)
				y.f <- unlist(y.f)
				weights <- unlist(weights)
			} else {
				x.s <- starts1
				x.e <- ends1
				x.i <- index1
				y.s <- starts2
				y.e <- ends2
				y.i <- index2
				x.f <- y.f <- rep(0L, length(starts1))
				weights <- subScore
			}
			
			# order by increasing sequence index in g1
			o <- order(x.i, x.s)
			x.s <- x.s[o]
			x.e <- x.e[o]
			x.i <- x.i[o]
			y.s <- y.s[o]
			y.e <- y.e[o]
			y.i <- y.i[o]
			x.f <- x.f[o]
			y.f <- y.f[o]
			weights <- weights[o]
			
			if (length(y.s) > 0) {
				d <- max(y.e) - (y.e + y.s)
			} else {
				d <- integer()
			}
			chains <- .Call("chainSegments",
				x.s,
				x.e,
				x.i,
				x.f,
				y.s + d,
				y.e + d,
				y.i,
				y.f,
				weights,
				sepCost,
				sepPower,
				gapCost,
				gapPower,
				shiftCost,
				codingCost,
				maxSep,
				maxGap,
				order(x.i, x.e) - 1L,
				minScore,
				maxW,
				allowOverlap,
				PACKAGE="DECIPHER")
			
			scores <- chains[[2]]
			chains <- chains[[1]]
			used <- unlist(chains,
				use.names=FALSE)
			for (i in seq_along(chains)) {
				chains[[i]] <- (count + 1L):(count + length(chains[[i]]))
				count <- count + length(chains[[i]])
			}
			results[g2, g1][[1]] <- list(chain=c(results[g2, g1][[1]]$chain,
					chains),
				score=c(results[g2, g1][[1]]$score,
					scores))
			
			result2 <- matrix(c(x.i[used],
					y.i[used],
					rep(1L, length(used)),
					x.e[used] - x.s[used] + 1L,
					x.s[used],
					y.e[used],
					x.f[used],
					y.f[used]),
				nrow=length(used),
				ncol=8,
				dimnames=list(NULL,
					c("index1",
						"index2",
						"strand",
						"width",
						"start1",
						"start2",
						"frame1",
						"frame2")))
			
			result <- rbind(result1, result2)
			results[g1, g2][[1]] <- result
			
			# determine ranges
			chains <- results[g2, g1][[1]]$chain
			starts <- unlist(lapply(chains, `[`, 1L))
			ends <- unlist(lapply(chains,
				function(x) x[length(x)]))
			starts1 <- result[starts, "start1"]
			ends1 <- result[ends, "start1"] + result[ends, "width"] - 1L
			strand <- result[starts, "strand"]
			starts2 <- ifelse(strand == 0,
				result[starts, "start2"],
				result[ends, "start2"] + 1L - result[ends, "width"])
			ends2 <- ifelse(strand == 0,
				result[ends, "start2"] + result[ends, "width"] - 1L,
				result[starts, "start2"])
			index1 <- result[starts, "index1"]
			index2 <- result[starts, "index2"]
			score <- results[g2, g1][[1]]$score
			results[g2, g1][[1]] <- matrix(c(index1,
					index2,
					strand,
					as.integer(score),
					starts1,
					starts2,
					ends1,
					ends2,
					starts,
					ends),
					ncol=10,
					dimnames=list(NULL,
						c("index1",
							"index2",
							"strand",
							"score",
							"start1",
							"start2",
							"end1",
							"end2",
							"first_hit",
							"last_hit")))
			
			# remove overlap
			if (length(strand) > 1 && !allowOverlap) {
				remove <- logical(length(strand))
				hits <- list()
				
				ss1 <- result[, "start1"]
				ee1 <- result[, "start1"] + result[, "width"] - 1L
				ss2 <- ifelse(result[, "strand"] == 0,
					result[, "start2"],
					result[, "start2"] + 1L - result[, "width"])
				ee2 <- ifelse(result[, "strand"] == 0,
					result[, "start2"] + result[, "width"] - 1L,
					result[, "start2"])
				
				o <- order(index1,
					index2,
					starts1,
					-ends1)
				last <- 1L
				for (i in 2:length(strand)) { # each block
					if (starts1[o[i]] <= ends1[o[last]] &&
						index1[o[i]] == index1[o[last]] &&
						index2[o[i]] == index2[o[last]]) {
						if (ends1[o[i]] <= ends1[o[last]]) {
							# completely overlapping
							chain_l <- chains[[o[last]]]
							over_l <- ee1[chain_l] >= starts1[o[i]] &
								ss1[chain_l] <= ends1[o[i]]
							if (any(over_l)) { # overlapping hits
								if (score[o[i]] >= score[o[last]]) {
									remove[o[last]] <- TRUE
									last <- i
								} else {
									remove[o[i]] <- TRUE
								}
							}
						} else {
							# remove partial overlap from last
							chain_l <- chains[[o[last]]]
							over_l <- which(ss1[chain_l] >= (starts1[o[i]] - 2)) # min anchor width of 2
							if (length(over_l) == length(chain_l)) {
								remove[o[last]] <- TRUE
								last <- i
								next
							} else if (length(over_l) > 0) {
								# remove completely overlapping hits
								hits[[length(hits) + 1L]] <- chain_l[over_l]
								chain_l <- chain_l[-over_l]
								chains[[o[last]]] <- chain_l
							}
							
							cl <- chain_l[length(chain_l)]
							overlap <- ee1[cl] - starts1[o[i]] + 1L
							if (overlap > 0) {
								# remove overlap from last hit
								results[g1, g2][[1]][cl, "width"] <- results[g1, g2][[1]][cl, "width"] - overlap
								ee1[cl] <- ee1[cl] - overlap
								if (strand[o[last]] == 0) {
									ee2[cl] <- ee2[cl] - overlap
								} else {
									ss2[cl] <- ss2[cl] + overlap
								}
							}
							ends1[o[last]] <- ee1[cl]
							results[g2, g1][[1]][o[last], "end1"] <- ee1[cl]
							if (strand[o[last]] == 0) {
								ends2[o[last]] <- ee2[cl]
								results[g2, g1][[1]][o[last], "end2"] <- ee2[cl]
							} else {
								starts2[o[last]] <- ss2[cl]
								results[g2, g1][[1]][o[last], "start2"] <- ss2[cl]
							}
							last <- i
						}
					} else { # non-overlapping
						last <- i
					}
				}
				
				o <- order(index2,
					index1,
					starts2,
					-ends2)
				w <- which(!remove[o])
				last <- w[1]
				for (i in w[-1]) {
					if (starts2[o[i]] <= ends2[o[last]] &&
						index1[o[i]] == index1[o[last]] &&
						index2[o[i]] == index2[o[last]]) {
						if (ends2[o[i]] <= ends2[o[last]]) {
							# completely overlapping
							chain_l <- chains[[o[last]]]
							over_l <- ee2[chain_l] >= starts2[o[i]] &
								ss2[chain_l] <= ends2[o[i]]
							if (any(over_l)) { # overlapping hits
								if (score[o[i]] >= score[o[last]]) {
									remove[o[last]] <- TRUE
									last <- i
								} else {
									remove[o[i]] <- TRUE
								}
							}
						} else {
							# remove partial overlap from last
							chain_l <- chains[[o[last]]]
							over_l <- which(ss2[chain_l] >= (starts2[o[i]] - 2)) # min anchor width of 2
							if (length(over_l) == length(chain_l)) {
								remove[o[last]] <- TRUE
								last <- i
								next
							} else if (length(over_l) > 0) {
								# remove completely overlapping hits
								hits[[length(hits) + 1L]] <- chain_l[over_l]
								chain_l <- chain_l[-over_l]
								chains[[o[last]]] <- chain_l
							}
							if (strand[o[last]] == 0) {
								cl <- chain_l[length(chain_l)]
							} else {
								cl <- chain_l[1]
							}
							overlap <- ee2[cl] - starts2[o[i]] + 1L
							if (overlap > 0) {
								# remove overlap from last hit
								results[g1, g2][[1]][cl, "width"] <- results[g1, g2][[1]][cl, "width"] - overlap
								if (strand[o[last]] == 0) {
									ee1[cl] <- ee1[cl] - overlap
									ee2[cl] <- ee2[cl] - overlap
								} else {
									ss1[cl] <- ss1[cl] + overlap
									ee2[cl] <- ee2[cl] - overlap
									results[g1, g2][[1]][cl, "start1"] <- results[g1, g2][[1]][cl, "start1"] + overlap
									results[g1, g2][[1]][cl, "start2"] <- results[g1, g2][[1]][cl, "start2"] - overlap
								}
							}
							if (strand[o[last]] == 0) {
								results[g2, g1][[1]][o[last], "end1"] <- ee1[cl]
								results[g2, g1][[1]][o[last], "end2"] <- ee2[cl]
							} else {
								results[g2, g1][[1]][o[last], "start1"] <- ss1[cl]
								results[g2, g1][[1]][o[last], "end2"] <- ee2[cl]
							}
							last <- i
						}
					} else { # non-overlapping
						last <- i
					}
				}
				
				widths <- results[g2, g1][[1]][, "end1"] - results[g2, g1][[1]][, "start1"] + 1L
				remove <- which(remove |
					widths*log(size) < minScore) # max possible score is too low
				if (length(remove) > 0 || length(hits) > 0) {
					if (length(remove) > 0)
						results[g2, g1][[1]] <- results[g2, g1][[1]][-remove,, drop=FALSE]
					hits <- c(unlist(hits), unlist(chains[remove]))
					results[g1, g2][[1]] <- results[g1, g2][[1]][-hits,, drop=FALSE]
					if (length(remove) > 0)
						chains <- chains[-remove]
					count <- 0L
					for (i in seq_along(chains)) {
						results[g2, g1][[1]][i, "first_hit"] <- count + 1L
						count <- count + length(chains[[i]])
						results[g2, g1][[1]][i, "last_hit"] <- count
					}
				}
				
				index1 <- results[g2, g1][[1]][, "index1"]
				index2 <- results[g2, g1][[1]][, "index2"]
				starts1 <- results[g2, g1][[1]][, "start1"]
				starts2 <- results[g2, g1][[1]][, "start2"]
			}
			
			if (verbose) {
				its <- its + 0.5
				setTxtProgressBar(pBar, its/tot)
			}
		}
		
		store[[g1]] <- list()
	}
	
	if (!exists("w2",
		envir=environment(),
		inherits=FALSE)) {
		s2 <- SearchDB(dbConn,
			tblName=tblName,
			identifier=identifier[length(identifier)],
			type="DNAStringSet",
			removeGaps="all",
			processors=processors,
			verbose=FALSE)
		w2 <- width(s2)
	}
	results[l, l] <- list(widths=w2)
	
	for (i in seq_len(l)) {
		names(results[i, i][[1]]) <- dbGetQuery(dbConn,
			paste('select ',
				dbQuoteIdentifier(dbConn, 'description'),
				' from ',
				dbQuoteIdentifier(dbConn, tblName),
				' where ',
				dbQuoteIdentifier(dbConn, 'identifier'),
				' = "',
				identifier[i],
				'"',
				sep=''))$description
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	class(results) <- "Synteny"
	
	return(results)
}

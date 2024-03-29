AlignSynteny <- function(synteny,
	dbFile,
	tblName="Seqs",
	identifier="",
	processors=1,
	verbose=TRUE,
	...) {
	
	# error checking
	if (!is(synteny, "Synteny"))
		stop("synteny must be an object of class 'Synteny'.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (is.character(identifier)) {
		if (identifier[1] == "") {
			identifier <- rownames(synteny)
			index <- seq_along(identifier)
		} else {
			index <- match(identifier,
				rownames(synteny))
			if (any(is.na(index)))
				stop("identifier(s) not present in synteny.")
		}
	} else if (is.numeric(identifier)) {
		if (any(is.na(identifier)))
			stop("identifier cannot contain NA values.")
		if (any(identifier < 1))
			stop("identifier not positive.")
		if (any(identifier > nrow(synteny)))
			stop("identifer not present in synteny.")
		index <- identifier
		identifier <- rownames(synteny)[identifier]
	} else {
		stop("identifier must be character or numeric.")
	}
	l <- length(identifier)
	if (l < 2)
		stop("At least two identifiers are required.")
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
	
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	}
	
	o <- order(index)
	count <- 0L
	num <- (l*l - l)/2L
	results <- vector(mode="list",
		length=num)
	ns <- character(num)
	for (i in 1:(l - 1)) {
		index1 <- index[o[i]]
		seq1 <- SearchDB(dbConn,
			tblName=tblName,
			identifier=identifier[o[i]],
			type="DNAStringSet",
			removeGaps="all",
			processors=processors,
			verbose=FALSE)
		if (length(seq1) == 0 ||
			!all(width(seq1) == synteny[index1, index1][[1]]))
			stop("dbFile is mismatched with synteny object.")
		
		for (j in (i + 1):l) {
			index2 <- index[o[j]]
			
			seq2 <- SearchDB(dbConn,
				tblName=tblName,
				identifier=identifier[o[j]],
				type="DNAStringSet",
				removeGaps="all",
				processors=processors,
				verbose=FALSE)
			if (length(seq2) == 0 ||
				!all(width(seq2) == synteny[index2, index2][[1]]))
				stop("dbFile is mismatched with synteny object.")
			
			s <- synteny[index2, index1][[1]]
			starts <- s[, "first_hit"]
			ends <- s[, "last_hit"]
			
			s1 <- s[, "start1"]
			e1 <- s[, "end1"]
			i1 <- s[, "index1"]
			i1 <- lapply(seq_along(seq1),
				function(x) {
					which(i1 == x)
				})
			r1 <- lapply(i1,
				function(x) {
					IRanges(s1[x],
						e1[x])
				})
			seg1 <- extractAt(seq1,
				as(r1, "IRangesList"))
			seg1 <- unlist(seg1)
			i1 <- unlist(i1)
			seg1 <- seg1[order(i1)]
			names(seg1) <- rep(identifier[o[i]],
				length(seg1))
			
			s2 <- s[, "start2"]
			e2 <- s[, "end2"]
			i2 <- s[, "index2"]
			i2 <- lapply(seq_along(seq2),
				function(x) {
					which(i2 == x)
				})
			r2 <- lapply(i2,
				function(x) {
					IRanges(s2[x],
						e2[x])
				})
			seg2 <- extractAt(seq2,
				as(r2, "IRangesList"))
			seg2 <- unlist(seg2)
			i2 <- unlist(i2)
			seg2 <- seg2[order(i2)]
			w <- which(s[, "strand"] == 1)
			if (length(w) > 0)
				seg2[w] <- reverseComplement(seg2[w])
			names(seg2) <- rep(identifier[o[j]],
				length(seg2))
			
			result <- vector(mode="list",
				length=length(starts))
			if (verbose) {
				tot <- mapply(function(x, y) y - x + 1L,
					s[, "first_hit"],
					s[, "last_hit"])
				tot <- cumsum(tot)
			}
			S <- synteny[index1, index2][[1]]
			for (k in seq_len(dim(s)[1])) {
				chain <- s[k, "first_hit"]:s[k, "last_hit"]
				
				s1 <- S[chain, "start1"]
				s1 <- s1 - s[k, "start1"] + 1L
				e1 <- s1 + S[chain, "width"] - 1L
				
				if (s[k, "strand"] == 0) {
					s2 <- S[chain, "start2"]
					s2 <- s2 - s[k, "start2"] + 1L
					e2 <- s2 + S[chain, "width"] - 1L
				} else {
					e2 <- S[chain, "start2"] - S[chain, "width"] + 1L
					s2 <- S[chain, "start2"]
					e2 <- s[k, "end2"] - e2 + 1L
					s2 <- s[k, "end2"] - s2 + 1L
				}
				
				anchors <- matrix(c(s1, e1, s2, e2),
					nrow=4,
					byrow=TRUE)
				result[[k]] <- AlignProfiles(seg1[k],
					seg2[k],
					anchor=anchors,
					processors=processors,
					...)
				
				if (verbose)
					setTxtProgressBar(pBar,
						(count + tot[k]/tot[length(tot)])/num)
			}
			
			result <- relist(do.call(base::c,
					result),
				result)
			
			names(result) <- seq_len(dim(s)[1])
			
			count <- count + 1L
			results[[count]] <- result
			ns[count] <- paste(identifier[o[i]],
				identifier[o[j]],
				sep="/")
		}
	}
	names(results) <- ns
	
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
	
	return(results)
}

AlignTranslation <- function(myXStringSet,
	sense="+",
	direction="5' to 3'",
	readingFrame=NA,
	type=class(myXStringSet),
	geneticCode=GENETIC_CODE,
	...) {
	
	# error checking
	if (!is(myXStringSet, "DNAStringSet") && !is(myXStringSet, "RNAStringSet"))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	if (min(width(myXStringSet)) < 2)
		stop("All sequences in myXStringSet must be at least two nucleotides long.")
	if (length(myXStringSet) < 2)
		stop("At least two sequences are required.")
	if (sense != "+" && sense != "-")
		stop('sense must be either "+" or "-".')
	if (direction != "5' to 3'" && direction != "3' to 5'")
		stop('direction must be either "5\' to 3\'" or "3\' to 5\'".')
	if (!is.na(readingFrame) && !is.numeric(readingFrame))
		stop('readingFrame must be a numeric.')
	if (length(readingFrame) != 1 && length(readingFrame) != length(myXStringSet))
		stop('readingFrame must be a single numeric or the length of myXStringSet.')
	if (any(!is.na(readingFrame) & (readingFrame > 3 | readingFrame < 1)))
		stop('readingFrame must be either NA, 1, 2, or 3.')
	if (any(!is.na(readingFrame) & floor(readingFrame) != readingFrame))
		stop('readingFrame must be either NA, 1, 2, or 3.')
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') must be removed before alignment.")
	a <- vcountPattern("+", myXStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') must be removed before alignment.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') must be removed before alignment.")
	if (length(type) == 0)
		stop("No type specified.")
	if (length(type) > 1)
		stop("Only one type may be specified.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet", "both")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (type == 1) {
		if (!is(myXStringSet, "DNAStringSet"))
			stop("type cannot be 'RNAStringSet' when myXStringSet is a DNAStringSet.")
	} else if (type == 2) {
		if (!is(myXStringSet, "RNAStringSet"))
			stop("type cannot be 'DNAStringSet' when myXStringSet is a RNAStringSet.")
	}
	if (is.list(geneticCode)) {
		if (length(geneticCode) != length(myXStringSet))
			stop("The list geneticCode must have one item per sequence.")
		mapping <- selfmatch(geneticCode)
		mapping <- tapply(seq_along(mapping), mapping, c)
		group <- as.integer(names(mapping))
	} else { # all identifiers use the same geneticCode
		geneticCode <- list(geneticCode)
		mapping <- list(seq_along(myXStringSet))
		group <- 1L
	}
	
	if (sense == "-")
		myXStringSet <- reverseComplement(myXStringSet)
	if (direction == "3' to 5'")
		myXStringSet <- reverse(myXStringSet)
	if (length(readingFrame) == 1)
		readingFrame <- rep(readingFrame, length(myXStringSet))
	
	index <- c(0:2, 0) # circle back around to the first reading frame
	AA <- AAStringSet(rep("", length(myXStringSet)))
	for (i in seq_along(index)) {
		w <- which(((readingFrame - 1) == index[i] & # (specified reading frame AND
			i != length(index)) | # not the last possible reading frame) OR
			is.na(readingFrame)) # reading frame is unspecified
		if (length(w) == 0)
			next
		
		start <- index[i] + 1
		for (j in seq_along(group)) {
			W <- w[w %in% mapping[[j]]]
			if (length(W) == 0)
				next
			
			end <- width(myXStringSet)[W]
			offset <- end - start + 1
			end <- end - offset %% 3
			end <- ifelse(end < start - 1,
				start - 1,
				end)
			AA[W] <- translate(subseq(myXStringSet[W],
					start,
					end),
				genetic.code=geneticCode[[group[j]]],
				if.fuzzy.codon="solve",
				no.init.codon=TRUE)
		}
		end <- width(myXStringSet)[w]
		offset <- end - start + 1
		end <- end - offset %% 3
		end <- ifelse(end < start - 1,
			start - 1,
			end)
		
		# mask missing positions
		AA[w] <- xscat(ifelse(start > 1, "+", ""),
			AA[w],
			ifelse(end < width(myXStringSet)[w], "+", ""))
		
		a <- vcountPattern("*", AA[w])
		lastResidue <- substring(AA[w],
			width(AA)[w],
			width(AA)[w])
		readingFrame[w] <- ifelse(a == 0 | # no stops OR
			(a == 1 & lastResidue == "*") | # one stop at end OR
			i == length(index), # already checked all reading frames
			index[i] + 1,
			readingFrame[w])
	}
	
	names(AA) <- names(myXStringSet)
	AA <- AlignSeqs(myXStringSet=AA, ...)
	if (type > 2) {
		if (type == 4) {
			results <- list(NULL, AA)
		} else {
			return(AA)
		}
	}
	
	# correct for shifts introduced by "+"
	readingFrame <- ifelse(readingFrame == 1,
		1,
		ifelse(readingFrame == 2,
			-1,
			0))
	
	gaps <- vmatchPattern("-", AA)
	starts <- list()
	maxRF <- max(readingFrame)
	Ls <- numeric(length(myXStringSet))
	for (i in seq_along(myXStringSet)) {
		start <- start(gaps[[i]])
		if (length(start) > 0) {
			w <- which(start + (length(start) - 1):0 == width(AA)[i])
			if (length(w) > 0)
				length(start) <- length(start) - length(w)
			start <- start*3 - 3 + readingFrame[i]
			start <- sort(c(start, start + 1, start + 2))
			start <- start - 0:(length(start) - 1)
			w <- which(start == (readingFrame[i]))
			if (length(w) > 0)
				start[w] <- 1
		}
		start <- c(rep(1, maxRF - readingFrame[i]), start)
		starts[[i]] <- start
		Ls[i] <- length(start) + width(myXStringSet)[i]
	}
	
	# add trailing gaps
	maxWidth <- max(Ls)
	for (i in seq_along(myXStringSet)) {
		starts[[i]] <- c(starts[[i]],
			rep(width(myXStringSet)[i] + 1, maxWidth - Ls[i]))
	}
	
	# remove leading common gaps
	start <- min(unlist(lapply(starts, function(x) return(length(which(x == 1))))))
	if (start > 0)
		starts <- lapply(starts, function(x) return(x[-(1:start)]))
	
	myXStringSet <- replaceAt(myXStringSet,
		starts,
		"-")
	if (sense == "-")
		myXStringSet <- reverseComplement(myXStringSet)
	if (direction == "3' to 5'")
		myXStringSet <- reverse(myXStringSet)
	
	if (type == 4) {
		results[[1]] <- myXStringSet
		return(results)
	} else {
		return(myXStringSet)
	}
}

ExtractGenes <- function(x,
	myDNAStringSet,
	type="DNAStringSet",
	...) {
	
	if (!is(x, "Genes"))
		stop("x must be an object of class 'Genes'.")
	if (!is(myDNAStringSet, "DNAStringSet"))
		stop("myDNAStringSet must be a DNAStringSet.")
	TYPES <- c("DNAStringSet", "RNAStringSet", "AAStringSet")
	if (length(type) == 0)
		stop("No type specified.")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (type == 3 && any(1/x[, "Gene"] < 0))
		stop("Only open reading frames in x can be translated when type is 'AAStringSet'.")
	if (!all(attr(x, "widths") == width(myDNAStringSet)))
		stop("Mismatch between x and myDNAStringSet.")
	
	y <- rep(list(IRanges()),
		length(myDNAStringSet))
	t <- tapply(seq_len(nrow(x)),
		x[, "Index"],
		c)
	for (i in seq_along(t)) {
		y[[as.numeric(names(t)[i])]] <- IRanges(x[t[[i]], "Begin"],
			x[t[[i]], "End"])
	}
	y <- IRangesList(y)
	z <- unlist(extractAt(myDNAStringSet, y))
	w <-  which(x[, "Strand"] == 1)
	if (length(w) > 0)
		z[w] <- reverseComplement(z[w])
	
	if (type == 2L) {
		z <- RNAStringSet(z)
	} else if (type == 3L) {
		z <- translate(z,
			genetic.code=attr(x, "geneticCode"),
			...)
	}
	
	names(z) <- NULL
	
	return(z)
}

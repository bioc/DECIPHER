Disambiguate <- function(myXStringSet) {
	
	# error checking
	if (is(myXStringSet, "DNAStringSet")) {
		type <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		type <- 2L
	} else {
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	}
	
	myXStringSet <- .Call("expandAmbiguities",
		myXStringSet,
		ifelse(type == 1,
			"T",
			"U"),
		PACKAGE="DECIPHER")
	
	if (type == 1) {
		myXStringSetList <- DNAStringSetList(myXStringSet)
	} else {
		myXStringSetList <- RNAStringSetList(myXStringSet)
	}
	
	return(myXStringSetList)
}

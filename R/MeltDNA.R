MeltDNA <- function(myDNAStringSet,
	type="derivative",
	temps=50:100,
	ions=0.2) {
	
	# error checking
	if (is.character(myDNAStringSet)) {
		myDNAStringSet <- DNAStringSet(myDNAStringSet)
	} else if (!is(myDNAStringSet, "DNAStringSet")) {
		stop("myDNAStringSet must be a DNAStringSet.")
	}
	if (min(width(myDNAStringSet)) < 3)
		stop("All sequences in myDNAStringSet must be at least 3 nucleotides long.")
	if (length(type) == 0)
		stop("No type specified.")
	TYPES <- c("positional probabilities", "melt curves", "derivative curves")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(temps))
		stop("temps must be a numeric.")
	if (length(temps) == 0)
		stop("temps cannot be length 0.")
	if (!all(temps == cummax(temps)))
		stop("temps must be monotonically increasing.")
	if (length(unique(temps)) != length(temps))
		stop("temps cannot repeat.")
	if (type == 3 && length(temps) < 3)
		stop("At least three temperatures must be specified for a derivative curve.")
	if (!is.numeric(ions))
		step("ions must be a numeric.")
	if (ions < 0.01 || is.nan(ions))
		stop("Sodium equivilent concentration must be at least 0.01M.")
	
	u <- unique(myDNAStringSet)
	x <- .Call("meltPolymer",
		u,
		as.numeric(temps),
		ions,
		type,
		PACKAGE="DECIPHER")
	
	if (length(u) != length(myDNAStringSet)) {
		m <- match(myDNAStringSet, u)
		if (type == 1L) {
			x <- x[m]
		} else {
			x <- x[, m, drop=FALSE]
		}
	}
	
	return(x)
}

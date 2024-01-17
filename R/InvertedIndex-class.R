print.InvertedIndex <- function(x,
	...) {
	
	cat(paste("  An InvertedIndex built with:",
		"\n   * ",
		ifelse(length(x$alphabet) == 4L,
			"Nucleotide",
			"Amino acid"),
		" sequences: ",
		formatC(length(x$length), big.mark=",", format="f", digits=0L),
		"\n   * Total k-mers: ",
		formatC(length(x$index), big.mark=",", format="f", digits=0L),
		"\n   * Alphabet: ",
		paste(tapply(names(x$alphabet),
				x$alphabet,
				paste,
				collapse=""),
			collapse=", "),
		"\n   * K-mer size: ",
		x$k,
		"\n   * Step size: ",
		x$step,
		"\n",
		sep=""))
	
	invisible(x)
}

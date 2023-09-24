print.NonCoding <- function(x, ...) {
	if (!is(x, "NonCoding"))
		stop("x must be an object of class 'NonCoding'.")
	
	cat("NonCoding object with ",
		nrow(x$motifs),
		ifelse(nrow(x$motifs) == 1L,
			" motif, ",
			" motifs, "),
		nrow(x$hairpins),
		ifelse(nrow(x$hairpins) == 1L,
			" hairpin, and ",
			" hairpins, and "),
		attr(x, "K"),
		"-mer frequencies.\n",
		sep="")
	
	invisible(x)
}

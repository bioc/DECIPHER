`[.Genes` <- function(x, i, j, ...) {
	if (class(x)[1] != "Genes")
		stop("x must be an object of class 'Genes'.")
	a <- attributes(x)
	x <- unclass(x)
	if (missing(j)) {
		x <- x[i, j, drop=FALSE]
		class(x) <- "Genes"
		attr(x, "widths") <- a$widths
		attr(x, "geneticCode") <- a$geneticCode
		attr(x, "minGeneLength") <- a$minGeneLength
		attr(x, "allowEdges") <- a$allowEdges
	} else {
		x <- x[i, j, ...]
	}
	
	return(x)
}

print.Genes <- function(x, ...) {
	if (!is(x, "Genes"))
		stop("x must be an object of class 'Genes'.")
	N <- nrow(x)
	if (N == 0) {
		cat("Genes object of size 0.\n")
	} else {
		w1 <- which(x[, "Gene"] > 0)
		n1 <- length(w1)
		l1 <- x[w1, "End"] - x[w1, "Begin"] + 1L
		w2 <- which(1/x[, "Gene"] == Inf) # disregard negative zero
		n2 <- length(w2)
		l2 <- x[w2, "End"] - x[w2, "Begin"] + 1L
		w3 <- which(x[, "Gene"] < 0)
		n3 <- length(w3)
		l3 <- x[w3, "End"] - x[w3, "Begin"] + 1L
		
		cat("Genes object of size",
				formatC(N,
					big.mark=",",
					format="fg"),
			"specifying:")
		if (n1 == 1) {
			cat("\n1 protein coding gene of ",
				formatC(l1,
					big.mark=",",
					format="fg"),
				" nucleotides.",
				sep="")
		} else if (n1 > 1) {
			cat("\n",
				formatC(n1,
					big.mark=",",
					format="fg"),
				" protein coding genes from ",
				formatC(min(l1),
					big.mark=",",
					format="fg"),
				" to ",
				formatC(max(l1),
					big.mark=",",
					format="fg"),
				" nucleotides.",
				sep="")
		}
		if (n2 == 1) {
			cat("\n1 open reading frame of ",
				formatC(l2,
					big.mark=",",
					format="fg"),
				" nucleotides.",
				sep="")
		} else if (n2 > 1) {
			cat("\n",
				formatC(n2,
					big.mark=",",
					format="fg"),
				" open reading frames from ",
				formatC(min(l2),
					big.mark=",",
					format="fg"),
				" to ",
				formatC(max(l2),
					big.mark=",",
					format="fg"),
				" nucleotides.",
				sep="")
		}
		if (n3 == 1) {
			cat("\n1 non-coding RNA of ",
				formatC(l3,
					big.mark=",",
					format="fg"),
				" nucleotides.",
				sep="")
		} else if (n3 > 1) {
			cat("\n",
				formatC(n3,
					big.mark=",",
					format="fg"),
				" non-coding RNAs from ",
				formatC(min(l3),
					big.mark=",",
					format="fg"),
				" to ",
				formatC(max(l3),
					big.mark=",",
					format="fg"),
				" nucleotides.",
				sep="")
		}
		
		cat("\n\n")
		if (ncol(x) > 6) {
			print(cbind(as.data.frame(round(head(x[, 1:5, drop=FALSE]), 2)),
				`...`="...",
				as.data.frame(round(head(x[, ncol(x), drop=FALSE]), 2))))
		} else {
			print(as.data.frame(round(head(x[, seq_len(ncol(x)), drop=FALSE]), 2)))
		}
		
		if (N == 7) {
			cat("... with 1 more row.\n")
		} else if (N > 6) {
			cat("... with ",
				formatC(N - 6,
					big.mark=",",
					format="fg"),
				" more rows.\n",
				sep="")
		}
	}
	
	invisible(x)
}

plot.Genes <- function(x,
	xlim=c(1, 1e4),
	ylim=c(-100, 500),
	interact=FALSE,
	colorBy="Strand",
	colorRamp=colorRampPalette(c("darkblue", "darkred")),
	colorGenes="green4",
	...) {
	
	if (class(x)[1] != "Genes")
		stop("x must be an object of class 'Genes'.")
	if (is.character(colorBy)) {
		colorBy <- pmatch(colorBy[1], colnames(x))
		if (is.na(colorBy))
			stop("Invalid colorBy.")
		if (colorBy == -1)
			stop("Ambiguous colorBy.")
	} else {
		stop("colorBy must be a character string.")
	}
	if (typeof(colorRamp) != "closure")
		stop("colorRamp must be a function.")
	if (!is.na(colorGenes) &&
		typeof(colorGenes) != "character")
		stop("colorGenes must be NA or a character string.")
	if (!is.numeric(xlim) ||
		length(xlim) != 2L)
		stop("xlim must be a numeric vector of length 2.")
	if (min(xlim) < 1 ||
		max(xlim) > sum(attr(x, "widths")))
		stop("xlim must be within the bounds of the sequences specified in x.")
	if (!is.numeric(ylim) || length(ylim) != 2L)
		stop("ylim must be a numeric vector of length 2.")
	
	ws <- attr(x, "widths")
	ws <- cumsum(ws) - ws
	x[, "Begin"] <- x[, "Begin"] + ws[x[, "Index"]]
	x[, "End"] <- x[, "End"] + ws[x[, "Index"]]
	ws <- ws[-1L]
	
	var <- x[, colorBy]
	var <- .bincode(var,
		seq(min(var, na.rm=TRUE),
			max(var, na.rm=TRUE),
			length.out=101),
		include.lowest=TRUE)
	colors <- colorRamp(100)
	colors <- colors[var]
	
	genes <- which(x[, "Gene"] != 0)
	pos_genes <- which(x[, "Gene"] != 0 & x[, "Strand"] == 0)
	neg_genes <- which(x[, "Gene"] != 0 & x[, "Strand"] == 1)
	
	start <- xlim[1L]
	offset <- diff(xlim) + 1
	
	while(start > 0 &&
		start < max(x[, "End"] + offset)) {
		plot(NA,
			xlim=c(start, start + offset - 1L),
			ylim=ylim,
			xlab="Cumulative genome position",
			ylab="Total score")
		abline(h=0, v=ws, lwd=2, lty=3)
		segments(x[, "Begin"],
			x[, "TotalScore"],
			x[, "End"],
			x[, "TotalScore"],
			col=colors)
		segments(x[genes, "Begin"],
			x[genes, "TotalScore"],
			x[genes, "End"],
			x[genes, "TotalScore"],
			col=colorGenes,
			lwd=2)
		abline(v=x[pos_genes, "Begin"])
		abline(v=x[pos_genes, "End"],
			lty=5)
		abline(v=x[neg_genes, "Begin"],
			lty=5)
		abline(v=x[neg_genes, "End"])
		
		if (!interact)
			break
		
		l <- locator(n=1)
		if (length(l$x) != 1)
			break
		if (l$x > start + offset/2) {
			start <- start + offset - 1L
		} else {
			start <- start - offset + 1L
		}
	}
}

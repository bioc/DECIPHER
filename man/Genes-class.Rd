\name{Genes}
\alias{Genes-class}% the class
\alias{[.Genes}
\alias{print.Genes}
\alias{plot.Genes}
\title{
Genes objects and accessors
}
\description{
Gene prediction consist of delimiting the boundaries of regions that function as genes within a genome.  Class \code{Genes} provides objects and functions for storing the boundaries of genes and associated information resulting from gene prediction.
}
\usage{
\method{plot}{Genes}(x,
    xlim = c(1, 1e4),
    ylim = c(-100, 500),
    interact = FALSE,
    colorBy="Strand",
    colorRamp=colorRampPalette(c("darkblue", "darkred")),
    colorGenes="green4",
    \dots)

\method{print}{Genes}(x, \dots)

\method{[}{Genes}(x, i, j, \dots)
}
\arguments{
  \item{x}{
An object of class \code{Genes}.
}
  \item{xlim}{
Numeric vector of length 2 specifying the x-axis limits for plotting.
}
  \item{ylim}{
Numeric vector of length 2 specifying the y-axis limits for plotting.
}
  \item{interact}{
Logical determining whether the plot is interactive.  If \code{TRUE}, clicking the plot on the right or left side will scroll one frame in that direction.  To end interaction, either right-click, press the escape key, or press the stop button depending on the graphics device in use.
}
  \item{colorBy}{
Character string indicating the name of the column in \code{x} that should be used for coloring.  Unambiguous abbreviations are also permitted.
}
  \item{colorRamp}{
A function that will return \code{n} colors when given a number \code{n}.  Examples are \code{rainbow}, \code{heat.colors}, \code{terrain.colors}, \code{cm.colors}, or (the default) \code{colorRampPalette}.
}
  \item{colorGenes}{
Character string specifying the color of genes, or \code{NA} to color genes according to \code{colorBy}.
}
  \item{i}{
Numeric or character vector of row indices to extract from \code{x}.
}
  \item{j}{
Numeric or character vector of column indices to extract from \code{x}.  If \code{j} is missing, all columns are included and the returned object will also belong to class \code{Genes}.
}
  \item{\dots}{
Other optional parameters.
}
}
\details{
Objects of class \code{Genes} are stored as numeric matrices containing information pertaining to gene predictions.  The matrix columns include the index (\code{"Index"}) of the corresponding sequence in the original genome, the strand (\code{"Strand"}) where the gene is located (either \code{"+"} (\code{0}) or \code{"-"} (\code{1}), the beginning (\code{"Begin"}) and ending (\code{"End"}) positions of the gene, scores acquired during prediction, and whether (\code{!= 0}) or not (\code{0}) the region was predicted to be a gene.  Note that the start of the gene is at the beginning position when the strand is \code{"+"} and end when the strand is \code{"-"}.  By convention, rows with negative values in the \code{"Gene"} column represent non-coding RNAs and rows with positive values represent protein coding genes.

The \code{plot} method will show the total score of each prediction along the genome.  This is most useful when displaying the result of setting \code{allScores} to \code{TRUE} in \code{FindGenes}.  Here, possible genes on either strand will be shown (by default), with the predicted genes highlighted.  The beginning (solid) and ending (dashed) positions are denoted by vertical lines.  Note that the x-axis is cumulative genome position, and changes between genome sequences indices are demarcated by dashed vertical lines.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{ExtractGenes}}, \code{\link{FindGenes}}, \code{\link{WriteGenes}}
}
\examples{
# import a test genome
fas <- system.file("extdata",
	"Chlamydia_trachomatis_NC_000117.fas.gz",
	package="DECIPHER")
genome <- readDNAStringSet(fas)

x <- FindGenes(genome, allScores=TRUE)
x
head(unclass(x)) # the underlying structure

plot(x) # default coloring by "Strand"
# color by RBS score (blue is weak/low, red is strong/high)
plot(x, colorBy="RibosomeBindingSiteScore", colorGenes=NA)
# color by fraction of times a gene was chosen
plot(x, colorBy="FractionReps", colorGenes=NA)
# color by which codon model was selected for each ORF
plot(x, colorBy="CodonModel", xlim=c(1, 3e4))

# example of splitting a Genes object by "Strand"
tapply(seq_len(nrow(x)), x[, "Strand"], function(i) x[i])
}

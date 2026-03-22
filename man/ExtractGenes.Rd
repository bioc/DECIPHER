\name{ExtractGenes}
\alias{ExtractGenes}
\title{
Extract Predicted Genes from a Genome
}
\description{
Extracts predicted genes from the genome used for prediction.
}
\usage{
ExtractGenes(x,
             myDNAStringSet,
             type = "DNAStringSet",
             ...)
}
\arguments{
  \item{x}{
An object of class \code{Genes}.
}
  \item{myDNAStringSet}{
The \code{DNAStringSet} object used in generating \code{x}.
}
  \item{type}{
The class of sequences to return.  This should be (an unambiguous abbreviation of) one of \code{"AAStringSet"}, \code{"DNAStringSet"} (the default), or \code{"RNAStringSet"}.
}
  \item{\dots}{
Other parameters passed directly to \code{translate}.
}
}
\details{
Extracts a set of gene predictions as either DNA, mRNA, or proteins.
}
\value{
An \code{"AAStringSet"}, \code{"DNAStringSet"}, or \code{"RNAStringSet"} determined by \code{type}.
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\seealso{
\code{\link{FindGenes}}, \code{\link{Genes-class}}, \code{\link{WriteGenes}}
}
\examples{
# import a test genome
fas <- system.file("extdata",
	"Chlamydia_trachomatis_NC_000117.fas.gz",
	package="DECIPHER")
genome <- readDNAStringSet(fas)

x <- FindGenes(genome)
genes <- ExtractGenes(x, genome)
proteins <- ExtractGenes(x, genome, type="AAStringSet")
}

\name{InvertedIndex}
\alias{InvertedIndex-class}% the class
\alias{print.InvertedIndex}
\title{
InvertedIndex objects
}
\description{
\code{InvertedIndex} objects store k-mer locations and indexes in a set of sequences.
}
\usage{
\method{print}{InvertedIndex}(x,
     \dots)
}
\arguments{
  \item{x}{
An object of class \code{InvertedIndex}.
}
  \item{\dots}{
Other optional parameters.
}
}
\details{
Objects of class \code{InvertedIndex} are stored as a list.  The function \code{IndexSeqs} returns an object of class \code{InvertedIndex}. The information stored in an \code{InvertedIndex} can be displayed with \code{print}.
}
\references{
ES Wright (2024) "Fast and Flexible Search for Homologous Biological Sequences with DECIPHER v3". The R Journal, \bold{16(2)}, 191-200.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{IndexSeqs}}, \code{\link{SearchIndex}}
}
\examples{
# import target sequences
fas <- system.file("extdata", "PlanctobacteriaNamedGenes.fas.gz", package="DECIPHER")
seqs <- readAAStringSet(fas)

# build an inverted index
index <- IndexSeqs(seqs, K=6L)
index # print the index
}

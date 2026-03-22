\name{NonCoding}
\alias{NonCoding-class}% the class
\alias{print.NonCoding}
\title{
NonCoding Objects and Methods
}
\description{
Non-coding RNAs can be represented by their conserved sequence motifs, secondary structure, and k-mer frequencies.  Class \code{NonCoding} provides objects and functions for representing non-coding RNAs.
}
\usage{
\method{print}{NonCoding}(x, \dots)
}
\arguments{
  \item{x}{
An object of class \code{NonCoding}.
}
  \item{\dots}{
Other optional parameters.
}
}
\details{
Objects of class \code{NonCoding} are stored as lists containing a compact representation of a family of non-coding RNAs.  The first list component is a matrix of sequence motifs that identify the non-coding RNAs, the second is a matrix of hairpin loops that are conserved across the family, the third is a list of k-mer frequencies derived from representative sequences, and the fourth is a vector of log-odds scores for sequence lengths.  An optional fifth list component denotes the log-odds scores for dependencies among patterns.  Patterns are defined by their distance to either end of the non-coding RNA, which helps to identify the boundaries of the non-coding RNA in a genome.
}
\references{
Wright, E. S. (2021). FindNonCoding: rapid and simple detection of non-coding RNAs in genomes. Bioinformatics. https://doi.org/10.1093/bioinformatics/btab708
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\seealso{
\code{\link{LearnNonCoding}}, \code{\link{FindNonCoding}}
}
\examples{
data(NonCodingRNA_Bacteria)
x <- NonCodingRNA_Bacteria
print(x)
class(x)
attributes(x[[1]])
x[[1]] # the first non-coding RNA
x[[1]][["motifs"]] # sequence motifs
x[[1]][["hairpins"]] # hairpin loops
head(x[[1]][["kmers"]]) # k-mer frequencies
}

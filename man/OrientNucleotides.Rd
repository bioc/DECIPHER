\name{OrientNucleotides}
\alias{OrientNucleotides}
\title{
Orient Nucleotide Sequences
}
\description{
Orients nucleotide sequences to match the directionality and complementarity of specified reference sequences.
}
\usage{
OrientNucleotides(myXStringSet,
                  reference = which.max(width(myXStringSet)),
                  type = "sequences",
                  orientation = "both",
                  threshold = 0.05,
                  verbose = TRUE,
                  processors = 1)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} of unaligned sequences.
}
  \item{reference}{
The index of reference sequences with the same (desired) orientation.  By default the first sequence with maximum width will be used.
}
  \item{type}{
Character string indicating the type of results desired.  This should be (an abbreviation of) either \code{"sequences"}, \code{"orientations"}, or \code{"both"}.
}
  \item{orientation}{
Character string(s) indicating the allowed reorientation(s) of non-reference sequences.  This should be (an abbreviation of) either \code{"all"}, \code{"reverse"}, \code{"complement"}, and/or \code{"both"} (for reverse complement).
}
  \item{threshold}{
Numeric giving the decrease in k-mer distance required to adopt the alternative orientation.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
}
\details{
Biological sequences can sometimes have inconsistent orientation that interferes with their analysis.  \code{OrientNucleotides} will reorient sequences by changing their directionality and/or complementarity to match specified \code{reference} sequences in the same set, which are assume to have the same orientation.  The process works by finding the k-mer distance between the reference sequence(s) and each allowed \code{orientation} of the sequences.  Alternative orientations that lessen the distance by at least \code{threshold} are adopted.  Note that this procedure requires a moderately similar \code{reference} sequence be available for each sequence that needs to be reoriented.  Sequences for which a corresponding reference is unavailable will most likely be left alone because alternative orientations will not pass the \code{threshold}.  For this reason, it is recommended to specify several diverse sequences as \code{reference}s.
}
\value{
\code{OrientNucleotides} can return two types of results:  the relative orientations of sequences and/or the reoriented sequences.  If \code{type} is \code{"sequences"} (the default) then the reoriented sequences are returned.  If \code{type} is \code{"orientations"} then a character vector is returned that specifies whether sequences were reversed (\code{"r"}), complemented (\code{"c"}), reversed complemented (\code{"rc"}), or in the same orientation (\code{""}) as the reference sequences (marked by \code{NA}).  If \code{type} is \code{"both"} then the output is a list with the first component containing the \code{"orientations"} and the second component containing the \code{"sequences"}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{CorrectFrameshifts}}
}
\examples{
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
dna <- RemoveGaps(dna)
DNA <- dna # 175 sequences

# reorient subsamples of the first 169 sequences
s <- sample(169, 30)
DNA[s] <- reverseComplement(dna[s])
s <- sample(169, 30)
DNA[s] <- reverse(dna[s])
s <- sample(169, 30)
DNA[s] <- complement(dna[s])

DNA <- OrientNucleotides(DNA, reference=170:175)
DNA==dna # all were correctly reoriented
}

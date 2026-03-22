\name{RemoveGaps}
\alias{RemoveGaps}
\title{
Remove Gap Characters in Sequences
}
\description{
Removes gaps ("-" or "." characters) in a set of sequences, either deleting all gaps or only those shared by all sequences in the set.
}
\usage{
RemoveGaps(myXStringSet,
           removeGaps = "all",
           includeMask = FALSE,
           processors = 1)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object containing sequences.
}
  \item{removeGaps}{
Determines how gaps ("-" or "." characters) are removed in the sequences.  This should be (an unambiguous abbreviation of) one of \code{"none"}, \code{"all"} or \code{"common"}.
}
  \item{includeMask}{
Logical specifying whether to consider the mask character ("+") as a gap.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
}
\details{
The \code{removeGaps} argument controls which gaps are removed in \code{myXStringSet}.  Setting \code{removeGaps} to \code{"all"} will remove all gaps in the input sequences, whereas setting \code{removeGaps} to \code{"common"} will remove only gaps that exist in the same position in every sequence.  Therefore, the latter method will leave gaps in place that are not shared by every sequence, requiring that the sequences in \code{myXStringSet} all have the same width (i.e., be aligned).  Setting \code{removeGaps} to \code{"none"} will simply return \code{myXStringSet} unaltered.
}
\value{
An \code{XStringSet} of the same type as \code{myXStringSet}.
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\seealso{
\code{\link{AlignSeqs}}
}
\examples{
dna <- DNAStringSet(c("ACT-G-", "AC--G-"))
dna
RemoveGaps(dna, "all")
RemoveGaps(dna, "common")
}

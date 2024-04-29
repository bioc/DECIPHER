\name{BLOSUM}
\alias{BLOSUM}
\docType{data}
\title{
BLOSUM Amino Acid Substitution Matrices
}
\description{
The BLOSUM amino acid substitution matrices defined by Henikoff, S., & Henikoff, J. (1992).
}
\usage{data("BLOSUM")}
\format{
  The format is:
 num [1:24, 1:24, 1:15] 2.4 -0.6 0 0 -1.8 0.6 0 0 -1.2 0 ...
 - attr(*, "dimnames")=List of 3
  ..$ : chr [1:24] "A" "R" "N" "D" ...
  ..$ : chr [1:24] "A" "R" "N" "D" ...
  ..$ : chr [1:15] "30" "35" "40" "45" ...
}
\details{
Substitution matrix values represent the log-odds of observing an aligned pair of amino acids versus the likelihood of finding the pair by chance.  The PFASUM substitution matrices are stored as an array named by each sub-matrix's similarity threshold.  (See examples section below.)  In all cases values are in units of third-bits (\eqn{log(odds\ ratio)*3/log(2)}).
}
\source{
Henikoff, S., & Henikoff, J. (1992). Amino Acid Substitution Matrices from Protein Blocks. \emph{PNAS}, \bold{89(22)}, 10915-10919.
}
\examples{
data(BLOSUM)
BLOSUM62 <- BLOSUM[,, "62"] # the BLOSUM62 matrix
BLOSUM62["A", "R"] # score for A/R pairing

data(PFASUM)
plot(PFASUM[1:20, 1:20, "62"], BLOSUM62[1:20, 1:20])
abline(a=0, b=1)
}
\keyword{datasets}

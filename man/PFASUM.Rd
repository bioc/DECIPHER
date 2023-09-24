\name{PFASUM}
\alias{PFASUM}
\docType{data}
\title{
PFASUM Amino Acid Substitution Matrices
}
\description{
The PFASUM amino acid substitution matrices defined by Keul, F., et al. (2017).
}
\usage{data("PFASUM")}
\format{
  The format is:
 num [1:25, 1:25, 1:90] 0.9492 -1.7337 0.2764 1.8153 0.0364 ...
 - attr(*, "dimnames")=List of 3
  ..$ : chr [1:25] "A" "R" "N" "D" ...
  ..$ : chr [1:25] "A" "R" "N" "D" ...
  ..$ : chr [1:90] "11" "12" "13" "14" ...
}
\details{
Substitution matrix values represent the log-odds of observing an aligned pair of amino acids versus the likelihood of finding the pair by chance.  The PFASUM substitution matrices are stored as an array named by each sub-matrix's similarity threshold.  (See examples section below.)  In all cases values are in units of third-bits (\eqn{log(odds\ ratio)*3/log(2)}).
}
\source{
Keul, F., et al. (2017). PFASUM: a substitution matrix from Pfam structural alignments. \emph{BMC Bioinformatics}, \bold{18(1)}, 293.
}
\examples{
data(PFASUM)
PFASUM62 <- PFASUM[,, "62"] # the PFASUM62 matrix
PFASUM62["A", "R"] # score for A/R pairing

data(BLOSUM62)
plot(BLOSUM62[1:20, 1:20], PFASUM62[1:20, 1:20])
abline(a=0, b=1)
}
\keyword{datasets}

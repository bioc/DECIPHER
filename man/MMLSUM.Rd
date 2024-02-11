\name{MMLSUM}
\alias{MMLSUM}
\docType{data}
\title{
MMLSUM Amino Acid Substitution Matrices
}
\description{
The MMLSUM amino acid substitution matrices defined by Sumanaweera, D., et al. (2022).
}
\usage{data("MMLSUM")}
\format{
  The format is:
 num [1:21, 1:21, 1:9] 10.31 -10.46 -11.95 -13.59 -8.33 ...
 - attr(*, "dimnames")=List of 3
  ..$ : chr [1:21] "A" "R" "N" "D" ...
  ..$ : chr [1:21] "A" "R" "N" "D" ...
  ..$ : chr [1:9] "10" "20" "30" "40" ...
}
\details{
Substitution matrix values represent the log-odds of observing an aligned pair of amino acids versus the likelihood of finding the pair by chance.  The MMLSUM substitution matrices are stored as an array named by each sub-matrix's similarity threshold.  (See examples section below.)  In all cases values are in units of third-bits (\eqn{log(odds\ ratio)*3/log(2)}).
}
\source{
Sumanaweera, D., et al. (2022). Bridging the gaps in statistical models of protein alignment. \emph{Bioinformatics}, \bold{38(Supplement_1)}, i228-i237.
}
\examples{
data(MMLSUM)
MMLSUM60 <- MMLSUM[,, "60"] # the MMLSUM60 matrix
MMLSUM60["A", "R"] # score for A/R pairing

data(BLOSUM62)
BLOSUM62 <- BLOSUM62*3/2 # convert from half-bits to third-bits
plot(BLOSUM62[1:20, 1:20], MMLSUM60[1:20, 1:20])
abline(a=0, b=1)
}
\keyword{datasets}

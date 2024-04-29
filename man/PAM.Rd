\name{PAM}
\alias{PAM}
\docType{data}
\title{
PAM Amino Acid Substitution Matrices
}
\description{
The PAM amino acid substitution matrices defined by Dayhoff, M., et al. (1978).
}
\usage{data("PAM")}
\format{
  The format is:
 num [1:24, 1:24, 1:50] 10.5 -15 -10.5 -9 -15 -10.5 -7.5 -6 -16.5 -12 ...
 - attr(*, "dimnames")=List of 3
  ..$ : chr [1:24] "A" "R" "N" "D" ...
  ..$ : chr [1:24] "A" "R" "N" "D" ...
  ..$ : chr [1:50] "10" "20" "30" "40" ...
}
\details{
Substitution matrix values represent the log-odds of observing an aligned pair of amino acids versus the likelihood of finding the pair by chance.  The PAM substitution matrices are stored as an array named by each sub-matrix's percentage accepted mutations.  (See examples section below.)  In all cases values are in units of third-bits (\eqn{log(odds\ ratio)*3/log(2)}).
}
\source{
Dayhoff, M., et al. (1978). A model of evolutionary change in proteins. \emph{Atlas of protein sequence and structure}, \bold{5}, 345-358.
}
\examples{
data(PAM)
PAM250 <- PAM[,, "250"] # the PAM250 matrix
PAM250["A", "R"] # score for A/R pairing

data(BLOSUM)
plot(BLOSUM[1:20, 1:20, "62"], PAM250[1:20, 1:20])
abline(a=0, b=1)
}
\keyword{datasets}

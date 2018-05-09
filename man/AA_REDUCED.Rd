\name{AA_REDUCED}
\alias{AA_REDUCED}

\title{
Reduced amino acid alphabets
}
\description{
The \code{AA_REDUCED} list contains reductions of the standard amino acid alphabet (\code{AA_STANDARD}).
}
\usage{
AA_REDUCED
}
\details{
Reduced amino alphabets can sometimes improve sensitivity and specificity of finding homologous matches between amino acid sequences.  The first twelve \code{AA_REDUCED} alphabets were optimized for finding synteny between genomic sequences.  The remaining alphabets are from a review of published amino acid alphabets (Solis, 2015).
}
\references{
Solis, A. (2015). Amino acid alphabet reduction preserves fold information contained in contact interactions in proteins. \emph{Proteins: Structure, Function, and Genetics}, \bold{83(12)}, 2198-2216.
}
\seealso{
\code{\link{FindSynteny}}
}
\examples{
str(AA_REDUCED)
AA_REDUCED[[1]]
}
\keyword{data}

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
Reduced amino acid alphabets can sometimes improve sensitivity and specificity of finding homologous matches between amino acid sequences.  The first 12 \code{AA_REDUCED} alphabets were optimized for finding synteny between genomic sequences.  The next 113 alphabets are from a review of published amino acid alphabets (Solis, 2015).  The following 17 alphabets were optimized for amino acid classification.  The subsequent 18 alphabets are progressive mergers based on average similarities in \code{PFASUM}.  The following 25 alphabets were optimized for protein search.  The final 8 alphabets were optimized for clustering.
}
\references{
Solis, A. (2015). Amino acid alphabet reduction preserves fold information contained in contact interactions in proteins. \emph{Proteins: Structure, Function, and Genetics}, \bold{83(12)}, 2198-2216.
}
\seealso{
\code{\link{FindSynteny}}, \code{\link{LearnTaxa}}, \code{\link{PFASUM}}
}
\examples{
str(AA_REDUCED)
AA_REDUCED[[1]]
}
\keyword{data}

\name{WriteDendrogram}
\alias{WriteDendrogram}
\title{
Write a Dendrogram to Newick Format
}
\description{
Writes a dendrogram object to a file in Newick (also known as New Hampshire) parenthetic format.
}
\usage{
WriteDendrogram(x,
                file = "",
                quote = "'",
                space = " ",
                internalLabels = TRUE,
                digits = 10,
                append = FALSE)
}
\arguments{
  \item{x}{
An object of class \code{dendrogram}.
}
  \item{file}{
A connection or a character string naming the file path where the tree should be exported.  If "" (the default), the tree is printed to the standard output connection, the console unless redirected by sink.
}
  \item{quote}{
A single character used to quote labels, or an empty character string (i.e., \code{""}) to avoid quoting labels.
}
  \item{space}{
A single character (e.g., \code{"_"}) used to replace spaces in labels, or a space (i.e., \code{" "}) to leave spaces intact.
}
  \item{internalLabels}{
Logical indicating whether to write any ``edgetext'' preceding a node as an internal node label.
}
  \item{digits}{
The maximum number of digits to print for edge lengths.
}
  \item{append}{
Logical indicating whether to append to an existing \code{file}.  Only applicable if \code{file} is a character string.  If \code{FALSE} (the default), then the file is overwritten.
}
}
\details{
\code{WriteDendrogram} will write a dendrogram object to a \code{file} in standard Newick format.  Note that special characters (commas, square brackets, colons, semi-colons, and parentheses) present in leaf labels will likely cause a broken Newick file unless \code{quote} is a single or double quotation mark (the default).
}
\value{
\code{NULL}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{TreeLine}}, \code{\link{ReadDendrogram}}
}
\examples{
dists <- matrix(c(0, 10, 20, 10, 0, 5, 20, 5, 0),
    nrow=3,
    dimnames=list(c("dog", "elephant", "horse")))
dend <- TreeLine(myDistMatrix=dists, method="NJ")
WriteDendrogram(dend)
}

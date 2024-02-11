\name{IdLengths}
\alias{IdLengths}
\title{
Determine the Number of Characters in Each Sequence
}
\description{
Counts the number of standard and non-standard characters in each sequence.
}
\usage{
IdLengths(dbFile,
          tblName = "Seqs",
          type = "DNAStringSet",
          add2tbl = FALSE,
          batchSize = 10000,
          processors = 1,
          verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A database connection object or a character string specifying the path to a SQLite database file.
}
  \item{tblName}{
Character string specifying the table where the sequences are located.
}
  \item{type}{
The type of \code{XStringSet} being processed.  This should be (an abbreviation of) one of \code{"AAStringSet"}, \code{"DNAStringSet"}, or \code{"RNAStringSet"}.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{batchSize}{
Integer specifying the number of sequences to process at a time.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{IdLengths} is designed to efficiently determine the number of standard and non-standard characters in every sequence within a database.  Standard and non-standard characters are defined with respect to the \code{type} of the sequences.  For DNA and RNA sequences there are four standard characters and 11 non-standard characters (i.e., ambiguity codes).  For amino acid sequences there are 20 standard and seven non-standard characters (including stops).  Gap (``-''), missing (``.''), and mask (``+'') characters count toward the \code{width} but not the number of standard or non-standard characters.
}
\value{
A \code{data.frame} with the number of \code{standard} characters, \code{nonstandard} characters, and \code{width} of each sequence.  The \code{row.names} of the \code{data.frame} correspond to the "row_names" in the \code{tblName} of the \code{dbFile}.
}
\references{
ES Wright (2016) "Using DECIPHER v2.0 to Analyze Big Biological Sequence Data in R". The R Journal, \bold{8(1)}, 352-359.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{Add2DB}}
}
\examples{
if (require("RSQLite", quietly=TRUE)) {
	db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
	l <- IdLengths(db)
	head(l)
}
}

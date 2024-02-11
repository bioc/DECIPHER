\name{FormGroups}
\alias{FormGroups}
\title{
Forms Groups By Rank
}
\description{
Agglomerates sequences into groups within a specified size range based on taxonomic rank.
}
\usage{
FormGroups(dbFile,
           tblName = "Seqs",
           goalSize = 50,
           minGroupSize = 25,
           maxGroupSize = 5000,
           includeNames = FALSE,
           add2tbl = FALSE,
           verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A database connection object or a character string specifying the path to a SQLite database file.
}
  \item{tblName}{
Character string specifying the table where the taxonomic rank (i.e., ``organism'') information is located.
}
  \item{goalSize}{
Number of sequences required in each group to stop adding more sequences.
}
  \item{minGroupSize}{
Minimum number of sequences in each group required to stop trying to recombine with a larger group.
}
  \item{maxGroupSize}{
Maximum number of sequences in each group allowed to continue agglomeration.
}
  \item{includeNames}{
Logical indicating whether to include the formal scientific name in the group name.
}
  \item{add2tbl}{
Logical or a character string specifying the table name in which to add the result.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{FormGroups} uses the ``organism'' field in the \code{dbFile} table to group sequences with similar taxonomic rank.  Taxonomic rank information must be present in the \code{tblName}, such as that created by default when importing sequences from a GenBank formatted file.

Organism information contains the formal scientific name on the first line, followed by the taxonomic lineage on subsequent lines.  When \code{includeNames} is \code{TRUE} the formal scientific name is appended to the end of the group name, otherwise only the taxonomic lineage is used as the group name.

The algorithm ascends the taxonomic tree, agglomerating taxa into groups until the \code{goalSize} is reached.  If the group size is below \code{minGroupSize} then further agglomeration is attempted with a larger group.  If additional agglomeration results in a group larger than \code{maxGroupSize} then the agglomeration is undone so that the group is smaller.  Setting \code{minGroupSize} to \code{goalSize} avoids the creation of polyphyletic groups.  Note that this approach may often result in paraphyletic groups.
}
\value{
A \code{data.frame} with the \code{organism} and corresponding group name as \code{identifier}.  Note that quotes are stripped from group names to prevent problems that they may cause.  The \code{origin} gives the \code{organism} preceding the \code{identifier}.  The \code{count} denotes number of sequences corresponding to each \code{organism}.  If \code{add2tbl} is not \code{FALSE} then the ``identifier'' and ``origin'' columns are updated in \code{dbFile}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{IdentifyByRank}}
}
\examples{
if (require("RSQLite", quietly=TRUE)) {
	db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
	g <- FormGroups(db, goalSize=10, minGroupSize=5, maxGroupSize=20)
	head(g)
	tapply(g$count, g$identifier, sum)
}
}

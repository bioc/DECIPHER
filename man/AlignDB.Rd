\name{AlignDB}
\alias{AlignDB}
\title{
Align Two Sets of Aligned Sequences In A Sequence Database
}
\description{
Merges the two separate sequence alignments in a database.  The aligned sequences must have separate identifiers in the same table or be located in different database tables.
}
\usage{
AlignDB(dbFile,
        tblName = "DNA",
        identifier = "",
        type = "DNAStringSet",
        add2tbl = "DNA",
        batchSize = 10000,
        perfectMatch = NULL,
        misMatch = NULL,
        gapOpening = NULL,
        gapExtension = NULL,
        terminalGap = -1,
        substitutionMatrix = NULL,
        processors = NULL,
        verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A SQLite connection object or a character string specifying the path to the database file.
}
  \item{tblName}{
Character string specifying the table(s) where the sequences are located.  If two \code{tblName}s are provided then the sequences in both tables will be aligned.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If "" then all identifiers are selected.  If two \code{identifier}s are provided then the set of sequences matching each \code{identifier} will be aligned.
}
  \item{type}{
The type of \code{XStringSet} being processed.  This should be (an unambiguous abbreviation of) one of \code{"AAStringSet"}, \code{"DNAStringSet"}, or \code{"RNAStringSet"}.
}
  \item{add2tbl}{
Character string specifying the table name in which to add the aligned sequences.
}
  \item{batchSize}{
Integer specifying the number of sequences to process at a time.
}
  \item{perfectMatch}{
Numeric giving the reward for aligning two matching nucleotides in the alignment, or \code{NULL} to determine the value based on input type (DNA/RNA/AA).
}
  \item{misMatch}{
Numeric giving the cost for aligning two mismatched nucleotides in the alignment, or \code{NULL} to determine the value based on input type (DNA/RNA/AA).
}
  \item{gapOpening}{
Numeric giving the cost for opening a gap in the alignment, or \code{NULL} to determine the value based on input type (DNA/RNA/AA).
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap in the alignment, or \code{NULL} to determine the value based on input type (DNA/RNA/AA).
}
  \item{terminalGap}{
Numeric giving the cost for allowing leading and trailing gaps in the alignment.  Either two numbers, the first for leading gaps and the second for trailing gaps, or a single number for both.
}
  \item{substitutionMatrix}{
Either a substitution matrix representing the substitution scores for an alignment or the name of the amino acid substitution matrix to use in alignment.  The latter may be one of the following:  ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250''.  The default (NULL) will use the \code{perfectMatch} and \code{misMatch} penalties for DNA/RNA or ``BLOSUM62'' for AA.  (See examples section below.)
}
  \item{processors}{
The number of processors to use, or \code{NULL} (the default) for all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Sometimes it is useful to align two large sets of sequences, where each set of sequences is already aligned but the two sets are not aligned to each other.  \code{AlignDB} first builds a profile of each sequence set in increments of \code{batchSize} so that the entire sequence set is not required to fit in memory.  Next the two profiles are aligned using dynamic programming.  Finally, the new alignment is applied to all the sequences as they are incrementally added to the \code{add2tbl}.

Two \code{identifier}s or \code{tblName}s must be provided, indicating the two sets of sequences to align.  The sequences corresponding to the first \code{identifier} and \code{tblName} will be aligned to those of the second \code{identifier} or \code{tblName}.  The aligned sequences are added to \code{add2tbl} under a new identifier formed from the concatenation of the two \code{identifier}s or \code{tblName}s.  (See examples section below.)
}
\value{
Returns the number of newly aligned sequences added to the database.
}
\author{
Erik Wright \email{DECIPHER@cae.wisc.edu}
}
\seealso{
\code{\link{AlignProfiles}}, \code{\link{AlignSeqs}}, \code{\link{AlignTranslation}}
}
\examples{
gen <- system.file("extdata", "Bacteria_175seqs.gen", package="DECIPHER")
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")

# Align two tables and place result into a third
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(gen, "GenBank", dbConn, "Seqs1", tblName="Set1")
Seqs2DB(fas, "FASTA", dbConn, "Seqs2", tblName="Set2")
AlignDB(dbConn, tblName=c("Set1", "Set2"), add2tbl="AlignedSets")
l <- IdLengths(dbConn, "AlignedSets", add2tbl=TRUE)
BrowseDB(dbConn, tblName="AlignedSets") # all sequences have the same width
dbDisconnect(dbConn)

# Align two identifiers and place the result in the same table
dbConn <- dbConnect(SQLite(), ":memory:")
Seqs2DB(gen, "GenBank", dbConn, "Seqs1")
Seqs2DB(fas, "FASTA", dbConn, "Seqs2")
AlignDB(dbConn, identifier=c("Seqs1", "Seqs2"))
l <- IdLengths(dbConn, add2tbl=TRUE)
BrowseDB(dbConn) # note the sequences with a new identifier
dbDisconnect(dbConn)
}
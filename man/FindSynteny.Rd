\name{FindSynteny}
\alias{FindSynteny}
\title{
Finds Synteny in a Sequence Database
}
\description{
Finds syntenic blocks between groups of sequences in a database.
}
\usage{
FindSynteny(dbFile,
            tblName = "Seqs",
            identifier = "",
            useFrames = TRUE,
            alphabet = AA_REDUCED[[172]],
            geneticCode = GENETIC_CODE,
            sepCost = -3,
            sepPower = 0.5,
            gapCost = -12,
            gapPower = 0.5,
            shiftCost = 0,
            codingCost = 0,
            maxSep = 150,
            maxGap = 15,
            minScore = 100,
            N = 10,
            dropScore = -6,
            maskRepeats = TRUE,
            maskLCRs = TRUE,
            allowOverlap = FALSE,
            storage = 0.5,
            processors = 1,
            verbose = TRUE)
}
\arguments{
  \item{dbFile}{
A database connection object or a character string specifying the path to a SQLite database file.  The database should contain DNA sequences, typically with a distinct \code{identifier} for sequences belonging to each genome or chromosome.
}
  \item{tblName}{
Character string specifying the table where the sequences are located.
}
  \item{identifier}{
Optional character string used to narrow the search results to those matching a specific identifier.  If \code{""} then all identifiers are selected.  Repeated identifiers will find synteny between a sequence and itself, while blocking identical positions from matching in both sequences.
}
  \item{useFrames}{
Logical specifying whether to use 6-frame amino acid translations to help find more distant hits.  Using the \code{alphabet} is helpful when the genome is largely composed of coding DNA.  If \code{FALSE} then faster but less sensitive to distant homology.
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between sequences.  A non-reduced amino acid alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.
}
  \item{geneticCode}{
Either a character vector giving the genetic code to use in translation, or a list containing one genetic code for each identifier.  If a list is provided then it must be named by the corresponding identifiers in the database.
}
  \item{sepCost}{
Cost per nucleotide separation between hits to apply when chaining hits into blocks.
}
  \item{sepPower}{
Positive numeric specifying the power applied to the separation between hits before multiplying by \code{sepCost}.
}
  \item{gapCost}{
Cost for gaps between hits to apply when chaining hits into blocks.
}
  \item{gapPower}{
Positive numeric specifying the power applied to the number of gaps between hits before multiplying by \code{gapCost}.
}
  \item{shiftCost}{
Cost for shifting between different reading frames when chaining reduced amino acid hits into blocks.
}
  \item{codingCost}{
Cost for switching between coding and non-coding hits when chaining hits into blocks.
}
  \item{maxSep}{
Maximal separation (in nucleotides) between hits in the same block.
}
  \item{maxGap}{
The maximum number of gaps between hits in the same block.
}
  \item{minScore}{
The minimum score required for a chain of hits to become a block.  Higher values of \code{minScore} are less likely to yield false positives.
}
  \item{N}{
Numeric indicating the approximate number of k-mers that can be randomly selected before one is found by chance on average.  For example, the default value of \code{10} will set k-mer length such that every 10th k-mer is expected to have a match by chance.
}
  \item{dropScore}{
The change from maximal score required to stop extending blocks.
}
  \item{maskRepeats}{
Logical specifying whether to mask repeats when searching for hits.
}
  \item{maskLCRs}{
Logical indicating whether to mask low complexity regions when searching for hits.
}
  \item{allowOverlap}{
Logical specifying whether to permit blocks to overlap on the same sequence.
}
  \item{storage}{
Excess gigabytes available to store objects so that they do not need to be recomputed in later steps.  This should be a number between zero and a (modest) fraction of the available system memory.  Note that more than \code{storage} gigabytes may be required, but will not be stored for later reuse.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Long nucleotide sequences, such as genomes, are often not collinear or may be composed of many smaller segments (e.g., contigs).  \code{FindSynteny} searches for ``hits'' between sequences that can be chained into collinear ``blocks'' of synteny.  Hits are defined as k-mer exact nucleotide matches or k-mer matches in a reduced amino acid alphabet (if \code{useFrames} is \code{TRUE}).  Hits are chained into blocks as long as they are:  (1) within the same sequence, (2) within \code{maxSep} and \code{maxGap} distance, and (3) help maintain the score above \code{minScore}.  Blocks are extended from their first and last hit until their score drops below \code{dropScore} from the maximum that was reached.  This process results in a set of hits and blocks stored in an object of class ``Synteny''.
}
\value{
An object of class ``Synteny''.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\note{
\code{FindSynteny} is intended to be used on sets of sequences with up to ~200 million nucleotides total per identifier.  For this reason, better performance can sometimes be achieved by assigning a unique identifier to each chromosome belonging to a large genome.
}
\seealso{
\code{\link{AlignSynteny}}, \code{\link{Synteny-class}}
}
\examples{
if (require("RSQLite", quietly=TRUE)) {
	db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
	synteny <- FindSynteny(db)
	synteny
	pairs(synteny) # scatterplot matrix
}
}

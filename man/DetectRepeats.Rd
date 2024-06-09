\name{DetectRepeats}
\alias{DetectRepeats}
\title{
Detect Repeats in a Sequence
}
\description{
Detects approximate copies of sequence patterns that likely arose from duplication events and therefore share a common ancestor.
}
\usage{
DetectRepeats(myXStringSet,
              type = "tandem",
              minScore = 8,
              allScores = FALSE,
              maxCopies = 1000,
              maxPeriod = 1000,
              maxFailures = 3,
              maxShifts = 5,
              alphabet = AA_REDUCED[[17]],
              useEmpirical = TRUE,
              correctBackground = TRUE,
              processors = 1,
              verbose = TRUE,
              \dots)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of unaligned sequences.
}
  \item{type}{
Character string indicating the type of repeats to detect.  This should be (an abbreviation of) one of \code{"tandem"}, \code{"interspersed"}, or \code{"both"}.  Only \code{"tandem"} is possible when \code{myXStringSet} is an \code{AAStringSet}.  (See details section below.)
}
  \item{minScore}{
Numeric giving the minimum score of repeats in \code{myXStringSet} to report.
}
  \item{allScores}{
Logical specifying whether all repeats should be returned (\code{TRUE}) or only the top scoring repeat when there are multiple overlapping matches in the same region.
}
  \item{maxCopies}{
Numeric defining the maximum copy number of tandem repeat.  Since alignment complexity is quadratic in the number of repeat copies, setting a limit on the repeat copy number prevents very long repeats from becoming rate limiting.
}
  \item{maxPeriod}{
Numeric indicating the maximum periodicity of tandem repeats to consider.  Interspersed repeats will only be detected that are at least \code{maxPeriod} nucleotides apart.
}
  \item{maxFailures}{
Numeric determining the maximum number of failing attempts to extend a repeat that are permitted.  Numbers greater than zero may increase accuracy at the expense of speed, with decreasing marginal returns as \code{maxFailures} gets higher.
}
  \item{maxShifts}{
Numeric determining the maximum number of failing attempts to shift a repeat left or right that are permitted.  Numbers greater than zero may increase accuracy at the expense of speed, with decreasing marginal returns as \code{maxShifts} gets higher.
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between sequences.  A non-reduced amino acid alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.  Only applicable if \code{myXStringSet} is an \code{AAStringSet}.
}
  \item{useEmpirical}{
Logical specifying whether to use empirical signals of structurally-determined tandem repeats when scoring. Empirical signals include the number of repeats, their periodicity, and their amino acid composition when \code{myXStringSet} is an \code{AAStringSet}.
}
  \item{correctBackground}{
Logical controlling whether to correct the substitution matrix for the background distribution of letter frequencies on a per sequence basis.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{FindSynteny}} if \code{type} is \code{"interspersed"} or \code{"both"}.
}
}
\details{
Many sequences are composed of a substantial fraction of repetitive sequence.  Two main forms of repetition are tandem repeats and interspersed repeats, which can be caused by duplication events followed by divergence.  Detecting duplications is challenging because of variability in repeat length and composition due to evolution.  The significance of a repeat can be quantified by its time since divergence from a common ancestor.  \code{DetectRepeats} uses a seed-and-extend approach to identify candidate repeats, and tests whether a set of repeats is statistically significant using a background-corrected substitution matrix and gap (i.e., insertion and deletion) penalties.  A higher score implies the repeats are more conserved and, therefore, are more likely to have diverged within a finite amount of time from a common ancestor.  When \code{myXStringSet} is an \code{AAStringSet}, similarity includes agreement among predicted secondary structures.

Two possible \code{type}s of repeats are detectable:
* \code{type} is "tandem" (the default) or "both"
Contiguous approximate copies of a nucleotide or amino acid sequence.  First, repeated k-mers are identified along the length of the input sequence(s).  Once a k-mer seed is identified, repeated attempts are made to extend the repeat left and right, as well as optimize the beginning and ending positions.
* \code{type} is "interspersed" or "both"
Dispersed approximate copies of a nucleotide sequence on the same strand or opposite strands.  These are identified with \code{\link{FindSynteny}}, aligned with \code{\link{AlignSynteny}}, and then scored using the same statistical framework as tandem repeats.

The highest scoring repeat in each region is returned, unless \code{allScores} is \code{TRUE}, in which case overlapping repeats are permitted in the result.
}
\value{
If \code{type} is \code{"tandem"}, a \code{data.frame} giving the \code{"Index"} of the sequence in \code{myXStringSet}, \code{"Begin"} and \code{"End"} positions of tandem repeats, \code{"Left"} and \code{"Right"} positions of each repeat, and its \code{"Score"}.

If \code{type} is \code{"interspersed"}, a \code{data.frame} similar to the matrix in the lower diagonal of Synteny objects (see \code{\link{Synteny-class}}).

If \code{type} is \code{"both"}, a \code{list} with the above two elements.
}
\references{
Schaper, E., et al. (2012). Repeat or not repeat?-Statistical validation of tandem repeat prediction in genomic sequences. Nucleic Acids Research, \bold{40(20)}, 10005-17.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{ScoreAlignment}}
}
\examples{
fas <- system.file("extdata", "Human_huntingtin_cds.fas.gz", package="DECIPHER")
dna <- readDNAStringSet(fas)

x <- DetectRepeats(dna)
x

# number of tandem repeats
lengths(x[, "Left"])

# average periodicity of tandem repeats
per <- mapply(function(a, b) b - a + 1,
	x[, "Left"],
	x[, "Right"],
	SIMPLIFY=FALSE)
sapply(per, mean)

# extract a tandem repeat
i <- 1
reps <- extractAt(dna[[x[i, "Index"]]],
	IRanges(x[[i, "Left"]], x[[i, "Right"]]))
reps
reps <- AlignSeqs(reps, verbose=FALSE) # align the repeats
reps
BrowseSeqs(reps)

# highlight tandem repeats in the sequence
colors <- c("deeppink", "deepskyblue")
colors <- lapply(colors, function(x) col2rgb(x)/255)
cols <- vector("list", length(dna))
for (i in seq_along(cols)) {
	cols[[i]] <- matrix(0, nrow=3, ncol=width(dna)[i])
	for (j in which(x[, "Index"] == i)) {
		left <- x[[j, "Left"]]
		right <- x[[j, "Right"]]
		n <- 0
		for (k in seq_along(left)) {
			r <- left[k]:right[k]
			n <- n + 1
			if (n > length(colors))
				n <- 1
			cols[[i]][, r] <- colors[[n]]
		}
	}
}
BrowseSeqs(dna, patterns=cols)

# find interspersed (rather than tandem) repeats
data(yeastSEQCHR1)
chr1 <- DNAStringSet(yeastSEQCHR1)

if (require("RSQLite", quietly=TRUE)) {
	z <- DetectRepeats(chr1, type="interspersed")
	z
}
}

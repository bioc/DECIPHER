\name{AlignPairs}
\alias{AlignPairs}
\title{
Align pairs of sequences
}
\description{
Aligns pairs of sequences globally or locally using anchored adaptive banding.
}
\usage{
AlignPairs(pattern,
           subject,
           pairs = NULL,
           type = "values",
           perfectMatch = 1,
           misMatch = -2,
           gapOpening = -16,
           gapExtension = -1.2,
           substitutionMatrix = NULL,
           bandWidth = 50,
           dropScore = -100,
           processors = 1,
           verbose = TRUE)
}
\arguments{
  \item{pattern}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of (unaligned) sequences.
}
  \item{subject}{
A \code{XStringSet} object of (unaligned) sequences matching the type of the \code{pattern}.
}
  \item{pairs}{
Either \code{NULL} or a \code{data.frame} with \code{Pattern} and \code{Subject} indices to align and, optionally, a \code{Position} column containing matrices of anchor ranges.  If \code{Position} is provided then local alignment is performed around the anchors, unless virtual anchors are provided to force global alignment with terminal gap penalties.  If \code{pairs} is \code{NULL} (the default) then global alignment of respective \code{pattern} and \code{subject} sequences is performed without terminal gap penalties.  (See examples section below.)
}
  \item{type}{
Character string indicating the type of output desired.  This should be (an abbreviation of) one of \code{"values"}, \code{"sequences"}, or \code{"both"}.  (See value section below.)
}
  \item{perfectMatch}{
Numeric giving the reward for aligning matching nucleotides, which is used in the absence of a \code{substitutionMatrix} when the \code{pattern} is a \code{DNAStringSet} or \code{RNAStringSet}.
}
  \item{misMatch}{
Numeric determining the cost for aligning mismatched nucleotides, which is used in the absence of a \code{substitutionMatrix} when the \code{pattern} is a \code{DNAStringSet} or \code{RNAStringSet}.
}
  \item{gapOpening}{
Numeric giving the cost for opening a gap in the alignment.
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap in the alignment.
}
  \item{substitutionMatrix}{
Either a substitution matrix representing the substitution scores for an alignment (in third-bits) or the name of the amino acid substitution matrix to use in alignment.  The default (NULL) will use the \code{perfectMatch} and \code{misMatch} penalties for DNA/RNA or \code{PFASUM50} for AA.
}
  \item{bandWidth}{
Integer determining the number of positions included in the adaptive band, which should be at least as large as the largest expected insertion or deletion (i.e., gap).  Smaller values will accelerate alignment, potentially at the expense of accuracy.
}
  \item{dropScore}{
Numeric giving the decrease in score required to stop extending the region to the left or right of flanking anchors when performing local alignment.  Lower values find longer alignments at the expense of speed.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Performs pairwise alignment of \code{pattern} and \code{subject} sequences, either in their respective pairs or as specified by \code{pairs}.  Uses adaptive banding and (optionally) anchoring to accelerate the alignment.  Unlike \code{AlignProfiles}, \code{AlignPairs} can perform local alignment around anchor positions to align local regions of the \code{pattern} and/or \code{subject} sequences.  In the absence of anchors, \code{AlignPairs} performs global alignment without terminal gap penalties or with terminal gap penalties when virtual anchors are provided immediately outside the bounds of the sequences.  (See examples section below.)

\code{AlignPairs} is designed to maximize speed, and provides slightly less accuracy than using \code{AlignProfiles} for pairwise alignment.  Adaptive banding is applied with a fixed \code{bandWidth} to reduce memory consumption, rather than the dynamic band width used by \code{AlignProfiles}.  There are a few other differences: \code{AlignPairs} applies affine rather than Zipfian gap penalties, there is no option to incorporate secondary structures, scores are not standardized by length, gap penalties are not modulated around specific residues, and any anchors must be supplied in \code{pairs} rather than determined automatically.  For very dissimilar sequences, it is preferable to use \code{AlignTranslation} (best for coding sequences), \code{AlignSeqs} (best for nucleotide/protein sequences), or \code{AlignProfiles}.

\code{AlignPairs} does not directly output the pairwise alignments.  Instead, it outputs statistics about the alignment and the position(s) of gaps in each sequence.  This makes alignment more efficient because no sequences are copied.  For many applications only the percent identity or number of gaps is needed, which can be calculated directly from the returned \code{data.frame}. However, the aligned sequences can also easily be obtained from the output if desired.  (See examples section below.)
}
\value{
If \code{type} is \code{"values"} (the default), a \code{data.frame} is returned with one alignment per input \code{pattern} or row of \code{pairs} if not \code{NULL}.  Columns are defined as the sequences' index in \code{pattern} (\code{Pattern}), start position in the \code{pattern} sequence (\code{PatternStart}), end position in the \code{pattern} sequence (\code{PatternEnd}), index in the \code{subject} sequence (\code{Subject}), start position in the \code{subject} sequence (\code{SubjectStart}), end position in the \code{subject} sequence (\code{SubjectEnd}), number of matching positions in the alignment (\code{Matches}), number of mismatched positions in the alignment (\code{Mismatches}), total number of positions in the alignment (\code{AlignmentLength}), alignment score (\code{Score}), and the position/length of gaps in the \code{pattern} and \code{subject} (i.e., \code{PatternGapPosition}, \code{PatternGapLength}, \code{SubjectGapPosition}, & \code{SubjectGapLength}).

If \code{type} is \code{"sequences"}, a \code{list} containing two components: the aligned \code{pattern} and \code{subject}.

If \code{type} is \code{"both"}, a \code{list} with three components: the values, aligned \code{pattern}, and aligned \code{subject}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignProfiles}}, \code{\link{IndexSeqs}}, \code{\link{SearchIndex}}
}
\examples{
# import target sequences and build an inverted index
fas <- system.file("extdata", "PlanctobacteriaNamedGenes.fas.gz", package="DECIPHER")
target <- readAAStringSet(fas)
index <- IndexSeqs(target, K=6L)
index

# import query sequences and search the index
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
query <- translate(dna)
hits <- SearchIndex(query, index, scoreOnly=FALSE)
head(hits)
nrow(hits) # number of query/target pairs

# locally pairwise align the query/target pairs representing each hit
aligned <- AlignPairs(query, target, hits) # local alignment around hits
head(aligned)

# plot the percent identity of each alignment versus score
PIDs <- aligned$Matches/aligned$AlignmentLength
plot(PIDs, hits$Score) # versus the hit score
plot(PIDs, aligned$Score) # versus the alignment score

# plot the number of pattern gaps versus subject gaps per alignment
subjectGaps <- sapply(aligned$SubjectGapLength, sum)
patternGaps <- sapply(aligned$PatternGapLength, sum)
plot(subjectGaps, patternGaps, pch=16, col="#00000011")
abline(a=0, b=1) # identity line (y = x)

# extract the pairwise aligned regions
alignments <- AlignPairs(query, target, hits, type="sequences")
BrowseSeqs(c(alignments[[1]][1], alignments[[2]][1])) # view the first pair

# perform global pairwise alignment by creating virtual endpoint anchors
hits$Position <- mapply(function(x, y, z)
		cbind(matrix(0L, 4), x, matrix(c(y, y, z, z), 4)),
	hits$Position,
	width(query)[hits$Pattern] + 1L,
	width(target)[hits$Subject] + 1L,
	SIMPLIFY=FALSE)
aligned <- AlignPairs(query, target, hits) # penalizes terminal gaps
head(aligned) # all alignments now span start to end of the sequences

# perform global pairwise alignment of sequences without anchors (approach 1)
pattern <- query[rep(1, length(query))] # first sequence
subject <- query
aligned1 <- AlignPairs(pattern, subject) # no terminal gap penalties
head(aligned1)
# perform global pairwise alignment of sequences without anchors (approach 2)
pairs <- data.frame(Pattern=1L, Subject=seq_along(query))
aligned2 <- AlignPairs(query, query, pairs) # no terminal gap penalties
head(aligned2) # note the Pattern column is always 1 (first sequence)
}

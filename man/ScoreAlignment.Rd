\name{ScoreAlignment}
\alias{ScoreAlignment}
\title{
Score a Multiple Sequence Alignment
}
\description{
Calculates a score for a multiple sequence alignment based on either sum-of-pairs or sum-of-adjacent-pairs scoring.
}
\usage{
ScoreAlignment(myXStringSet,
               method = "pairs",
               type = "sum",
               perfectMatch = 1,
               misMatch = 0,
               gapOpening = -7.5,
               gapExtension = -0.6,
               substitutionMatrix = NULL,
               structures = NULL,
               structureMatrix = NULL,
               includeTerminalGaps = FALSE,
               weight = 1)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of aligned sequences.
}
  \item{method}{
Character string indicating the method of scoring.  This should be (an abbreviation of) one of \code{"pairs"} for sum-of-pairs or \code{"adjacent"} for sum-of-adjacent-pairs.  (See details section below.)
}
  \item{type}{
Character string giving the type of result.  This should be (an abbreviation of) one of \code{"sum"} for the total score or \code{"scores"} for a vector with one score per site (column) in \code{myXStringSet}.
}
  \item{perfectMatch}{
Numeric giving the reward for aligning two matching nucleotides in the alignment.  Only used for \code{DNAStringSet} or \code{RNAStringSet} inputs.
}
  \item{misMatch}{
Numeric giving the cost for aligning two mismatched nucleotides in the alignment.  Only used for \code{DNAStringSet} or \code{RNAStringSet} inputs.
}
  \item{gapOpening}{
Numeric giving the cost for opening and closing a gap in the alignment.
}
  \item{gapExtension}{
Numeric giving the cost for extending an open gap in the alignment.
}
  \item{substitutionMatrix}{
Either a substitution matrix representing the substitution scores for an alignment (in third-bits) or the name of the amino acid substitution matrix to use in alignment.  The latter may be one of the following:  ``BLOSUM45'', ``BLOSUM50'', ``BLOSUM62'', ``BLOSUM80'', ``BLOSUM100'', ``PAM30'', ``PAM40'', ``PAM70'', ``PAM120'', ``PAM250'', ``MIQS'', ``MMLSUM10'' to ``MMLSUM90'' (by tens), or ``PFASUM11'' to ``PFASUM100'' (by ones).  The default (NULL) will use the \code{perfectMatch} and \code{misMatch} penalties for DNA/RNA or \code{PFASUM50} for AA.
}
  \item{structures}{
Either \code{NULL} (the default) or a list of matrices with one list element per sequence in \code{myXStringSet}.
}
  \item{structureMatrix}{
A structure matrix if \code{structures} are supplied, or \code{NULL} otherwise.
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" or "." characters on each end of the sequence) into the calculation of score.
}
  \item{weight}{
A numeric vector of weights for each sequence, or a single number implying equal weights.
}
}
\details{
Sum-of-pairs scoring is the standard way to judge whether a set of sequences are homologous.  \code{ScoreAlignment} calculates the sum-of-pairs score for \code{myXStringSet} when \code{method} is \code{"pairs"}.  This score can also be used to compare among different alignments of the same sequences.  If \code{method} is \code{"adjacent"} then the sum-of-adjacent-pairs scores is calculated, where each sequence is compared to the next sequence.  Hence, the input order of sequences in \code{myXStringSet} matters when \code{method} is \code{"adjacent"}.

Both scores are linearly related to the number of sequences in the alignment and the number of sites in the alignment.  Therefore, it is possible to normalize the score by dividing by the \code{width} and \code{length} (minus \code{1}) of the \code{myXStringSet}.
}
\value{
A single \code{numeric} score.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AlignSeqs}}, \code{\link{PFASUM}}
}
\examples{
# small example
x <- DNAStringSet(c("C-G", "CTG", "C-G", "CTG"))
ScoreAlignment(x, method="pairs", gapOpening=-1) # +3 -1 +3 = 5
ScoreAlignment(x, method="adjacent", gapOpening=-1) # +3 -3 +3 = 3

# DNA alignment with the defaults
fas <- system.file("extdata", "Streptomyces_ITS_aligned.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
dna # input alignment
ScoreAlignment(dna, method="pairs")
ScoreAlignment(dna, method="adjacent")

# provide a DNA substitution matrix for greater discerning power
sub <- matrix(c(1.5, -2.134, -0.739, -1.298,
		-2.134, 1.832, -2.462, 0.2,
		-0.739, -2.462, 1.522, -2.062,
		-1.298, 0.2, -2.062, 1.275),
	nrow=4,
	dimnames=list(DNA_BASES, DNA_BASES))
ScoreAlignment(dna, substitutionMatrix=sub)

# use structures with an amino acid alignment
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
aa <- AlignTranslation(dna, type="AAStringSet")
structureMatrix <- matrix(c(0.187, -0.8, -0.873,
		-0.8, 0.561, -0.979,
		-0.873, -0.979, 0.221),
	nrow=3,
	dimnames=list(c("H", "E", "C"), c("H", "E", "C")))
ScoreAlignment(aa,
	structures=PredictHEC(aa, type="probabilities"),
	structureMatrix=structureMatrix)
}

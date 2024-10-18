\name{DistanceMatrix}
\alias{DistanceMatrix}
\title{
Calculate the Distances Between Sequences
}
\description{
Calculates a distance matrix for an \code{XStringSet}.  Each element of the distance matrix corresponds to the dissimilarity between two sequences in the \code{XStringSet}.
}
\usage{
DistanceMatrix(myXStringSet,
               method = "overlap",
               type = "matrix",
               includeTerminalGaps = FALSE,
               penalizeGapLetterMatches = TRUE,
               minCoverage = 0,
               correction = "none",
               processors = 1,
               verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet} object of aligned sequences.
}
  \item{type}{
Character string indicating the type of output desired.  This should be either \code{"matrix"} or \code{"dist"}.  (See value section below.)
}
  \item{method}{
Character string determining the region in which distance is calculated.  This should be (an unambiguous abbreviation of) one of \code{"overlap"}, \code{"shortest"}, or \code{"longest"}.  The default \code{method} (\code{"overlap"}) calculates distance in the overlapping region between terminal gaps when \code{includeTerminalGaps} is \code{FALSE} and the entire alignment otherwise.  Setting \code{method} to \code{"shortest"} or \code{"longest"} will use the region between the start and end of the shortest or longest sequence, respectively, for each pairwise distance.  The \code{method} is only relevant when \code{includeTerminalGaps} is \code{TRUE}.
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" or "." characters on each end of the sequence) into the calculation of distance.
}
  \item{penalizeGapLetterMatches}{
Logical specifying whether or not to consider gap-to-letter matches as mismatches.  If \code{FALSE}, then gap-to-letter matches are not included in the total length used to calculate distance.  The default (\code{TRUE}) is to penalize gap-to-letter mismatches the same as letter-to-letter mismatches.  If \code{NA} then gap-to-letter mismatches are only penalized once per insertion or deletion, i.e., when changing to gap-to-letter or letter-to-gap.
}
  \item{minCoverage}{
Numeric giving the minimum fraction of sequence positions (not gap or mask) that must be overlapping in each pair.  If positive then coverage is relative to the shortest sequence.  If negative then coverage is relative to both sequences.  Sequences failing to meet \code{minCoverage} will be assigned \code{NA} distances.  Note that completely non-overlapping sequences are always given \code{NA} distances, regardless of \code{minCoverage}, unless \code{includeTerminalGaps} is \code{TRUE} (i.e., distance = 100\%).
}
  \item{correction}{
The evolutionary model used for distance correction.  This should be (an abbreviation of) either \code{"none"}, \code{"Jukes-Cantor"} (i.e., \code{"JC69"}), \code{"F81"}, or \code{"Poisson"}.  For \code{"F81"} letter frequencies are derived from \code{myXStringSet}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The uncorrected (\code{correction = "none"}) distance matrix represents the hamming distance between each of the sequences in \code{myXStringSet}.  Ambiguity can be represented using the characters of the \code{IUPAC_CODE_MAP} for \code{DNAStringSet} and \code{RNAStringSet} inputs, or using the \code{AMINO_ACID_CODE} for an \code{AAStringSet} input.  For example, the distance between an 'N' and any other nucleotide base is zero.  The letters B (N or D), J (I or L), Z (Q or E), and X (any letter) are degenerate in the \code{AMINO_ACID_CODE}.

If \code{includeTerminalGaps = FALSE} then terminal gaps ("-" or "." characters) are not included in sequence length.  This can be faster since only the positions common to each pair of sequences are compared.  Sequences with no overlapping region in the alignment are given a value of \code{NA}, unless \code{includeTerminalGaps = TRUE}, in which case distance is 100\%.  Masked characters (\code{"+"}) in either sequence are not considered in distance.

Penalizing gap-to-letter mismatches specifies whether to penalize these special mismatch types and include them in the total length when calculating distance.  Both "-" and "." characters are interpreted as gaps.  The default behavior is to calculate distance as the fraction of positions that differ across the region of the alignment shared by both sequences (not including gap-to-gap matches).

Three correction factors are available to transform distances into an expected number of changes per site.  All can be applied to nucleotides or amino acids, but two are intended for nucleotides and one for amino acids.  The two nucleotide correction factors are \code{"JC69"} and \code{"F81"}, which are described in \code{\link{MODELS}}.  Both transform raw distance (\eqn{d}) by \deqn{-E * \log \left( 1 - d / E \right)}, where \deqn{E = 1 - \sum_{i \in sym} f_i^2} and (\eqn{f}) is the relative frequency of each symbol (\eqn{sym}).  In the \code{"JC69"} model symbols are assumed to have equal frequency, whereas in the \code{"F81"} model the symbol frequencies are empirically derived from the input \code{myXStringSet}.  The correction factor intended for amino acids is \code{"Poisson"}, which assumes the number of substitutions per site follows the Poisson distribution.  The raw distance (\eqn{d}) is transformed by \deqn{- \log \left( 1 - d \right)}.  Note that gaps are treated as an additional symbol when \code{penalizeGapLetterMatches} is \code{TRUE}, and all correction factors will result in an infinite distance when the sequences are 100\% different (i.e., \eqn{d = 1}).
}
\value{
If \code{type} is \code{"matrix"}, a symmetric matrix where each element is the distance between the sequences referenced by the respective row and column.  The \code{dimnames} of the matrix correspond to the \code{names} of the \code{XStringSet}.

If \code{type} is \code{"dist"}, an object of \code{class} \code{"dist"} that contains one triangle of the distance matrix as a vector.  Since the distance matrix is symmetric, storing only one triangle is more memory efficient.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{TreeLine}}
}
\examples{
# example of using the defaults:
dna <- DNAStringSet(c("ACTG", "ACCG"))
dna
DistanceMatrix(dna)

# changing the output type to "dist":
d <- DistanceMatrix(dna, type="dist")
d
length(d) # minimal memory space required
m <- as.matrix(d)
length(m) # more memory space required

# supplying an AAStringSet
aa <- AAStringSet(c("ASYK", "ATYK", "CTWN"))
aa
DistanceMatrix(aa)

# correct for multiple substitutions per site:
DistanceMatrix(aa, correction="Poisson")

# defaults compare intersection of internal ranges:
dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
dna
d <- DistanceMatrix(dna)
d
# d[1,2] is 1 base in 4 = 0.25

# compare union of internal positions, without terminal gaps:
dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
dna
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE)
d
# d[1,2] is now 2 bases in 5 = 0.40

# gap ("-") and unknown (".") characters are interchangeable:
dna <- DNAStringSet(c("ANGCT.", ".ACCT-"))
dna
d <- DistanceMatrix(dna, includeTerminalGaps=TRUE)
d
# d[1,2] is still 2 bases in 5 = 0.40

# compare different methods for calculating distance:
dna <- DNAStringSet(c("--ACTG", "TGAGT-"))
dna
DistanceMatrix(dna, method="overlap") # 1/3
DistanceMatrix(dna, method="shortest",
               includeTerminalGaps=FALSE) # 1/3
DistanceMatrix(dna, method="shortest",
               includeTerminalGaps=TRUE) # 2/4
DistanceMatrix(dna, method="shortest",
               includeTerminalGaps=TRUE,
               penalizeGapLetterMatches=FALSE) # 1/3
DistanceMatrix(dna, method="longest",
               includeTerminalGaps=FALSE) # 1/3
DistanceMatrix(dna, method="longest",
               includeTerminalGaps=TRUE) # 3/5
DistanceMatrix(dna, method="longest",
               includeTerminalGaps=TRUE,
               penalizeGapLetterMatches=FALSE) # 1/3
DistanceMatrix(dna, method="overlap",
               minCoverage=1) # NA (insufficient overlap)
DistanceMatrix(dna, method="overlap",
               minCoverage=0.75) # 3/4 sites covered in shorter
DistanceMatrix(dna, method="overlap",
               minCoverage=-0.75) # 3/5 sites covered in longer

# neither internal nor external gap/gap matches are considered:
dna <- DNAStringSet(c("--A-CTA", "-AG-C--"))
dna
DistanceMatrix(dna) # 1/2
DistanceMatrix(dna, includeTerminalGaps=TRUE) # 4/5
DistanceMatrix(dna, includeTerminalGaps=TRUE,
               method="shortest") # 2/3
DistanceMatrix(dna, includeTerminalGaps=TRUE,
               method="longest") # 3/4
}

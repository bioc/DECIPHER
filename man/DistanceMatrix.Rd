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
               penalizeGapLetterMatches = FALSE,
               minCoverage = 0,
               correction = NA,
               substitutionMatrix = NULL,
               frequencies = NULL,
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
Logical specifying whether to correct distance for gap-to-letter mismatches.  The default (\code{FALSE}) ignores gaps paired with letters, and \code{TRUE} considers them as mismatches.  If \code{NA} then gap-to-letter mismatches are only penalized once per run of gaps (i.e., insertion or deletion).  If \code{correction} is not \code{NA} or \code{frequencies} is not \code{NULL}, the evolutionary distance due to insertion and deletion (Tajima, 1984) is added to the standard distance when \code{penalizeGapLetterMatches} is not \code{FALSE}.
}
  \item{minCoverage}{
Numeric giving the minimum fraction of sequence positions (not gap or mask) that must be overlapping in each pair.  If positive then coverage is relative to the shortest sequence.  If negative then coverage is relative to both sequences.  Sequences failing to meet \code{minCoverage} will be assigned \code{NA} distances.  Note that completely non-overlapping sequences are always given \code{NA} distances, regardless of \code{minCoverage}, unless \code{includeTerminalGaps} is \code{TRUE} and \code{penalizeGapLetterMatches} is not \code{FALSE} (i.e., distance = 100\%).
}
  \item{correction}{
The evolutionary model used for distance correction.  This should be either \code{NA} or (an abbreviation of) either \code{"JC69"}, \code{"K80"}, \code{"F81+F"}, \code{"HKY85+F"}, \code{"T92+F"}, \code{"TN93+F"}, \code{"Poisson"}, or one of the protein \code{MODELS}.  The default (\code{NA}) will return uncorrected (i.e., Hamming) distances.  Models with fixed frequencies (\code{+F}) are empirically derived from \code{myXStringSet} (e.g., \code{"F81+F"} or \code{"WAG+F"}) unless provided as \code{frequencies}.  (See details section below.)
}
  \item{substitutionMatrix}{
A symmetric matrix, \code{dist} object, or a single character string specifying a substitution matrix.  Interpreted as costs if \code{frequencies} is \code{NULL} (the default), and interpreted as relative substitution rates otherwise.  If a built-in \code{substitutionMatrix} is named (e.g., \code{"BLOSUM62"}), the negative of the matrix is used to obtain distances rather than similarities.  Only used if \code{correction} is \code{NA}.  (See examples section below.)
}
  \item{frequencies}{
Numeric vector containing relative letter \code{frequencies} or \code{NULL} (the default) if unnecessary.  Must be named and in the same order as \code{substitutionMatrix}.  Only considered when \code{substitutionMatrix} is specified, in which case values are interpreted as stationary \code{frequencies} and \code{substitutionMatrix} must contain the corresponding substitution when \code{correction} is \code{NA} (the default).
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{DistanceMatrix} computes pairwise distances between all sequences in \code{myXStringSet}.  There are three main modes in which \code{DistanceMatrix} is typically used: (1) to obtain Hamming distances with or without gap penalties, (2) to obtain phylogenetic distances by correcting for multiple substitutions per site based on a model of evolution, and (3) to obtain normalized distances based on a \code{substitutionMatrix} containing transition costs.  For (1), distances are in units of differences/site, and the key parameters to consider are \code{penalizeGapLetterMatches} and \code{method} if \code{includeTerminalGaps} is \code{TRUE}.  For (2), distances are in units of substitutions per site, and the key parameters are \code{correction} or \code{substitutionMatrix} and \code{frequencies}.  For (3), distances are in the same units as the \code{substitutionMatrix} normalized by the number of shared sites, and \code{substitutionMatrix} is the main parameter to consider.  These parameters are described in more detail below and several examples are provided in the examples section.

The defaults compute an uncorrected distance matrix representing the Hamming distance between each of the sequences in \code{myXStringSet}, also known as the p-distance.  That is, matches are assigned a cost of \code{0} and mismatches are given a cost of \code{1}, and distances is computed by dividing the number of mismatches by the number of matches plus mismatches.  Ambiguity can be represented using the characters of the \code{IUPAC_CODE_MAP} for \code{DNAStringSet} and \code{RNAStringSet} inputs, or using the \code{AMINO_ACID_CODE} for an \code{AAStringSet} input.  For example, the distance between an 'N' and any other nucleotide base is zero.  The letters B (N or D), J (I or L), Z (Q or E), and X (any letter) are degenerate in the \code{AMINO_ACID_CODE}.

If \code{includeTerminalGaps = FALSE} then terminal gaps ("-" or "." characters) are not included in pairwise comparisons.  This can be faster since only the positions common to each pair of sequences are compared.  Sequences with no overlapping region in the alignment are given a value of \code{NA}, unless \code{includeTerminalGaps = TRUE}, in which case distance is 100\%.  Both "-" and "." characters are interpreted as gaps, and masked characters (\code{"+"}) in either sequence are not considered in distance.  The default behavior is to calculate distance as the fraction of positions that differ across the region of the alignment shared by both sequences (not including gaps).  Gap-to-letter matches are treated as mismatches if \code{penalizeGapLetterMatches} is \code{TRUE}.  Only the first gap per run of gaps is penalized if \code{penalizeGapLetterMatches} is \code{NA}, which represents an event based model where insertions and deletions are considered as single events analogous to individual substitutions (also known as a block model).

Multiple \code{correction} factors are available to transform distances into an expected number of changes per site.  Three models are generic, although two are intended for nucleotides (\code{"JC69"} and \code{"F81+F"}) and one for amino acids (\code{"Poisson"}).  There are four models specifically for nucleotides (\code{"K80"}, \code{"HKY85+F"}, \code{"T92+F"}, and \code{"TN93+F"}) that introduce different transition (versus transversion) rates.  The analytical formulas for \code{"JC69"} and \code{"Poisson"} give the exact distance, while the others use approximation formulas to circumvent parameter estimation (McGuire, 1999; Tamura, 1993).  In order of decreasing exactness they are: \code{"K80"}, \code{"F81+F"}, \code{"HKY85+F"}, \code{"TN93+F"}, and \code{"T92+F"}.  Models with fixed frequencies (\code{+F}) derive empirical frequencies from the input \code{myXStringSet}, unless they are given as \code{frequencies}.

The remaining \code{correction}s are defined by \code{MODELS} of amino acid evolution that provide both a matrix of substitution rates and stationary \code{frequencies}.  These models can be modified with empirical frequencies (\code{+F}), or \code{frequencies} can be given separately.  It is also feasible to specify the \code{substitutionMatrix} (rates) and \code{frequencies} separately when \code{correction} is \code{NA} (the default).  A maximum likelihood distance is returned based on fitting \deqn{P = e^(Q*d)}, where \eqn{P} is the observed substitution probabilities and \eqn{Q} is derived from \code{substitutionMatrix} such that \deqn{Q*frequencies = 0}.  Here, \eqn{d} represents twice the estimated number of substitutions per sites that occurred since a pair of sequences most recent common ancestor, i.e., their phylogenetic distance.

If \code{penalizeGapLetterMatches} is not \code{FALSE} when using a model of evolution, the component of distance due to insertion and deletion is added to the distance (Tajima, 1984).  That is, the traditional evolutionary distance due to multiple substitutions is added to \deqn{-2 * \log \left( n_xy / \left( n_x + n_y \right)^0.5 \right)}, where \eqn{n_xy} is the number of shared sites, \eqn{n_x} is the number of sites in one sequence and \eqn{n_y} is the number in the other.  This partly corrects for multiple insertion and deletion events in a related manner to correcting for multiple substitutions per site.  The resulting distances can be considered to have units of changes per site, since substitutions, insertions, and deletions all contribute to the distance.

If a \code{substitutionMatrix} is given without \code{frequencies}, the matrix is interpreted as the cost for each type of substitution.  It is also possible to specify \code{substitutionMatrix} as a character string (e.g., \code{"BLOSUM62"}), and the negative of the built-in substitution matrix is used to obtain distances rather than similarities.  This is useful for obtaining pairwise substitution scores according to a cost matrix.
}
\value{
If \code{type} is \code{"matrix"}, a symmetric matrix where each element is the distance between the sequences referenced by the respective row and column.  The \code{dimnames} of the matrix correspond to the \code{names} of the \code{XStringSet}.

If \code{type} is \code{"dist"}, an object of \code{class} \code{"dist"} that contains the lower triangle of the distance matrix as a vector.  Since the distance matrix is symmetric, storing only one triangle is more memory efficient.
}
\references{
McGuire G., \emph{et al.} (1999) Improved error bounds for genetic distances from DNA sequences. \emph{Biometrics}, \bold{55(4)}, 1064-1070.

Tajima, F. and Nei, M. (1984) Estimation of evolutionary distance between nucleotide sequences. \emph{Molecular Biology and Evolution}, \bold{1(3)}, 269-285.

Tamura, K. and Nei, M. (1993) Estimation of the number of nucleotide substitutions in the control region of mitrochondrial DNA in humans and chimpanzees. \emph{Molecular Biology and Evolution}, \bold{10(3)}, 512-526.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{MODELS}}, \code{\link{Treeline}}
}
\examples{
# example of using the defaults
dna <- DNAStringSet(c("ACTG", "ACCG"))
dna
DistanceMatrix(dna)

# changing the output type to "dist"
d <- DistanceMatrix(dna, type="dist")
d
length(d) # minimal memory space required
m <- as.matrix(d)
length(m) # more memory space required

# supplying an AAStringSet
aa <- AAStringSet(c("ASYK", "ATYK", "CTWN"))
aa
DistanceMatrix(aa)

# correct for multiple substitutions per site
DistanceMatrix(aa, correction="Poisson")

# defaults compare intersection of internal ranges
dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
dna
d <- DistanceMatrix(dna)
d
# d[1,2] is 1 base in 4 = 0.25

# compare union of internal positions, without terminal gaps
dna <- DNAStringSet(c("ANGCT-", "-ACCT-"))
dna
d <- DistanceMatrix(dna,
                    includeTerminalGaps=TRUE,
                    penalizeGapLetterMatches=TRUE)
d
# d[1,2] is now 2 bases in 5 = 0.40

# gap ("-") and unknown (".") characters are interchangeable
dna <- DNAStringSet(c("ANGCT.", ".ACCT-"))
dna
d <- DistanceMatrix(dna,
                    includeTerminalGaps=TRUE,
                    penalizeGapLetterMatches=TRUE)
d
# d[1,2] is still 2 bases in 5 = 0.40

# compare different methods for calculating distance
dna <- DNAStringSet(c("--ACTG", "TGAGT-"))
dna
DistanceMatrix(dna, method="overlap") # 1/3
DistanceMatrix(dna, method="shortest",
               includeTerminalGaps=FALSE) # 1/3
DistanceMatrix(dna, method="shortest",
               includeTerminalGaps=TRUE,
               penalizeGapLetterMatches=TRUE) # 2/4
DistanceMatrix(dna, method="shortest",
               includeTerminalGaps=TRUE) # 1/3
DistanceMatrix(dna, method="longest",
               includeTerminalGaps=FALSE) # 1/3
DistanceMatrix(dna, method="longest",
               includeTerminalGaps=TRUE,
               penalizeGapLetterMatches=TRUE) # 3/5
DistanceMatrix(dna, method="longest",
               includeTerminalGaps=TRUE) # 1/3
DistanceMatrix(dna, method="overlap",
               minCoverage=1) # NA (insufficient overlap)
DistanceMatrix(dna, method="overlap",
               minCoverage=0.75) # 3/4 sites covered in shorter
DistanceMatrix(dna, method="overlap",
               minCoverage=-0.75) # 3/5 sites covered in longer

# neither internal nor external gap/gap matches are considered
dna <- DNAStringSet(c("--A-CTA", "-AG-C--"))
dna
DistanceMatrix(dna) # 1/2
DistanceMatrix(dna,
               includeTerminalGaps=TRUE,
               penalizeGapLetterMatches=TRUE) # 4/5
DistanceMatrix(dna,
               includeTerminalGaps=TRUE,
               penalizeGapLetterMatches=TRUE,
               method="shortest") # 2/3
DistanceMatrix(dna,
               includeTerminalGaps=TRUE,
               penalizeGapLetterMatches=TRUE,
               method="longest") # 3/4

# examples using real sequences
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
ali <- AlignTranslation(dna, type="both", verbose=FALSE)
DNA <- ali[[1L]]
AA <- ali[[2L]]

# Hamming distances with different gap penalties
d_NT_noGap <- DistanceMatrix(DNA, type="dist") # penalizeGapLetter=FALSE
d_NT_allGap <- DistanceMatrix(DNA, type="dist", penalizeGapLetter=TRUE)
d_NT_oneGap <- DistanceMatrix(DNA, type="dist", penalizeGapLetter=NA)
pairs(data.frame(d_NT_noGap, d_NT_allGap, d_NT_oneGap), pch=46, col="#00000011")

# correcting distances for multiple substitutions per site
Q <- 1 - diag(4)
rownames(Q) <- colnames(Q) <- DNA_BASES
freqs <- setNames(rep(0.25, 4), DNA_BASES)
d_NT_JC <- DistanceMatrix(DNA, type="dist", correction="JC")
d_NT_even <- DistanceMatrix(DNA, type="dist", substitutionMatrix=Q, frequencies=freqs)
pairs(data.frame(d_NT_noGap, d_NT_JC, d_NT_even), pch=46, col="#00000011")

# specifying a rate matrix with uneven transition/transversion rates
Q["C", "T"] <- Q["T", "C"] <- Q["A", "G"] <- Q["G", "A"] <- 2
d_NT_uneven <- DistanceMatrix(DNA, type="dist", substitutionMatrix=Q, frequencies=freqs)
pairs(data.frame(d_NT_noGap, d_NT_uneven, d_NT_even), pch=46, col="#00000011")

# comparing amino acid and nucleotide distances
d_AA_noGap <- DistanceMatrix(AA, type="dist") # penalizeGapLetter=FALSE
d_AA_WAG <- DistanceMatrix(AA, type="dist", correction="WAG")
pairs(data.frame(d_NT_noGap, d_AA_noGap, d_AA_WAG), pch=46, col="#00000011")

# using a model with letter frequencies derived from the input sequences
d_NT_F81_F <- DistanceMatrix(DNA, type="dist", correction="F81+F")
d_AA_WAG_F <- DistanceMatrix(AA, type="dist", correction="WAG+F")
pairs(data.frame(d_NT_F81_F, d_AA_WAG, d_AA_WAG_F), pch=46, col="#00000011")

# choosing a model with different transition/transversion rates
d_NT_K80 <- DistanceMatrix(DNA, type="dist", correction="K80")
d_NT_TN93_F <- DistanceMatrix(DNA, type="dist", correction="TN93+F")
pairs(data.frame(d_NT_F81_F, d_NT_K80, d_NT_TN93_F), pch=46, col="#00000011")

# incorporating the indel component of distance
d_NT_JC_indel <- DistanceMatrix(DNA, type="dist", correction="JC", penalize=TRUE)
d_NT_JC_indelblock <- DistanceMatrix(DNA, type="dist", correction="JC", penalize=NA)
pairs(data.frame(d_NT_JC, d_NT_JC_indel, d_NT_JC_indelblock), pch=46, col="#00000011")

# specifying a substitutionMatrix (without frequencies)
d_AA_BLOSUM <- DistanceMatrix(AA, type="dist", substitutionMatrix="BLOSUM62")
d_AA_PFASUM <- DistanceMatrix(AA, type="dist", substitutionMatrix="PFASUM50")
pairs(data.frame(d_AA_WAG, d_AA_BLOSUM, d_AA_PFASUM), pch=46, col="#00000011")
}

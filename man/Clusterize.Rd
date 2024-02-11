\name{Clusterize}
\alias{Clusterize}
\title{
Cluster Sequences By Distance
}
\description{
Groups the sequences into approximate clusters of similarity.
}
\usage{
Clusterize(myXStringSet,
           cutoff = 0,
           method = "overlap",
           includeTerminalGaps = FALSE,
           penalizeGapLetterMatches = NA,
           minCoverage = 0.5,
           maxPhase1 = 2e4,
           maxPhase2 = 2e3,
           maxPhase3 = 2e3,
           maxAlignments = 200,
           rareKmers = 50,
           probability = 0.99,
           invertCenters = FALSE,
           singleLinkage = FALSE,
           alphabet = AA_REDUCED[[186]],
           processors = 1,
           verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
The (unaligned) \code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet} to cluster.
}
  \item{cutoff}{
A vector of maximum distances (approximately) separating sequences in the same cluster (i.e., 1 - similarities).  Multiple cutoffs may be provided in ascending or descending order.  (See details section below.)
}
  \item{method}{
Character string determining the region in which distance is calculated.  This should be (an unambiguous abbreviation of) one of \code{"overlap"}, \code{"shortest"}, or \code{"longest"}.  The default \code{method} (\code{"overlap"}) calculates distance from the overlapping region between terminal gaps when \code{includeTerminalGaps} is \code{FALSE} and the entire alignment otherwise.  Setting \code{method} to \code{"shortest"} or \code{"longest"} will use the region between the start and end of the shortest or longest sequence, respectively, for each pairwise distance.  The \code{method} is only applicable when \code{includeTerminalGaps} is \code{TRUE}.
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" characters) in the pairwise alignments into the calculation of distance.
}
  \item{penalizeGapLetterMatches}{
Logical specifying whether or not to consider gap-to-letter matches as mismatches.  If \code{FALSE}, then gap-to-letter matches are not included in the total length used to calculate distance, and if \code{TRUE} all gaps-to-letter pairs are considered mismatches.  The default (\code{NA}) is to penalize gap-to-letter mismatches once per insertion or deletion, which treats runs of gaps (i.e., indels) as equivalent to a single mismatch.
}
  \item{minCoverage}{
Numeric giving the minimum fraction of sequence positions that must be overlapping for a sequence to be clustered with the cluster representative.  If positive then coverage is calculated relative to the sequence being clustered (i.e., the shorter sequence).  If negative then coverage is computed relative to the cluster representative (i.e., the longest sequence in the cluster).
}
  \item{maxPhase1}{
An integer specifying the maximum number of related sequences to consider in the initial partitioning of the sequences.
}
  \item{maxPhase2}{
An integer giving the maximum number of replicates to perform when sorting sequences based on their k-mer similarity.
}
  \item{maxPhase3}{
An integer determining the number of comparisons per sequence to perform when attempting to find cluster centers.
}
  \item{maxAlignments}{
An integer designating the maximum number of alignments to perform when attempting to assign a sequence to an existing cluster.
}
  \item{rareKmers}{
An integer setting the number of rare k-mers to record per sequence. Larger values require more memory but may improve accuracy with diminishing returns.
}
  \item{probability}{
Numeric between 0 and 1 (exclusive) defining the approximate probability of clustering sequences that are exactly \code{cutoff} distant. Typically near, but always less than, \code{1}. Lower values result in faster clustering at the expense of effectiveness.
}
  \item{invertCenters}{
Logical controlling whether the cluster center is inverted (i.e., multiplied by \code{-1}), which allows the centers to be determined from the results.  The default (\code{FALSE}) only returns positive cluster numbers.  If \code{TRUE}, the center sequence(s) of each cluster are negative.
}
  \item{singleLinkage}{
Logical specifying whether to perform single-linkage clustering. The default (\code{FALSE}) only establishes linkage to the cluster center. Single-linkage clustering creates broader clusters that may better correspond to natural groups depending on the application.
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between sequences.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{Clusterize} groups the input sequences into approximate clusters using a heuristic algorithm with linear time and memory complexity.  In phase 1, the sequences are partitioned into groups of similarity.  In phase 2, the sequences are ordered by k-mer similarity by relatedness sorting.  In phase 3, the sequences are iteratively clustered in this order by their similarity to surrounding sequences in the sorting.  That is, the first sequence becomes the representative of cluster #1.  If the second sequence is within \code{cutoff} distance then it is added to the cluster, otherwise it becomes a new cluster representative.  The remaining sequences are matched to cluster representatives in a similar fashion until all sequences belong to a cluster.  In the majority of cases, this process results in clusters with members separated by less than \code{cutoff} distance, and all cluster members must be within \code{cutoff} distance of their cluster representative.

The calculation of distance can be controlled in multiple ways, with each parameterization of distance having advantages and disadvantages.  By default, distance is the fraction of positions that are different, including gaps, within the overlapping region in a pairwise alignment.  The defaults will handle partial-length sequences well, but also cluster sequences with high similarity between their opposite ends.  For this reason, it is important to set \code{minimumCoverage} such that distances are based off of considerable overlap between sequences in the pairwise alignment.  The default (\code{0.5}) requires sequences to share at least 50\% of their positions with the cluster representative.  This distance parameterization works well, but there are reasonable alternatives.

If \code{penalizeGapLetterMatches} is \code{FALSE}, the distance will exclude gap regions.  If \code{includeTerminalGaps} is \code{TRUE}, the calculation of distance will use the entire (global) alignment.  If \code{method} is \code{"shortest"} and \code{includeTerminalGaps} is set to \code{TRUE}, then the distance is calculated for the region encompassed by the shorter sequence in each pair, which is the common definition of distance used by other clustering programs.  This common definition of distance will sometimes separate partly overlapping sequences, which is why it is not the default.  Another option is to set \code{minCoverage} to a negative fraction (e.g., \code{-0.5}).  This requires sequences to have substantial overlap with the cluster representative, which is the longest sequence in the cluster.

The algorithm requires time proportional to the number of input sequences in \code{myXStringSet}.  The phase 1, up to \code{maxPhase1} sequences sharing a k-mer are tabulated while partitioning each sequence.  In phase 2, the sequences are compared with up to \code{maxPhase2} passes that each take linear time.  Ordering of the sequences is performed in linear time using radix sorting.  In phase 3, each sequence is compared with up to \code{maxPhase3} previous cluster representatives of sequences sharing \code{rareKmers} or nearby sequences in the relatedness ordering.  This is possible because the sequences are sorted by relatedness, such that more recent cluster representatives are more similar.  Hence, the complete algorithm scales in linear time asymptotically and returns clusters of sequences within \code{cutoff} distance of their center sequence.

Multiple cutoffs can be provided in sorted order, which saves time because phases 1 and 2 only need to be performed once.  If the \code{cutoff}s are provided in \emph{descending} order then clustering at each new value of \code{cutoff} is continued within the prior \code{cutoff}'s clusters.  In this way clusters at lower values of \code{cutoff} are completely contained within their ``umbrella'' clusters at higher values of \code{cutoff}.  This slightly accelerates the clustering process, because each subsequent group is only clustered within the previous group.  If multiple cutoffs are provided in \emph{ascending} order then clustering at each level of \code{cutoff} is independent of the prior level.

Note, the clustering algorithm is stochastic.  Hence, clusters can vary from run-to-run unless the random number seed is set for repeatability (i.e., with \code{set.seed}).  Also, \code{invertCenters} can be used to determine the center sequence of each cluster from the output.  Since identical sequences will always be assigned the same cluster numbers, it is possible for more than one input sequence in \code{myXStringSet} to be assigned as the center of a cluster if they are identical.
}
\value{
A data.frame is returned with dimensions \eqn{N*M}, where each one of \eqn{N} sequences is assigned to a cluster at the \eqn{M}-level of cutoff.  The row.names of the data.frame correspond to the \emph{names} of \code{myXStingSet}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AA_REDUCED}}, \code{\link{DistanceMatrix}}, \code{\link{TreeLine}}
}
\examples{
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
aa <- translate(dna)

# typical usage (e.g., clustering at >= 90 percent similarity)
clusters <- Clusterize(aa, cutoff=0.1) # set processors = NULL for max speed
head(clusters)

# typical usage (e.g., obtaining cluster representatives)
clusters <- Clusterize(aa, cutoff=0.1, invertCenters=TRUE)
aa[clusters[[1]] < 0]

# cluster each cutoff within the previous cluster (nested groups)
clusters <- Clusterize(aa, cutoff=seq(0.7, 0, -0.1))
head(clusters)
apply(clusters, 2, max) # number of clusters per cutoff

# cluster each cutoff independently (possibly fewer clusters per cutoff)
clusters <- Clusterize(aa, cutoff=seq(0, 0.7, 0.1))
head(clusters)
apply(clusters, 2, max) # number of clusters per cutoff

# make cluster center(s) negative for tracking
clusters <- Clusterize(aa, cutoff=0.5, invertCenters=TRUE)
head(clusters)
clusters[clusters$cluster < 0,, drop=FALSE]
unique(aa[clusters$cluster < 0]) # unique cluster centers
apply(clusters, 2, function(x) max(abs(x))) # number of clusters

# cluster nucleotide sequences
clusters <- Clusterize(dna, cutoff=0.2, invertCenters=TRUE)
head(clusters)
apply(clusters, 2, function(x) max(abs(x))) # number of clusters
}

\name{Treeline}
\alias{Treeline}
\alias{TreeLine}
\title{
Construct a Phylogenetic Tree
}
\description{
Builds a phylogenetic tree from a set of sequences or distance matrix.
}
\usage{
Treeline(myXStringSet = NULL,
         myDistMatrix = NULL,
         method = "ME",
         type = "dendrogram",
         cutoff = -Inf,
         showPlot = FALSE,
         standardDeviation = 0.2,
         fracRandomNNIs = 0.4,
         goalPercent = NA,
         tolerance = 5e-5,
         minIterations = 40,
         maxIterations = 400,
         maxTime = Inf,
         root = 0,
         collapse = -1,
         reconstruct = FALSE,
         costMatrix = NULL,
         model = MODELS,
         informationCriterion = "AICc",
         quadrature = FALSE,
         processors = 1,
         verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet}.  Required if (a) not supplying \code{myDistMatrix}, (b) \code{method} is \code{"ML"}, (c) \code{method} is \code{"MP"}, or (d) \code{reconstruct} is not \code{FALSE}.
}
  \item{myDistMatrix}{
A symmetric \eqn{N} x \eqn{N} distance matrix with the values of dissimilarity between \eqn{N} sequences or an object of class \code{'dist'}.  If \code{NULL} (the default), \code{myDistMatrix} will be automatically determined from \code{myXStringSet} based on the distance between sequences using the \code{DistanceMatrix} function and the \code{model} (if specified).  If not \code{NULL}, \code{myDistMatrix} is used to construct the initial tree, even when \code{method} is \code{"ML"} or \code{"MP"}.  If \code{method} is \code{"ME"}, \code{myDistMatrix} determines the branch lengths and, therefore, overall tree length that is optimized.  Missing values (i.e., \code{NA}) in \code{myDistMatrix} are imputed using the ultrametric method described by Makarenkov and Lapointe (2004).
}
  \item{method}{
The phylogenetic method to be used.  This should be (an abbreviation of) one of \code{"ME"} for (balanced) minimum evolution (the default), \code{"ML"} for maximum likelihood, \code{"MP"} for maximum parsimony, \code{"NJ"} for neighbor joining, \code{"complete"} for complete-linkage, \code{"single"} for single-linkage, \code{"UPGMA"} for ``unweighted pair group method with arithmetic mean'', or \code{"WPGMA"} for ``weighted pair group method with arithmetic mean''.  (See details section below.)
}
  \item{type}{
Character string indicating the type of output desired.  This should be (an abbreviation of) one of \code{"dendrogram"} (the default), \code{"clusters"}, or \code{"both"}.  (See value section below.)
}
  \item{cutoff}{
A vector with the maximum edge length separating the sequences in the same cluster.  A negative value (the default) will prevent clustering.  Multiple cutoffs may be provided in ascending or descending order.  (See details section below.)
}
  \item{showPlot}{
Logical specifying whether or not to plot the resulting dendrogram.
}
  \item{standardDeviation}{
Numeric determining the extent to which cophenetic distances are perturbed prior to constructing the initial candidate tree in each iteration.  Only applicable if \code{method} is \code{"ME"}, \code{"ML"}, or \code{"MP"}.
}
  \item{fracRandomNNIs}{
Numeric giving the fraction of stochastic nearest neighbor interchanges to perform when perturbing the tree to develop a new candidate tree in each iteration after the first.  Only applicable if \code{method} is \code{"ME"}, \code{"ML"}, or \code{"MP"}.
}
  \item{goalPercent}{
Numeric providing the target percent difference in score relative to the best observed tree after optimizing the candidate tree through nearest neighbor interchanges (NNIs).  If \code{goalPercent} is not \code{NA} (the default), \code{fracRandomNNIs} is iteratively adjusted to reach \code{goalPercent} relative score.  Only applicable if \code{method} is \code{"ME"}, \code{"ML"}, or \code{"MP"}.
}
  \item{tolerance}{
Numeric determining the relative convergence tolerance.  Iterating will cease when the relative score has changed by less than \code{tolerance} for \code{minIterations} and is expected to change by less than \code{tolerance} per iteration.  Only applicable if \code{method} is \code{"ME"}, \code{"ML"}, or \code{"MP"}.
}
  \item{minIterations}{
Integer indicating the minimum number of iterations of optimization to perform.  Lower values will result in faster convergence but may insufficiently explore tree space.  Only applicable if \code{method} is \code{"ME"}, \code{"ML"}, or \code{"MP"}.
}
  \item{maxIterations}{
Integer indicating the maximum number iterations of optimization to perform.  More iterations will potentially better search tree space at the expense of added runtime.  Only applicable if \code{method} is \code{"ME"}, \code{"ML"}, or \code{"MP"}.
}
  \item{maxTime}{
Numeric giving the maximum number of hours the algorithm is allowed to run before returning a result.  Once \code{maxTime} is reached, the algorithm will proceed at the next available opportunity, even if \code{minIterations} is unmet.  Only applicable if \code{method} is \code{"ME"}, \code{"ML"}, or \code{"MP"}.  Note that setting a time limit may prevent reproducibility when using a random number seed.
}
  \item{root}{
Integer specifying the index of the outgroup sequence or \code{0} (the default) to midpoint root the dendrogram.
}
  \item{collapse}{
Numeric controlling which internal edges of the tree are removed by collapsing their nodes.  If \code{collapse} is zero then nodes at the same height will be collapsed to a single node, resulting in a multifurcating tree.  When \code{collapse} is greater than zero, nodes that are within \code{collapse} difference in height are made into a single node.  A value of \code{collapse} less than zero (the default) will ensure that the \code{dendrogram} is purely bifurcating.  Note that \code{collapse} has no effect on cluster numbers or \code{cutoff}.
}
  \item{reconstruct}{
Logical or numeric determining whether to perform ancestral state reconstruction when \code{myXStringSet} is specified.  If \code{TRUE}, ancestral character states are determined at internal nodes of the \code{dendrogram} and provided as a \code{"state"} attribute.  Ancestral states are determined as the most parsimonious state according to the \code{"costMatrix"}, unless \code{method} is \code{"ML"}, in which case ancestral states are determined as those with the highest (marginal) likelihood.  A numeric value between zero and one (exclusive) can be provided when \code{method} is \code{"ML"} to require that fraction of the state's likelihood to be greater than that of all alternative states, otherwise a more ambiguous degeneracy code is used.  Only applicable if \code{type} is \code{"dendrogram"} (the default) or \code{"both"}.
}
  \item{costMatrix}{
Either \code{NULL} (the default) or a symmetric matrix setting the penalties used in Sankoff parsimony.  The default (\code{NULL}) will apply a cost of \code{1} for standard character state changes and \code{0} otherwise (i.e., equivalent to Fitch parsimony).  If a \code{matrix} then the states are taken from its row or column names.  Only applicable if \code{method} is \code{"MP"}, or \code{reconstruct} is \code{TRUE} and \code{method} is not \code{"ML"}.
}
  \item{model}{
One or more of the available \code{MODELS} of evolution provided as a character vector or list with components `Protein' and/or `Nucleotide' if \code{method} is \code{"ML"}.  Automatic model selection based on the \code{informationCriterion} will be performed for \code{"ML"} if more than one \code{model} is provided.  For \code{"ME"}, if \code{myDistMatrix} is \code{NULL}, a single model of evolution can be specified for generating the distance matrix from \code{myXStringSet}.
}
  \item{informationCriterion}{
Character string specifying which information criterion to use in automatic model selection when \code{method} is \code{"ML"}.  Must be (an abbreviation of) either \code{"AICc"} or \code{"BIC"}.  The best \code{model} is automatically chosen based on the \code{informationCriterion} calculated from the likelihood and the sample size (defined as the average number of sites for sequences in \code{myXStringSet}).  Only applicable if \code{method} is \code{"ML"}.
}
  \item{quadrature}{
Logical determining whether to use the Laguerre quadrature or equal-sized bins when discretizing the rate distribution across sites when \code{method} is \code{"ML"}.  The default (\code{FALSE}) is to use equal-sized bins for direct comparison among likelihoods computed by other programs, although the Laguerre quadrature offers a modest improvement in accuracy for the same number of rates (Felsenstein, 2001).
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.  Note, the number of processors in some steps is automatically selected between \code{1} and \code{processors} to optimize performance.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{Treeline} builds a phylogenetic tree using either \code{myXStringSet} and/or \code{myDistMatrix}.  The output is either a \code{dendrogram} and/or \code{data.frame} containing cluster numbers.

Multiple \code{method}s of tree building are supported:

(1) Minimum evolution: \code{ME} iteratively minimizes the (balanced) tree length according to a distance matrix, as described by Pauplin (2000).  According to Gonnet (2012) and Spirin \emph{et al.} (2024), the \code{ME} criterion performs best overall on benchmarks constructed from real genes and is therefore the default \code{method}.  \code{ME} is also the fastest of the optimized \code{method}s.  Branch lengths are determined from \code{myDistMatrix} or, if missing (\code{NULL}), calculated from the hamming distances between sequences in \code{myXStringSet} unless a \code{model} is specified.

(2) Maximum likelihood: \code{ML} iteratively maximizes the likelihood of the tree and any free \code{model} parameters given aligned sequences (\code{myXStringSet}).  One or more \code{MODELS} of sequence evolution must be specified, of which the best \code{model} is automatically selected based on the \code{informationCriterion}.  The \code{ML} criterion performs well on benchmarks constructed from real genes when given sufficiently long and variable alignments.  Branch lengths are in units of substitution per site.

(3) Maximum parsimony: \code{MP} iteratively maximizes the parsimony of a tree given aligned sequences (\code{myXStringSet}) and a \code{costMatrix}.  The default cost matrix is binary, corresponding to Fitch (1971) parsimony.  The \code{costMatrix} often has a large influence over the result, and more biologically reasonable cost matrices will likely improve the resulting tree.  See the examples below for non-uniform cost matrices.  Here, branch lengths are in units of the average cost per site.

(4) Neighbor-joining: \code{NJ} uses the Neighbor-Joining method proposed by Saitou and Nei (1987), which creates an approximate minimum evolution tree from a distance matrix (\code{myDistMatrix}).  The \code{NJ} criterion is a greedy heuristic approach to minimum evolution, and can be considered a less accurate alternative to \code{ME}.  Notably, \code{NJ} is faster than \code{ME} for moderate numbers of sequences, since it only computes a single tree, but is slower than \code{ME} for large numbers of sequences due to its cubic time complexity.

(5) Ultrametric: The method \code{complete} assigns clusters using complete-linkage so that sequences in the same cluster are no more than \code{cutoff} distance apart.  The method \code{single} assigns clusters using single-linkage so that sequences in the same cluster are within \code{cutoff} of at least one other sequence in the same cluster.  \code{UPGMA} and \code{WPGMA} assign clusters using average-linkage which is a compromise between the sensitivity of complete-linkage clustering to outliers and the tendency of single-linkage clustering to connect distant relatives that are not closely related.  \code{UPGMA} produces an unweighted tree, where each leaf contributes equally to the average edge lengths, whereas \code{WPGMA} produces a weighted result.

For \code{method}s \code{"ME"}, \code{"ML"}, and \code{"MP"}, candidate trees are iteratively optimized through rounds of nearest neighbor interchanges (``climbs'') followed by fusion of remaining differences to the best observed tree (``grafts'').  Candidate trees are generated with a heuristic variant of neighbor joining after perturbing the best observed tree's cophenetic distance matrix, followed by stochastic NNIs if \code{fracRandomNNIs} is greater than \code{0}.  The value of \code{fracRandomNNIs} is adaptively varied to reach \code{goalPercent} relative score on average after optimization via climbs, unless \code{goalPercent} is \code{NA} (the default).  This process results in gradual improvement until reaching a best tree, which is returned if insufficient score improvement is made for \code{minIterations}, unless \code{maxIterations} or \code{maxTime} is reached.

The objective of optimization is to find the global optimum, but it is necessary to set a random seed for exact reproducibility since optimization is stochastic and may only find a local optimum.  Setting \code{maxTime} may prevent reproducibility if iteration terminates prior to reaching convergence or \code{maxIterations}.  Also, not all math operations (e.g., logarithm) have standard implementations across platforms, so reproducibility is not guaranteed on different machines.  The results are not necessarily significant even if they are reproducible, and it possible to use bootstrapping to gauge support for different partitions.  When \code{method} is \code{"ML"}, aBayes support values are provided that are a reasonable representation of support.

The returned dendrogram has information stored in its attributes, which can be accessed with the \code{attributes} or \code{attr} functions.  For maximum likelihood trees, the edges have a \code{"probability"} attribute representing the aBayes support probability (Anisimova, 2011).  If \code{reconstruct} is not \code{FALSE}, each edge of the tree will have a \code{"state"} attribute representing the node's predicted ancestral state.  Also, when \code{reconstruct} is not \code{FALSE}, maximum likelihood trees will provide the likelihoods at each site and other \code{method}s will provide a state transition matrix.

When non-negative \code{cutoff}(s) are supplied, \code{Treeline} will assign clusters based on edge lengths in the tree.  Multiple cutoffs may be provided in sorted order.  If the \code{cutoff}s are provided in \emph{descending} order then clustering at each new value of \code{cutoff} is continued within the prior \code{cutoff}'s clusters.  In this way clusters at lower values of \code{cutoff} are completely contained within their umbrella clusters at higher values of \code{cutoff}.  This is useful for defining taxonomy, where groups need to be hierarchically nested.  If multiple cutoffs are provided in \emph{ascending} order then clustering at each level of \code{cutoff} is independent of the prior level.
}
\value{
If \code{type} is \code{"dendrogram"} (the default), then a tree of class \code{dendrogram} is returned with attributes including pertinent information.
If \code{type} is \code{"clusters"}, then a \code{data.frame} is returned with dimensions \eqn{N*M}, where each one of \eqn{N} sequences is assigned to a cluster at the \eqn{M}-level of \code{cutoff}.
If \code{type} is \code{"both"} then a list is returned containing both the \code{"clusters"} and \code{"dendrogram"} outputs.
}
\note{
Cophenetic distance between leaves of the \code{dendrogram} is often defined as the sum of branch lengths separating the leaves, also known as the patristic distance.  This is the typical phylogenetic interpretation but different than that for dendrograms produced by \code{hclust} where leaves are merged at a height equal to their cophenetic distance.  Hence, always use \code{Cophenetic} (rather than \code{cophenetic}) to compute cophenetic distances if the length between leaves is desired.
}
\references{
Anisimova M., \emph{et al.} (2011) Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. \emph{Syst Biol.}, \bold{60(5)}, 685-99.

Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a maximum likelihood approach. \emph{Journal of Molecular Evolution}, \bold{17(6)}, 368-376.

Felsenstein J. (2001) Taking variation of evolutionary rates between sites into account in inferring phylogenies. \emph{Journal of molecular evolution}, \bold{53}, 447-455.

Fitch, W. M. (1971) Toward defining the course of evolution: minimum change for a specified tree topology. \emph{Systematic Zoology}, \bold{20}:406-416.

Gonnet, G. H. (2012) Surprising results on phylogenetic tree building methods based on molecular sequences. \emph{BMC Bioinformatics}, \bold{13}, 148.

Makarenkov V., and Lapointe, F. (2004) A weighted least-squares approach for inferring phylogenies from incomplete distance matrices. \emph{Bioinformatics}, \bold{20(13)}, 2113-2121.

Pauplin, Y. (2000) Direct Calculation of a Tree Length Using a Distance Matrix. \emph{J Mol Evol}, \bold{51(1)}, 41-47.

Saitou, N. and Nei, M. (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. \emph{Molecular Biology and Evolution}, \bold{4(4)}, 406-425.

Sankoff, D. (1975) Minimal mutation trees of sequences. \emph{SIAM Journal of Applied Math}, \bold{28}.

Spirin, S., \emph{et al.} (2024) PhyloBench: A Benchmark for Evaluating Phylogenetic Programs. \emph{Molecular Biology and Evolution}, \bold{41(6)}, 1-11.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{Clusterize}}, \code{\link{Cophenetic}}, \code{\link{DistanceMatrix}}, \code{\link{MapCharacters}}, \code{\link{MODELS}}, \code{\link{ReadDendrogram}}, \code{\link{WriteDendrogram}}
}
\examples{
# using the matrix from the original paper by Saitou and Nei (1987)
m <- matrix(0,8,8) # only the lower triangle is used
m[2:8,1] <- c(7, 8, 11, 13, 16, 13, 17)
m[3:8,2] <- c(5, 8, 10, 13, 10, 14)
m[4:8,3] <- c(5, 7, 10, 7, 11)
m[5:8,4] <- c(8, 11, 8, 12)
m[6:8,5] <- c(5, 6, 10)
m[7:8,6] <- c(9, 13)
m[8,7] <- 8

# returns an object of class "dendrogram"
tree <- Treeline(myDistMatrix=m, cutoff=10, method="NJ", showPlot=TRUE)

# example of specifying multiple cutoffs
clusters <- Treeline(myDistMatrix=m, method="UPGMA", type="clusters", cutoff=c(2,6,10,20))
head(clusters)

# example of creating a complete-linkage tree from an alignment
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
alignments <- AlignTranslation(dna, verbose=FALSE, type="both")
dna <- alignments[[1]]
aa <- alignments[[2]]
dna # input alignment
d <- DistanceMatrix(dna, type="dist") # returns an object of class 'dist'
complete <- Treeline(myDistMatrix=d, method="complete", cutoff=0.05, showPlot=TRUE)

# example of minimum evolution (ME) tree optimization (the default)
treeME <- Treeline(dna, processors=1L) # the recommended way to build trees
plot(treeME)

# compare the distance matrix to the cophenetic (patristic) distance matrix
plot(d, Cophenetic(treeME),
	xlab="Pairwise distance", ylab="Patristic distance",
	asp=1, pch=46, col="#00000033")
abline(a=0, b=1)

# example of maximum parsimony (MP) tree optimization
costs <- matrix(c(0, 2, 1, 2, 2, 0, 2, 1, 1, 2, 0, 2, 2, 1, 2, 0), 4)
dimnames(costs) <- list(DNA_BASES, DNA_BASES)
costs # the cost matrix
treeMP <- Treeline(dna, method="MP", reconstruct=TRUE, costMatrix=costs, maxTime=0.001)

# display ancestral states on each edge
start <- 50 # starting position in alignment
end <- 52 # ending position in alignment
tree <- dendrapply(treeMP,
	function(x) {
		attr(x, "edgetext") <- substring(attr(x, "state"),
			start,
			end)
		x
	})
plot(tree,
	edgePar=list(p.col=NA,
		p.border=NA,
		t.col=c("#AA3355", "#33FFFF")),
	edge.root=FALSE,
	leaflab="none")

# example of maximum likelihood tree optimization
treeML <- Treeline(head(dna, 10), method="ML", reconstruct=TRUE, maxTime=0.001)

# example of accessing and using the attributes
attributes(treeML) # show all attributes at a node
plot(dendrapply(treeML,
	function(x) {
		s <- attr(x, "probability")
		if (!is.null(s))
			attr(x, "edgetext") <- formatC(as.numeric(s), digits=2, format="f")
		attr(x, "edgePar") <- list(p.col=NA, p.border=NA, t.col="#CC55AA", t.cex=0.7)
		x
	}),
	horiz=TRUE)

# construct a tanglegram comparing two trees back-to-back
layout(matrix(1:2, ncol=2))
tree1 <- treeME # tree on left
tree2 <- treeMP # tree on right
tree1 <- reorder(tree1, unlist(tree1))
tree2 <- reorder(tree2, unlist(tree2))
layout(matrix(1:2, nrow=1))
par(mai=c(0.5, 0, 0.5, 0.5)) # add space on right
plot(tree1,
	main="First tree",
	horiz=TRUE, leaflab="none")
par(mai=c(0.5, 0.5, 0.5, 0)) # add space on left
plot(tree2,
	main="Second tree",
	horiz=TRUE, leaflab="none",
	xlim=c(0, attr(tree2, "height")))
segments(par("usr")[1] - 1.05*diff(grconvertX(0:1, 'inches', 'user')),
	match(unlist(tree2), unlist(tree1)),
	-0.05*diff(grconvertX(0:1, 'inches', 'user')),
	seq_len(attr(tree2, "members")),
	xpd=NA)

# example of supplying an amino acid cost matrix for MP
codons <- getGeneticCode("1") # the standard genetic code
codons <- tapply(names(codons), codons, c)
costs <- outer(codons,
	codons,
	function(x, y)
		mapply(function(x, y) mean(adist(x, y)), x, y))
costs # minimum number of nucleotide changes to switch amino acids
treeMP_AA1 <- Treeline(aa, method="MP", costMatrix=costs, maxTime=0.001)
plot(Cophenetic(treeMP_AA1), Cophenetic(treeMP))

# alternative approach to obtaining an amino acid cost matrix
data(BLOSUM)
subM <- BLOSUM[AA_STANDARD, AA_STANDARD, "62"] # BLOSUM62 matrix
subM <- diag(subM) - subM # standardize to diagonal = 0
subM <- (subM + t(subM))/2 # make symmetric
treeMP_AA2 <- Treeline(aa, method="MP", costMatrix=subM, maxTime=0.001)
plot(Cophenetic(treeMP_AA1), Cophenetic(treeMP_AA2)) # different lengths
}

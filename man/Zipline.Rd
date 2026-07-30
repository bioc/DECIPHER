\name{Zipline}
\alias{Zipline}
\title{
Summarize Multiple Phylogenetic Trees
}
\description{
Generates a summary dendrogram or distance matrix from a set of phylogenetic trees using linear regression of transformed distances.
}
\usage{
Zipline(x,
        type = "dendrogram",
        distance = "length",
        weights = 1,
        bootstraps = 0,
        power = 0.5,
        verbose = TRUE,
        ...)
}
\arguments{
  \item{x}{
A list containing one \code{dendrogram} object per list element.
}
  \item{type}{
Character string indicating the type of output desired.  This should be one of \code{"dendrogram"} (the default), \code{"dist"}, or \code{"both"}.  (See value section below.)
}
  \item{distance}{
The type of distance measure(s) to use.  This should be \code{"length"} (the default), \code{"edges"}, and/or named edge attributes.  (See the \code{\link{Cophenetic}} manual page for more information.)
}
  \item{weights}{
Either a single number (the default) or a numeric vector with one weight corresponding to each \code{dendrogram} in \code{x}.  Specifying \code{weights} as proportional to trustworthiness (e.g., number of aligned sites originally used to construct each dendrogram) may modestly improve summary tree accuracy.
}
  \item{bootstraps}{
The number of bootstrap replicates desired or zero (the default) for no replicates.  Bootstrapping entails resampling (with replacement) the set of trees in \code{x} and recording the fraction of times an edge partitioning the same sets of leaves was observed across replicates.
}
  \item{power}{
A single number determining the transformed space in which regression is performed.  The default (\code{0.5}) will perform regression in square root space, while setting \code{power} to \code{1} would perform traditional regression.  Smaller values will diminish the influence of large outliers at the expense of losing some of their information.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Additional arguments to be passed directly to \code{\link{Treeline}} for generating the output \code{dendrogram}.
}
}
\details{
The genealogies of gene trees are often discordant with each other and the species tree due to biological processes and technical difficulties.  Given sufficient gene trees, it is feasible to reconstruct the species tree that best summarizes all gene trees.  This is the purpose of \code{Zipline}, which is named for its objective of zippy (fast) and accurate summary tree inference via linear regression to connect different trees.

\code{Zipline} works by first regressing gene tree distances onto the average inter-leaf distances to calculate the relative rate of evolution for each gene tree.  Next, the function regresses inter-leaf distances onto the per tree rate to derive the best estimate of distances between leaves.  By default, this distance matrix is used to construct a summary tree with \code{Treeline}.

Regressions are performed in square root space (by default) to mitigate the effect of outliers on estimated slopes.  The use of regression makes the \code{Zipline} algorithm robust to error and effortlessly capable of handling missing data (e.g., absent genes).  \code{Zipline} also handles duplicated labels in the same tree (e.g., multi-copy genes) by collapsing to the minimum distance.  In this manner, all input trees do not need to have the same set of (leaf) labels.  Regression also permits the incorporation of \code{weights} to bias to the distance estimate to towards more trustworthy trees in \code{x}.

Multiple measures of distance are supported.  Historically the internode distance (\code{"edges"}) was commonly used for species tree inference, but recently path length (\code{"length"}) and other measures have been shown to perform better.  Dendrograms that include branch supports (e.g., \code{"edgetext"}) may offer higher accuracy and make hybrid distance measures possible.  (See examples section below.)

Output distances (and edge lengths) are in the same units as the distance measure.  That is, if \code{"length"} is specified then output distances are in the units of branch length (e.g., substitutions per site).  Other distance measures, such as the sum of edge support values, are less easily interpretable although they remain in the expected units.  Note that all trees in \code{x} are expected to share the same units.

Bootstrapping is incorporated at the level of sampling gene trees within \code{x} (i.e., multi-locus bootstrapping).  As with most bootstrap procedures, this requires a large enough sample set (i.e., \code{x}) to adequately represent the data distribution.  Bootstrap supports (as fractions) are provided via ``edgetext'' attributes on internal branches of the tree.
}
\value{
If \code{type} is \code{"dendrogram"} (the default), then a summary tree of class \code{dendrogram} is returned with ``edgetext'' attributes if \code{bootstraps} is positive.
If \code{type} is \code{"dist"}, then a summary distance matrix of class \code{dist} is returned.
If \code{type} is \code{"both"} then a list is returned containing both the \code{"dendrogram"} and \code{"dist"} outputs.
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\seealso{
\code{\link{Cophenetic}}, \code{\link{ReadDendrogram}}, \code{\link{Treeline}}, \code{\link{WriteDendrogram}}
}
\examples{
# import a set of single-copy gene trees from plant genomes
newick <- system.file("extdata", "Plant_gene_trees.newick.gz", package="DECIPHER")
trees <- ReadDendrogram(newick)
length(trees) # number of dendrograms
plot(trees[[1]], main="First gene tree")

# using the default distance (path length)
species_tree <- Zipline(trees)
plot(species_tree, main="Species tree under default settings")

# using the internode distance
species_tree <- Zipline(trees, distance="edges")
plot(species_tree, main="Species tree from internode distances")

# using branch supports as distance
species_tree <- Zipline(trees, distance="edgetext")
plot(species_tree, main="Species tree from bootstrap supports")

# normalizing bootstrap supports to zero to one
trees <- lapply(trees,
	dendrapply,
	function(x) {
		b <- as.numeric(attr(x, "edgetext"))
		if (length(b) == 0L) {
			attr(x, "bootstrap") <- 0
		} else {
			attr(x, "bootstrap") <- b/100
		}
		x
	})

# using a hybrid of branch confidence times length
species_tree <- Zipline(trees, distance=c("bootstrap", "length"))
plot(species_tree, main="Species tree from branch length*support")
}

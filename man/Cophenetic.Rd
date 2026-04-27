\name{Cophenetic}
\alias{Cophenetic}
\title{
Compute cophenetic distances on dendrogram objects
}
\description{
Calculates the matrix of cophenetic distances represented by a dendrogram object.
}
\usage{
Cophenetic(x,
           distance = "length")
}
\arguments{
  \item{x}{
A dendrogram object.
}
  \item{distance}{
The type of distance measure(s) to use.  This should be \code{"length"} (the default), \code{"edges"}, and/or named edge attributes.  Abbreviations are not permitted except in named edge attributes.  (See details section below.)
}
}
\details{
The cophenetic distance between two observations is defined by default as their separation on a dendrogram, also known as the patristic distance.  This function differs from the \code{cophenetic} function in that it does not assume the tree is ultrametric and outputs (by default) the branch length separating pairs of observations (leaves) rather than the height of their merger.

The default distance is the path \code{length} between leaves of the dendrogram.  Specifying \code{"edges"} will sum the number of edges separating leaves, while disregarding the root node of a bifurcating tree and any other edges without splits.

Alternatively, it is possible to customize the measure of distance by specifying named edge attributes.  The values associated with each attribute are converted to numbers and summed without \code{NA} values if containing more than one number.  If the attribute is missing then it is considered to be zero.

If multiple \code{distance} measures are given, these are multiplied together (where they exist) along each edge to define more complex distance measures, such as the \code{"length"} times support probability of an edge.  In this case, one \code{distance} must be \code{"length"} or \code{"edges"}, where the latter has no effect on distance (i.e., multiply by one).
}
\value{
An object of class 'dist'.
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\seealso{
\code{\link{Treeline}}
}
\examples{
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
d1 <- DistanceMatrix(dna, type="dist", correction="JC69", verbose=FALSE)
dend <- Treeline(myDistMatrix=d1)
d2 <- Cophenetic(dend)
cor(d1, d2)

# examples of inter-node and inter-edge distances
d3 <- Cophenetic(dend, distance="edges")
cor(d2, d3) # perfect correlation
range(d3) # 2+ edges between leaves

# example of using an attribute as the distance
tree <- Treeline(dna,
	method="ML", # ML trees get aBayes support values
	model="JC69",
	maxTime=0) # set a stringent time cutoff for speed
attr(tree[[1]], "probability") # support value for first edge
d4 <- Cophenetic(tree, "probability")
cor(d1, d4)
# rescale aBayes (1/3 to 1) supports between 0 and 1
tree_rescaled <- dendrapply(tree,
	function(x) {
		if (!is.leaf(x))
			attr(x, "probability") <- (3*attr(x, "p") - 1)/2
		x
	})
d5 <- Cophenetic(tree_rescaled, "probability")
cor(d4, d5)

# example of specifying multiple distance measures
d6 <- Cophenetic(tree, c("length", "probability"))
cor(d4, d6)
}

\name{DateTree}
\alias{DateTree}
\title{
Time Calibrate a Dendrogram
}
\description{
Uses the relative rate framework to time calibrate a phylogenetic tree with known tip dates or root node age.
}
\usage{
DateTree(x,
         time = 1,
         penalty = 0,
         epsilon = 1e-06)
}
\arguments{
  \item{x}{
An object of class \code{dendrogram} representing a rooted and strictly bifurcating phylogenetic tree with branch lengths reflecting evolutionary distance.
}
  \item{time}{
Either a single number representing the root age, where all the tips are at assumed to be at the same \code{time} (i.e., an isochronous tree), or a named vector containing sampling dates for two or more leaves, where the root age is inferred (i.e., a heterochronous tree).
}
  \item{penalty}{
Either zero (the default) to infer an unconstrained relative rate at each internal node (as in the original algorithm) or a positive regularization \code{penalty}.  (See details section below.)
}
  \item{epsilon}{
Numeric giving the minimum assumed branch length, since branch lengths must be non-zero in the relative rate framework.
}
}
\details{
The relative rate framework (Miura \emph{et al.}, 2020) is used to date the input phylogenetic tree (\code{x}) given the known \code{time} of the root (i.e., node dating) or leaves (i.e., tip dating).  If one value is provided for \code{time}, then the root of the output dendrogram will be at zero height and all the leaves are assumed to occur at \code{time}.  If multiple (named) \code{time} values are provided, then the corresponding leaves are at these heights and the root height is inferred with linear regression.

A single \code{time} can be provided when all leaves are sampled at effectively the same time relative to the (presumably much older) root.  This is useful for dating phylogenetic trees over long evolutionary timescales, where only the root age is known and all tips are assumed to be at approximately the same age.  The default value (i.e., \code{time = 1}) will calibrate the tree into units of relative time (i.e., unit time from zero to one).

Alternatively, \code{time} can contain multiple (named) values providing the known sampling dates for leaf labels.  The root time is then inferred by linear regression between root-to-tip distance and the known sampling \code{time}, which requires that the samples differ across the time scale.  This regression is critical for accurate time calibration towards the root, and a strong phylogenetic signal is desirable (\eqn{R^2 > 0.3}).  The resulting slope represents the rate of evolutionary divergence under a strictly clock-like model.

A relaxed molecular clock is used to make the output tree ultrametric, such that branches represent chronological time rather than evolutionary distance.  The method works by projecting internal node distances from the tips toward the root (Miura \emph{et al.}, 2020), which is based on the assumption that the relative rates of branches must be equal where they intersect.  In this way, a rate parameter is estimated for each internal node of the tree.  Unlike many alternative methods, the relative rate framework does not make prior assumptions about how rates are distributed.

Since the number of inferred rate parameters is high, a regularization \code{penalty} can be introduced to encourage more clock-like rates.  Under this alternative model, deviations from a strict molecular clock (\eqn{rate = 1}) are penalized such that \eqn{\text{rate}_{\text{new}} = \text{rate}_{\text{original}}^{1/(1 + \text{penalty})}}.  Thus, the rate is unaltered when \code{penalty} is zero (the default), corresponding to the unconstrained local rate variation in the original model.  Lower penalties place more trust in the branch lengths, whereas higher penalties encourage more clock-like rates.  Results on simulated benchmarks suggest that a \code{penalty} value of \code{1} generally yields slightly more accurate results than the default value of zero.

Since the underlying algorithm is relative fast, confidence intervals can be established by repeated sampling of the branch lengths and \code{time} values from appropriate distributions.  This approach requires additional information about the expected variance in branch lengths and \code{time}, as well as any variability in tree topology.  It is possible to combine time calibration with bootstrap or jackknife subsampling during tree inference to obtain joint error estimates of node ages.  See the examples section below for a demonstration of using variability in sampling \code{time} to quantify uncertainty in root age.

Several considerations are important for obtaining accurate time trees (Rieux & Balloux, 2016).  First, past demographic changes, recombination, or natural selection can interfere with date estimates.  Second, root location is central to time projection, since the root polarizes the direction of evolutionary change.  When possible, outgroup rooting is preferred, as root placement cannot be disentangled from varying rates of evolution.  Third, branch lengths of the input dendrogram (\code{x}) are expected to be in units proportional to evolutionary change (e.g., substitutions per site).  Branch length or topological inaccuracies in the input dendrogram (\code{x}) will hinder accurate time calibration.
}
\value{
A time calibrated dendrogram (also known as a timetree or chronogram) where the root is at a lower height (further in the past) than the tips.  This directionality is reversed from the input dendrogram (\code{x}), because distance on the output tree represents chronological time that elapsed since the root rather than accumulated evolutionary distance from the tips.
}
\references{
Miura S., \emph{et al.} (2020) A new method for inferring timetrees from temporally sampled molecular sequences. \emph{PLOS Comp. Biol.}, \bold{16(1)}, e1007046.

Rieux A. and Balloux F. (2016) Inferences from tip-calibrated phylogenies: a review and a practical guide. \emph{Molecular Ecology}, \bold{25}, 1911-24.
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\note{
The root of dated trees is always older (lower in value) than the leaves.  This time axis is flipped from that of conventional dendrograms, and can be offset from zero depending on the values of \code{time}.  Therefore, it is necessary to specify axis limits whenever plotting dated dendrograms.
}
\seealso{
\code{\link{Zipline}}, \code{\link{Treeline}}, \code{\link{ReadDendrogram}}, \code{\link{WriteDendrogram}}
}
\examples{
# read a Newick formatted tree of SARS-CoV-2 relationships
fpath <- system.file("extdata", "SARS-CoV-2.newick", package="DECIPHER")
tree <- ReadDendrogram(fpath)
plot(tree)

# parse sampling dates from labels
times <- labels(tree)
times <- setNames(sapply(strsplit(times, " "), tail, n=1),
	times)
times <- as.Date(times)
times <- times - as.Date("2020-01-01")
times <- setNames(as.numeric(times), names(times))
times <- 2020 + times/365
head(times) # sampling year

# evaluate phylogenetic strength
dists <- attr(tree, "height") - rapply(tree, attr, which="height")
plot(times,
	dists,
	xlab="Sampling time",
	ylab="Root-to-tip distance")
fit <- lm(dists~times)
abline(fit, h=0)
summary(fit)$adj.r.squared # R^2 > 0.3 is preferable

# convert from distances to times
timetree <- DateTree(tree, times)

# need to adjust y-axis limits for plotting
lim <- range(attr(timetree, "height"),
	rapply(timetree, attr, which="height"))
plot(timetree, ylim=rev(lim))

# display sampling time as leaf label
relabeled <- dendrapply(timetree,
	function(x) {
		if (is.leaf(x))
			attr(x, "label") <- format(times[attr(x, "label")])
		x
	})
plot(relabeled, ylim=rev(lim))

# root height represents origin time of SARS-CoV-2
root <- attr(timetree, "height") # mid-November 2019

# measure variability induced by uncertainty in sampling time
est <- numeric(100) # replicate measurements
for (i in seq_along(est)) {
	delta <- sample(c(-1, 0, 1), length(times), replace=TRUE)
	alt_times <- times + delta/365 # shift by +/- 1 day
	est[i] <- attr(DateTree(tree, time=alt_times), "height")
}
hist(est, xlab="Estimated uncertainty in root time")
abline(v=root, lwd=2, lty=2) # unperturbed estimate

# timetree for samples sharing an ancient ancestor
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
x <- Treeline(dna, model="TN93+F")
plot(x, leaflab="none") # not ultrametric
y <- DateTree(x)
plot(y, leaflab="none", ylim=c(1, 0))
}

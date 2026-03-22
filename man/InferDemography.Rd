\name{InferDemography}
\alias{InferDemography}
\title{
Infer Demographic History from Allele Frequencies
}
\description{
Fits a population genetics model to an input site frequency spectrum (SFS) or one derived from aligned sequences.  Returns estimated effective population sizes at time intervals in the past.
}
\usage{
InferDemography(x,
                readingFrame = NA,
                mu = 1e-08,
                ploidy = 1,
                informationCriterion = "BIC",
                showPlot = FALSE,
                verbose = TRUE)
}
\arguments{
  \item{x}{
A numeric vector of folded site frequencies (i.e., SFS[0], SFS[1], SFS[2], ..., SFS[\code{n}/2]), or a \code{DNAStringSet} or \code{RNAStringSet} of \code{n} aligned sequences from which to derive the folded SFS.
}
  \item{readingFrame}{
Either \code{NA} to use all sites, or the \code{readingFrame} (i.e., \code{1}, \code{2}, or \code{3}) to limit analysis to only third codon positions.  Only applicable if \code{x} is a \code{XStringSet}.
}
  \item{mu}{
A numeric giving the mutation rate per site per generation, which linearly affects the scaling of inferred demographic variables.
}
  \item{ploidy}{
An integer providing the ploidy (e.g., \code{1} for haploid or \code{2} for diploid organisms).
}
  \item{informationCriterion}{
Character string specifying which information criterion (i.e., \code{"AIC"} or \code{"BIC"}) to use in determining the number of time intervals, or an integer specifying the number of intervals in the output.
}
  \item{showPlot}{
Logical specifying whether to show the fitted site frequency spectrum and inferred changes in effective population size (i.e., stairway plot).
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The site frequency spectrum (SFS) represents the distribution of allele frequencies within a population, which is sensitive to natural selection and changes in population size (Johri, et al., 2022).  Using sites that are presumed to be neutrally evolving, it is possible to reconstruct historical changes in population size.  This approach assumes frequent recombination, migration, and an unstructured population.  Effective population sizes are inferred as a step function corresponding to different levels of the coalescent, with transitions times scaled to generations in the past.

\code{InferDemography} implements an extensive, but inexhaustive, procedure for optimizing transition points based on the likelihood equations described in Lynch, et al. (2020).  The results provide an estimate of population expansions and contractions that reasonably fit simulations of the Wright-Fisher model given sufficient data.  Axis scaling is dependent on the mutation rate (\code{mu}) and \code{ploidy}, which affect the quantitative (but not qualitative) results.  Error bounds can be determined through bootstrapping by repeated sampling of sequences in \code{x} with replacement or by Poisson sampling frequency classes in \code{x}.

Increasing the number of sampled input species and neutral sites permits more accurate inference, although all sequences must originate from the same population.  For this reason, the observed SFS that is output can be aggregated across multiple gene alignments of the same organisms and then used as input to increase the total number of polymorphisms.  Notably, some organisms' life histories may poorly fit the standard Kingman coalescent used here (Freund, et al., 2023), and it is important to always verify a close match between the observed SFS and estimated SFS.

This approach requires multiple sequences sampled from a population with minimal divergence, since sites with more than two different alleles are discarded under the infinite sites model.  The results can be used to compare demographics of genes within a species or contrast different populations' histories.  In general, the trajectory of population expansions and contractions can provide insights into the balance between drift and selective forces previously acting on a population or gene.
}
\value{
A named numeric vector with the following elements:
(1) \code{Intervals} - number of different effective population sizes\cr
(2) \code{LogLikelihood} - fitted model's log-likelihood\cr
(3) \code{Time} - estimated time boundaries of each interval in units of generations ago\cr
(4) \code{Ne} - estimated effective population size during each interval\cr
(5) \code{Observed} - the observed (or input) folded site frequency spectrum\cr
(6) \code{Estimated} - the estimated (fitted) folded site frequency spectrum

Different \code{Time} and \code{Ne} estimates are named by their merger level in the coalescent (\code{n} to \code{2}).
}
\note{
It is possible to exclude specific frequency classes in \code{x} by specifying them as \code{NA}.  For example, if singletons (i.e., SFS[1]) are deemed unreliable due to sequencing error, it is feasible to exclude them by setting \code{x[2]} to \code{NA}, where \code{x} is the input numeric vector of folded site frequencies.
}
\references{
Freund, F., et al. (2023). Interpreting the pervasive observation of U-shaped Site Frequency Spectra. PLOS Genetics, \bold{19(3)}, e1010677.

Johri, P., et al. (2022). Recommendations for improving statistical inference in population genomics. PLOS Biology, \bold{20(5)}, e3001669.

Lynch, M., et al. (2020). Inference of Historical Population-Size Changes with Allele-Frequency Data. G3, \bold{10(1)}, 211-223.
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\seealso{
\code{\link{InferRecombination}}, \code{\link{InferSelection}}

Run \code{vignette("PopulationGenetics", package = "DECIPHER")} to see a related vignette.
}
\examples{
# example of providing a folded SFS from the supplement of Lynch, et al. (2020)
SFS <- c(9526, 3998, 2315, 1487, 1075, 873, 689, 570, 567, 627, 547, 605, 460, 573, 236)
SFS <- c(1e7 - sum(SFS), SFS) # add monomorphic sites as first value
InferDemography(SFS, mu=1.2e-8, ploidy=2, show=TRUE)

# resample the folded SFS to test for stability
SFS2 <- rbinom(length(SFS), sum(SFS), SFS/sum(SFS))
InferDemography(SFS2, mu=1.2e-8, ploidy=2, show=TRUE)

# example of providing sequences from a (haploid) species
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
dna <- dna[startsWith(names(dna), "Helicobacter pylori")]
DNA <- AlignTranslation(dna) # align the translation then reverse translate
DNA
InferDemography(DNA, readingFrame=1, mu=1e-9, ploidy=1, show=TRUE)
}

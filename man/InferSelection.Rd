\name{InferSelection}
\alias{InferSelection}
\title{
Infer Codon Selection on Protein Coding Sequences
}
\description{
Infers the magnitude and direction of natural selection by fitting a population genetics model to aligned protein coding sequences.  Returns estimated selection parameters, including the Ka/Ks ratio (omega) by region.
}
\usage{
InferSelection(myDNAStringSet,
               readingFrame = 1L,
               windowSize = NA,
               tolerance = 5e-05,
               geneticCode = GENETIC_CODE,
               showPlot = FALSE,
               verbose = TRUE)
}
\arguments{
  \item{myDNAStringSet}{
A \code{DNAStringSet} of aligned sequences.
}
  \item{readingFrame}{
A numeric vector giving the starting position of the first codon in the alignment, either \code{1} (the default), \code{2}, or \code{3}.
}
  \item{windowSize}{
Either \code{NA} to estimate the average Ka/Ks ratio (omega) across all codons, or an integer giving the size (in codons) of non-overlapping windows in which to estimate individual Ka/Ks ratios (omegas).  For example, setting the \code{windowSize} to \code{1} will estimate an omega for every codon.
}
  \item{tolerance}{
Numeric determining the relative convergence tolerance.  Optimization will cease when the relative likelihood has changed by less than \code{tolerance}.
}
  \item{geneticCode}{
A character vector giving the genetic code in the same format as \code{GENETIC_CODE} (the default).
}
  \item{showPlot}{
Logical specifying whether to show the estimated Ka/Ks ratio(s) (omega) along the alignment.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The Ka/Ks ratio, also known as omega or dN/dS, is a measure of the magnitude and direction of natural selection operating on protein coding genes.  It represents the rate of non-synonymous versus synonymous codon substitutions, and a bias in this rate is indicative of selection.  Ratios of \code{1} are expected under neutral evolution, but ratios less than \code{1} are typically observed due to negative (purifying) selection acting to maintain the protein sequence.  In contrast, ratios greater than \code{1} are of particular interest, as they may represent the presence of positive (Darwinian) selection for protein changes.  Ratios are rarely greater than \code{1} on average for entire coding sequences, so it is common to rank genes by the fraction (or total number) of codons under positive selection (Moutinho, et al., 2023).

\code{InferSelection} fits a three parameter NY98 substitution rate matrix (Nielsen & Yang, 1998) to the observed distribution of codons at sites in an alignment using a population genetics method (Wilson, D & CRyPTIC Consortium, 2020).  This approach derives estimates of kappa (transition/transversion ratio), theta (population scaled mutation rate), and omega (Ka/Ks ratio).  The Ka/Ks ratio (omega) can be estimated for the entire alignment or non-overlapping windows of \code{windowSize} codons.  This estimation approach effectively considers omega as a mutational bias in the rate of non-synonymous versus synonymous changes.  Codon frequencies are derived from nucleotide frequencies in the input sequences and fitted parameters are determined by maximum likelihood estimation.

The method used by \code{InferSelection} requires a low mutation rate approximation, so it is expected to work best on a single population from the same species.  It is noteworthy that omega is known to be closer to \code{1} than expected when measured within populations, rather than between populations where mutations are fixed (Kryazhimskiy & Plotkin, 2008).  This prevents omega from being accurately transformed into (and interpreted as) a population-scaled selection coefficient.  Notwithstanding this limitation, the advantage of this approximation is that a phylogenetic tree is not required, and the method can easily scale to large numbers of input sequences.  Many sequences are typically needed for accurate estimates since variation is relatively rare within populations.

The population genetics model employed here also assumes independence between sites.  Nevertheless, simulations show the method is robust to model violations when there is a lack of recombination (Wilson, D & CRyPTIC Consortium, 2020).  The independence assumption confers the advantages of scalability, simple handling of missing data (e.g., gaps), and no requirement for haplotype information.  Sequences are assumed to be randomly sampled from an unstructured population with constant population size over time.  The population sample is very important, as biased sampling could result in correspondingly biased inferences about selection.

Statistical significance is obtained from the likelihood ratio test based on the chi-squared distribution with \code{1} degree of freedom, comparing the fitted value of omega to \code{1} when \code{windowSize} is \code{NA} (the default) or each omega to the average omega when there are multiple windows.  The likelihood ratio test is conservative but powerful when testing for selection (Anisimova, et al., 2001).  Alternative estimates of statistical significance can be obtained by bootstrapping the input sequences when the number of sequences is sufficiently large.  More significant p-values are displayed in a brighter green color when \code{showPlot} is \code{TRUE}.
}
\value{
A named numeric vector with the following elements:
(1) \code{LogLikelihood} - fitted model's log-likelihood\cr
(2) \code{theta} - expected substitution rate in units of 2*ploidy*Ne generations\cr
(3) \code{kappa} - estimated transition to transversion ratio\cr
(4) \code{omega} - estimated Ka/Ks (dN/dS) ratio(s) per window\cr
(5) \code{pvalue} - approximate p-value per window\cr
}
\references{
Anisimova, M., et al. (2001). Accuracy and power of the likelihood ratio test in detecting adaptive molecular evolution. Molecular Biology and Evolution, \bold{18(8)}, 1585-1592.

Kryazhimskiy, S. & Plotkin, J. (2008). The population genetics of dN/dS. PLoS Genetics, \bold{4(12)}, e1000304.

Moutinho, A., et al. (2019). Variation of the adaptive substitution rate between species and within genomes. Evolutionary Ecology, \bold{34(3)}, 315-338.

Nielsen, R. & Yang, Z. (1998). Likelihood models for detecting positively selected amino acid sites and applications to the HIV-1 envelope gene. Genetics, \bold{148(3)}, 929-936.

Wilson, D. & CRyPTIC Consortium. (2020) GenomegaMap: Within-Species Genome-Wide dN/dS Estimation from over 10,000 Genomes. Molecular Biology and Evolution, \bold{37(8)}, 2450-2460.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{InferDemography}}, \code{\link{InferRecombination}}

Run \code{vignette("PopulationGenetics", package = "DECIPHER")} to see a related vignette.
}
\examples{
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
DNA <- readDNAStringSet(fas)
DNA <- DNA[startsWith(names(DNA), "Helicobacter")] # subset to species
DNA <- AlignTranslation(DNA)

InferSelection(DNA, windowSize=NA, showPlot=TRUE)
# note: set windowSize=1 to estimate omega per codon
}

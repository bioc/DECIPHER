\name{InferRecombination}
\alias{InferRecombination}
\title{
Fit Recombination Parameters to Correlation Profiles
}
\description{
Derives a correlation profile from aligned sequences and fits a coalescent-based model of recombination.  Returns recombination parameters along with the fitted correlation profile.
}
\usage{
InferRecombination(x,
                  readingFrame = NA,
                  position = 1:3,
                  N = 249,
                  showPlot = FALSE,
                  verbose = TRUE)
}
\arguments{
  \item{x}{
A \code{DNAStringSet} or \code{RNAStringSet} of aligned sequences, or a list or XStringSetList containing multiple \code{XStringSet} objects.
}
  \item{readingFrame}{
Either \code{NA}, a numeric vector giving a single reading frame for all alignments, or an individual reading frame for each alignment in \code{x}.  The \code{readingFrame} can be either \code{1}, \code{2}, \code{3} to signify codons beginning on the first, second, and third nucleotide position, respectively. If \code{NA} (the default), analysis is performed without reference to a reading frame.  (See details section below.)
}
  \item{position}{
Numeric vector containing the codon \code{position}(s) to be analyzed, or \code{NA} to analyze all positions together.  Typically, the third (\code{3}) codon position is analyzed for coding sequences and all positions (\code{NA}) for non-coding sequences.  (See details section below.)
}
  \item{N}{
Numeric giving the maximum number of positions away from the reference (initial) position to include in the correlation profile.
}
  \item{showPlot}{
Logical specifying whether or not to plot the correlation profile(s) and fitted curves.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Recombination accelerates the adaptive evolution of many organisms.  The transfer of genetic fragments leaves behind a signature in the form of decaying autocorrelation among substitutions.  This signal can be visualized with a correlation profile that shows the probability of a genetic difference at increasing distances away from another genetic difference.  Flat correlation profiles are evidence for an absence of recombination since the last common ancestor.  In contrast, more ``L'' shaped correlation profiles result from greater degrees of recombination (Steinberg, 2023).

\code{InferRecombination} fits a coalescent-based model to a correlation profile derived from one or more sets of aligned sequences (\code{x}).  This model assumes the sequences correspond to a population whose genealogical structure is characterized by a single average coalescence rate.  This assumption is automatically true for a pair of sequences, and may hold true within (e.g., species) or between clusters of sequences (Steinberg, 2022).  The model can be applied to any organisms belonging to a focal population recombining genetic fragments from an external pool.

The population genetics model consist of three free variables that are fitted to the correlation profile (Lin & Kussell, 2019).  These variables can be used to infer six other parameters of interest corresponding to the sample (\code{x}) or (external) pool in which they recombine.  (See value section below.)  This approach is fast, accurate, and does not rely on phylogenetic reconstruction.  Confidence intervals can be determined by bootstrapping the input alignment(s).

It is possible to indicate the \code{readingFrame} and codon \code{position}(s) that should be assessed.  Typically, the analysis is performed on the third \code{position} of synonymous codons to mitigate the influence of adaptive substitutions.  If the \code{readingFrame} is \code{NA}, positions are analyze with respect to the (initial) substitutions, and a positional pattern may emerge in due to variability in substitution rates at each codon position.  It is also possible to include all positions (\code{NA}), such as when analyzing non-coding RNAs.

The output consists of a matrix of measured values and inferred parameters.  The recombined coverage is 0\% when the sequences have evolved independently since their last common ancestor and 100\% when the (input) sample has combined all of its nucleotides with the (external) pool.  The ratio of recombination rate to mutation rate can be interpreted as the ratio of substitutions due to recombination relative to point mutation.  The ratio is greater than \code{1} when recombination contributes more to genetic diversity than point mutation, and vice versa when less than \code{1}.
}
\value{
A \code{matrix} with named rows for each value and a column for each position.  The meaning of each parameter is described in Lin & Kussell (2019).

The first three rows are the fitted parameters:
(1) fragment - mean recombined genetic fragment size (in base pairs)
(2) Theta_sample - mutational divergence within the (input) sample
(3) Phi_sample - recombinational divergence within the (input) sample

The next six rows are derived from the fitted parameters:
(4) Theta_pool - mutational divergence within the (external) pool
(5) Phi_pool - recombinational divergence within the (external) pool
(6) ratio - the relative rate of recombination to mutation (also known as gamma/mu or r/m)
(7) coverage - proportion of sites in the (input) sample whose diversity originated from the (external) pool
(8) d_pool - pairwise diversity of the (external) pool (i.e., probability that any two sequences in the pool differ at a site)
(9) d_clonal - pairwise diversity of the (input) sample due to accumulated substitutions (i.e., probability any two sequences in the sample differ at a site due to clonal mutation)

The remaining rows contain fixed values:
(10) d_sample - measured pairwise diversity due to both recombination and mutation (i.e., probability two sequences differ at a site)
(11) Position - distance away from the initial substitution (in base pairs)
(12) Profile - measured probability of a difference at each position (i.e., the measured correlation profile)
(13) Fitted - estimated probability of a difference at each position (i.e., the fitted correlation profile)
}
\references{
Lin, M. & Kussell, E. (2019). Inferring bacterial recombination rates from large-scale sequencing datasets. Nature Methods, \bold{16(2)}, 199-204.

Steinberg, A., et al. (2022). Core genes can have higher recombination rates than accessory genes within global microbial populations. eLife, \bold{11}, e78533.

Steinberg, A., et al. (2023). Correlated substitutions reveal SARS-like coronaviruses recombine frequently with a diverse set of structured gene pools. PNAS, \bold{120(5)}, e2206945119.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\examples{
# example for an alignment of coding sequences
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
DNA <- readDNAStringSet(fas)
DNA <- AlignTranslation(DNA)
ans <- InferRecombination(DNA, position=3, readingFrame=1L, showPlot=TRUE)
head(ans, n=10)

# example for an alignment of non-coding sequences
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
DNA <- readDNAStringSet(fas)
ans <- InferRecombination(DNA, position=NA, showPlot=TRUE)
head(ans, n=10)

if (require("RSQLite", quietly=TRUE)) {
	# example of inferring recombination among genomes
	db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
	synteny <- FindSynteny(db, minScore=50)
	DNA <- AlignSynteny(synteny, db)
	DNA <- do.call(c, unname(DNA))
	ans <- InferRecombination(DNA, showPlot=TRUE)
	head(ans, n=10)
}
}

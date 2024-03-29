\name{CalculateEfficiencyPCR}
\alias{CalculateEfficiencyPCR}
\title{
Predict Amplification Efficiency of Primer Sequences
}
\description{
Calculates the amplification efficiency of primers from their hybridization efficiency and elongation efficiency at the target site.
}
\usage{
CalculateEfficiencyPCR(primer,
                       target,
                       temp,
                       P,
                       ions,
                       batchSize = 1000,
                       taqEfficiency = TRUE,
                       maxDistance = 0.4,
                       maxGaps = 2,
                       processors = 1)
}
\arguments{
  \item{primer}{
A \code{DNAStringSet} object or character vector with unaligned primer sequences in 5' to 3' orientation.
}
  \item{target}{
A \code{DNAStringSet} object or character vector with unaligned target or non-target sequences in 5' to 3' orientation.
}
  \item{temp}{
Numeric specifying the annealing temperature used in the PCR reaction.
}
  \item{P}{
Numeric giving the molar concentration of primers in the reaction.
}
  \item{ions}{
Numeric giving the molar sodium equivalent ionic concentration.  Values may range between 0.01M and 1M.
}
  \item{batchSize}{
Integer specifying the number of primers to simulate hybridization per batch.  See the Description section below.
}
  \item{taqEfficiency}{
Logical determining whether to make use of elongation efficiency and maxDistance to increase predictive accuracy for \emph{Taq} DNA Polymerase amplifying primers with mismatches near the 3' terminus.  Note that this should be set to FALSE if using a high-fidelity polymerase with 3' to 5' exonuclease activity.
}
  \item{maxDistance}{
Numeric specifying the maximal fraction of mismatched base pairings on a rolling basis beginning from the 3' end of the primer.  Only used if \code{taqEfficiency} is \code{TRUE}.
}
  \item{maxGaps}{
Integer specifying the maximum number of insertions or deletions (indels) in the primer/target alignment.  Only used if \code{taqEfficiency} is \code{TRUE}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
}
\details{
Amplification of pairwise \code{primer}/\code{target} pairs is simulated \emph{in silico}.  A complex model of hybridization is employed that takes into account the side reactions resulting from probe-folding, target-folding, and primer-dimer formation.  The resulting hybridization efficiency is multiplied by the elongation efficiency to predict the overall efficiency of amplification.

Free energy is obtained from system calls to OligoArrayAux, which must be properly installed (see the Notes section below).  Primer/target pairs are sent to OligoArrayAux in batches of \code{batchSize}, which prevents systems calls from being too many characters.  Note that OligoArrayAux does not support degeneracy codes (non-base letters), although they are accepted without error.  Any sequences with ambiguity should be expanded into multiple permutations with \code{\link{Disambiguate}} before input.
}
\value{
A vector of predicted efficiencies for amplifying each \code{primer}/\code{target} pair of sequences.
}
\references{
ES Wright et al. (2013) "Exploiting Extension Bias in PCR to Improve Primer Specificity in Ensembles of Nearly Identical DNA Templates." Environmental Microbiology, doi:10.1111/1462-2920.12259.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\note{
The program OligoArrayAux (\url{http://www.unafold.org/Dinamelt/software/oligoarrayaux.php}) must be installed in a location accessible by the system.  For example, the following code should print the installed OligoArrayAux version when executed from the R console:

\code{system("hybrid-min -V")}
}
\seealso{
\code{\link{AmplifyDNA}}, \code{\link{DesignPrimers}}, \code{\link{DesignSignatures}}
}
\examples{
primers <- c("AAAAACGGGGAGCGGGGGG", "AAAAACTCAACCCGAGGAGCGCGT")
targets <- reverseComplement(DNAStringSet(primers))
# not run (must have OligoArrayAux installed first):
\dontrun{CalculateEfficiencyPCR(primers, targets, temp=75, P=4e-7, ions=0.225)}
}

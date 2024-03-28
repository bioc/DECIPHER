\name{IndexSeqs}
\alias{IndexSeqs}
\title{
Build an inverted index
}
\description{
Builds an inverted index from a set of amino acid, DNA, or RNA sequences.
}
\usage{
IndexSeqs(subject,
          K,
          sensitivity,
          percentIdentity,
          patternLength,
          step = 1,
          alphabet = AA_REDUCED[[171]],
          maskRepeats = TRUE,
          maskLCRs = TRUE,
          maskNumerous = TRUE,
          batchSize = 1e+07,
          processors = 1,
          verbose = TRUE)
}
\arguments{
  \item{subject}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of target (unaligned) sequences to use as the subject of the search.
}
  \item{K}{
Integer providing the k-mer length.  Typical values are \code{5} to \code{7} for amino acids (in the default \code{alphabet}) and \code{8} to \code{12} for nucleotides.
}
  \item{sensitivity}{
Numeric giving the goal search sensitivity, which is used in the absence of specifying \code{K}.  Typically near, but always less than, \code{1}.
}
  \item{percentIdentity}{
Numeric identifying the goal percent identity for the given \code{sensitivity}, which is used in the absence of specifying \code{K}.  Values closer to \code{100} allow for larger \code{K} and, thereby, faster searches.
}
  \item{patternLength}{
Integer setting the expected (minimum) length of query sequences, which is used in the absence of specifying \code{K}.
}
  \item{step}{
Integer determining the number of positions between the start of adjacent k-mers.  Must be between \code{1} and \code{K}.  Larger values reduce the memory required for the inverted index at the expense of search sensitivity.
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between protein sequences.  A non-reduced amino acid alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.  Note that choice of \code{alphabet} affects the optimal value of \code{K}, with smaller alphabets requiring larger \code{K} and vice versa.
}
  \item{maskRepeats}{
Logical specifying whether to mask repeats in the \code{subject}.  (See details section below.)
}
  \item{maskLCRs}{
Logical indicating whether to mask low complexity regions in the \code{subject}.  (See details section below.)
}
  \item{maskNumerous}{
Logical indicating whether to mask frequent k-mers in the \code{subject}, or a positive numeric specifying the degree of masking to apply.  (See details section below.)
}
  \item{batchSize}{
Integer defining the number of sequences to process in a batch.  Smaller values reduce the function's memory footprint, potentially at the cost of increased runtime.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
An \code{InvertedIndex} object functions much like a the index at the back of a book, except with the goal of finding homologous sequence regions.  It primarily contains the locations of length \code{K} subsequences (i.e., k-mers) in \code{subject}.  Only the set of unmasked k-mers separated every \code{step} positions are stored.  In general, lower values of \code{K} and, especially, \code{step} are preferable for improving search sensitivity and specificity.  If an appropriate value for \code{K} is unknown, it is possible to automatically calculate \code{K} by providing a goal search \code{sensitivity} for sequences of \code{patternLength} positions with a given \code{percentIdentity} to a target sequence.  The fastest matching value of \code{K} will be automatically selected, or an error will be returned when the constraints cannot be met.

Masking repeats and low complexity regions helps to avoid spurious hits that are not related to sequence homology (i.e., common descent), while masking extremely frequent k-mers can improve search speed and decrease the size of the inverted index.  If \code{maskRepeats} is \code{TRUE} (the default), masking is first applied to tandem repeats, while maintaining the first copy of the repeat.  Next, if \code{maskLCRs} is \code{TRUE} (the default), any remaining low complexity regions are masked that are at least as long as the k-mer length (i.e., \code{K}).  Finally, if \code{maskNumerous} is \code{TRUE} (the default), extremely frequent k-mers that remain are masked on a per sequence basis.  It is also possible to control the degree of masking frequent k-mers by supplying a positive number for \code{maskNumerous}, representing the -log(p-value) of observing that many k-mers by chance.  Giving larger positive values for \code{maskNumerous} will mask fewer k-mers.

Proper masking is critical to balance search speed and sensitivity.  When mapping reads to a \code{subject} genome it is often best to leave \code{maskNumerous} as \code{TRUE} but set \code{maskRepeats} and \code{maskLCRs} to \code{FALSE}.  This is because masking whole regions may preclude some reads from mapping, whereas masking the most frequent k-mers is sometimes necessary for speed.  In contrast, when the \code{subject} is a set of shorter sequences (e.g., proteins) there are rarely extremely frequent k-mers in a single sequence, so \code{maskNumerous} has less (or no) effect.  In this case, keeping \code{maskRepeats} and \code{maskLCRs} as \code{TRUE} will prevent k-mer matches between sequences that are not related by common descent.  This is because some biological mechanisms can generate similar sequence patterns in unrelated sequences, which would otherwise mislead the statistical model of homology based on the likelihood of shared k-mers.
}
\value{
An object of class \code{InvertedIndex}, which is stored as a list with components:
  \item{k}{
The input or automatically calculated value for \code{K}.
}
  \item{step}{
The input value of \code{step}.
}
  \item{alphabet}{
Numbered letter groups present in the sequence alphabet.
}
  \item{frequency}{
Numeric frequencies of each letter grouping.
}
  \item{count}{
The number of times every possible (unmasked) k-mer was observed.
}
  \item{length}{
An integer giving the number of (unmasked) k-mers in each sequence in \code{subject}.
}
  \item{location}{
The location of each k-mer within every sequence in \code{subject}, ordered by k-mer.
}
  \item{index}{
The index of each k-mer within every sequence in \code{subject}, ordered by k-mer.
}
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{SearchIndex}}, \code{\link{AlignPairs}}
}
\examples{
# import target sequences
fas <- system.file("extdata", "PlanctobacteriaNamedGenes.fas.gz", package="DECIPHER")
seqs <- readAAStringSet(fas)

# build an inverted index, specifying K
index <- IndexSeqs(seqs, K=6L)
index # K = 6

# alternatively, determine K automatically
index <- IndexSeqs(seqs, sensitivity=0.99, percentIdentity=60, patternLength=300)
index # K = 3
}

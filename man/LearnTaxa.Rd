\name{LearnTaxa}
\alias{LearnTaxa}
\title{
Train a Classifier for Assigning Taxonomy
}
\description{
Trains a classifier based on a reference taxonomy containing sequence representatives assigned to taxonomic groups.
}
\usage{
LearnTaxa(train,
          taxonomy,
          rank = NULL,
          K = floor(log(100*quantile(width(train), 0.99), 4)),
          minFraction = 0.01,
          maxFraction = 0.06,
          maxIterations = 10,
          multiplier = 100,
          maxChildren = 200,
          verbose = TRUE)
}
\arguments{
  \item{train}{
A \code{DNAStringSet} or \code{RNAStringSet} of unaligned sequences.
}
  \item{taxonomy}{
Character string providing the reference taxonomic assignment for each sequence in \code{train}.  Taxonomic ranks are separated by semicolons (``;'') beginning with ``Root''.
}
  \item{rank}{
Optionally, a \code{data.frame} with 5 named columns giving the ``Index'' (i.e., 0 to the number of unique taxa), ``Name'' (i.e., taxon name), ``Parent'' (i.e., ``Index'' of the parent taxon), ``Level'' (i.e., integer rank level), and ``Rank'' (e.g., ``genus'') of each taxonomic rank.  This information is often provided in a separate ``taxid'' file along with publicly available training sequence sets.
}
  \item{K}{
Integer specifying the k-mer size.  The default value of \code{K} is such that k-mer matches between sequences are rarely (~1\%) found by chance.
}
  \item{minFraction}{
Numeric giving the minimum fraction of k-mers to sample during the initial tree descent phase of the classification algorithm.  (See details section below.)
}
  \item{maxFraction}{
Numeric giving the maximum fraction of k-mers to sample during the initial tree descent phase of the classification algorithm.  (See details section below.)
}
  \item{maxIterations}{
Integer specifying the maximum number of iterations to attempt re-classification of a training sequence before declaring it a ``problem sequence''.  (See details section below.)
}
  \item{multiplier}{
Numeric indicating the degree to which individual sequences have control over the fraction of k-mers sampled at any edge during the initial tree descent phase of the classification algorithm.  (See details section below.)
}
  \item{maxChildren}{
Integer giving the maximum number of child taxa of any taxon at which to consider further descending the taxonomic tree.  A value of \code{1} will prevent use of the tree descent algorithm altogether.  Lower values may decrease classification speed, but result in output objects that require less memory.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Learning about the training data is a two part process consisting of (i) forming a taxonomic tree and then (ii) ensuring that the \code{train}ing sequences can be correctly reclassified.  The latter step relies on reclassifying the sequences in \code{train} by descending the taxonomic tree, a process termed ``tree descent''.  Ultimately, the goal of tree descent is to quickly and accurately narrow the selection of groups where a sequence may belong.  During the learning process, tree descent is tuned so that it performs well when classifying new sequences.

The process of training the classifier first involves learning the taxonomic tree spanning all of the reference sequences in \code{train}.  Typically, reference taxonomic classifications are provided by an authoritative source, oftentimes along with a ``taxid'' file containing taxonomic \code{rank} information.  The taxonomic tree may contain any number of levels (e.g., Root, Phylum, Class, Order, Family, Genus) as long as they are hierarchically nested and always begin with ``Root''.

The second phase of training the classifier, tree descent, involves learning the optimal set of k-mers for discerning between the different sub-groups under each edge.  Here a fraction of the k-mers with the greatest discerning power are matched to a training sequence, and this process is repeated with 100 random subsamples to decide on the set of possible taxonomic groups to which a training sequence may belong.

The learning process works by attempting to correctly re-classify each training sequence in the taxonomy.  Initially, \code{maxFraction} of informative k-mers are repeatedly sampled at each edge during tree descent.  Training sequences that are incorrectly classified at an edge will lower the fraction of k-mers that are sampled by an amount that is proportional to \code{multiplier}.  As the fraction of sampled k-mers decreases, the tree descent process terminates at higher rank levels.

A major advantage of tree descent is that it both speeds up the classification process and indicates where the training set likely contains mislabeled sequences or incorrectly-placed taxonomic groups.  Training sequences that are not correctly classified within \code{maxIterations} are marked as ``problem sequences'', because it is likely that they are mislabeled.  If enough sequences have difficulty being correctly classified at an edge that the fraction drops below \code{minFraction}, then the edge is recorded as a ``problem group''.

The final result is an object that can be used for classification with \code{\link{IdTaxa}}, as well as information about \code{train} that could be used to help correct any errors in the taxonomy.
}
\value{
An object of class \code{Taxa} and subclass {Train}, which is stored as a list with components:
  \item{taxonomy}{
A character vector containing all possible groups in the taxonomy.
}
  \item{taxa}{
A character vector containing the basal taxon in each taxonomy.
}
  \item{ranks}{
A character vector of rank names for each taxon, or \code{NULL} if \code{rank} information was not supplied.
}
  \item{levels}{
Integer giving the rank level of each taxon.
}
  \item{children}{
A list containing the index of all children in the taxonomy for each taxon.
}
  \item{parents}{
An integer providing the index of the parent for each taxon.
}
  \item{fraction}{
A numeric between \code{minFraction} and \code{maxFraction} that represents the learned fraction of informative k-mers to sample for each taxon during the initial tree descent phase of the classification algorithm.  Problem groups are marked by a fraction of \code{NA}.
}
  \item{sequences}{
List containing the integer indices of sequences in \code{train} belonging to each taxon.
}
  \item{kmers}{
List containing the unique sorted k-mers (converted to integers) belonging to each sequence in \code{train}.
}
  \item{crossIndex}{
Integer indicating the index in taxonomy of each sequence's taxonomic label.
}
  \item{K}{
The value of \code{K} provided as input.
}
  \item{IDFweights}{
Numeric vector of length \code{4^K} providing the inverse document frequency weight for each k-mer.
}
  \item{decisionKmers}{
List of informative k-mers and their associated relative frequencies for each internal edge in the taxonomy.
}
  \item{problemSequences}{
A \code{data.frame} providing the ``Index'', ``Expected'' label, and ``Predicted'' taxon for sequences that could not be correctly classified during the initial tree descent phase of the algorithm.
}
  \item{problemGroups}{
Character vector containing any taxonomic groups that repeatedly had problems with correctly re-classifying sequences in \code{train} during the initial tree descent phase of the classification algorithm.  Problem groups likely indicate that a number of the sequences (or an entire group of sequences) assigned to the problem group are incorrectly placed in the taxonomic tree.
}
}
\references{
Coming soon!
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{IdTaxa}}, \code{\link{Taxa-class}}
}
\examples{
# import training sequences
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

# parse the headers to obtain a taxonomy
s <- strsplit(names(dna), " ")
genus <- sapply(s, `[`, 1)
species <- sapply(s, `[`, 2)
taxonomy <- paste("Root", genus, species, sep="; ")
head(taxonomy)

# train the classifier
trainingSet <- LearnTaxa(dna, taxonomy)
trainingSet

# view information about the classifier
plot(trainingSet)
}
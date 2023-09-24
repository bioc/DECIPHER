\name{PredictDBN}
\alias{PredictDBN}
\title{
Predict RNA Secondary Structure in Dot-Bracket Notation
}
\description{
Predicts a consensus RNA secondary structure from a multiple sequence alignment using mutual information.
}
\usage{
PredictDBN(myXStringSet,
           type = "states",
           minOccupancy = 0.4,
           impact = c(1, 1.2, 0.4, -1),
           avgProdCorr = 1,
           slope = 2,
           shift = 1.3,
           threshold = 0.4,
           pseudoknots = 1,
           weight = NA,
           useFreeEnergy = TRUE,
           deltaGrules = NULL,
           processors = 1,
           verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} object containing aligned sequences.
}
  \item{type}{
Character string indicating the type of results desired.  This should be (an unambiguous abbreviation of) one of \code{"states"}, \code{"pairs"}, \code{"evidence"}, \code{"scores"}, \code{"structures"}, or \code{"search"}.  (See value section below.)
}
  \item{minOccupancy}{
Numeric specifying the minimum occupancy (1 - fraction of gaps) required to include a column of the alignment in the prediction.
}
  \item{impact}{
A vector with four elements giving the weights of A/U, G/C, G/U, and other pairings, respectively.  The last element of \code{impact} is the penalty for pairings that are inconsistent with two positions being paired (e.g., A/- or A/C).
}
  \item{avgProdCorr}{
Numeric specifying the weight of the average product correction (APC) term, as described in Buslje et al. (2009).
}
  \item{slope}{
Numeric giving the slope of the sigmoid used to convert mutual information values to scores ranging from zero to one.
}
  \item{shift}{
Numeric giving the relative shift of the sigmoid used to convert mutual information values to scores ranging from zero to one.
}
  \item{threshold}{
Numeric specifying the score threshold at which to consider positions for pairing.  Only applicable if \code{type} is \code{"states"} or \code{"pairs"}.
}
  \item{pseudoknots}{
Integer indicating the maximum order of pseudoknots that are acceptable.  A value of \code{0} will prevent pseudoknots in the structure, whereas \code{1} (the default) will attempt to find first-order psuedoknots.  Only used if \code{type} is \code{"states"} or \code{"pairs"}.
}
  \item{weight}{
Either a numeric vector of weights for each sequence, a single number implying equal weights, or \code{NA} (the default) to automatically calculate sequence weights based on \code{myXStringSet}.
}
  \item{useFreeEnergy}{
Logical determining whether RNA free energy predictions should be incorporated along with mutual information into the secondary structure prediction.
}
  \item{deltaGrules}{
Free energy rules for all possible base pairings in quadruplets.  If NULL, defaults to pseudoenergies (\code{\link{deltaGrulesRNA}}).  Only applicable if \code{useFreeEnergies} is \code{TRUE}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{PredictDBN} employs an extension of the method described by Freyhult et al. (2005) for determining a consensus RNA secondary structure.  It uses the mutual information (\eqn{H}) measure to find covarying positions in a multiple sequence alignment.  The original method is modified by the addition of different weights for each type of base pairing and each input sequence.  The formula for mutual information between positions \eqn{i} and \eqn{j} then becomes:

\deqn{H(i,j) = \sum_{XY \in bp}^{} \left( impact(XY) \cdot f_{i,j}(XY) \cdot \log_2 \left( \frac{f_{i,j}(XY)}{f_{i}(X) \cdot f_{j}(Y)} \right) \right)}

where, \eqn{bp} denotes the base pairings A/U, C/G, and G/U; \code{impact} is their weight; \eqn{f} is the frequency of single bases or pairs weighted by the corresponding \code{weight} of each sequence.

A penalty is then added for bases that are inconsistent with pairing:

\deqn{H_{mod}(i,j) = H(i,j) + \sum_{XY \notin bp}^{} \Big( impact(XY) \cdot f_{i,j}(XY) \Big)}

Next an average product correction (Buslje et al., 2009) is applied to the matrix \eqn{H}:

\deqn{H_{APC}(i,j) = H_{mod}(i,j) - avgProdCorr \cdot \frac{\overline{H_{mod}(i,.)} \cdot \overline{H_{mod}(.,j)}}{\overline{H_{mod}(.,.)}}}

The mutual information values are then rescaled between \code{0} and \code{1} by applying a sigmoidal transformation, which is controlled by \code{shift} and \code{slope}:

\deqn{H_{final}(i,j) = \left( 1 + \exp \left( slope \cdot log_e \left( \frac{H_{APC}(i,j)}{shift \cdot H_{APC}[n]} \right) \right) \right)^{-1}}

where, \eqn{n} is the number of positions having \code{minOccupancy} divided by two (i.e., the maximum possible number of paired positions) and \eqn{H_{APC}[n]} denotes the \eqn{n^{th}} highest value in the matrix \eqn{H_{APC}}.

If \code{useFreeEnergies} is \code{TRUE}, mutual information is supplemented with a probabalistic model of folding based on \code{deltaGrules}.  That is, palindromes in each sequence are ranked by their free energy, and converted to probabilities of base pairing by assuming an exponential distribution of free energies.  This tends to improve predictive accuracy when the aligned sequences are insufficiently diverse for considerable evidence of compensatory mutations.

If \code{type} is \code{"states"} or \code{"pairs"}, the secondary structure is determined using a variant of the Nussinov algorithm similar to that described by Venkatachalam et al. (2014).  Pairings with a score below \code{threshold} are not considered during the traceback.  If \code{psuedoknots} is greater than \code{0}, paired positions are removed from consideration and the method is applied again to find pseudoknots.

In practice the secondary structure prediction is most accurate when the input alignment is of high quality, contains a wide diversity of sequences, the number of sequences is large, no regions are completely conserved across all sequences, and most of the sequences span the entire alignment (i.e., there are few partial/incomplete sequences).
}
\value{
If \code{type} is \code{"states"} (the default), then the output is a character vector with the predicted secondary structure assignment for each position in \code{myXStringSet}.  Standard dot-bracket notation (DBN) is used, where ``.'' signifies an unpaired position, ``('' and ``)'' a paired position, and successive ``[]'', ``\{\}'', and ``<>'' indicate increasing order pseudoknots.  Columns below \code{minOccupancy} are denoted by the ``-'' character to indicate that they contained too many gaps to be included in the consensus structure.

If \code{type} is \code{"pairs"}, then a matrix is returned with one row for each base pairing and three columns giving the positions of the paired bases and their \code{pseudoknot} order.

If \code{type} is \code{"evidence"}, then a matrix is returned with one row for each base pairing and three columns giving the positions of the paired bases and their respective scores (greater than or equal to \code{threshold}).  This differs from \code{type} \code{"pairs"} in that \code{"evidence"} does not perform a traceback.  Therefore, it is possible to have conflicting evidence where a single base has evidence for pairing with multiple other bases.

If \code{type} is \code{"scores"}, then a matrix of three rows is returned, where the values in a column represent the maximum score for a state in each position.  Columns sum to \code{1} if the position was above \code{minOccupancy} and \code{0} otherwise.

If \code{type} is \code{"structures"}, then the output is a list with one element for each sequence in \code{myXStringSet}.  Each list element contains a matrix of dimension 3 (each state) by the number of nucleotides in the sequence.  Columns of the matrix sum to zero where the nucleotide was located in a position that was below \code{minOccupancy}.  Otherwise, positions are considered paired if they are consistent with pairing (i.e., A/U, C/G, or G/U) in the consensus secondary structure.

When \code{type} is \code{"search"} the results are similar to \code{"structures"}, but an attempt is made to find additional secondary structure beyond positions exhibiting covariation.  First, anchors are identified as pairs of covarying positions with their score above \code{threshold}.  Next, the regions between anchors are searched for previously unidentified stem loops.  Finally, any helices are assigned a score according to their length, i.e. one minus the probability of finding that many consecutive pairs within the anchor boundaries by chance.  Hence, output \code{type} \code{"search"} will find secondary structure outside of the consensus structure shared by most sequences, and can identify secondary structure in conserved alignment regions.
}
\references{
Buslje, C., et al. (2009). Correction for phylogeny, small number of observations and data redundancy improves the identification of coevolving amino acid pairs using mutual information. Bioinformatics, \bold{25(9)}, 1125-1131.

Freyhult, E., et al. (2005). Predicting RNA Structure Using Mutual Information. Applied Bioinformatics, \bold{4(1)}, 53-59.

Venkatachalam, B., et al. (2014). Faster algorithms for RNA-folding using the Four-Russians method. Algorithms for Molecular Biology : AMB, \bold{9(1)}, 1-12.

Wright, E. S. (2020). RNAconTest: comparing tools for noncoding RNA multiple sequence alignment based on structural consistency. RNA 2020, 26, 531-540.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{PredictHEC}}
}
\examples{
# load the example non-coding RNA sequences
db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
rna <- SearchDB(db, type="RNAStringSet")

# predict the secondary structure in dot-bracket notation (dbn)
p <- PredictDBN(rna, "states") # predict the secondary structure in dbn
p # pairs are denoted by (), and (optionally) pseudoknots by [], {}, and <>

# convert the dot-bracket notation into pairs of positions within the alignment
p <- PredictDBN(rna, "pairs") # paired positions in the alignment
head(p) # matrix giving the pairs and their pseudoknot order (when > 0)

# plot an arc diagram with the base pairings
plot(NA, xlim=c(0, 1), ylim=c(0, 1),
	xaxs="i", yaxs="i",
	xlab="Alignment position", ylab="",
	bty="n", xaxt="n", yaxt="n")
ticks <- pretty(seq_len(width(rna)[1]))
axis(1, ticks/width(rna)[1], ticks)
rs <- c(seq(0, pi, len=100), NA)
r <- (p[, 2] - p[, 1] + 1)/width(rna)[1]/2
r <- rep(r, each=101)
x <- (p[, 1] + p[, 2])/2/width(rna)[1]
x <- rep(x, each=101) + r*cos(rs)
y <- r*sin(rs)/max(r, na.rm=TRUE)
lines(x, y, xpd=TRUE)

# show all available evidence of base pairing
p <- PredictDBN(rna, "evidence") # all pairs with scores >= threshold
head(p) # matrix giving the pairs and their scores

# determine the score at every alignment position
p <- PredictDBN(rna, "scores") # score in the alignment
p["(", 122] # score for left-pairing at alignment position 122
p[")", 260] # score for right-pairing at alignment position 260

# find the scores individually for every sequence in the alignment
p <- PredictDBN(rna, "structures") # scores per sequence
p[[1]][, 1] # the scores for the first position in the first sequence
p[[2]][, 10] # the scores for the tenth position in the second sequence
# these positional scores can be used as shades of red, green, and blue:
BrowseSeqs(rna, patterns=p) # red = unpaired, green = left-pairing, blue = right
# positions in black are not part of the consensus secondary structure

# search for additional secondary structure between the consensus pairs
p <- PredictDBN(rna, "search") # scores per sequence after searching
BrowseSeqs(rna, patterns=p) # red = unpaired, green = left-pairing, blue = right
# note that "search" identified many more pairings than "structures"
}
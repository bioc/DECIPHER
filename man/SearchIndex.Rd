\name{SearchIndex}
\alias{SearchIndex}
\title{
Search an inverted index
}
\description{
Searches an inverted index for homologous sequences.
}
\usage{
SearchIndex(myXStringSet,
            invertedIndex,
            type = "one",
            minScore = NA,
            scoreOnly = FALSE,
            sepCost = -0.4,
            gapCost = -2.5,
            maskRepeats = TRUE,
            maskLCRs = TRUE,
            processors = 1,
            verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of query (unaligned) sequences to use as the pattern that will be queried.
}
  \item{invertedIndex}{
An object of class \code{InvertedIndex} compatible with the class of \code{myXStringSet}.
}
  \item{type}{
Character string determining which search hits are returned.  This should be (an unambiguous abbreviation of) one of \code{"all"}, \code{"one"}, or \code{"top"}.  The default \code{type} (\code{"one"}) will return the best hit per target sequence in the \code{invertedIndex}.  Setting \code{type} to \code{"all"} or \code{"top"} will return all hits or only a single top scoring hit, respectively.
}
  \item{minScore}{
Numeric specifying the minimum score of hits to return.  The default (\code{NA}) will automatically determine \code{minScore} from the size of the \code{invertedIndex} to only return statistically significant hits, such that the number of expected hits with at least that \code{Score} by chance is \code{exp(minScore - Score)}.  Therefore, the default \code{minScore} will only allow hits to be reported that are expected to occur once (i.e., \code{exp(0)}) by chance per query sequence search.
}
  \item{scoreOnly}{
Logical determining whether to return only the hits and their scores or also the \code{Position} of k-mer hits.
}
  \item{sepCost}{
Numeric giving the penalty applied to sequence positions separating neighboring k-mer hits.
}
  \item{gapCost}{
Numeric providing the penalty applied to the minimum number of implied inserted or deleted positions (i.e., gaps) separating neighboring k-mer hits.
}
  \item{maskRepeats}{
Logical specifying whether to mask repeats when searching for hits.
}
  \item{maskLCRs}{
Logical indicating whether to mask low complexity regions when searching for hits.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The \code{invertedIndex} is searched for all umasked k-mers shared by \code{myXStringSet}, and the set of statistically significant matches is returned.  By default, \code{SearchIndex} returns the highest scoring match per target sequence (above \code{minScore}), but it is possible to set \code{type} to output \code{"all"} matches or only the \code{"top"} scoring match per query sequence in \code{myXStringSet}.  A \code{data.frame} is returned with (by default) or without the \code{Position}(s) of matches, depending on the value of \code{scoreOnly}.  The resulting \code{Score} approximates the negative log probability of observing the set of k-mer matches by chance between the query and target sequence.  For example, a \code{Score} of 7 should be found every approximately 1,097 (i.e., \code{exp(7)} searches of (random) target sequences by chance.  Since typically more than one target sequence is searched, the default \code{minScore} (\code{NA}) is set such that hits are only returned when they are expected fewer than once per query search in a database the size of \code{invertedIndex}, thereby correcting for multiple target searches.
}
\value{
A \code{data.frame} is returned with dimensions with columns \code{Pattern}, \code{Subject}, \code{Score}, and (optionally) \code{Position}.  The \code{Pattern} is the index of the sequence in \code{myXStringSet} and the \code{Subject} is the index of the sequence in the set used to build the \code{invertedIndex}.  Each row contains a hit with \code{Score} meeting the \code{minScore}.  If \code{scoreOnly} is \code{FALSE} (the default), the \code{Position} column contains a list of matrices with four rows: start/end positions of k-mer hits in the \code{Pattern} and start/end positions of k-mer hits in the \code{Subject}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{IndexSeqs}}, \code{\link{AlignPairs}}
}
\examples{
# import target sequences
fas <- system.file("extdata", "PlanctobacteriaNamedGenes.fas.gz", package="DECIPHER")
target <- readAAStringSet(fas)

# build an inverted index
index <- IndexSeqs(target, K=6L)
index

# import query sequences
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
query <- translate(dna)

# search the index, using the defaults
hits <- SearchIndex(query, index)
head(hits)
dim(hits) # number of hits

# search the index, requesting only the top hits
tophits <- SearchIndex(query, index, type="top")
head(tophits)
dim(tophits) # number of hits
tophits$Position[[1]] # query/target k-mer positions supporting first hit

# search the index, requesting only the score for all hits
allhits <- SearchIndex(query, index, type="all", scoreOnly=TRUE)
head(allhits)
dim(allhits) # number of hits
}

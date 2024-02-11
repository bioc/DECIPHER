\name{SearchIndex}
\alias{SearchIndex}
\title{
Search an inverted index
}
\description{
Searches an inverted index for homologous sequences.
}
\usage{
SearchIndex(pattern,
            invertedIndex,
            subject=NULL,
            type = "one",
            minScore = NA,
            scoreOnly = FALSE,
            sepCost = -0.4,
            gapCost = -2.5,
            maskRepeats = TRUE,
            maskLCRs = TRUE,
            dropScore = -10,
            processors = 1,
            verbose = TRUE)
}
\arguments{
  \item{pattern}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of query (unaligned) sequences to use as the pattern that will be queried.
}
  \item{invertedIndex}{
An object of class \code{InvertedIndex} compatible with the class of \code{pattern}.
}
  \item{subject}{
The \code{XStringSet} object used to build the \code{invertedIndex} or \code{NULL} (the default) to skip the extension of k-mer matches.  (See details section below.)
}
  \item{type}{
Character string determining which search hits are returned.  This should be (an unambiguous abbreviation of) one of \code{"all"}, \code{"one"}, or \code{"top"}.  The default \code{type} (\code{"one"}) will return the best hit per target sequence in the \code{invertedIndex}.  Setting \code{type} to \code{"all"} or \code{"top"} will return all hits or only a single top scoring hit, respectively.
}
  \item{minScore}{
Numeric specifying the minimum score of hits to return.  The default (\code{NA}) will automatically determine \code{minScore} from the size of the \code{invertedIndex} to correct for multiple testing (i.e., searching more than one target sequence).
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
  \item{dropScore}{
Numeric giving the decrease from maximum score required to stop extending k-mer matches.  Only applicable when \code{subject} is not \code{NULL}.  Values closer to zero will increase speed, potentially at the expense of sensitivity.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The \code{invertedIndex} is searched for all umasked k-mers shared by \code{pattern}, and the set of matches meeting the \code{minScore} is returned.  By default, \code{SearchIndex} returns the highest scoring match per target sequence (above \code{minScore}), but it is possible to set \code{type} to output \code{"all"} matches or only the \code{"top"} scoring match per query sequence in \code{pattern}.  A \code{data.frame} is returned with (by default) or without the \code{Position}(s) of matches, depending on the value of \code{scoreOnly}.

If the set of \code{subject} sequences is provided (i.e., not \code{NULL}), then k-mer matches are extended to increase search sensitivity.  Extension proceeds to the left and right of each k-mer match until another match is encountered or the score falls below \code{dropScore}.  This can decrease search speed, depending on \code{dropScore}, but may help to find more distant matches.  The \code{Score} of any hits is defined by their log-odds regardless of whether \code{subject} is provided.
}
\value{
A \code{data.frame} is returned with dimensions with columns \code{Pattern}, \code{Subject}, \code{Score}, and (optionally) \code{Position}.  The \code{Pattern} is the index of the sequence in \code{pattern} and the \code{Subject} is the index of the sequence in the set used to build the \code{invertedIndex}.  Each row contains a hit with \code{Score} meeting the \code{minScore}.  If \code{scoreOnly} is \code{FALSE} (the default), the \code{Position} column contains a list of matrices with four rows: start/end positions of k-mer hits in the \code{Pattern} and start/end positions of k-mer hits in the \code{Subject}.
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
index <- IndexSeqs(target, K=6)
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

# include the target sequences to improve sensitivity (but slower)
allhits <- SearchIndex(query, index, target, type="all", scoreOnly=TRUE)
head(allhits)
dim(allhits) # number of hits
}

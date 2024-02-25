\name{MapCharacters}
\alias{MapCharacters}
\title{
Map Changes in Ancestral Character States
}
\description{
Maps character changes on a phylogenetic tree containing reconstructed ancestral states.
}
\usage{
MapCharacters(x,
              refPositions = seq_len(nchar(attr(x, "state")[1])),
              labelEdges = FALSE,
              type = "dendrogram",
              chars = LETTERS,
              ignoreIndels = TRUE)
}
\arguments{
  \item{x}{
An object of class \code{dendrogram} with \code{"state"} attributes for each node.
}
  \item{refPositions}{
Numeric vector of reference positions in the original sequence alignment.  Only changes at \code{refPositions} are reported, and state changes are labeled according to their position in \code{refPositions}.
}
  \item{labelEdges}{
Logical determining whether to label edges with the number of changes along each edge.
}
  \item{type}{
Character string indicating the type of output desired.  This should be (an abbreviation of) one of \code{"dendrogram"}, \code{"table"}, or \code{"both"}.  (See value section below.)
}
  \item{chars}{
Character vector specifying the characters to consider in state changes at each site.  The default (\code{LETTERS}) is to consider any upper case letter. Alternatively, \code{chars} could be \code{AA_STANDARD}, \code{DNA_BASES}, or \code{RNA_BASES}.
}
  \item{ignoreIndels}{
Logical specifying whether to report insertions and deletions (indels).  If \code{TRUE} (the default), only substitutions of one state with another are reported.
}
}
\details{
Ancestral state reconstruction affords the ability to identify character changes that occurred along edges of a rooted phylogenetic tree.  Character changes are reported according to their index in \code{refPositions}.  If \code{ignoreIndels} is \code{FALSE}, adjacent insertions and deletions are merged into single changes occurring at their first position.  The table of changes can be used to identify parallel, convergent, and divergent mutations.
}
\value{
If \code{type} is \code{"dendrogram"} (the default) then the original \code{dendrogram} \code{x} is returned with the addition of \code{"change"} attributes on every edge except the root.  If \code{type} is \code{"table"} then a sorted \code{table} of character changes is returned with the most frequent parallel changes at the beginning.  If \code{type} is \code{"both"} then a \code{list} of length 2 is provided containing both the \code{dendrogram} and \code{table}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{TreeLine}}
}
\examples{
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

d <- DistanceMatrix(dna, type="dist", correction="JC")
tree <- TreeLine(myDistMatrix=d,
                 method="NJ",
                 type="dendrogram",
                 myXStringSet=dna,
                 reconstruct=TRUE)

out <- MapCharacters(tree,
                     labelEdges=TRUE,
                     type="both",
                     chars=DNA_BASES)

# plot the tree with defaults
tree <- out[[1]]
plot(tree, horiz=TRUE) # edges show number of changes

# color edges by number of changes
maxC <- 200 # changes at maximum of color spectrum
colors <- colorRampPalette(c("black", "darkgreen", "green"))(maxC)
colorEdges <- function(x) {
   num <- attr(x, "edgetext") + 1
   if (length(num)==0)
       return(x)
   if (num > maxC)
       num <- maxC
   attr(x, "edgePar") <- list(col=colors[num])
   attr(x, "edgetext") <- NULL
   return(x)
}
colorfulTree <- dendrapply(tree, colorEdges)
plot(colorfulTree, horiz=TRUE, leaflab="none")

# look at parallel changes (X->Y)
parallel <- out[[2]]
head(parallel) # parallel changes

# look at convergent changes (*->Y)
convergent <- gsub(".*?([0-9]+.*)", "\\\\1", names(parallel))
convergent <- tapply(parallel, convergent, sum)
convergent <- sort(convergent, decreasing=TRUE)
head(convergent)

# look at divergent changes (X->*)
divergent <- gsub("(.*[0-9]+).*", "\\\\1", names(parallel))
divergent <- tapply(parallel, divergent, sum)
divergent <- sort(divergent, decreasing=TRUE)
head(divergent)

# plot number of changes by position
changes <- gsub(".*?([0-9]+).*", "\\\\1", names(parallel))
changes <- tapply(parallel, changes, sum)
plot(as.numeric(names(changes)),
     changes,
     xlab="Position",
     ylab="Total independent changes")

# count cases of potential compensatory mutations
compensatory <- dendrapply(tree,
    function(x) {
        change <- attr(x, "change")
        pos <- as.numeric(gsub(".*?([0-9]+).*", "\\\\1", change))
        e <- expand.grid(seq_along(pos), seq_along(pos))
        e <- e[pos[e[, 1]] < pos[e[, 2]],]
        list(paste(change[e[, 1]], change[e[, 2]], sep=" & "))
    })
compensatory <- unlist(compensatory)
u <- unique(compensatory)
m <- match(compensatory, u)
m <- tabulate(m, length(u))
compensatory <- sort(setNames(m, u), decreasing=TRUE)
head(compensatory) # ranked list of concurrent mutations
}

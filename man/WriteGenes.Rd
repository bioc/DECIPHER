\name{WriteGenes}
\alias{WriteGenes}
\title{
Write Genes to a File
}
\description{
Writes predicted genes to a file in GenBank (gbk) or general feature format (gff).
}
\usage{
WriteGenes(x,
           file = "",
           format = "gbk",
           append = FALSE)
}
\arguments{
  \item{x}{
An object of class \code{Genes}.
}
  \item{file}{
A connection or a character string naming the file path where the tree should be exported.  If "" (the default), the tree is printed to the standard output connection, the console unless redirected by sink.
}
  \item{format}{
Character specifying \code{"gbk"} or \code{"gff"} output format.
}
  \item{append}{
Logical indicating whether to append to an existing \code{file}.  Only applicable if \code{file} is a character string.  If \code{FALSE} (the default), then the file is overwritten.
}
}
\details{
\code{WriteGenes} will write a \code{"Genes"} object to a GenBank (if \code{format} is \code{"gbk"}) or general feature format (if \code{format} is \code{"gff"}) file.
}
\value{
\code{NULL}.
}
\author{
Erik Wright \email{Erik.Wright@bcm.edu}
}
\seealso{
\code{\link{ExtractGenes}}, \code{\link{FindGenes}}, \code{\link{Genes-class}}
}
\examples{
# import a test genome
fas <- system.file("extdata",
	"Chlamydia_trachomatis_NC_000117.fas.gz",
	package="DECIPHER")
genome <- readDNAStringSet(fas)

x <- FindGenes(genome)
WriteGenes(x[1:10,], format="gbk")
WriteGenes(x[1:10,], format="gff")
}

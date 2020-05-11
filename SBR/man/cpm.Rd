\name{cpm}
\alias{cpm}
\title{ Counts Per Million reads (CPM) from counts}
\description{
  This function maps read count matrics  to CPM matrices
}
\usage{
cpm(counts)
}
\arguments{
  \item{counts}{ Matrix of counts.}
}
\value{ A matrix of CPM values.}
\references{See https://support.bioconductor.org/p/91218/ and https://www.biostars.org/p/335187/.}
\author{ Tom Radivoyevitch}

\examples{
##Toy example from this video https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
library(SBR) 
(X=matrix(c(10,20,5,0,12,25,8,0,30,60,15,1),ncol=3))
cpm(X)/1e5
}
\keyword{ math }% at least one, from doc/KEYWORDS

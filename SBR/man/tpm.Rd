\name{tpm}
\alias{tpm}
\title{ Transcripts per million from counts}
\description{
  This function maps read count matrics and gene lengths to TPM matrices
}
\usage{
tpm(counts, lengths)
}
\arguments{
  \item{counts}{ Matrix of counts.}
  \item{lengths}{ Vector of gene lengths.}
}
\value{ A matrix of TPM values.}
\references{See https://support.bioconductor.org/p/91218/ and https://www.biostars.org/p/335187/.}
\author{ Tom Radivoyevitch}

\examples{
##Toy example from this video https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
library(SBR) 
(X=matrix(c(10,20,5,0,12,25,8,0,30,60,15,1),ncol=3))
(lens=c(2,4,1,10))
tpm(X,lens)/1e5
}
\keyword{ math }% at least one, from doc/KEYWORDS

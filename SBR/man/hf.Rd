\name{hf}
\alias{hf}
\docType{data}
\title{Header for Ferrets (hf)}
\description{
  Information about ferret genes in GSE147507.  
}
\usage{data(hf)}
\format{
  A data frame with the following columns. Only rows with non-blank Symbols are included.
  \describe{
    \item{\code{id}}{Gene id, e.g. ENSMPUG00000000007}
    \item{\code{Symbol}}{Gene symbol (short acronym). }
    \item{\code{Length}}{Length of gene. Used by tpm().}
    \item{\code{desc}}{Short description. }
  }
}



\keyword{datasets}

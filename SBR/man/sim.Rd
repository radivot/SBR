\name{sim}
\alias{sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a model of S3 class SBMLR}
\description{
  This function simulates a model given report times and optional modulators. 
  It uses \code{lsoda} of the deSolve package.
}
\usage{
sim(model, times, modulator=NULL,X0=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ The S3 model object to be simulated. Initial conditions are passed through this object.}
  \item{times}{ The sequence of time points to be sampled and provided as rows of the output matrix.}
  \item{modulator}{ Null if there are no modulators (default), 
  a vector of numbers if there are steady state Vmax modulators,  
  and a list of interpolating functions if there are time course Vmax modulators. }
  \item{X0}{ Override model initial conditions in simulations, particularly piece-wise perturbation simulations.} 
  \item{...}{ To pass extra args such as event data frames to deSolve.}
}
\details{
  This is a wrapper for ode. }
\value{ The data frame output that comes out of \code{ode}.}
\references{For the folate cycle example given below: Morrison PF, Allegra CJ: Folate cycle kinetics in human 
  breast cancer cells. JBiolChem 1989, 264(18):10552-10566.}
\author{ Tom Radivoyevitch}
\note{ Rules are implemented through time varying boundary conditions updated at each time point 
as a side effect within the (now internal) function \code{fderiv}. }

\examples{
##---- The following perturbs PRPP from 5 to 50 uM in Curto et al.'s model.
library(SBR) 
curto=readSBMLR(file.path(system.file(package="SBR"), "models/curto.r"))  
(dPRPP10 <- data.frame(var = "PRPP", time = 0, value = 10,method = "mult"))
(out=sim(curto,times=seq(-20,70,1),events = list(data = dPRPP10) ) )
plot(out,which=c("PRPP","den","IMP","HX","Gua","aprt","XMP","Xa","UA"))
}
\keyword{ math }% at least one, from doc/KEYWORDS

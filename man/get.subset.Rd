\name{get.subset}
\alias{get.subset}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Return the subset of controlling species.
}

\description{
Obtain the smallest subset that contains the species with the smallest alpha-th quantile with at least pre-specified certainty.
}

\usage{

get.subset(pis, Pstar) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{pis}{
A vector of length K; k-th element must be a posterior probability that species k has the smallest alpha-th quantile.
}
\item{Pstar}{
The pre-specified posterior probability of certainty. 
The selected subset is the smallest one containing the species with the smallest alpha-th quantile 
with at least \code{Pstar} posterior probability.
}
}


\value{
\item{Size}{The size of the selected subset}
\item{Subset}{The index of the species selected in the subset}
\item{sumPisSubset}{The posterior probability that the selected subset contains the species with the smallest alpha-th quantile}
}

\references{
Yumi Kondo, James V. Zidek, Carolyn Taylor, and Constance van Eeden
"BAYESIAN SUBSET SELECTION PROCEDURES WITH AN APPLICATION TO LUMBER STRENGTH PROPERTIES"
}

\author{
Yumi Kondo
}


\seealso{

  \code{\link{CDF.eBayes}},
  \code{\link{MCMC.DPw}},
  \code{\link{ASTM}},
  \code{\link{get.subset}},
  \code{\link{useSamp}}
}


\examples{

pis <- c(0.1,0.2,0.3,0.4)
Pstar <- 0.7
get.subset(pis, Pstar) 
}

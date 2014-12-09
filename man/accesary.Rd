\name{useSamp}
\alias{useSamp}
\alias{simdata}
\alias{Nuniq}
\alias{muinv}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Discrad burn-in from the MCMC output and thin the remaining output
}

\description{
Return a vector that indicates which MCMC output to keep.
}

\usage{

useSamp(thin, burnin, B) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{


  \item{thin}{
MCMC samples are thinned at every \code{thin} iterations. 
}
\item{burnin}{
first \code{burnin} MCMC samples are discarded.
}

  \item{B}{
 The original MCMC sample size.
}
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


B <- 5500
burnin <- 500
thin <- 2
useSamp(thin, burnin, B) 


}

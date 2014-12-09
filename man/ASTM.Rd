\name{ASTM}
\alias{ASTM}
\alias{ASTMTL}
\alias{ASTMMed}
\alias{TL}

\title{
ASTM D1990 procedures for finding the subset of controlling species.
}

\description{
ASTM D1990 procedures vary for different characteristic values (5th percentile, median).
However, they are all based on sequential hypothesis testings, and the idea is to combine all "weak" species whose samples are statistically indistinguishable with respect to the quantile of relevance as the subset of controlling species.
See reference for the details.
}

\usage{
ASTM(data, alpha = c(0.05, 0.5), Sig = 0.01, detail = FALSE) 
}

%- maybe also 'usage' for other objects documented here.


\arguments{
\item{data}{
data.frame containing variables named \code{resp} and \code{S}.
\code{data$resp} must contains the strength measures.
\code{data$S} must contain the labels to indicate to which species the strengths belong.
}
\item{alpha}{
The characteristic values. It must be either 0.05 (5th percentile) or 0.5 (median).
With 0.05, the sequential hypothesis test is based on the chi-square procedure.
With 0.5, the sequential hypothesis test is based on the nonparametric analysis of variance that test the equality of medians.
}
\item{Sig}{
The significant level. It should be 0.01.  
}
\item{detail}{Internal use only}
}




\details{


See references.}

\value{
\item{subset.out}{Selected subset of controlling species}
}


\references{
ASTM Standard D1990 (2007). Standard practice for establishing allowable properties for
visually-graded dimension lumber from in-grade tests of full-size specimens. Technical
Report DOI: 10.1520/D1990-07, ASTM International, DOI: 10.1520/D1990-07.

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
## generate strength dataset from 5 species
## Tk ~ weibull(shape=k*0.5scale=5)
m <- 100
K <- 5

resp <- NULL
for (ik in 1 : K) resp <- c(resp,rweibull(m,shape=4,scale=ik*0.5))
species <- factor(rep(1:K,each=m))
dat <- data.frame(resp=resp,S=species)
alpha <- 0.05
## The selected subset
ASTM(data=dat, alpha = 0.05) 
}

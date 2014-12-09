\name{CDF.eBayes}
\alias{CDF.eBayes}
\alias{nu}
\alias{dweib}
\alias{sprob}
\alias{H}
\alias{pr.discrete}
\alias{F}
\alias{F_inv}
\alias{beliefToW}
\title{An empirical Bayes approach to the subset selection problem.}

 \description{
For each species' strength, our empirical Bayes approach models the CDF of the strengths with Dirichlet Process (DP) with the Weibull base prior.In this model, the precision parameter, and shape and scale parameters of the Weibull are fixed hyperparameters.
The shape and scale parameters are estimated via the empirical Bayes procedure:
They are replaced with the maximum likelihood estimates with Weibull likelihood on strengths.
The precision parameter must be pre-specified.

Under this model, the posterior probability that species k has the smallest alpha-th quantile of strength is estimated in two steps:

STEP 1: Estimate shape and scale parameters via the empirical Bayes procedure for each species, and compute the posterior CDF of the alpha-th quantile of strengths at pre-specified grids for each species. This step is done by \code{CDF.eBayes}

STEP 2: Given the output from \code{cdf} \code{sprob} computes the posterior probability that species k has the smallest alpha-th quantile for each species via the quasi-Monte Carlo method.
}

 \usage{
 CDF.eBayes(data, c0s, c0.belief.G0 = NULL, ts = NULL, length.out = 1e+06, 
            qs = c(1e-10, 1 - 1e-10), alpha = 0.05, C = TRUE) 
 sprob(outputCDF, N.MC, ExcludeS = NULL)
 }


 \arguments{

--- CDF.eBayes ---

 \item{data}{
 A list of size K; each item of the list contains a vector of strength measures from a single species.
 }

 \item{c0s}{
 A prior for precision parameter \emph{v} of the Dirichlet Process.
 If it is a vector, it must have length K, and each element contains the prior precision parameter for each species.
 If it is a scalar, the same precision parameter is used for all species.
 }

\item{c0.belief.G0}{
 A scalar value between 0 and 1.
 Instead of specifying c0s, one can specify the weight on the base prior relative to the weight on empirical CDF of the dataset. 
}

\item{ts}{ 
 A vector containing equally spaced grid of points where the posterior CDFs of alpha-th quantiles are evaluated.
 If \code{ts} is not specified then \code{qs} and \code{length.out} must be specified.
}

\item{length.out}{
  The number of points to evaluate posterior CDF of the alpha-th quantiles.
  Necessary if ts is not specified.
} 

\item{qs}{
A vector of length 2 containing the values between 0 and 1.

If ts is not specified, then the grids are specified as follows:

[1] compute \code{qs[1]} quantile of the posterior CDF of the alpha-th quantile for each species.
The minimum of such values across species is selected as the smallest value at which the posterior CDF of alpha-th quantile is evaluated. 

[2] compute \code{qs[2]} quantile of the posterior CDF of the alpha-th quantile for each species.
The maximum of such values across species is selected as the largest value at which the posterior CDF of alpha-th quantile is evaluated. 

Finally, \code{length.out} many equally spaced points between the points selected from [1] and [2] as grids at which posterior CDF of alpha-th quantiles are evaluated.
	
} 

\item{alpha}{ The size of the quantile of interest for e.g. 0.05.}	
\item{C}{Internal use only}



--- sprob ---
\item{outputCDF}{
The first element of the list output of a function i.e., \code{out.CDF.eBayes$outputCDF}
where \code{out.CDF.eBayes} is an output of \code{CDF.eBayes}.
}

\item{N.MC}{
The number of quasi-Monte Carlo samples to approximate the posterior probability that species k has the smalles alpha-th quantile.
}

\item{ExcludeS}{Internal use only}

 }

 \value{

--- CDF.eBayes ---

 \item{outputCDF}{
  \emph{k+1} by \emph{M} matrix.
  The first row contains the grid of points at which the posterior CDF of the alpha-th quantiles are evaluated. i.e., \code{ts}
 }

\item{MLE}{
  The maximum likelihood estimates of scale and shape parameters for Weibull base prior of each species.
}

 --- sprob ---
\item{ALL}{
The estimated posterior probability that species k has the smalles alpha-th quantile for each species.
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
 \dontrun{





## generate strength dataset from 5 species
## Tk ~ weibull(shape=0.5*k,scale=5)
m <- 100
K <- 5
datlist <- rep(list(NULL),K)
names(datlist) <- paste("species",1:K,sep="")
xlim <- c(0,5)
par(mfrow=c(K,1),mar=c(2,2,1,1))
for (ik in 1 : K){
  datlist[[ik]] <- rweibull(m,shape=4,scale=ik*0.5)
  hist(datlist[[ik]],xlim=xlim,
       border="white",col="gray",
       breaks=seq(xlim[1],xlim[2]*2,0.1),
       main=names(datlist)[ik])
}

## ----------------------------------------------------------------------- ##
## STEP 1: compute the posterior CDF of the alpha-th quantile of strengths ##
## ----------------------------------------------------------------------- ##
## 50% weight on Weibull base prior 
c0.belief.G0 <- 0.5
length.out <- 1E+6 ## the number of grids of points
alpha <- 0.05 
nonCDFs <- CDF.eBayes(data=datlist, c0.belief.G0=c0.belief.G0,
                      alpha = alpha,length.out=length.out)


## -------------------------------------------------------------------------------------------- ##
## STEP 2: compute the posterior probability that species k has the smallest alpha-th quantile ##
## -------------------------------------------------------------------------------------------- ##
N.MC <- 5000
re <- sprob(N.MC=N.MC,outputCDF=nonCDFs$outputCDF)

## The estimated posterior probability that species k has the smallest alpha-th percentile
re$All

## Find the smallest subset that contains the smallest alpha-th percentile
## with posterior probability at least Pstar
Pstar <- 0.99
get.subset(pis=re$All,Pstar=Pstar)
 }
 }

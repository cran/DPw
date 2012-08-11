\name{DPfitAllUnif}
\alias{DPfitAllUnif}
\alias{simdata}
\alias{get.subset.p}
\alias{getden}
\alias{dWEI2}
\alias{rWEI2}
\alias{muinv}
\alias{plotigamma}
\alias{plotKwei}
\alias{plotgamma}
\alias{cdf}
\alias{nu}
\alias{dweib}
\alias{sprob}
\alias{H}
\alias{getM}
\alias{useSamp}
\alias{F}
\alias{F_inv}
\alias{Nuniq}
\alias{ASTM}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Fit the Weibull Dirichlet process mixture.
}

\description{

The purpose of the functions implemented in this package is to find the subset containing the weakest species with some acceptable certainty;
Suppose we have sampled some quality measures T_{k,j} (>=0) from K species (k=1...K, j=1,...,m_j= sample size of jth species),
and wish to find the species with the smallest alpha quantile (say, 5th quantile) of T (namely the weakest species). 
Due to the smalle sample size (m_k) or large number of species to be considered (K), we might not be able to select a single species as a weakest species with acceptably large certainty. Instead of selecting one specie with a low certainty, the subset selection proceeds to select a subset of species which contains the weakest specie with predetermined certainty. Our subset selection procedures are semiparametric to avoid unrealistic assumptions.

Our subset selection procedure first (1) fits a nonparametric model to estimate the density of each \emph{T_k};fit a DP mixture of weibull distributions to each \emph{T_k}, separately. This step is done by the function \code{DPfitAllUnif}.

Our MCMC sampler is designed to return multiple copies of the posterior density of \emph{T_k}.  
Hence (2) the second step is to invert the posterior densities to obtain the posterior sample of the alpha quantiles.
This step is done by the function \code{getdens}.

The steps (1) and (2) are performed for each species separately then the final step (3) evaluates the posterior probability of specie k to be the weakest one by estimating the probability by Monte Carlo method.
This step is done by the function \code{get.subset.p}

Once the probabilities \emph{pi_k} that a specie k to be the weakest are computed, we order the species in the decreasing manner of the estimated probabilities \emph{pi_k} then choose the first few species such that the sum of \emph{pi_k} of selected species exceeds the predetermined certainty.

See the examples to learn the use of these functions.
See reference for more details of the subset selection procedure.

}

\usage{

DPfitAllUnif(ts,B=20000,
             aphi=2,    
             loc.phi=10, 
             aGam=2,
             loc.Gam=10,
             minD=0.01,
             maxD=5,
             burnin=ceiling(B/4),
             M = NULL,
             printFreq = B, 
             epsilonM=0.01
             )


}

%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{ts}{
 A vector of length sample size, containing a single specie dataset.
}

  \item{B}{
 A positive integer for MCMC sample size.
}


  \item{aphi}{
The hyperparameter for \emph{phi}.
}

  \item{loc.phi}{
The hyperparameter for \emph{phi}.
}


\item{aGam}{
 The hyperparameter for \emph{gam}. 
}

\item{loc.Gam}{
The hyperparameter for \emph{gam}. 
}

  \item{minD}{
The hyperparameter for \emph{V}.
The minimum value of \emph{V}.
}
  \item{maxD}{
The hyperparameter for \emph{V}.
The maximum value of \emph{V}.
}

\item{burnin}{
%
A scalar for burnin.
%
The \emph{L} samples of size \code{B}-\code{burnin} from the posterior predictive distribution are drawn.
}
\item{M}{
The turncation value.
If it is \code{NULL}, the function select the \code{M} that assigns the probability less than \emph{epsilonM} on the fintal probability given the maximum value of \code{v}.
}

\item{printFreq}{
Integer. The Gibbs iteration is counted and printed at every \code{printFreq}th iteration.
If (\code{printFreq > B}), then no printing statement is made.
}
\item{epsilonM}{
Necessary if M=NULL
}


}

\details{

Given a sampled vector T, the functions \code{DPfitAllUnif} and \code{getden} essentially fit the infinite mixture of weibull distributions, and estimate the density of T. At pre-specified grid of points t, the interval estimate of density is obtained.
By inverting the estimated density, the function also estimate the distribution of alpha quantile of T.


The function \code{DPfitAllUnif} sample from the posterior distribution of the hierarchical model of DP mixture of weibull distribution.
The hierarchy of the model is;

 \emph{t_j} | \emph{beta_j} \emph{Klambda_j} ~ Weibull( shape=\emph{beta_j},scale=\emph{lambda_j}), j= 1,...,N

 \emph{beta_j}, \emph{lambda_j} | \emph{G} ~ \emph{G}, j= 1,...,N

 \emph{G} ~ DP(\emph{V,G0})
 where \emph{G0} =  unif(\emph{beta}; \code{min}=0,\code{max}=\emph{phi}) unif(\emph{lambda}; \code{min}=0,\code{max}=\emph{Gam})

 \emph{gam}     ~ pareto(\code{shape}=\emph{aphi},\code{location}=\emph{loc.phi}) 

 \emph{phi}     ~ pareto(\code{shape}=\emph{aGam},\code{location}=\emph{loc.Gam}) 

 \emph{V}       ~ unif(\emph{beta}; \code{min}=\emph{minD},\code{max}=\emph{maxD}) 

Given the output of \code{DPfitAllUnif}, function \code{getden} returns the interval estimate of the density of \emph{t_j} evaluated at the grid of points \emph{xs}, and the distribution of the \emph{alpha}th quantile of \emph{t_j}.


Given \emph{T_k} from more than one specie, the output of \code{getden} can be used for selecting the subset which contains the weakest specie with some predetermined certainty \emph{P*}.
The function \code{get.subset.p} 


}

\value{



\item{S}{
A B by N matrix, containing the posterior cluster labels of each observation.
} 

\item{betas}{
A B by N matrix, containing the posterior betas of each observation.
}

\item{lambdas}{
A B by N matrix, containing the posterior Klambdas of each observation.
}

\item{betas_uniq}{
A matrix of B by M, containing the unique values of the M sampled beta at each draw.
}

\item{lambdas_uniq}{
A matrix of B by M, containing the unique values of the M sampled lambdas at each draw.
}

\item{D}{
A vector of length B, containing the sampled precision parameter.
}

\item{phi}{
A vector of length B.
}

\item{gam}{
A vector of length B.
}

\item{can}{The proposal variance of the Metropolis hasting algorithms for beta, lambda and D}

\item{AR}{The acceptance rates of the Metropolis hasting algorithms for beta, lambda and D}


}
\references{
Kondo, Yumi and Zidek, James V. "Bayesian nonparametric subset selection procedures with Weibull components"
}

\author{
Yumi Kondo
}

% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{

%}

\examples{


\dontrun{



## generate sample of size 360 from mixture of two weibull distribution
ts1 <- c(rweibull(80,shape=4.726,scale=5.8), rweibull(20,shape=6.4,scale=10.803))
ts2 <- c(rweibull(20,shape=4.726,scale=5.8), rweibull(80,shape=6.4,scale=10.803))

## parameters for the 
B <- 15000
lts <- 10000
xs <- seq(0.01,20,length.out=lts)
burnin <- 5000
thin <- 5
useSample <- useSamp(thin=thin,burnin=burnin,B=B)

## Sample from the posterior distribution given ts1
re1 <- DPfitAllUnif(ts=ts1,B=B,burnin=burnin)

## Sample from the posterior distribution given ts2
re2 <- DPfitAllUnif(ts=ts2,B=B,burnin=burnin)
  
##====== compute the density curve =========
## posterior estimates of the density of T given ts1
dens1 <- getden(xs=xs,
               weightS=re1$weightS[useSample,],
               Klambdas=re1$lambdas_uni[useSample,],
               betas=re1$betas_uni[useSample,],
               alpha = 0.05   
               )
mat1 <- apply(dens1[[1]],2,quantile,prob=c(0.025,0.5,0.975))
qq1 <- quantile(dens1$quantiles,prob=c(0.025,0.5,0.975))
mat1[2,] <- colMeans(dens1[[1]])


## posterior estimates of the density of T given ts2
dens2 <- getden(xs=xs,
               weightS=re2$weightS[useSample,],
               Klambdas=re2$lambdas_uni[useSample,],
               betas=re2$betas_uni[useSample,],
               alpha = 0.05   
               )
mat2 <- apply(dens2[[1]],2,quantile,prob=c(0.025,0.5,0.975))
qq2 <- quantile(dens2$quantiles,prob=c(0.025,0.5,0.975))
mat2[2,] <- colMeans(dens2[[1]])


#### ============= Assess the result ===================
par(mar=c(3,4,1,1),mfrow=c(1,2))
hist(ts1,
     breaks=seq(0,max(ts1)*1.3,0.5),
     prob=TRUE,
     ylim=c(0,0.5),
     xlab="",main="",
     col="gray"
     )
points(xs,mat1[2,],type="l",col="blue")
points(xs,mat1[3,],type="l",lty=2,col="blue")
points(xs,mat1[1,],type="l",lty=2,col="blue")
abline(v=qq1[2],lwd=2,col="blue")
abline(v=qq1[1],lty=2,lwd=2,col="blue")
abline(v=qq1[3],lty=2,lwd=2,col="blue")

## For estimated density t2
hist(ts2,
     breaks=seq(0,max(ts1)*1.3,0.5),
     prob=TRUE,
     ylim=c(0,0.5),
     xlab="",main="",
     col="gray"
     )
points(xs,mat2[2,],type="l",col="blue")
points(xs,mat2[3,],type="l",lty=2,col="blue")
points(xs,mat2[1,],type="l",lty=2,col="blue")
abline(v=qq2[2],lwd=2,col="blue")
abline(v=qq2[1],lty=2,lwd=2,col="blue")
abline(v=qq2[3],lty=2,lwd=2,col="blue")

## Check the posterior probability of each specie being the weakest
quanMat <- rbind(dens1$quantiles,dens2$quantiles)
(postProb <- get.subset.p(quanMat))
}


}


% \name{sprob}
% \alias{cdf}
% %\alias{quantcdf}
% \alias{nu}
% \alias{dweib}
% \alias{sprob}
% \alias{H}

% %- Also NEED an '\alias' for EACH other topic documented here.
% \title{
% An empirical Bayes approach to the subset selection problem.
% }
% \description{

% The probability measure of the response \emph{T_j} is directly modeled by Dirichlet process with its baseline distribution for the Dirichlet process prior.


% }
% \usage{
% cdf(data,c0=ncol(data),ts=seq(0,15,0.001),alpha=0.05)
% sprob(outputCDF,N.MC=5000)


% }
% %- maybe also 'usage' for other objects documented here.

% \arguments{

% \item{data}{
% A \emph{k} by n matrix, k is the number of species and n is the number of observations (m_k in manuscript).
% This code currently only accept the balanced design, where m_k=m_k' for all k. 
% }

% \item{c0}{
% A scalar, precision parameter \emph{v} of DP.
% The default value is the sample size.
% }

% \item{ts}{ 
% A vector of length M containing  equally spaced values that define the subintervals into which \emph{(0,U)} is partitioned.
% \code{U} is a specified upper limit to the support of the distributions of the \code{alpha}th quantiles for the \emph{K} species assumed to be equal.  
% For example if a resolution of \emph{0.001} were chosen the vector ts=\code{seq(0.001,10,0.001)} would have length \code{M}=10,000.

% %% \code{U}=10.00 has been chosen based on parameter values in JEG.
% }


% \item{alpha}{ The size of the quantile of interest for e.g. 0.05.}
% \item{CDF}{The outputCDF of function cdf.}

% \item{N.MC}{The number of Monte Carlo samples.}


% }

% \value{

% \item{outputCDF}{
% \emph{k+1} by \emph{M} 
% }

% }


% \seealso{
% \code{\link{DPfit}} 
% }

% \references{
% Kondo, Yumi and Zidek, James V. "Bayesian nonparametric subset selection procedures with Weibull components"
% }

% \author{
% Yumi Kondo
% }

% \examples{
% \dontrun{

% x <- simdata(simM=1,iseed=sample(10000,1));
% re <- rep(NA,7);
% c0 <- 10
% alpha <- 0.05
% ts <- seq(0.001,10,0.001)
% M <- length(ts)
% outputcdf <- cdf(ts=ts,x=x,c0=c0,alpha=alpha);
% for (s in 1 : nrow(x))    re[s] <- sprob(M=M,s=s,CDF=outputcdf,x=x,c0=c0,alpha=alpha)

% }
% }

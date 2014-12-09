\name{MCMC.DPw}
\alias{MCMC.DPw}
\alias{MCMC.DPw.dp}
\alias{getden}
\alias{sprob.semi}
\alias{getM}



%- Also NEED an '\alias' for EACH other topic documented here.

\title{
A semiparametric Bayesian approach to the subset selection problem
}

\description{

Our semiparametric Bayesian subset selection procedure has 3 main steps:

STEP 1: obtain the posterior samples from a DP mixture of Weibulls model for strengths of each species;
This step is done by the function \code{MCMC.DPw}.

STEP 2: based on the posterior samples from STEP 1, derive the posterior samples of strengths' density and the posterior samples of strengths' aalpha-th quantiles. Given the output of \code{MCMC.DPw}, this step is done by the function \code{getden}.

(STEPs 1 and 2 are performed for each species separately.)

STEP 3: evaluates the posterior probability that specie k has the smallest alpha-th quantile by comparing the posterior samples of alpha-th quantiles from each species. Given the output of \code{getden}, this step is done by the function \code{sprob.semi}


See the examples to learn the use of these functions.
See reference for more details of the subset selection procedure.


}



\usage{


MCMC.DPw(ts, species, B, aphi = 2, loc.phi = 10, aGam = 2, loc.Gam = 10, 
         minD = 0.01, maxD = 5, burnin = ceiling(B/4), M = NULL, printFreq = B, 
         epsilonM = 0.01, model="indp")

getden(xs, weightS, Klambdas, betas, alpha = 0.05, densi = TRUE, 
       para = TRUE, Kwei = FALSE) 

sprob.semi(s.etas, model = c("indp", "dep")) 

}

%- maybe also 'usage' for other objects documented here.
\arguments{

--- MCMC.DPw ---
  \item{ts}{
 A vector containing strength measures from K species group.
}
\item{species}{
A index vector of length the same as \code{ts}, indicating which species each entry of \code{ts} belongs.
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
\item{model}{internal use only}



--- getden ---
\item{xs}{The grid of points at which the density is evaluated.} 
\item{weightS}{ 
B by M matrix. Each row contains M stick-breaking probabilities sampled via MCMC.
It is usually out.DP$weightS where out.DP is the output of MCMC.DPw.
}

\item{Klambdas}{
B by M matrix. Each row contains M unique values of lambda (scale parameter of Weibull) sampled via MCMC.
It is usually out.DP$lambdas_uni where out.DP is the output of MCMC.DPw.
}

\item{betas}{
B by M matrix. Each row contains M unique values of beta (shape parameter of Weibull) sampled via MCMC.
It is usually out.DP$betas_uni where out.DP is the output of MCMC.DPw.
}

\item{alpha}{quantile of interest}
\item{densi}{
If TRUE then the posterior density estimates of the strength is evaluated at every grid of \code{xs}
} 
\item{para}{Internal use only} 
\item{Kwei}{Internal use only}

--- sprob.semi ---
\item{s.etas}{ 
K by B matrix, each row must contains the posterior samples of alpha-th quantiles from a single species. 
}
}

\details{

We assume that the strengths are sampled independently across species given population parameters, and the population parameters are independent across species priori. The strengths of each species is modelled separately. 

For each species, our DP mixture of Weibull model has following hierarchy:

 \emph{t_j} | \emph{beta_j}\emph{lambda_j} ~ Weibull( \code{shape}=\emph{beta_j},\code{scale}=\emph{lambda_j}), j= 1,...,N

 \emph{beta_j}, \emph{lambda_j} | \emph{G} ~ \emph{G}, j= 1,...,N

 \emph{G} ~ DP(\emph{V,G0})
 where \emph{G0} =  unif(\emph{beta}; \code{min}=0,\code{max}=\emph{phi}) unif(\emph{lambda}; \code{min}=0,\code{max}=\emph{gamma})

 \emph{gamma}     ~ pareto(\code{shape}=\emph{a_phi},\code{location}=\emph{loc.phi}) 

 \emph{phi}     ~ pareto(\code{shape}=\emph{a_gamma},\code{location}=\emph{loc.gamma}) 

 \emph{V}       ~ unif(\emph{beta}; \code{min}=\emph{minD},\code{max}=\emph{maxD}) 



}

\value{

--- MCMC.DPw ---

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




--- getden ---
\item{densities}{B by length(xs) matrix containing \code{B} posterior samples of density of strength at specified grid of points (\code{xs}). }

\item{quantiles}{A vector of length B containing \code{B} posterior samples of alpha-th quantile of the strength. }

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





## generate sample of size 360 from mixture of two weibull distribution
ts1 <- c(rweibull(80,shape=4.726,scale=5.8), rweibull(20,shape=6.4,scale=10.803))
ts2 <- c(rweibull(20,shape=4.726,scale=5.8), rweibull(80,shape=6.4,scale=10.803))
ts.both <- c(ts1,ts2)
species.label <- c(rep(1,length(ts1)),rep(2,length(ts2)))
K <- 2 ## The number of species


## --------------------------- ##
## STEP 1: obtain MCMC samples ##
## --------------------------- ##
## parameters for the MCMC samples
B <- 15000
burnin <- 5000
thin <- 5
useSample <- useSamp(thin=thin,burnin=burnin,B=B)
reMCMC <- MCMC.DPw(ts=ts.both,species=species.label,B=B,burnin=burnin)

## ---------------------------------------------------------------------- ##
## STEP 2: obtain the posterior distribution of strength's 5th percentile ##
## ---------------------------------------------------------------------- ##
alpha <- 0.05  
## Grid of points at which the posterior densities are evaluated
lts <- 10000
xs <- seq(0.01,20,length.out=lts)
dens <- mat <- qq <- list()
for (ik in 1 : K){
  mcmc <- reMCMC[[ik]]
  dens[[ik]] <- getden(xs=xs,
                       weightS=mcmc$weightS[useSample,],
                       Klambdas=mcmc$lambdas_uni[useSample,],
                       betas=mcmc$betas_uni[useSample,],
                       alpha = alpha
                       )
  mat1 <- apply(dens[[ik]]$densities ,2,quantile,prob=c(0.025,0.5,0.975))
  qq1 <- quantile(dens[[ik]]$quantiles,prob=c(0.025,0.5,0.975))
  mat1[2,] <- colMeans(dens[[ik]]$densities)
  qq[[ik]] <- qq1
  mat[[ik]] <- mat1
}
names(dens) <- names(mat) <- names(qq) <- paste("species",1:K,sep="")

##  -- The posterior density estimates and the 95\% credible intervals for 5th quantile --
par(mar=c(3,4,1,1),mfrow=c(1,K))
for (ik in 1 : K){
  hist(ts.both[species.label==ik],
       breaks=seq(0,max(ts1)*1.3,0.5),
       prob=TRUE,
       ylim=c(0,0.5),
       xlab="",main="",
       col="gray",border="white"
       )
  points(xs,mat[[ik]][2,],type="l",col="red")
  points(xs,mat[[ik]][3,],type="l",lty=2,col="blue")
  points(xs,mat[[ik]][1,],type="l",lty=2,col="blue")
  abline(v=qq[[ik]][2],lwd=2,col="red")
  abline(v=qq[[ik]][1],lty=2,lwd=2,col="blue")
  abline(v=qq[[ik]][3],lty=2,lwd=2,col="blue")
}

## ------------------------------------------------------------------ ##
## Step 3: the posterior probability of each specie being the weakest ##
## ------------------------------------------------------------------ ##
quanMat <- rbind(dens$species1$quantiles,dens$species2$quantiles)
(postProb <- sprob.semi(s.etas=quanMat))


## Find the smallest subset that contains the smallest alpha-th percentile
## with posterior probability at least Pstar
Pstar <- 0.99
get.subset(pis=postProb,Pstar=Pstar)
}


}


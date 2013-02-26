###=== function for simulation study ====

F = function(x,w,u,s) sum( w*plnorm(x,meanlog=u,sdlog=s) )
Fwei = function(x,w,u,s) sum( w*pweibull(x,shape=u,scale=s) )

F_inv = function(p,w,u,s,br=c(-1000,1000),wei=FALSE)
{
  if (wei) G = function(x) Fwei(x,w,u,s) - p
  else G = function(x) F(x,w,u,s) - p
  return( uniroot(G,br)$root ) 
}


ASTM <- function(matTs,lowerLimit=0.05,Sig=0.01,confidence=0.75)
  {
    ## The current code works only for the same sample size
    ## matTs: K by mk matrix containing MOR This code only accept the senario when mk=102 for all k or mk=347 for all k
    ## The tolerance limit value for the combined grouping

    m <- ncol(matTs)
    K <- nrow(matTs)
    ## The order statistics of the common sample sizes are computed in advance
    if (length(matTs)==(360*6)){ TLCall <- sort(matTs)[101]
    }else if (length(matTs)==(360*7)) {TLCall <- sort(matTs)[119] ##else  if (length(matTs)==(360*8)) TLCall <- sort(matTs)[136]
    }else if (length(matTs)==(100*6)) {TLCall <- sort(matTs)[26]
    }else if (length(matTs)==(100*7)) {TLCall <- sort(matTs)[31] ##else if (length(matTs)==(100*8)) TLCall <- sort(matTs)[36]
    }else{
      for (k in 1 : m) {
        ## Find 5 % lower tolerance Limit which is the largest order statistics such that P(X_(k) < eta_{0.05}) > confidence
        ## cat("\n k",k," ",1-pbinom(k-1,size=length(matTs),prob=lowerLimit) )
        if( 1-pbinom(k-1,size=length(matTs),prob=lowerLimit)  < confidence)
          {
            TLCall <- sort(matTs)[k-1]
            break;
          }
      }
    }
    ## determine the number of pieces in each species group below/above the group tolerance limit value
    Conting <- table(c(matTs < TLCall),rep(1:K,m))
    ## conduct a chi-square test to determine if the percent of pieces below the group value is statistically significant
    ## for each species in the group.
    re <- chisq.test(Conting,correct = FALSE)
   
    if (re$p.value < Sig)
      {
        ## If the test is significant at the 0.01 level, begin with a subgroup consisting of the
        ## two species with highest percent of pieces below the group value.
        cLessAlpha <- Conting[rownames(Conting)=="TRUE",]
        orderSp <- order(-cLessAlpha)
        subset.out <- orderSp[1]
        for (iSpe in 2 : K)
          {
            pickSp <- orderSp[1:iSpe] ## pick the iSpe species with highest percent of pieces below the group value
            p <- chisq.test(Conting[,pickSp],correct = FALSE)$p.value ## test with the subset of contingency table
            if (p < Sig) break;
            subset.out <- pickSp
          }
        ## Use the chi-square test to determine if the percent of pieces below the group value are comparable.
        
      }else{
        subset.out <- 1:K
      }
    return(subset.out)
  }

D1tr <- function (y, x = 1) 
{
  ## first order numerical derivative  
    n <- length(y)
    if (length(x) == 1) 
        c(y[2] - y[1], 0.5 * (y[-(1:2)] - y[-((n - 1):n)]), y[n] - y[n - 1])/x
    else {
        if (n != length(x)) 
            stop("lengths of 'x' & 'y' must equal")
        if (is.unsorted(x)) 
            stop("'x' must be sorted !")
        c(y[2] - y[1],0.5 * (y[-(1:2)] - y[-((n - 1):n)]),y[n] - y[n - 1])/
          c(x[2] - x[1], 0.5 * (x[-(1:2)] - x[-((n - 1):n)]),x[n] - x[n - 1])
    }
}

weiMLE <- function(Ts)
  {
    dweib <- function(par,x)
      {
        ## density of the weibull distribution 
        return(- sum(dweibull(x,
                              shape=exp(par[1]),
                              scale=exp(par[2]),log=TRUE)) )
      }
    hh <- optim(par=c(log(1),log(1)), ## initial value of dweib,
                fn=dweib,x=Ts)$par
    temp <- exp(hh)
    names(temp) <- c("shape:beta","scale:lambda")
    return(temp)
  }





simdata <- function(dType,iseed=1,n=360,output=0,xs=seq(0.001,20,0.001))
  {
    if (output >= 3) n = 100000
    else set.seed(iseed)
    ## output = 0: generate a simulation dataset
    ## output = 1: return the true values of the 5^th quantiles
    ## output = 2: return the density value of the true distribution. dType==1 only
    ## output = 3: return good estimates of lambda (note Klambda) and beta to develop uninformative prior
    ## output = 4: return good estimates of lambda (note Klambda) and beta to develop uninformative prior
    if (dType==1){
      
      shape <- 3.85
      scale <- 7.32
      
      if (output == 0 || output == 3 || output == 4)
        {
          
          x1 <- rweibull(n,shape=shape,scale=scale)
          
        }else if (output==1){
          
          x1 <- qweibull(0.05,shape=shape,scale=scale)
          
        }else if (output==2)
          
          x1 <- dweibull(xs,shape=shape,scale=scale)
      ## qlnorm(0.05,meanlog=m_,sdlog=sdlog_)
      ## = 3.38

      m_<- 1.85; sdlog_ <- 0.31
      if (output == 0 || output == 3 || output == 4)
        {
          
          x2 <- rlnorm(n,meanlog=m_,sdlog=sdlog_)
          
        }else if (output==1){

          x2 <- qlnorm(0.05,meanlog=m_,sdlog=sdlog_)

        }else if (output==2)
          x2 <- dlnorm(xs,meanlog=m_,sdlog=sdlog_)
      ## qlnorm(0.05,meanlog=m_,sdlog=sdlog_)
      ## = 6.665939
      
      shape <- 5.19
      scale <- 7.27
      if (output == 0 || output == 3 || output == 4)
        {
          x3 <- rweibull(n,shape=shape,scale=scale)

        }else if(output==1){
          
          x3 <- qweibull(0.05,shape=shape,scale=scale)

        }else if (output==2)
          x3 <- dweibull(xs,shape=shape,scale=scale)
      ## qlnorm(0.05,meanlog=m_,sdlog=sdlog_)
      ## = 6.827983

      probs_<- c(0.67,0.33); m_ <- c(1.98,1.74); sdlog_ <- c(0.17,0.23)
      if (output == 0 || output == 3 || output == 4)
        {
          n1 <- sum(sample(1:2,n,prob=probs_,replace=TRUE)==1)
          x4 <- c(rlnorm(n1,meanlog=m_[1],sdlog=sdlog_[1]),rlnorm(n-n1,meanlog=m_[2],sdlog=sdlog_[2]))
      
        }else if(output==1){
          
          x4 <- F_inv(p=0.05,w=probs_,u=m_,s=sdlog_)

        }else if (output==2)
          x4 <- probs_[1]*dlnorm(xs,meanlog=m_[1],sdlog=sdlog_[1])+
          probs_[2]*dlnorm(xs,meanlog=m_[2],sdlog=sdlog_[2])
      ## F_inv(p=0.05,w=probs_,u=m_,s=sdlog_)
      ## = 7.077309

      probs_<- c(0.79,0.21); shapes <- c(5.43,12.01); scales <- c(7.64,6.19)
      if (output == 0 || output == 3 || output == 4)
        {
          n1 <- sum(sample(1:2,n,prob=probs_,replace=TRUE)==1)
          x5 <- c(rweibull(n1,shape=shapes[1],scale=scales[1]),
                  rweibull(n-n1,shape=shapes[2],scale=scales[2]))
          ## F_inv(p=0.05,w=probs_,u=m_,s=sdlog_)
          ## = 5.66328
        }else if(output==1){
          
          x5 <- F_inv(p=0.05,w=probs_,u=shapes,s=scales,wei=TRUE)

        }else if (output==2)
          x5 <- probs_[1]*dweibull(xs,shape=shapes[1],scale=scales[1])+
          probs_[2]*dweibull(xs,shape=shapes[2],scale=scales[2])

      probs_<- c(0.74,0.26); shapes <- c(5.49,15.81); scales <- c(7.60,5.98)
      if (output == 0 || output == 3 || output == 4)
        {
          n1 <- sum(sample(1:2,n,prob=probs_,replace=TRUE)==1)
          x6 <- c(rweibull(n1,shape=shapes[1],scale=scales[1]),
                  rweibull(n-n1,shape=shapes[2],scale=scales[2]))
          ## F_inv(p=0.05,w=probs_,u=m_,s=sdlog_)
          ## = 5.66328
        }else if(output==1){
          
          x6 <- F_inv(p=0.05,w=probs_,u=shapes,s=scales,wei=TRUE)
          
        }else if (output==2)
          
          x6 <- probs_[1]*dweibull(xs,shape=shapes[1],scale=scales[1])+
            probs_[2]*dweibull(xs,shape=shapes[2],scale=scales[2])
      
      probs_<- c(0.98,0.02); m_ <- c(1.90,1.25); sdlog_ <- c(0.19,0.10)
      if (output == 0 || output == 3 || output == 4)
        {
          n1 <- sum(sample(1:2,n,prob=probs_,replace=TRUE)==1)
          x7 <- c(rlnorm(n1,meanlog=m_[1],sdlog=sdlog_[1]),rlnorm(n-n1,meanlog=m_[2],sdlog=sdlog_[2]))
      
        }else if(output==1){
          
          x7 <- F_inv(p=0.05,w=probs_,u=m_,s=sdlog_)

        }else if (output==2)
          x7 <- probs_[1]*dlnorm(xs,meanlog=m_[1],sdlog=sdlog_[1])+
          probs_[2]*dlnorm(xs,meanlog=m_[2],sdlog=sdlog_[2])

      ## m_<- 2.08 ; sdlog_ <- 0.31
      ## if (output == 0 || output == 3 || output == 4)
      ##   {
          
      ##     x8 <- rlnorm(n,meanlog=m_,sdlog=sdlog_)
          
      ##   }else if (output==1){

      ##     x8 <- qlnorm(0.05,meanlog=m_,sdlog=sdlog_)

      ##   }else if (output==2)
      ##     x8 <- dlnorm(xs,meanlog=m_,sdlog=sdlog_)
      ## ##
      qlnorm(0.05,meanlog=m_,sdlog=sdlog_)
      ## = 6.665939
      
      x <- rbind(x1,x2,x3,x4,x5,x6,x7)
      
      if (output==1) x <- c(x)

      
    }else if(dType==2){
      
      ## scales <- seq(10.803-0.6*6,10.803,length.out=7)
      scales <- c(10.803,11.560,11.906,12.858,12.988,15.194,15.651)
      ## The distribution of the weakest specie is set as the same as the estimated dist from JEG
      ## These values are selected to make the percentage gaps between the species' 5th quantiles to be the same as
      ## our confidential dataset.
      if (output==0 || output == 3 || output == 4){
        x <- NULL
        for ( i in 1 : 7) 
          x <- rbind(x, rweibull(n, shape=4.726, scale = scales[i]))

      }else if(output==1){

        x <- rep(NA,7) ;
        for (i in 1 : 7)  x[i] <- qweibull(0.05,shape=4.726,scale=scales[i])                 

      }else if (output==2){

        x <- matrix(NA,nrow=7,ncol=length(xs))
        for ( i in 1 : 7) x[i,] <- dweibull(xs,shape=4.726, scale = scales[i])

      }
      ## for (i in 1 : 7) { h <- qweibull(0.05,shape=4.726,scale=scales[i]);cat(" i:", h)}
      ## i: 5.762343 i: 6.029044 i: 6.295745 i: 6.562446 i: 6.829147 i: 7.095848 i: 7.362549
      ## the difference of about 0.3 return interesting simulation result 
    }else if (dType==3){

      scales <-c(rep(10.803,2),rep(11.906,5))
      
      if (output==0 || output == 3 || output == 4){
        x <- NULL
        for ( i in 1 : 7)
          x <- rbind(x,rweibull(n, shape=4.726, scale = scales[i]))
        
      }else if(output==1){
        
        x <- rep(NA,7)
        for ( i in 1 : 7) x[i] <- qweibull(0.05,shape=4.726, scale = scales[i])
        
      }else if (output==2){
        
        x <- matrix(NA,nrow=7,ncol=length(xs))
        for ( i in 1 : 7) x[i,] <- dweibull(xs,shape=4.726, scale = scales[i])

      }
    }else if (dType==4){

       scales <- c(5.8,10.803)

       ps <- seq(0.9,0,length.out=7)
       x <- NULL
       for (iK in 1 : 7)
         {
           p <- ps[iK]
           if (output == 0 || output == 3 || output == 4)
             {
               n1 <- sum(sample(1:2,n,prob=c(p,1-p),replace=TRUE)==1)
               x <- rbind(
                          x,
                          c(rweibull(n1,shape=4.726,scale=scales[1]),
                            rweibull(n-n1,shape=4.726,scale=scales[2]))
                          )
             }else if(output==1){
               
               x <-rbind(x, F_inv(p=0.05,w=c(p,1-p),u=rep(4.726,2),s=scales,wei=TRUE))
               
             }else if (output==2)
               x <- rbind(x,
                          p*dweibull(xs,shape=4.726,scale=scales[1])+(1-p)*dweibull(xs,shape=4.726,scale=scales[2])
                          )
         }

    }
    ## uninformative prior 
    if (output == 3) x <- weiMLE(c(x))
    else if (output == 4){
      ## informative prior 
      temp <- matrix(NA,nrow=7,ncol=2)
      for (iK in 1 : 7) temp[iK,] <- weiMLE(x[iK,])
      x <- temp
    }
    return(x)
  }

## getInfoPrior <- function(res)
##   {
##     ## gam
##     B <- length(res$gam)
##     useSample <-useSamp(thin=10,burnin=B/3,B=B)
##     mean_gam <- mean(res$gam[useSample])
##     sd_gam <- sd(res$gam[useSample])
    
##     bGam <- mean_gam/sd_gam
##     aGam <- mean_gam*bGam

##     ## dum
##     max_dum <- max(res$dum[useSample])

##     ## phi
##     mean_phi <- mean(res$phi[useSample])
##     median_phi <- median(res$phi[useSample])
##     G <- function(aphi) log(mean_phi)+log(aphi-1)-log(aphi)+log(2)/aphi-log(median_phi)
##     aphi <- uniroot(G,c(1.001,100))$root
##     bphi <- mean_phi*(aphi-1)/aphi
    
##     return ( list(max_dum=max_dum,aGam=aGam,bGam=bGam,aphi=aphi,bphi=bphi) ) 
##   }


ChoosePrior <- function(lambda,beta,aGam=5,aphi=5,
                     delta.min = 1.001,
                     ## delta = 1.001
                     delta.max = 1.001001, 
                     beta.min=0.1)
{
  ## Klambda: the MLE of Klambda or good guess of Klambda 
  ## beta: the MLE of beta or good guess of beta
  ## the hierarchy model:
  
  ## Kwei(scale=Klambda,shape=beta)

  ## Klambda ~ igamma(shape=delta,scale=gam)
  ## delta ~ unif(delta.min,delta.max) where delta.min > 1 and given
  ## gam ~ gamma(shape=aGam,rate=bGam)
  
  ## beta ~ unif(beta.min,phi)
  ## phi ~ pareto(shape=aphi,location=bphi)

  ## The prior specifications procedures:
  ## goal: specify aphi, bphi, max_dum (i.e. delta.max)
  ## idea: given a good guess of shape (beta*) and scale (Klambda*) parameters of a single two-parameter kottas weibull distribution,
  ## select the hyperparameters to satisfy:
  ## E(Klambda) = E(gam/(delta-1)) = E(gam)E(1/(delta-1)) = Klambda* (where delta > 1 ) and
  ## E(beta) = E(beta.min+phi)/2 = beta.min/2+E(phi)/2 = beta*

  ## ::max_dum::
  ## As there are more equations than unknown parameters, we accept an ad-hoc criterio that E(1/(delta-1)) = E1.delm1 = 0.8,
  ## which allows us to select an unique max_dum such that log(max_dum-1)-log(delta.min-1)= 0.1(max_dum-delta.min)
  ## 
  ## ::aGam,bGam::
  ## As gam is from a gamma distribution, the prior has E(gam)/sqrt(a_G)=SD(gam)
  ## So the selection of a_G determines the variance of gam with respect to E(gam).
  ## We allows add-hoc criterion that a_G = 5. Then b_gam = E1.delm1*aGam/Klambda*
  ##
  ## ::aphi,bphi::
  ## we set aphi = 5 then select bphi to satisfy beta.min/2+E(phi)/2 = beta*;
  ## bphi = (aphi-1)*(2*beta-beta.min)/aphi
  lam2Klam <- function(beta,lambda) lambda^beta 
  root.delta.max <- function(delta.max,delta.min,E1.delm1)
    {
      log(abs(delta.max-1))-log(abs(delta.min-1))-E1.delm1*(delta.max-delta.min)
    }
  
  Klambda <- lam2Klam(lambda=lambda, beta= beta)
  cat("Klam",Klambda)

  delta <- (delta.min+delta.max)/2
  bGam <- 1/(delta-1)*aGam/Klambda

  bphi <- (aphi-1)*(2*beta-beta.min)/aphi

  par = list(lambda=lambda,
    beta=beta,
    aGam=aGam,aphi=aphi,
    delta.min=delta.min,
    beta.min=beta.min)
  
  return(
         list(aphi=aphi,
              bphi=bphi,
              max_dum=delta.max,
              aGam=aGam,
              bGam=bGam,
              par=par
              )
         )
}








## ====

plotprior <- function(b.phi,b.gam,
                       lambda.domain=seq(0.01,80,0.5),
                       gam.domain=seq(0,150,0.1),
                       a.gam=1, ## selected to 1 by Kottas
                       beta.domain=seq(0.01,30,0.1),
                       phi.domain=seq(0,30,0.01),
                       a.phi = 2, ## selected by Kottas
                       dum=2)## selected by Kottas
  {
    qpareto <- function (p, location, shape) 
      {
        if (any(p <= 0) || any(p >= 1)) 
          stop("argument 'p' must be between 0 and 1")
        ans = location/(1 - p)^(1/shape)
        ans[location <= 0] = NaN
        ans[shape <= 0] = NaN
        ans
      }
    
    dpareto <- function (x, location, shape, log = FALSE) 
      {
        if (!is.logical(log.arg <- log)) 
          stop("bad input for argument 'log'")
        rm(log)
        L = max(length(x), length(location), length(shape))
        x = rep(x, length.out = L)
        location = rep(location, length.out = L)
        shape = rep(shape, length.out = L)
        logdensity = rep(log(0), length.out = L)
        xok = (x > location)
        logdensity[xok] = log(shape[xok]) + shape[xok] * log(location[xok]) - 
          (shape[xok] + 1) * log(x[xok])
        if (log.arg) 
          logdensity
        else exp(logdensity)
      }

    densigamma <- function (x, alpha, beta) 
      {
        if (alpha > 0 & beta > 0 & all(x > 0)) 
          (beta^alpha)/gamma(alpha) * x^(-alpha - 1) * exp(-beta/x)
        else stop("densigamma: invalid parameters\n")
      }
    ## Assess the values of Hyper parameters based on prior distributions##
    par(mfrow=c(2,2),mar=c(4,2,3,1))
    ## (1) gam
    ## lambda ~ IG(dum,gam) scale
    ## I want the support of the lambda to contain the interval (0, 30)
    ## prior dist'n for gam for the selected hyper-parameter
    ## comment: gam was originally specified to be 25 
    ## so domain of gam should contain 25.
    dens.gam.domain <- dgamma(gam.domain,shape=a.gam,rate=b.gam)
    plot(gam.domain,
         dens.gam.domain,
         type="l",
         xlab=paste("The prior density of gamma \n with b.gamma=",b.gam)
         )
    ## add selected percentiles
    pick.prob <- c(0.25,0.5,0.75)
    abline(v=qgamma(p=pick.prob,shape=a.gam,rate=b.gam),
           col="red",lty=2:4)
    text(qgamma(p=pick.prob,shape=a.gam,rate=b.gam),rep(0.02,3),
         pick.prob) 
    plot(0,0,ylab="",type="n",
         xlim=range(lambda.domain),ylim=c(0,0.05),
         xlab=paste("The prior density of lambda \n with the selected percentile if gamma"))
    for ( i in 1:3)
      {
        points(lambda.domain,
               densigamma(lambda.domain,dum,
                          qgamma(p=pick.prob,shape=a.gam,rate=b.gam)[i]),
               type="l", lty=(2:4)[i]
               )
      }
    legend("topleft",legend=paste(pick.prob),lty=2:4)

    ## (2) phi
    ##  beta ~ Uni(beta|0,phi) shape
    ## prior dist'n for phi for the selected hyper-parameter
    ## comment: phi was originally specified to be 10
    ## so domain of phi should contain 10.
    plot(phi.domain,
         dpareto(phi.domain,shape=a.phi,location=b.phi),
         type="l",xlab=paste("The prior density of phi with b.phi=",b.phi))
    ## add selected percentiles
    pick.prob <- c(0.25,0.5,0.75)
    abline(v=qpareto(p=pick.prob,shape=a.phi,location=b.phi),
           col="red",lty=2:4)
    text(qpareto(p=pick.prob,shape=a.phi,location=b.phi),rep(0.02,3),
         pick.prob)
    plot(0,0,ylab="",type="n",
         xlim=range(beta.domain),ylim=c(0,0.3),
         xlab=paste("The prior density of beta \n with the selected percentile of phi"))
    for ( i in 1:3)
      {
        points(beta.domain,
               dunif(beta.domain,min=0,
                     max=qpareto(p=pick.prob,shape=a.phi,location=b.phi)[i]),
               type="l", lty=(2:4)[i]
               )
      }
    legend("topleft",legend=paste(pick.prob),lty=2:4)

    
  }

plotigamma <- function(shape,scale,poi=FALSE,xs=seq(0,5000,0.001))
{
  ys <-( scale^shape)/gamma(shape) * xs^(-shape - 1) * exp(- scale/xs)
  if (poi)
    points(xs,ys,type="l")
  else
    plot(xs,ys,type="l")
}

plotKwei <- function(beta,Klambda,poi=FALSE)
{
  ts <- seq(0,20,0.001)
  ys <- dweibull(ts,shape=beta,scale=muinv(Klambda,beta))
  if (poi)
    points(ts,ys,type="l")
  else
    plot(ts,ys,type="l")
}

plotgamma<-function(shape,scale,xmin=0,xmax=5,main="",poi=FALSE)
  {
    ## plotting the density of gamma dist'n with given shape and scale
    xs <-seq(xmin,xmax,length.out=1000)
    ys <- dgamma(xs,shape=shape,scale=scale)
    if (!poi)
      plot(xs,ys,main=main,type="l")
    else
      points(xs,ys,type="l")
  }

Nuniq <-function(test,sampleindex){
  tt <- apply(test$S[sampleindex,],1,unique)
  return( unlist(lapply(tt,length)))
}

## E(K|N,D) digamma(D+N)-digamma(D) ~= D*log(1+N/D)
K_N.D <- function(D,N=200) D*log(1+N/D)

## E(K|N) = E(E(K|N,D)) where D ~ gamma(shape=aD,scale1/rD)
## E(L|N) depends on aD and rD and N
K_N <- function(aD,rD,N)
  {
    f <- function(D,N,rD,aD) log(1+N/D)*exp(-rD*D)*D^aD
    k <- integrate(f=f,N=N,rD=rD,aD=aD,lower=0.00001,upper=Inf)
    return(k$value*rD^aD/gamma(aD))
  }


K_NisK <- function(aD,rD,N,K)
  {
    return(K_N(aD,rD,N)-K)
  }

getM <- function(epsilonM,a_D,r_D,prob)
  {
    ## returns the appropriate truncation number M
    Dbig <- qgamma(prob,shape=a_D,scale=1/r_D)
    return(round(1+ log(epsilonM)/(log(Dbig/(1+Dbig))),0))
  }


useSamp <- function(thin,burnin,B)
  {
    ind <- (1:(B/thin))*thin
    return (ind[burnin < ind])
  }

rWEI2 <- function (n, Klambda = 1, beta = 1) 
{

    if (any(Klambda <= 0)) 
        stop(paste("Klambda Klambdast be positive", "\n", ""))
    if (any(beta <= 0)) 
        stop(paste("beta must be positive", "\n", ""))
    if (any(n <= 0)) 
        stop(paste("n must be a positive integer", "\n", ""))
    r <- rweibull(n, scale = muinv(Klambda, beta), shape = beta)
    r
}



qWEI2<-
function (p, Klambda = 1, beta = 1) 
{
    if (any(Klambda <= 0)) 
        stop(paste("Klambda Klambdast be positive", "\n", ""))
    if (any(beta <= 0)) 
        stop(paste("beta must be positive", "\n", ""))
    if (any(p <= 0)) 
        stop(paste("p must be a positive btw 0 and 1", "\n", ""))
    r <-qweibull(p, scale = muinv(Klambda, beta), shape = beta)
    r
}



dWEI2<-function (x, Klambda = 1, beta = 1, log = FALSE) 
  {

    if (any(Klambda <= 0)) 
      stop(paste("Klambda must be positive", "\n", ""))
    if (any(beta <= 0)) 
      stop(paste("beta must be positive", "\n", ""))
    if (any(x < 0)) 
      stop(paste("x must be positive", "\n", ""))
    fy <- dweibull(x, scale = muinv(Klambda, beta), shape = beta, 
                   log = log)
    fy
  }

## density function for Kottas' weibull distribution
muinv <- function(Klambda, beta)
  {
    lambda <- Klambda^(1/beta)
      lambda
  }


DPfit <- function(ts,
		  B=20000,
		  bphi=8,
		  bGam=0.01, ## rate of gamma
		  a_D=NULL,
                  ib_D=0.9,
                  wantK=2,
		  aphi=2,    ## shape of pareto
		  aGam=1,    ## shape of gamma
		  max_dum=1,
                  burnin=ceiling(B/3),
                  M = NULL,
                  printFreq = B, ##B means no printing printing 5 times during the Gibbs run
                  epsilonM=0.01,
                  prob=0.9
                  )
{
  ## This function is designed for ts which ranges between 0 and 10
  if (is.null(a_D)) a_D <- uniroot(f=K_NisK,interval=c(0.0001,10),rD=ib_D,N=length(ts),K=wantK)$root
  if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=ib_D,prob=prob)
  re <- .Call("gibbs", as.numeric(ts), as.integer(B), as.integer(M),
              as.numeric(a_D), 
              as.numeric(ib_D), as.numeric(aphi), as.numeric(bphi), as.numeric(aGam),
              as.numeric(bGam), as.numeric(max_dum),
              as.integer(burnin),as.integer(printFreq),package="DPw"
              )
  N <- length(ts)
  for (i in 1 : 3 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
  for (i in 4 : 7 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
  names(re) <- c("S","betas","Klambdas",
                 "betas_uniq","Klambdas_uniq","vs","weightS",
                 "D","phi","gam","dum",
                 "accept_can")
  names(re$accept_can) <- c("beta_accept","beta_can","dum_accept","dum_can")
                   

  re$para = list(ts=ts,B=B,a_D=a_D,ib_D=ib_D,
    aphi=aphi,bphi=bphi,aGam=aGam,bGam=bGam,max_dum=max_dum,
    burnin = burnin, M = M,pringFreq= printFreq)
  
  return(re)
}




DPfitUnif <- function(ts,
                      B=20000,
                      aphi=2,    ## shape of pareto
                      bphi=8,
                      a_D=0.01, ## D ~ unif (a_D,ib_D)
                      ib_D=10,
                      max_dum=100,
                      burnin=ceiling(B/3),
                      M = NULL,
                      printFreq = B, ##B means no printing printing 5 times during the Gibbs run
                      epsilonM=0.01,
                      prob=0.9
                      )
{
  ## This function is designed for ts which ranges between 0 and 10
  if (is.null(M)) M  <- 1 + log(epsilonM)/log(ib_D/(ib_D+1)) ## vlog((v + m)/v).
  re <- .Call("gibbsUnif",
              as.numeric(ts), as.integer(B), as.integer(M),
              as.numeric(a_D), 
              as.numeric(ib_D), as.numeric(aphi), as.numeric(bphi),
              ##as.numeric(aGam),
              ##as.numeric(bGam),
              as.numeric(max_dum),
              as.integer(burnin),as.integer(printFreq),package="DPw"
              )
  N <- length(ts)
  for (i in 1 : 3 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
  for (i in 4 : 7 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
  names(re) <- c("S","betas","Klambdas",
                 "betas_uniq","Klambdas_uniq","vs","weightS",
                 "D","phi","gam","dum",
                 "AR","can")
  names(re$AR) <-  names(re$can) <- c("beta","dum","gam","D")
                   

  re$para = list(ts=ts,B=B,a_D=a_D,ib_D=ib_D,
    aphi=aphi,bphi=bphi,max_dum=max_dum,
    burnin = burnin, M = M,pringFreq= printFreq)
  
  return(re)
}

DPfitAllUnif <- function(ts,
                         B=20000,
                         aphi=2,    ## shape of pareto
                         loc.phi=10,   ## location 
                         aGam=2,
                         loc.Gam=10, ## location
                         minD=0.01, ## D ~ unif (minD,maxD)
                         maxD=5,
                         burnin=ceiling(B/4),
                         M = NULL,
                         printFreq = B, ##B means no printing printing 5 times during the Gibbs run
                         epsilonM=0.01
                         )
{
  ## This function is designed for ts which ranges between 0 and 10
  if (is.null(M)) M  <- round(1 + log(epsilonM)/log(maxD/(maxD+1))) ## vlog((v + m)/v).
  cat("\n M",M," B",B,"\n")
  re <- .Call("gibbsAllUnif",
              as.numeric(ts), as.integer(B), as.integer(M),
              as.numeric(minD), 
              as.numeric(maxD), as.numeric(aphi), as.numeric(loc.phi),
              as.numeric(aGam),
              as.numeric(loc.Gam),
              ##as.numeric(max_dum),
              as.integer(burnin),as.integer(printFreq),package="DPw"
              )
  N <- length(ts)
  for (i in 1 : 3 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
  for (i in 4 : 7 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
  names(re) <- c("S","betas","lambdas",
                 "betas_uniq","lambdas_uniq","vs","weightS",
                 "D","phi","gam",##"dum",
                 "AR","can")
  names(re$AR) <-  names(re$can) <- c("beta","lambda","D")
                   
  re$S <- re$S + 1
  re$para = list(
    ts=ts,B=B,minD=minD,maxD=maxD,
    aphi=aphi,loc.phi=loc.phi,
    aGam=aGam,loc.Gam=loc.Gam,
    burnin = burnin, M = M,pringFreq= printFreq)
  
  return(re)
}




getden <- function(xs,weightS,Klambdas,betas,alpha=0.05,densi=TRUE,para=TRUE,Kwei=FALSE)
  {
    if (is.vector(weightS)||is.vector(Klambdas)||is.vector(betas))
      {
        weightS <- matrix(weightS,nrow=1)
        Klambdas <- matrix(Klambdas,nrow=1)
        betas <- matrix(betas,nrow=1)
      }
    M <- ncol(weightS)
    B <- nrow(weightS)

    if (Kwei)
      {
        re <- .Call("getDens",
                    as.numeric(xs),
                    as.numeric(c(t(weightS))),
                    as.numeric(c(t(Klambdas))),
                    as.numeric(c(t(betas))),
                    as.integer(M),as.integer(B),
                    as.numeric(alpha),as.integer(densi),
                    package="DPw"
                    )
      }else{
        ## standard weibull
        re <- .Call("getDens_weib",
                    as.numeric(xs),
                    as.numeric(c(t(weightS))),
                    as.numeric(c(t(Klambdas))),
                    as.numeric(c(t(betas))),
                    as.integer(M),as.integer(B),
                    as.numeric(alpha),as.integer(densi),
                    package="DPw"
                    )
      }
    names(re) <- c("densities","quantiles","NbeforeConv")
    if (para) re$para <- list(alpha=alpha,densi=densi,xs=xs)
    re[[1]] <- matrix(re[[1]],nrow=B,ncol=length(xs),byrow=TRUE)
    return(re)
  }


get.subset.p <-function(s.etas) ## K by B-Burnin
  {
    ## objective: compute the Prob(Fs^{-1}(alpha) < Fk^{-1}(alpha)|data) for all k != s
    ##            for s = 1,..., # species
    B <- ncol(s.etas)
    K <- nrow(s.etas)
    counts <- rep(NA,K)
    for (k in 1 : K)
      {
        count <- 0
        ss <- s.etas[k,,drop=FALSE]
        s.etas.wo.s <- s.etas[-k,,drop=FALSE]
        for ( i.B in 1 : B)
          {
            s <- ss[i.B]
            count <- count + prod(rowSums(s < s.etas.wo.s ))
          }
        counts[k] <-count
      }
    return(counts/(B^K))
    ## return a vector of length K, containing estimates of probability that
    ## alpha^th quantile of k^th species is the smallest among all species.
  }

### ================ nonparametric method ================


## ====
## The engine of the analysis is the CDF calculator for supplies CDF values
## for all points in the input vector t and each species i =1,...,k:

cdf <- function(data,c0,ts=seq(0,15,0.001),alpha=0.05)
  {
    M<- length(ts)
    k <- nrow(data)
    n <- ncol(data)
    ## ts:    A vector of length m
    ## data:     A k by n matrix, k is the number of species and
    ##        n is the number of repeated measures m_i in manuscript
    ## c0:    A scalar v0 in the manuscript
    ## alpha: A scalar alpha^th quantile
    z <- matrix(NA,nrow=k,ncol=M)
    logMLEs <- matrix(NA,nrow=k,ncol=2)
    rownames(logMLEs) <- paste("specie",1:k,sep="")
    colnames(logMLEs) <- c("shape","scale")
    for (i in 1 : k)
      {
        ## get MLE of shape and scale parameters of weibull
        ## par = c(log(shape),log(scale))
        par <- logMLEs[i,] <- optim(par=c(log(1),log(1)), ## initial value of dweib,
                                 fn=dweib,x=data[i,])$par
        ## cat("\n shape: ",round(exp(par[1]),3)," scale:",round(exp(par[2]),3))
        for ( j in 1 : M)
          {
            d  <- nu(v=ts[j],x=data[i,],c0=c0,par=par)
            shp2 <- c0+n-d
            if (d<1e-20) d <- 1e-20
            if (shp2<1e-20) shp2 <- 1e-20
            z[i,j] <- 1 - pbeta(alpha,shape1=d,shape2=shp2)
          }
      }
    return(list(outputCDF=rbind(ts,z),MLEs=exp(logMLEs)))
    ## returns k+1 by M matrix: the 1st row is t and each successive row i represents
    ## the value of the posterior cdf { P(eta_{k,alpha} < t |data) }_{k,t}
    ## for the random alpha quantiles associated with species samples.
    ## If t = seq(0.01,10,0.01) s+1^th row of the output corresponds to:
    ## P( Fs^{-1}(alpha) < 0.01 | data) P( Fs^{-1}(alpha) < 0.02 | data) ... P( Fs^{-1}(alpha) < 10 | data)
  }


nu <- function(v,x,c0,par)
  {
    ## x: single specie dataset
    ## v: a point at which  P( Fs^{-1}(alpha) < v | data) is evaluated
    F0 = pweibull(q=v, shape=exp(par[1]), scale=exp(par[2]));
    return( c0*F0 + sum(x < v) # *w/w
           ) 
  }


dweib <- function(par,x)
  {
    ## density of the weibull distribution 
    return(- sum(dweibull(x,
                          shape=exp(par[1]),
                          scale=exp(par[2]),log=TRUE)) )
  }


## Once cdf has been run to yiled the matrix CDF, values of the alpha^th quantile
## for any specified species S can be randomly 


sprob <- function(outputCDF,N.MC)
  {
    ## probability calculator
    ## P(kappa = ispecie | data) \arrpox
    ## sum_{j=1}^M prod_{i \ne s} prob(eta^{(j)}_{s,alpha} < eta_{i,alpha}| data) / M

    ## Given outputCDF matrix, generate N.MC posterior samples quantiles via inverse transforming method
    ts <- outputCDF[1,]
    K <- nrow(outputCDF)-1
    re <- rep(NA,K)
    names(re) <- paste("specie",1:K,sep="")
    for (ispecie in 1 : K) ## compute the 
      {
        ## STEP 1: generate sample of size N.MC from uniform dist'n
        w <- outputCDF[ispecie+1,]
        p <- runif(N.MC)
        
        ## STEP 2:
        ## we use MC approximation 
        ## integral prod_{i \ne k} P(u < eta_{i alpha}|data) dP(eta_{i alpha} < u  | data)
        ## \approx sum_{j}^{N.MC} prod_{i \ne k} P(u^{j} < eta_{i alpha}|data) / N.MC
        ## \approx sum_{j}^{N.MC} prod_{i \ne k} ( 1 - P( eta_{i alpha} =< u_{j} |data)) / N.MC
        z <- 1:N.MC
        temp <- 0;
        for (i in 1 : N.MC)
          {
            ## t[w > p[i]] is the values of F; F(t) > p[i] the smallest of such is the p[i]^th quantile 
            sampledQ <- min(ts[w > p[i]]) ## Quantile F^{-1}(p) = inf {k; F(k) > p }
            temp = temp + H(sampledQ=sampledQ,ispecie=ispecie,outputCDF=outputCDF)
          }
        re[ispecie] <- (temp/N.MC)
        ## returns sum_j prod_{k !=s } P( Fs^{-1}(alpha)^(j) <= Fk^{-1}(alpha) | data) / N.MC
      }
    return(re)
  }

H <- function(sampledQ,ispecie,outputCDF,k1=nrow(outputCDF),m=ncol(outputCDF))
  {
    ## subroutine of sprob
    ## u: A scalar, random sample of the alpha^th quantile of specie s
    maxt <- outputCDF[1,m]
    J <- round(m*sampledQ/maxt,1) ## we can do this because ts are a grid of points lying in the same interval.
    ##   outputCDF[,J]
    ## = P(F_k^{-1}(alpha)< ts[J]|data) for each k=1,..,K
    ## = P(Fk^{-1}(alpha)< u* |data) where k != s where u* is the closest value to u among ts


    CDF1 <- 1 - outputCDF[,J]
    CDF1 <- CDF1[2:k1] ## A vector of length # specie
    return( prod(CDF1[ (1:(k1-1))!= ispecie]) ) ## k1-1 = # species
    ## returns prod_{k !=ispecie } P( u* <= Fk^{-1}(alpha) | data)
  }


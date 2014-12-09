
###=== function simulation study ====

F = function(x,w,u,s) sum( w*plnorm(x,meanlog=u,sdlog=s) )
Fwei = function(x,w,u,s) sum( w*pweibull(x,shape=u,scale=s) )

F_inv = function(p,w,u,s,br=c(-1000,1000),wei=FALSE)
{
  if (wei) G = function(x) Fwei(x,w,u,s) - p
  else G = function(x) F(x,w,u,s) - p
  return( uniroot(G,br)$root ) 
}




get.subset <- function(pis,Pstar)
  {
    sortpis <- sort(pis,decreasing = TRUE)
    rankpis <- rank(-pis)
    ## rank pis in increasing way; the position with the largest probability in pis receive 1 
    selectedSet <- 0
    for (countSize in 1 : length(rankpis))
      {
        ## selectedSet contains the sum of the probabilities of the selected species
        selectedSet <- selectedSet + sortpis[countSize]
        if (selectedSet > Pstar) break; 
      }
    Subset <- which(rankpis <= countSize)
    return(list(Size=countSize,
                Subset=Subset,
                pisSubset = pis[Subset],
                sumPisSubset = sum(pis[Subset])  ))
  }


## when iTrueSpecie is a vector (i.e., there are two species with the lowest quantile), the condition becomes TRUE if both 

TL <- function(ts,lowerLimit=0.05,confidence=0.75)
  {
    Nt <-length(ts)
    
    ## Find 5 % lower tolerance Limit which is the smallest order statistics such that
    ## Pr( F(X[k]) >= lowerLimit ) > confidence
    ## Notice that F(X[k]) is from the beta(t;shape1=k,shape2=n-k+1) (See p265 of Casella book)
    ## Therefore the algorithm is:
    ## increase k from 1,... until when
    ## Pr( F(X[k]) >= lowerLimit ) = 1 - Pr( F(X[k]) < lowerLimit ) = 1 - pbeta(lowerLimit,shape1=k,shape2=n-k+1)
    ## become greater than "confidence"
    order <- which(pbeta(lowerLimit,shape1=1:Nt,shape2=Nt-1:Nt+1) > confidence)
    order <- order[length(order)]
    TLCall <- sort(ts)[order]
        
      ## for (k in 1 : Nt) {
      ## cat("\n k",k," ",1-pbinom(k-1,size=length(matTs),prob=lowerLimit) )
      ## if( 1-pbinom(k-1,size=Nt,prob=lowerLimit)  < confidence)
      ##   {
      ##     TLCall <- sort(ts)[k-1]
      ##     break;
      ##   }
      ##}
    return(TLCall)
  }


ASTM <- function(data, alpha  = c(0.05,0.5),Sig=0.01,detail=FALSE )
  {
    ## LoweLimit and confidence do not have to be specified if alpha = 0.05
    if (length(alpha)> 1) alpha <- alpha[1]
    if (! alpha %in% c(0.05,0.5) ) stop("alpha must be 0.05 or 0.5")
    lowerLimit <- 0.05
    confidence=0.75
    if (alpha == 0.05){
      
      return(ASTMTL(data=data,
                    Sig=Sig,detail=detail))
    }else{
      return(ASTMMed(data=data,Sig=Sig,detail=detail))
    }
  }

ASTMMed <- function(data,Sig=0.01,detail=FALSE)
  {
    resp <- data$resp
    S <- factor(data$S)

    uniS <- levels(S)
    K <- length(uniS)
    r.resp <- rank(resp)
    relme <- lm(r.resp ~ S)
    reANOVA <- anova(relme)

    
    if (reANOVA$"Pr(>F)"[1] > Sig){
      subset.out <- uniS
      mat <- tukey <- NULL
    }else{
      ## m <- tapply(r.resp,S,mean)
      ## rank.m <- rank(m)

      sp.min.median <- names(which.min(tapply(resp,S,median)))
      
      a1 <- aov(r.resp ~ S)
      tukey <- TukeyHSD(x=a1, 'S', conf.level=1-Sig)
      posthoc <- tukey$S
      
     
      
      mat <- matrix(0,K,K)
      ## adjusted p-values of each pair difference
      ## p-value < Sig means significant difference
      mat[lower.tri(mat, diag = FALSE)] <- posthoc[,4] 
      mat <- mat + t(mat) ## symetric
      diag(mat) <- 1
      colnames(mat) <- rownames(mat) <- uniS
      ## Which group is in-significantly different from the smallest sample median group?
      subset.out <- whichInSig <- uniS[which(mat[ uniS ==sp.min.median,] > Sig)]
      ## reorder S in increasing the way that the sample median increases
      ## uniS.order <- uniS[order(m)]
      ## Among the species which are significantly different from the one with the smallest median,
      ## which species has the smallest median?
      ## min.which.Sig.smallest <- min(which( uniS.order %in% whichInSig == FALSE))
      ## subset.out <- uniS.order[1 : (min.which.Sig.smallest - 1)]
    }
    
    if (detail){
      return(list(subset.out=subset.out,TukeyHSD=tukey,
                  matPv=mat))
    }else return(subset.out)
  }



ASTMTL <- function(data,Sig=0.01,detail=FALSE)
  {
    ## data must contains resp and S
    resp <- data$resp
    S <- data$S
    ## The current code works only for the same sample size
    ## matTs: K by mk matrix containing MOR This code only accept the senario when mk=102 for all k or mk=347 for all k
    ## The tolerance limit value for the combined grouping
    lowerLimit <- 0.05
    confidence <- 0.75
    ## K <- length(uniS)
    ## ms <- rep(NA,K)
    ## for (iK in 1 : K )
    ##   {
    ##     ## Sample size of each species
    ##     ms[iK] <- sum(uniS[iK]==S)
    ##   }
    ## The order statistics of the common sample sizes are computed in advance
    if (length(resp)==(360*6)){
      TLCall <- sort(resp)[101]
    }else if (length(resp)==(360*7)) {
      TLCall <- sort(resp)[119] ##else  if (length(resp)==(360*8)) TLCall <- sort(resp)[136]
    }else if (length(resp)==(100*6)) {
      TLCall <- sort(resp)[26]
    }else if (length(resp)==(100*7)) {
      TLCall <- sort(resp)[31] ##else if (length(resp)==(100*8)) TLCall <- sort(resp)[36]
    }else{
      TLCall <- TL(ts=resp,lowerLimit=lowerLimit,confidence=confidence)
    }
    ## determine the number of pieces in each species group below/above the group tolerance limit value
    Counting <- table(c(resp < TLCall),S)
    uniS <- colnames(Counting)
    ms <- colSums(Counting)
    K <- ncol(Counting)
    ## conduct a chi-square test to determine if the percent of pieces below the group value is statistically significant
    ## for each species in the group.
    re <- chisq.test(Counting,correct = FALSE)
    
    if (re$p.value < Sig)
      {
        ## If the test is significant at the 0.01 level, begin with a subgroup consisting of the
        ## two species with highest percent of pieces below the group value.
        cLessAlpha <- Counting[rownames(Counting)=="TRUE",]
        orderSp <- order(-cLessAlpha/ms) ## order the species in terms of the precent of pieces below the TL
        subset.out <- orderSp[1]
        for (iSpe in 2 : K)
          {
            pickSp <- orderSp[1:iSpe] ## pick the iSpe species with highest percent of pieces below the group value
            p <- chisq.test(Counting[,pickSp],correct = FALSE)$p.value ## test with the subset of contingency table
            if (p < Sig) break;
            subset.out <- pickSp
          }
        ## Use the chi-square test to determine if the percent of pieces below the group value are comparable.
        
      }else{
        subset.out <- 1:K
      }
    subset.out <- uniS[subset.out]
    if (detail)return(list(subset.out=subset.out,Freq=Counting,RelFreq=scale(Counting,scale=ms,center=FALSE),strength=TLCall))
    else return(subset.out)
        
  }

## D1tr <- function (y, x = 1) 
## {
##   ## first order numerical derivative  
##     n <- length(y)
##     if (length(x) == 1) 
##         c(y[2] - y[1], 0.5 * (y[-(1:2)] - y[-((n - 1):n)]), y[n] - y[n - 1])/x
##     else {
##         if (n != length(x)) 
##             stop("lengths of 'x' & 'y' must equal")
##         if (is.unsorted(x)) 
##             stop("'x' must be sorted !")
##         c(y[2] - y[1],0.5 * (y[-(1:2)] - y[-((n - 1):n)]),y[n] - y[n - 1])/
##           c(x[2] - x[1], 0.5 * (x[-(1:2)] - x[-((n - 1):n)]),x[n] - x[n - 1])
##     }
## }

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
    if (output == 0 )   rownames(x) <- paste("x",1:7,sep="")
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



MCMC.DPw <- function(
                         ## each element contain the sample from one species
                         ts,
                         species,
                         B,
                         aphi=2
                         ,    ## shape of pareto
                         loc.phi=10
                         ,   ## location 
                         aGam=2
                         ,
                         loc.Gam=10
                         , ## location
                         minD=0.01
                         , ## D ~ unif (minD,maxD)
                         maxD=5
                         ,
                         burnin=ceiling(B/4)
                         ,
                         M = NULL
                         ,
                         printFreq = B
                         , ##B means no printing printing 5 times during the Gibbs run
                         epsilonM=0.01
                         ,
                         model="indp"
                         )
{
  if (length(ts) != length(species)) stop("!")
  if (!is.factor(species)) species <- factor(species)
  if (is.null(M)) M  <- round(1 + log(epsilonM)/log(maxD/(maxD+1))) ## vlog((v + m)/v).
  cat("\n M",M," B",B,"\n")


  if (model=="dep")
    {
      re <- MCMC.DPw.dp(ts=ts,species=species,
                            B=B, aphi=aphi,    ## shape of pareto
                            loc.phi=loc.phi,  
                            aGam=aGam,loc.Gam=loc.Gam, ## location
                            minD=minD, maxD=maxD,
                            burnin=burnin,M = M,
                            printFreq = printFreq)

    }else if (model=="indp"){

      uniS <- levels(factor(species))
      re <- list()
      for (ik in 1 : length( uniS ))
        {
          ## Semiparametric
          ts.each <- ts[species==uniS[ik]]
          re[[ik]] <- MCMC.DPw.indp(ts=ts.each, B=B, aphi=aphi,    ## shape of pareto
                                        loc.phi=loc.phi,  
                                        aGam=aGam,loc.Gam=loc.Gam, ## location
                                        minD=minD, maxD=maxD,
                                        burnin=burnin,M = M,
                                        printFreq = printFreq)
        }
      names(re) <- uniS

    }

  re$para <- list(ts=ts,B=B,minD=minD,maxD=maxD,
    aphi=aphi,loc.phi=loc.phi,
    aGam=aGam,loc.Gam=loc.Gam,
    burnin = burnin, M = M,printFreq= printFreq)
  return(re)
}

  
MCMC.DPw.dp <- function(
                            ## each element contain the sample from one species
                            ts,
                            species,
                            B,
                            aphi=2
                            ,    ## shape of pareto
                            loc.phi=10
                            ,   ## location 
                            aGam=2
                            ,
                            loc.Gam=10
                            , ## location
                            minD=0.01
                            , ## D ~ unif (minD,maxD)
                            maxD=5
                            ,
                            burnin=ceiling(B/4)
                            ,
                            M = NULL
                            ,
                            printFreq = B
                            )
  {
    
    uniS <- levels(species)
    K <- length(uniS)
    mks <- rep(NA,K)
    ts.rightOrder <- NULL
    for (ik in 1 : K)
      {
        pick <- species==uniS[ik]
        ts.rightOrder <- c(ts.rightOrder,ts[pick])
        mks[ik]<- sum(pick)
      }

    ## the sample size of each species
    m_tot <- sum(mks)
    
    re.orig <- .Call("gibbsAllUnifAllSample",
                     as.numeric(ts.rightOrder), 
                     as.integer(B), 
                     as.integer(M),
                     as.numeric(minD), 
                     as.numeric(maxD), as.numeric(aphi), as.numeric(loc.phi),
                     as.numeric(aGam),
                     as.numeric(loc.Gam),
                     as.integer(burnin),
                     as.integer(printFreq),
                     as.integer(K),
                     as.integer(mks),
                     as.integer(max(mks)))

    cumsum_mks0 <- c(0, cumsum(mks)) ## length K + 1


    ## Structure the data so that output is grouped by species
    re.spec <- rep(list(NULL),K) 
    names(re.spec) <- uniS
    
    ## re.orig.mat <- list()
    ## for (ire in 1 : 7)
    ##   {
    ##     if (ire <= 3)
    ##       re.orig.mat[[ire]] <- matrix(re.orig[[ire]],nrow=B,ncol=m_tot,byrow=TRUE)
    ##     if (ire == 1) re.orig.mat[[ire]] <- re.orig.mat[[ire]] + 1## cluster labels
    ##     if (ire > 3)
    ##       re.orig.mat[[ire]] <- matrix(re.orig[[ire]],B,M*K,byrow=TRUE)
    ##   }    
    for (ik in 1 : K)
      {
        eachSpec <- list()
        for (ire in 1 : 3 )
          {
            ## cat("\nik",ik,"ire",ire)
            ## print(Sys.time())
            pick <- rep(1:mks[ik],B) + rep(m_tot*(0:(B-1)),each=mks[ik])+cumsum_mks0[ik]
            eachSpec[[ire]] <- matrix(re.orig[[ire]][pick],byrow=TRUE,ncol=mks[ik])
            #eachSpec[[ire]] <- re.orig.mat[[ire]][, (cumsum_mks0[ik]+1) : cumsum_mks0[ik+1]]
          }
        for (ire in 4 : 7 )
          {
            ## cat("\nik",ik,"ire",ire)
            ## print(Sys.time())
            pick <- rep(1:M,B) + rep(K*M*(0:(B-1)),each=M) + M*(ik-1)
            eachSpec[[ire]] <- matrix(re.orig[[ire]][pick],byrow=TRUE,ncol=M)
            #eachSpec[[ire]] <- re.orig.mat[[ire]][, ( M*(ik-1) + 1) : ( M*ik ) ]
          }
        names(eachSpec) <- c("S","betas","lambdas","betas_uniq","lambdas_uniq","vs","weightS")
        eachSpec$S <- eachSpec$S + 1 ## cluster labels
        re.spec[[ik]] <- eachSpec
      }

    re <- list()
    for (ire in 8 : length(re.orig)) re[[ire-7]] <- re.orig[[ire]]
    names(re) <- c("D","phi","gam","AR","can")
    re$sp <- re.spec
    names(re$AR) <-  names(re$can) <- c("beta","lambda","D")
   
 
    return(re)
  }

MCMC.DPw.indp <- function(ts,
                              B=20000,
                              aphi=2,    ## shape of pareto
                              loc.phi=10,   ## location 
                              aGam=2,
                              loc.Gam=10, ## location
                              minD=0.01, ## D ~ unif (minD,maxD)
                              maxD=5,
                              burnin=ceiling(B/4),
                              M = NULL,
                              printFreq = B ##B means no printing printing 5 times during the Gibbs run
                              )
{
  ## This function is designed for ts which ranges between 0 and 10

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

## This used to be called "get.subset.p"
sprob.semi <-function(s.etas,model=c("indp","dep")) ## K by B-Burnin
  {
    ## objective: compute the Prob(Fs^{-1}(alpha) < Fk^{-1}(alpha)|data) for all k != s
    ##            for s = 1,..., # species
    if (length(model)>1) model <- model[1]
    B <- ncol(s.etas)
    K <- nrow(s.etas)
    counts <- rep(NA,K)
    for (k in 1 : K)
      {
        
        ss <- s.etas[k,]
        s.etas.wo.s <- s.etas[-k,,drop=FALSE]
        if (model=="indp")
          {
            count <- 0
            for ( i.B in 1 : B)
              {
                s <- ss[i.B]
                count <- count + prod(rowSums(s < s.etas.wo.s ))
              }
          }else if(model=="Dep"){
            ## for each
            aug.ss <- matrix(ss,nrow=K-1,ncol=B,byrow=TRUE) 
            count <- sum(colSums(s.etas.wo.s > aug.ss) == K-1)
            ##count <- sum(colSums(ss < s.etas.wo.s)==K-1)
              ## for (i.B in 1 : B)
              ##   {
              ##     etaK <- ss[i.B]
              ##     etaRest <- s.etas.wo.s[,i.B]
              ##     ## apply takes too long time 
              ##     ## count <- count + sum( apply(s < s.etas.wo.s,2,all ) );
              ##     ## count <- count + sum( s < colMins(s.etas.wo.s))
              ##     ## count <- count + sum(colSums( s < s.etas.wo.s )  == (K-1) )
              ##     ##if(i.B%%100 == 0) cat(i.B)
              ##   }
            }
        counts[k] <-count
      }
    if(model=="indp"){
      re <- counts/(B^K)
    }else if(model == "dep")
      re <- counts/B
    return(re)
    ## return a vector of length K, containing estimates of probability that
    ## alpha^th quantile of k^th species is the smallest among all species.
  }


## ====
## The engine of the analysis is the CDF calculator for supplies CDF values
## for all points in the input vector t and each species i =1,...,k:
beliefToW <- function(c0.belief.G0,n=1) c0.belief.G0*n/(1 - c0.belief.G0)
CDF.eBayes <- function(data,c0s,c0.belief.G0=NULL,ts=NULL,length.out=1E+6,qs=c(1E-10,1-1E-10),alpha=0.05,C=TRUE)
  {
    
    if (is.null(ts) & (length(length.out) !=1| length(qs) !=2))
      stop("If ts (the grid of points to evaluate the posterior CDF), then qs and length.out must be specified!")
    if (is.null(ts))
      {
        if (length(qs) != 2)
          stop("qs must be a vector of length 2 containing minimum and maximum quantiles to evaluate")
        
      }

    
    if (missing(c0s) & missing(c0.belief.G0))
      stop("One of c0 or c0.belief.G0 must be specified")
    if (is.matrix(data))
      {
        temp <- data
        data <- list()
        for (irow in 1 : nrow(temp))
          data[[irow]] <- temp[irow,]
        names(data) <- rownames(temp)
      }
    if (!is.list(data)) stop("!!!")
    K <- length(data)

    ms <- unlist(lapply(data,length))
    if (!is.null(c0.belief.G0))
      {
        if (0 < c0.belief.G0 & c0.belief.G0 <  1){
          c0s <- beliefToW(c0.belief.G0,n=ms)
        }else stop("c0.belief.G0 must be within [0,1)")
      }
    if (length(c0s)==1) c0s <- rep(c0s,K)
    ## ts:    A vector of length m
    ## data:     A K by n matrix, K is the number of species and
    ##        n is the number of repeated measures m_i in manuscript
    ## c0s:    A scalar v0 in the manuscript
    ## alpha: A scalar alpha^th quantile

    logMLEs <- matrix(NA,nrow=K,ncol=2)
    rownames(logMLEs) <- names(data)
    colnames(logMLEs) <- c("shape","scale")
    ## Determine the MLEs of the parameters of beta distribution
    tsbound <- matrix(NA,K,2)
    for (ik in 1 : K)
      {
        dat.single <- data[[ik]]
        n <- length(dat.single)
        ## get MLE of shape and scale parameters of weibull
        ## par = c(log(shape),log(scale))
        logMLEs[ik,] <- optim(par=c(log(1),log(1)), ## initial value of dweib,
                             fn=dweib,x=dat.single)$par
        if (is.null(ts)){
          for (iq in 1 : length(qs))
            tsbound[ik,iq] <- uniroot(petapost.root,interval=c(0,max(dat.single)),
                                     par=logMLEs[ik,],alpha=alpha,data=dat.single,
                                     c0=c0s[ik],q=qs[iq])$root
        }
      }
    if (is.null(ts)) ts <- seq( min(tsbound[,1]),max(tsbound[,2]),length.out=length.out)
    M <- length(ts)
    z <- matrix(NA,nrow=K,ncol=M)
    for (i in 1 : K)
      {
        dat.single <- sort(data[[i]])
        m <- ms[i]
        ## get MLE of shape and scale parameters of weibull
        ## par = c(log(shape),log(scale))
        par <- logMLEs[i,]

        c0.single <- c0s[i]
        ## cat("\n",i,"th species","prec.para v=",c0.single )
        ## cat("\n shape: ",round(exp(par[1]),3)," scale:",round(exp(par[2]),3))
        if (C){
          z[i,] <- .Call("cdf_sub",
                         as.numeric(ts),
                         as.numeric(exp(par[1])),
                         as.numeric(exp(par[2])),
                         as.numeric(alpha), as.numeric(dat.single),as.numeric(c0.single),package="DPw")[[1]]
        }else{
          for ( j in 1 : M)
            {
              z[i,j] <- petapost(t=ts[j],par=par,alpha=alpha,
                                 data=dat.single,c0=c0.single,n=m)
            }
        }
      }
    outputCDF <- rbind(ts,z)
    rownames(outputCDF) <-c("ts.grid",names(data))
    return(list(outputCDF=outputCDF, MLEs=exp(logMLEs)))
    ## returns k+1 by M matrix: the 1st row is t and each successive row i represents
    ## the value of the posterior cdf { P(eta_{k,alpha} < t |data) }_{k,t}
    ## for the random alpha quantiles associated with species samples.
    ## If t = seq(0.01,10,0.01) s+1^th row of the output corresponds to:
    ## P( Fs^{-1}(alpha) < 0.01 | data) P( Fs^{-1}(alpha) < 0.02 | data) ... P( Fs^{-1}(alpha) < 10 | data)
  }

petapost.root <- function(t,par,alpha,data,c0,q)
  {
    petapost(t=t,par=par,alpha=alpha,data=data,
             c0=c0,n=length(data)) - q
  }

petapost <- function(t,par=NULL,alpha,data,c0,n=length(data))
  {
    if (is.null(par))
      {
        par <- optim(par=c(log(1),log(1)), ## initial value of dweib,
                     fn=dweib,x=data)$par
      }
    d <- nu(v=t,x=data,c0=c0,par=par)
    shp2 <- c0 + n - d
    if (d<1e-20) d <- 1e-20
    if (shp2<1e-20) shp2 <- 1e-20
    return( 1 - pbeta(alpha,shape1=d,shape2=shp2) )
  }





nu <- function(v,x,c0,par)
  {
    ## x: single specie dataset
    ## v: a point at which  P( Fs^{-1}(alpha) < v | data) is evaluated
    F0 = pweibull(q=v, shape=exp(par[1]), scale=exp(par[2]));
    return( c0*F0 + sum(x <= v) # *w/w
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

## sp <- function(nonCDF,T,S,N.MC=1E+5,c0.belief.G0=0.1,alpha=0.05)
##   {

##     outputCDF <- nonCDF$outputCDF
##     ts.grid <- outputCDF[1,]
##     parMLE <- nonCDF$MLEs

##     uniS <- rownames(outputCDF)[-1]
##     K <- length(uniS)
    
##     msALL <- tapply(S,S,length)
##     max.m <- max(msALL)
##     TsOthersALL <- tapply(T,S,
##                           function(x){
##                             sort.x <- sort(x)
##                             if( length(x)==max.m) return(sort.x)
##                             else if ( length(x) < max.m) return(c(sort.x,rep(Inf,max.m-length(x))))
##                             else stop()
##                           })
##     c0sALL <- beliefToW(c0.belief.G0,n=msALL)
##     re <- rep(NA,K)
##     for (ik in 1 : K)
##       {
##         cdf <- outputCDF[ik+1,]
##         p <- runif(N.MC)
##         par <- parMLE[rownames(parMLE)!=uniS[ik],]
        
##         ms <- msALL[names(msALL)!= uniS[ik]]
##         c0s <- c0sALL[names(c0sALL)!= uniS[ik]]
##         cat("\n\n",ik,"c0s",c0s)
##         TsOthers <- TsOthersALL
##         TsOthers[[which(names(TsOthersALL)==uniS[ik])]] <- NULL
##         TsOthers <- unlist(TsOthers)
##         re[ik]  <- .Call("sprobC2",
##                          as.numeric(sort(p)),
##                          as.numeric(ts.grid),
##                          as.numeric(cdf),
##                          as.numeric(TsOthers),
##                          as.integer(ms),
##                          as.numeric(t(par)),
##                          as.numeric(c0s),
##                          as.numeric(alpha),
##                          as.integer(max.m)
##                          ,package="DPw"
##                          )[[1]]
##       }
##   }


sprob <- function(outputCDF,N.MC,ExcludeS=NULL)
  {
    ## probability calculator
    ## P(kappa = ik | data) \arrpox
    ## sum_{j=1}^M prod_{i \ne s} prob(eta^{(j)}_{s,alpha} < eta_{i,alpha}| data) / M

    ## Given outputCDF matrix, generate N.MC posterior samples quantiles via inverse transforming method
    ts <- outputCDF[1,]
    K <- nrow(outputCDF)-1
    re1 <- re2 <- re3 <- rep(NA,K)
    names(re1) <- names(re2) <- rownames(outputCDF[-1,]) 
    Exclude_orig_k <- which(ExcludeS == names(re1)) 
    if (length(Exclude_orig_k) == 0) Exclude_orig_k <- -100
    for (ik in 1 : K) ## compute the 
      {
        ## STEP 1: generate sample of size N.MC from uniform dist'n
        w <- outputCDF[ik+1,] 
        ## STEP 2:
        ## we use MC approximation 
        ## integral prod_{i \ne k} P(u < eta_{i alpha}|data) dP(eta_{i alpha} < u  | data)
        ## \approx sum_{j}^{N.MC} prod_{i \ne k} P(u^{j} < eta_{i alpha}|data) / N.MC
        ## \approx sum_{j}^{N.MC} prod_{i \ne k} ( 1 - P( eta_{i alpha} =< u_{j} |data)) / N.MC
        ##if (C){
        outputCDF1 <- c( t(outputCDF[-c(1,ik+1),]) )
        if (Exclude_orig_k > 0){

          if (ik < Exclude_orig_k){
            Exclude_k <- Exclude_orig_k - 2
            ## -1 because R -> C index difference
            ## -1 because ik < Exclude_orig_k
          }else if (ik == Exclude_orig_k){
            Exclude_k <- -100
          }else if (ik > Exclude_orig_k){
            Exclude_k <- Exclude_orig_k - 1
            ## -1 because ik < Exclude_orig_k
          }
          
        }else Exclude_k <- -100
        ## re[ik]  <- .Call("sprobC",
        ##                       as.numeric(sort(runif.sobol(n=N.MC, dimension=1))),
        ##                       as.numeric(ts),
        ##                       as.numeric(w),
        ##                       as.numeric(outputCDF1),package="DPw"
        ##                       )[[1]]
        re  <- .Call("sprobC2",
                     as.numeric(sort(runif.sobol(n=N.MC, dimension=1))),
                     as.numeric(ts),
                     as.numeric(w),
                     as.numeric(outputCDF1),
                     as.integer(Exclude_k),
                     package="DPw"
                     )
        re1[ik]  <- re[[1]]
        re2[ik]  <- re[[2]]
        re3[ik] <- re[[3]]
        ## }else{
        ##   temp <- 0;
        ##   cat("specie",ik)
        ##   for (i in 1 : N.MC)
        ##     {
        ##                                 #cat("i",i)
        ##       ## t[w > p[i]] is the values of F; F(t) > p[i] the smallest of such is the p[i]^th quantile 
        ##       sampledQ <- min(ts[w > p[i]]) ## Quantile F^{-1}(p) = inf {k; F(k) > p }
        ##       temp = temp + H(sampledQ=sampledQ,ik=ik,outputCDF=outputCDF)
        ##     }
        ##   re1[ik] <- (temp/N.MC)
        ## }
        ## returns sum_j prod_{k !=s } P( Fs^{-1}(alpha)^(j) <= Fk^{-1}(alpha) | data) / N.MC
      }
    
    if (is.null(ExcludeS)) return(list(All=re1,ExpEta=re3)) else  return(list(All=re1,Exclude=re2,ExpEta=re3))
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


pr.discrete <- function(Ts,S,par.mat,c0.belief.G0,uniS,alpha)
  {
    re <- list()
    for (ik in 1 : length(uniS))
      {
        ts <- Ts[ S == uniS[ik]]
        sort.ts <- c(0,sort(ts),Inf)
        Pr.Eta <- rep(NA,(length(sort.ts)-1 ) )
        c0 <-beliefToW(c0.belief.G0,n=length(ts)) 
        for (i in 1 : (length(sort.ts)-1 ) )
          {
            x1 <- sort.ts[i]
            x2 <- sort.ts[i+1]
            Pr.Eta[i] <- Pr.EtaEqualObs(x1=x1,x2=x2,
                                        par=par.mat[rownames(par.mat)== uniS[ik],],
                                        c0=c0,x=ts,alpha=alpha)
          }
        obsTs <- sort.ts[-length(sort.ts)]
        Pr.order <- rbind(obsT=obsTs[order(-Pr.Eta)],Pr=-sort(-Pr.Eta))
        re[[ik]] <- list(Pr=rbind(obsT=obsTs,Pr=Pr.Eta),
                         Pr.order=Pr.order,
                         allPr=sum(Pr.Eta))
      }
    names(re) <- uniS
    return(re)
  }
Pr.EtaEqualObs <- function(x1,x2,c0,x,par,m=length(x),alpha)
  {
    if (x2==Inf) return(0)
    ## x1 < x2 adjacent order statistics
    G0 <- pweibull(q=x2, shape=par[1], scale=par[2]);
    F2 <- sum(x <= x2)
    F1 <- sum(x <= x1)
    
    a2 <- c0*G0 + F2
    b2 <- c0 + m - a2

    a1 <- c0*G0 + F1
    b1 <- c0 + m - a1

    if (a1==0 | a2 == 0 | b1 == 0 | b2 == 0) return(0)
    ## Pr( eta <= X(r) | data ) - lim_{ t -> X(r)- } Pr( eta <= t | data )
    return(pbeta(alpha,a1,b1) - pbeta(alpha,a2,b2))
  }

### ===================





## g <- function(x,a,b,log=FALSE){
##   if (a < 0 | b < 0) stop()
##   temp <- (a-1)*log(x) + (b-1)*log(1-x)
##   if (log) return(temp) else return( exp(temp) )
## }
## betalog<- function(x,a,b)
##   {
##     temp <- g(x=x,a=a,b=b)*( log(x)-log(1-x))

##     ## log x/(1-x)
##     return( temp )
##   }

## der.beta <- function(t,a,b,c0,par)
## {
##   ## Derivative of beta function without c0*dweibull(t,shape=par[1],scale=par[2])
##   term1 <- beta(a,b)
##   ## poly-gamma function
##   term2 <- digamma(a) - digamma(b)

##   return(term1 * term2 )
## }


## der.ibeta <- function(t,a,b,c0,alpha,par)
## {
##   if ( a < 1)
##     {
##       term <- integrate(betalog,a=a,b=b,lower=alpha,upper=1,
##                         subdivisions = 10000L,stop.on.error = FALSE,
##                         rel.tol = .Machine$double.eps^2)$value
##       out <- der.beta(t=t,a=a,b=b,c0=c0,par=par) - term
##     }else{
##       out <- integrate(betalog,a=a,b=b,lower=0,upper=alpha,
##                        subdivisions = 10000L,stop.on.error = FALSE,
##                        rel.tol = .Machine$double.eps^2)$value
##     }
##   return(out)
## }


## ibeta <- function(x,a,b,lbeta.ab=NULL,log=FALSE)
## {
##   ## compute int_0^x t^a (1-t)^b dt
##   ## lbeta.ab must contain int_0^1 x^a (1-x)^b dx
##   if (is.null(lbeta.ab)) lbeta.ab <- lbeta(a,b)
##   ## incomplete beta
##   temp <- pbeta(x,shape1=a,shape2=b,log.p=TRUE) + lbeta.ab
##   if (log) return(temp) else return(exp(temp))
## }


## get.nu <- function(v,x,c0,par)
##   {
##     ## x: single specie dataset
##     ## v: a point at which  P( Fs^{-1}(alpha) < v | data) is evaluated
##     F0 <- pweibull(q=v, shape=par[1], scale=par[2]);
##     m.hatF <-  sum(x <= v)
##     #cat("\n sum(x<=v)",m.hatF)
##     return( c0*F0 + m.hatF # *w/w
##            ) 
##   }


## get.nu.fixedF <- function(v,c0,par,m.hatF)
##   {
##     ## x: single specie dataset
##     ## v: a point at which  P( Fs^{-1}(alpha) < v | data) is evaluated
##     F0 <- pweibull(q=v, shape=par[1], scale=par[2]);
##     #m.hatF <-  sum(x <= v)
##     #cat("\n sum(x<=v)",m.hatF)
##     return( c0*F0 + m.hatF # *w/w
##            ) 
##   }


## der.Eta <- function(t,c0,par,alpha,m.hatF,m)
## {
##   a <- get.nu.fixedF(v=t,c0=c0,par=par,m.hatF=m.hatF)
##   b <- c0 + m - a ## 0 <= b <= c0 + m
##   if (a <= 0 | b <= 0) return(0)


##   lbeta.ab <- lbeta(a,b)

##   #cat("\n t",t,"a",a,"b",b,"c0",c0,"m.hatF",m.hatF)

##   if (lbeta.ab == Inf){
##     return(0)
##   }


##   der.den <- der.beta(t=t,a=a,b=b,c0=c0,par=par)
##   der.num <- der.ibeta(t=t,a=a,b=b,c0=c0,alpha=alpha,par=par)

##   den <-  exp(lbeta.ab)
##   ## num <- ibeta(x=alpha,a=a,b=b,lbeta.ab=lbeta.ab,log=FALSE)
##   ##
##   ## term1 <- der.den*num/(den^2)
##   lnum <- ibeta(x=alpha,a=a,b=b,lbeta.ab=lbeta.ab,log=TRUE)
##   term1 <- der.den*exp(lnum-2*lbeta.ab)
##   term2 <- der.num/den
##   out <- term1 - term2
##   if (out < 0) out <- 0
##   const <- c0*dweibull(t,shape=par[1],scale=par[2])
##   return(const*out)
## }



## get.H <- function(t,all.ts,c0,par.mat,alpha,pickS,all.S,
##                   uniS=unique(all.S),ms=tapply(all.ts,all.S,length))
##   {
##     ## Hs(t) = prod_{k != pickS } Pr(t < eta | data)
##     ##       = prod_{k != pickS } pbeta(alpha, nu(t), v + m - nu(t))
##     if (t == 0)
##       {
##         ## if t = 0 then nut = 0 which is impossible
##         return(0)
##       }
##     lCDF <- 0
##     for (iS in 1 : length(uniS))
##       {
##         if (pickS == iS) next;

##         pick <- all.S == uniS[iS]
##         par <- par.mat[iS,]
##         ts <- all.ts[pick]
##         m <- ms[iS]
##         nut <- get.nu(v=t,x=ts,c0=c0,par=par)

##         if ( nut == 0 | c0 + m == nut)
##           {
##             lCDF <- -Inf
##             break
##           }else{
##             ltemp <- pbeta(alpha, shape1 = nut, shape2 = c0 + m - nut,log=TRUE)
##           }

##        ## cat("\n nut",nut,"c0+m-nut",c0+m-nut,"temp",temp)
##         lCDF <- lCDF + ltemp
##       }
##     return(exp(lCDF))
##   }



## intFun <- function(obs.etas,m.hatF,alpha,c0,all.ts,pickS,all.S,uniS,ms=tapply(all.S,all.S,length),par.mat)
##   {
##     derEtas <- rep(NA,length(obs.etas))
##     m <- ms[pickS]
##     par <- par.mat[pickS,]
##     for (it in 1 : length(obs.etas))
##       {
##         t <- obs.etas[it]
##         derE <-  der.Eta(t=t,m.hatF=m.hatF,c0=c0,
##                          par=par,
##                          alpha=alpha,m=m)
##         ## Cannot be negative!
##         #cat(" derEta",derE)
##         derEtas[it] <- derE*
##           get.H(t=t,all.ts=all.ts,c0=c0,
##                 par.mat=par.mat,alpha=alpha,pickS=pickS,
##                 uniS=uniS,ms=ms,all.S=all.S)
##       }
##     return(derEtas)
##   }

## all.ts <- dat1$MOR
## all.S <- dat1$Species.Group
## uniS <- unique(all.S)
## K <- length(uniS)
## par.mat <- matrix(NA,K,2)
## ms <- pis <- rep(NA,K)
## for (ik in 1 : K)
##   {
##     ts <- all.ts[all.S==uniS[ik]]
##     par.mat[ik,] <- exp(optim(par=c(log(1),log(1)), ## initial value of dweib,
##                               fn=dweib,x=ts)$par)
##     ms[ik] <- length(ts)
##   }
## Sys.time()
## c0 <- 300
## alpha <- 0.5
## for (ik in 1 : K)
##   {
##     cat("\n\n k",ik)
##     ts <- all.ts[all.S==uniS[ik]]
##     sort.ts <- c(0,sort(ts),Inf)
##     int1 <- int2 <- Pr.Eta <- rep(NA,length(sort.ts) - 1)
##     m <- ms[ik]
##     for (it in 1 : (length(sort.ts) - 1) )
##       {
##         x2 <- sort.ts[it+1]
##         x1 <- sort.ts[it]
##         int1[it] <- integrate(intFun,alpha=alpha,
##                               m.hatF = sum(ts <= (x1+x2)/2),
##                               c0 = c0,
##                               all.ts = all.ts,
##                               pickS = ik,
##                               par.mat = par.mat,
##                               all.S = all.S,
##                               uniS = uniS,
##                               ms = ms,
##                               lower = x1,
##                               upper = x2,
##                               subdivisions = 10000L,
##                               stop.on.error = FALSE,
##                               rel.tol = .Machine$double.eps^0.7)$value

##         Pr.Eta[it] <-  Pr.EtaEqualObs(x1=x1,x2=x2,
##                                       par=par.mat[ik,],
##                                       c0=c0,
##                                       x=ts,m=m)

##         int2[it] <- get.H(t=x2,all.ts=all.ts,
##                           c0=c0,
##                           par.mat=par.mat,
##                           alpha=alpha,
##                           pickS=ik,
##                           uniS=uniS,
##                           ms=ms,all.S=all.S)

##       }
##     pis[ik]<- sum(int1) + sum(int2*Pr.Eta)
##     cat("\n pi",pis[ik],"sumPr.Eta",sum(Pr.Eta))
##   }
## Sys.time()


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Linpack.h>


// ==== subroutines =====

double dWEI2(const double x,  const double beta, const double Klambda, const int logTF )
{
  // density function for Kottas' weibull distribution
  // res:   A vector of size 1, this function changes its value and contains the result
  // x:     A vector of size 1, containing the value to evaluate the Kweibull at.
  if (Klambda <= 0)
    {
      error("\v ERROR! Klambda must be positive");
    }
  if (beta <= 0)
    {
      error("\v ERROR! beta must be positive");
    }
  if (x < 0 )
    {
      error("\v ERROR! x must be positive");
    }
  return( dweibull(x,
		   beta, // shape
		   R_pow(Klambda,(1/beta)), //scale
		   logTF// log
		   )) ;
}

double pWEI2(const double x,  const double beta, const double Klambda, const int logTF )
{
  // density function for Kottas' weibull distribution
  // res:   A vector of size 1, this function changes its value and contains the result
  // x:     A vector of size 1, containing the value to evaluate the Kweibull at.
  if (Klambda <= 0)
    {
      error("\v ERROR! Klambda must be positive");
    }
  if (beta <= 0)
    {
      error("\v ERROR! beta must be positive");
    }
  if (x < 0 )
    {
      error("\v ERROR! x must be positive");
    }
  return( pweibull(x,
		   beta, // shape
		   R_pow(Klambda,(1/beta)), //scale
		   1, // Lower tail
		   logTF// log
		   )) ;
}

double dinvgamma(const double x, const double shape,const double scale, int logTF)
{
  if (shape < 0 || scale < 0 || x < 0)
    {
      //Rprintf("\v ERROR in dinvgamma!shape %f scale %f x %f",shape,scale,x);
      error("\v ERROR in dinvgamma!");
    }

  if (logTF)
    return shape*log(scale)-lgammafn(shape)-(shape+1.0)-scale/x;
  else
    return (R_pow(scale,shape))/gammafn(shape) * R_pow(x,(-shape - 1.0)) * exp(-scale/x);
}

int sample(double prob[], const int size)
{
  // randomly select the number between 0 to size-1 
  double p_tot=0.0;
  int ih,rN[size];
  for (ih = 0; ih < size; ih ++) p_tot += prob[ih];
  if (p_tot <0) error("At least one prob must be positive!");
  for (ih = 0; ih < size; ih ++) prob[ih] /= p_tot;

  rmultinom(1, prob, size, rN);

  for ( ih=0; ih < size; ih++)
    {
      if(rN[ih]==1) return(ih);
    }
}




// rpareto <-
// function (n, location, shape) 
// {
//     ans = location/runif(n)^(1/shape)
//     ans[location <= 0] = NaN
//     ans[shape <= 0] = NaN
//     ans
// }
// <environment: namespace:VGAM>
double rpareto(const double location, // location > 0
	       const double shape // shape >  0
	       )
{
  if (location <= 0)
    {
      //Rprintf("subroutine:rpareto, location must be positive %f",location);
      error("subroutine:rpareto, location must be positive");
    }
  if (shape <= 0)
    {
      //Rprintf("subroutine:rpareto, shape must be positive %f",shape);
      error("subroutine:rpareto, shape must be positive");
    }
  // generate a random number from pareto distribution
  return( location/R_pow(runif(0.0,1.0),1/shape));

}

//rigamma
//function (n, alpha, beta) 
//{
//    if (alpha > 0 & beta > 0) 
//        1/rgamma(n = n, alpha, beta) // alpha is scale, beta is rate of gamma
//    else stop("rigamma: invalid parameters\n")
//}
//<environment: namespace:pscl>
double rigamma( const double shape, const double scale)
{
  // If X ~ gamma(shape,scale)
  // then 1/X ~ inv-gamma(shape,1/scale)
  // NOTE:: sugar rgamma(n,shape,scale) but rgamma in R default is rgamma(n,shape,rate) 
  if (scale <= 0 ) // if scale =NaN then any logical statement is true so this trap catches the NaN.
    {
      //Rprintf("subroutine:rigamma: scale must be positive!! %f",scale);
      error("subroutine:rigamma: scale must be positive!! ");
    }
  if (shape <= 0 )
    {
      //Rprintf("subroutine:rigamma: shape must be positive!! %f",shape);
      error("subroutine:rigamma: shape must be positive!!");
    }
  return(1/rgamma(
		  shape,// shape
		  1/scale// scale
		  )
	 );
}





double max_VecScal(const double vec[],
		   const double scal,
		   const int M
		   )  
{
  
  double tempDbl1 = scal; // This may need to be initalised different base on your data
  for (int i = 0; i < M; i++)
    {
      if (tempDbl1 < vec[i]) tempDbl1 = vec[i]; // vec.at(i) = vec[i] 
    }
  return tempDbl1;
}


double sum_inv(const double Klambdas[],int M)
{
  double tempDbl1 = 0;
  for (int i=0; i < M;i++)
    {
      tempDbl1 += 1.0/Klambdas[i];
    }
  return tempDbl1;
}




// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP gibbs( SEXP inpts, SEXP inpB, SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inpbphi, SEXP inpaGam, SEXP inpbGam, SEXP inpmax_dum, SEXP inpburnin,  SEXP inpprintFreq) ;
// }



SEXP gibbs( SEXP inpts, SEXP inpB, SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inpbphi, SEXP inpaGam, SEXP inpbGam, SEXP inpmax_dum, SEXP inpburnin,  SEXP inpprintFreq){

  /* ==== :Model description: =====
     ti | beta_i Klambda_i ~ kW(lambdai,beta_i)
     beta_i, Klambda_i | G ~ G
     G ~ DP(D,G0)
     where G0 =  unif(beta; 0.0001,phi) * igamma(Klambda;shape=dum,scale=gam) 
     gam     ~  gamma(shape=aGam*,rate=bGam*)
     phi     ~  pareto(shape=aphi*,location=bphi*) 
     D       ~  gamma(shape=a_D*,scale=1/ib_D*)
     dum     ~  unif(1.001,max_gum)
     
     * are hyperparameters 
     Kottas recommended dum* = 2 (this code treat dum as random quantity), aGam* = 1, aphi* = 2, a_D = 2, and ib_D = 0.9  
     
     The weibull pdf is parametrized as:
     k_W(t|Klambda,beta) = beta*t^{beta-1}/lambda * exp(- t^{beta}/Klambda)
     　　　
     　　The relationship with the standard weibull Weibul(shape=beta,scale=lambda) is: lambda = Klambda^{1/beta}
  */
  
  // === Importing values === // 
  const double *ts=REAL(inpts); // array of length N
  const int N = length(inpts);
  const int B = INTEGER(inpB)[0]; 
  const int M = INTEGER(M_)[0]; // step 1
  const double a_D = REAL(inpa_D)[0];                   // step 2
  const double ib_D = REAL(inpib_D)[0]; 
  const double aphi = REAL(inpaphi)[0];               // step 3
  const double bphi = REAL(inpbphi)[0];               
  const double aGam = REAL(inpaGam)[0];               // step 4  
  const double bGam = REAL(inpbGam)[0]; 
  const double max_dum = REAL(inpmax_dum)[0];   // default value should be 2 by Kottas
  const int burnin = INTEGER(inpburnin)[0];
  const int printFreq = INTEGER(inpprintFreq)[0];

 
  double phi=0, gam=0, D=0, beta=0, qi0=0, candi=0,tempDbl1=0, w=0,i1=0, i2=0,MHrate=0.0,dum=(max_dum+1.0)/2.0; 
  double can[2]={2.0,2.0}, ct_tot[2]={1.0,1.0},ct_acc[2]={1.0,1.0};
  int Nuniq = 0,tempInt1 = 0, iN=0,ih=0,idh=0, S[N];
  double Klambdas[M],betas[M],vs[M],postprob[M],weightS[M]; // change the size 
  const double min_beta = 0.1;  
  const double min_dum = 1.001;
  // beta is supposed to be from uniform (0,max_gum), however if beta<0.01, dwei(shape=beta,scale=Klambda) becomes undefined
  // for Klambda > 20,000. As our model should favor the set (beta, Klambda) = about (3,20,000), we must restrict any possible values of beta
  // to accept Klambda = 20,000. For this reason, we turncate the left hand extreme value of beta
  // Create initial values of Klambdas, betas:
  D = rgamma(
	     a_D, // shape 
	     1/ib_D  // scale
	     );
    
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1.0,D);
  vs[M-1]=1.0;
  
  weightS[0] = vs[0];
  w = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      w *= (1.0-vs[ih-1]);
      weightS[ih] = vs[ih]*w;
    }


  gam = rgamma(aGam,1/bGam);

  phi = rpareto(bphi, // location > 0
		aphi  // shape >  0
		); // returns double


  for (iN = 0 ; iN < N; iN++ ) S[iN] = 0;

  for (ih = 0 ; ih < M ; ih ++)
    {
      Klambdas[ih] =  rigamma(dum,gam); // return double
      betas[ih] = runif(
			min_beta, // min
			phi // max
			);
    }

  // Storage



  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP,12)); // result is stored in res
  SEXP post_S = allocVector(INTSXP, B*N); 
  SET_VECTOR_ELT(res, 0, post_S); 
  SEXP post_betas = allocVector(REALSXP, B*N); 
  SET_VECTOR_ELT(res, 1, post_betas); 
  SEXP post_Klambdas = allocVector(REALSXP, B*N);
  SET_VECTOR_ELT(res, 2, post_Klambdas);
 
  SEXP post_betas_uni = allocVector(REALSXP, B*M); 
  SET_VECTOR_ELT(res, 3, post_betas_uni); 
  SEXP post_Klambdas_uni = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 4, post_Klambdas_uni); 
  SEXP post_vs = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 5, post_vs); 
  SEXP post_weightS = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 6, post_weightS); 

  SEXP post_D = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 7, post_D); 
  SEXP post_phi = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 8, post_phi); 
  SEXP post_gam = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 9, post_gam); 
  SEXP post_dum = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 10, post_dum); 

  SEXP post_ct_accept = allocVector(REALSXP,4); // ct_integ_tot:
  SET_VECTOR_ELT(res, 11, post_ct_accept);

  GetRNGstate();
  
  // Let's Roll! MCMC! 
  for (int iB = 0 ; iB < B ; iB++)
    {
      
      //  Update S, a vector of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of S has a categorical distribution with M categories. 
      //  if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
      for (iN = 0 ; iN < N ; iN++)
	{ 
	  // Given patient iN, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dWEI2(ts[iN],betas[ih],Klambdas[ih],0)*weightS[ih];
	  S[iN] = sample(postprob,M); 
	  
	}
      // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightS[0] = vs[0]
      // weightS[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      w = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  i1 = 0.0;    // # patients with S == ih
	  i2 = 0.0; // # patients with S > ih
	  for (iN = 0 ; iN < N; iN++ ) 
	    {
	      if ( S[iN]==ih ) i1++ ;
	      if ( S[iN] > ih ) i2++ ;
	    }
	  vs[ih] = rbeta(
			 1.0 + i1, 
			 D + i2 
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+i1+idh ,  D + (double)i2 + (double)ivec);
	  if (ih > 0) w *= (1.0-vs[ih-1]);
	  weightS[ih] = vs[ih]*w;
	}
      weightS[M-1] = w*(1-vs[M-2]);


      // Updating Klambdas[ih] and betas[ih]
      for (ih=0;ih < M ; ih ++)
	{
	  // Update Klambdas
	  i1 = 0.0;
	  i2 = 0.0;
	  for (iN = 0 ; iN < N ; iN ++) 
	    {
	      if (S[iN]==ih) 
		{
		  i1++; // = sum I (S[ih]==ih)
		  i2 += R_pow(ts[iN],betas[ih]);
		}
	    }
	  idh = 0.0;
	  while(!R_FINITE(dWEI2(6.0,betas[ih],Klambdas[ih],0)) || idh==0)  
	    {
	      // All the randomly generated Klambdas are rejected if they have infinite v1lue of the density function
	      // evaluated at 6.0.
	      Klambdas[ih]=rigamma(dum+i1,gam+i2);
	      if (Klambdas[ih] < 0.001) Klambdas[ih] = 0.001;
	      if (Klambdas[ih] > 1.0e+15)Klambdas[ih] = 1.0e+15;
	      idh++;
	      //if (idh== 3) Rprintf(" fail!");
	      if (idh == 5000 )
		{
		  //Rprintf("Klambdas[ih] %f betas[ih] %f The prior specification is unrealistic. MCMC does not converge.",
		  //Klambdas[ih],betas[ih]);
		  error("The prior specification is unrealistic. MCMC does not converge.");
		}
	    }
	  // Update betas M-H
	  candi = rnorm(betas[ih],can[0]);
	  if (weightS[ih] > 0.1) ct_tot[0]++;
	  // the algorithm requires to compute Klambda^beta hence Klambda and beta must not be too small
	  if (candi > min_beta && candi < phi)  
	    {
	      MHrate = 0.0;  // If there is no obs to h^th component, candi is always accepted if it lies within [0.0001, phi]
	      for (iN = 0 ; iN < N ; iN++)
		{
		  if (S[iN]==ih)
		    MHrate += dWEI2(ts[iN],candi,Klambdas[ih],1) - dWEI2(ts[iN],betas[ih],Klambdas[ih],1);
		}
	      // If the proposed beta does not return finite value of density of Weibull, such proposal is not accepted
	      if (runif(0.0,1.0) < exp(MHrate) && R_FINITE(dWEI2(6.0,candi,Klambdas[ih],0)) )
		{
		  if (weightS[ih] > 0.1) ct_acc[0]++;
		  betas[ih] = candi;
		}
	    }
	}
      // Update phi 
      /* Step 3: update phi
	 phi  ~ pareto(shape=aphi*,location=bphi*) 
	 => phi |  betas ~ pareto(shape=aphi+Nuniq,location=max(betas,bphi)) */
      phi = rpareto(max_VecScal(betas,bphi,M), //  const double location, // location > 0
		    M + aphi             //  const double shape // shape >  0
		    );

      /* Update gamma:
	 gam     ~ gamma(shape=aGam*,rate=bGam*)
	 => gam| lambda ~ gamma(shape=dum+) */
      
      gam = rgamma(
		   aGam+dum*M, // shape 
		   1/(bGam+sum_inv(Klambdas,M))  // scale  // sum_inv(Klambda) = sum(1/Klambda)
		   );
      
      /* Update dum:
	 dum ~ unif(0,max_dum) via M-H
       */
      MHrate =0.0;
      candi = rnorm(dum,can[1]);
      ct_tot[1]++;
      if (min_dum <= candi && candi <= max_dum )
	{
	  for (ih = 0; ih < M ; ih++)
	    MHrate += dinvgamma(Klambdas[ih], candi, gam,1) - dinvgamma(Klambdas[ih], dum, gam,1);
	  if ( runif(0,1.0) < exp(MHrate) )
	    {
	      dum = candi;
	      ct_acc[1]++;
	    }
	}
      // Update D
      // D ~ gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
      w = 0.0;
      for (ih=0;ih<(M-1);ih++) w += log(1.0-vs[ih]);
      D = rgamma( (a_D+ (double) M-1.0), (1.0/(-w+ib_D)) ); // N+w is the number of random var following DP given js
      if (ISNAN(D)){ 
	D = 1.0; 
	Rprintf("\v ISNAN(D): a_D %f ib_D %f sum log(1-vs[ih]) %f ",a_D,ib_D,w);
      }
      if( D < 0.01) D = 0.01; // if D is too small, rbeta(1,a,0) returns NaN for any a
      
      if (iB < burnin)
	{
	  for (ih = 0; ih < 2 ; ih++)
	    {
	      w = ct_acc[ih]/ct_tot[ih];
	      if (can[ih] > 0.001 && can[ih] < 15)
		{
		  if (w < 0.2) can[ih] *= 0.7;
		  else if (w > 0.5) can[ih] *=1.3;
		}
	      else if (can[ih] < 0.001) can[ih] = 0.001;
	      else if (can[ih] > 15) can[ih] = 15;
	    }
	}
      else if (iB==burnin)
	{
	  for (ih = 0; ih < 2 ; ih++)
	    {
	      ct_tot[ih] = 0.0;
	      ct_acc[ih] = 0.0;
	    }
	}

      // Store the result from this iteration
      for (iN = 0 ; iN < N ; iN++ )
	{
	  idh = iN + iB * N;
	  INTEGER(post_S)[idh] = S[iN];
	  REAL(post_betas)[idh] = betas[S[iN]];
	  REAL(post_Klambdas)[idh] = Klambdas[S[iN]]; 
	}
      for (ih = 0 ; ih < M; ih++ )
	{
	  idh = ih + iB * M;
	  REAL(post_betas_uni)[idh] = betas[ih];
	  REAL(post_Klambdas_uni)[idh] = Klambdas[ih]; 
	  REAL(post_weightS)[idh] = weightS[ih];
	  REAL(post_vs)[idh] = vs[ih]; 
	}
      REAL(post_D)[iB] = D;
      REAL(post_phi)[iB] = phi;
      REAL(post_gam)[iB] = gam;
      REAL(post_dum)[iB] = dum;
	
      
      if(iB % printFreq == 0) Rprintf(" %d iterations are done...   ", iB);
      
    }
  for (ih = 0 ; ih < 2 ; ih ++)
    {
      REAL(post_ct_accept)[2*ih] = ct_acc[ih]/ct_tot[ih];
      REAL(post_ct_accept)[2*ih+1] = can[ih];
    }
  PutRNGstate();
  UNPROTECT(1);
  return res;

}





// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP gibbsUnif( SEXP inpts, SEXP inpB, SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inpbphi, 
// 		  //SEXP inpaGam, SEXP inpbGam, 
// 		  SEXP inpmax_dum, SEXP inpburnin,  SEXP inpprintFreq) ;
// }



SEXP gibbsUnif( SEXP inpts, SEXP inpB, SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inpbphi, 
		//SEXP inpaGam, SEXP inpbGam, 
		SEXP inpmax_dum, SEXP inpburnin,  SEXP inpprintFreq){

  /* ==== :Model description: =====
     ti | beta_i Klambda_i ~ kW(lambdai,beta_i)
     beta_i, Klambda_i | G ~ G
     G ~ DP(D,G0)
     where G0 =  unif(beta; 0.0001,phi) * igamma(Klambda;shape=dum,scale=gam) 
     gam     ~  gamma(shape=aGam*,rate=bGam*)
     phi     ~  pareto(shape=aphi*,location=bphi*) 
     D       ~  gamma(shape=a_D*,scale=1/ib_D*)
     dum     ~  unif(1.001,max_gum)
     
     * are hyperparameters 
     Kottas recommended dum* = 2 (this code treat dum as random quantity), aGam* = 1, aphi* = 2, a_D = 2, and ib_D = 0.9  
     
     The weibull pdf is parametrized as:
     k_W(t|Klambda,beta) = beta*t^{beta-1}/lambda * exp(- t^{beta}/Klambda)
     　　　
     　　The relationship with the standard weibull Weibul(shape=beta,scale=lambda) is: lambda = Klambda^{1/beta}
  */
  
  // === Importing values === // 
  const double *ts=REAL(inpts); // array of length N
  const int N = length(inpts);
  const int B = INTEGER(inpB)[0]; 
  const int M = INTEGER(M_)[0]; // step 1
  const double a_D = REAL(inpa_D)[0];                   // step 2
  const double ib_D = REAL(inpib_D)[0]; 
  const double aphi = REAL(inpaphi)[0];               // step 3
  const double bphi = REAL(inpbphi)[0];               
  //const double aGam = REAL(inpaGam)[0];               // step 4  
  //const double bGam = REAL(inpbGam)[0]; 
  const double max_dum = REAL(inpmax_dum)[0];   // default value should be 2 by Kottas
  const int burnin = INTEGER(inpburnin)[0];
  const int printFreq = INTEGER(inpprintFreq)[0];

 
  double phi=0, gam=0, D=0, beta=0, qi0=0, candi=0,tempDbl1=0, w=0,i1=0, i2=0,MHrate=0.0,dum=(max_dum+1.0)/2.0; 
  const int Nacc = 4;
  double can[Nacc], ct_tot[Nacc],ct_acc[Nacc];
  can[0]=2.0;can[1]=2.0;can[2]=15.0;can[3]=2.0;
  ct_tot[0]=1.0;ct_tot[1]=1.0;ct_tot[2]=1.0;ct_tot[3]=1.0;
  ct_acc[0]=1.0;ct_acc[1]=1.0;ct_acc[2]=1.0;ct_acc[3]=1.0;

  int Nuniq = 0,tempInt1 = 0, iN=0,ih=0,idh=0, S[N];
  double Klambdas[M],betas[M],vs[M],postprob[M],weightS[M]; // change the size 
  const double min_beta = 0.1;  
  const double min_dum = 1.001;
  double maxCan[Nacc]; // beta, dum, gam, D
  maxCan[0]=15;maxCan[1]=15;maxCan[2]=R_PosInf;maxCan[3]=15;
  // beta is supposed to be from uniform (0,max_gum), however if beta<0.01, dwei(shape=beta,scale=Klambda) becomes undefined
  // for Klambda > 20,000. As our model should favor the set (beta, Klambda) = about (3,20,000), we must restrict any possible values of beta
  // to accept Klambda = 20,000. For this reason, we turncate the left hand extreme value of beta
  // Create initial values of Klambdas, betas:

  // Unlike gibbs code 
  // - D ~ unif(a_D,ib_D)
  // - gam ~ unif(aGam,bGam)
  D = runif(
	    a_D, // min 
	    ib_D  // max
	    );
    
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1.0,D);
  vs[M-1]=1.0;
  
  weightS[0] = vs[0];
  w = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      w *= (1.0-vs[ih-1]);
      weightS[ih] = vs[ih]*w;
    }


  gam = runif(1,1E+10);

  phi = rpareto(bphi, // location > 0
		aphi  // shape >  0
		); // returns double


  for (iN = 0 ; iN < N; iN++ ) S[iN] = 0;

  for (ih = 0 ; ih < M ; ih ++)
    {
      Klambdas[ih] =  rigamma(dum,gam); // return double
      betas[ih] = runif(
			min_beta, // min
			phi // max
			);
    }

  // Storage



  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP,13)); // result is stored in res
  SEXP post_S = allocVector(INTSXP, B*N); 
  SET_VECTOR_ELT(res, 0, post_S); 
  SEXP post_betas = allocVector(REALSXP, B*N); 
  SET_VECTOR_ELT(res, 1, post_betas); 
  SEXP post_Klambdas = allocVector(REALSXP, B*N);
  SET_VECTOR_ELT(res, 2, post_Klambdas);
 
  SEXP post_betas_uni = allocVector(REALSXP, B*M); 
  SET_VECTOR_ELT(res, 3, post_betas_uni); 
  SEXP post_Klambdas_uni = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 4, post_Klambdas_uni); 
  SEXP post_vs = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 5, post_vs); 
  SEXP post_weightS = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 6, post_weightS); 

  SEXP post_D = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 7, post_D); 
  SEXP post_phi = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 8, post_phi); 
  SEXP post_gam = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 9, post_gam); 
  SEXP post_dum = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 10, post_dum); 

  SEXP post_AR = allocVector(REALSXP,Nacc); // ct_integ_tot:
  SET_VECTOR_ELT(res, 11, post_AR);
  
  SEXP post_can = allocVector(REALSXP,Nacc); // ct_integ_tot:
  SET_VECTOR_ELT(res, 12, post_can);
  
  GetRNGstate();
  
  // Let's Roll! MCMC! 
  for (int iB = 0 ; iB < B ; iB++)
    {
      
      //  Update S, a vector of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of S has a categorical distribution with M categories. 
      //  if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
      for (iN = 0 ; iN < N ; iN++)
	{ 
	  // Given patient iN, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dWEI2(ts[iN],betas[ih],Klambdas[ih],0)*weightS[ih];
	  S[iN] = sample(postprob,M); 
	  
	}
      // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightS[0] = vs[0]
      // weightS[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      w = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  i1 = 0.0;    // # patients with S == ih
	  i2 = 0.0; // # patients with S > ih
	  for (iN = 0 ; iN < N; iN++ ) 
	    {
	      if ( S[iN]==ih ) i1++ ;
	      if ( S[iN] > ih ) i2++ ;
	    }
	  vs[ih] = rbeta(
			 1.0 + i1, 
			 D + i2 
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+i1+idh ,  D + (double)i2 + (double)ivec);
	  if (ih > 0) w *= (1.0-vs[ih-1]);
	  weightS[ih] = vs[ih]*w;
	}
      weightS[M-1] = w*(1-vs[M-2]);


      // Updating Klambdas[ih] and betas[ih]
      for (ih=0;ih < M ; ih ++)
	{
	  // Update Klambdas
	  i1 = 0.0;
	  i2 = 0.0;
	  for (iN = 0 ; iN < N ; iN ++) 
	    {
	      if (S[iN]==ih) 
		{
		  i1++; // = sum I (S[ih]==ih)
		  i2 += R_pow(ts[iN],betas[ih]);
		}
	    }
	  idh = 0.0;
	  while(!R_FINITE(dWEI2(6.0,betas[ih],Klambdas[ih],0)) || idh==0)  
	    {
	      // All the randomly generated Klambdas are rejected if they have infinite v1lue of the density function
	      // evaluated at 6.0.
	      Klambdas[ih]=rigamma(dum+i1,gam+i2);
	      if (Klambdas[ih] < 0.001) Klambdas[ih] = 0.001;
	      if (Klambdas[ih] > 1.0e+15)Klambdas[ih] = 1.0e+15;
	      idh++;
	      //if (idh== 3) Rprintf(" fail!");
	      if (idh == 5000 )
		{
		  //Rprintf("Klambdas[ih] %f betas[ih] %f The prior specification is unrealistic. MCMC does not converge.",
		  //Klambdas[ih],betas[ih]);
		  error("The prior specification is unrealistic. MCMC does not converge.");
		}
	    }
	  // Update betas M-H
	  candi = rnorm(betas[ih],can[0]);
	  if (weightS[ih] > 0.1) ct_tot[0]++;
	  // the algorithm requires to compute Klambda^beta hence Klambda and beta must not be too small
	  if (candi > min_beta && candi < phi)  
	    {
	      MHrate = 0.0;  // If there is no obs to h^th component, candi is always accepted if it lies within [0.0001, phi]
	      for (iN = 0 ; iN < N ; iN++)
		{
		  if (S[iN]==ih)
		    MHrate += dWEI2(ts[iN],candi,Klambdas[ih],1) - dWEI2(ts[iN],betas[ih],Klambdas[ih],1);
		}
	      // If the proposed beta does not return finite value of density of Weibull, such proposal is not accepted
	      if (runif(0.0,1.0) < exp(MHrate) && R_FINITE(dWEI2(6.0,candi,Klambdas[ih],0)) )
		{
		  if (weightS[ih] > 0.1) ct_acc[0]++;
		  betas[ih] = candi;
		}
	    }
	}
      // Update phi 
      /* Step 3: update phi
	 phi  ~ pareto(shape=aphi*,location=bphi*) 
	 => phi |  betas ~ pareto(shape=aphi+Nuniq,location=max(betas,bphi)) */
      phi = rpareto(max_VecScal(betas,bphi,M), //  const double location, // location > 0
		    M + aphi             //  const double shape // shape >  0
		    );

      /* Update gamma:
	 gamma ~ unif(aGam,bGam) via M-H
      */
      MHrate =0.0;
      candi = rnorm(gam,can[2]);
      ct_tot[2]++;
      if (1 <= candi) // any number is possible flat prior
	{
	  for (ih = 0; ih < M ; ih++)
	    MHrate += dinvgamma(Klambdas[ih], dum, candi,1) - dinvgamma(Klambdas[ih], dum, gam,1);
	  if ( runif(0,1.0) < exp(MHrate) )
	    {
	      gam = candi;
	      ct_acc[2]++;
	    }
	}

      
      /* Update dum:
	 dum ~ unif(0,max_dum) via M-H
       */
      MHrate =0.0;
      candi = rnorm(dum,can[1]);
      ct_tot[1]++;
      if (min_dum <= candi && candi <= max_dum )
	{
	  for (ih = 0; ih < M ; ih++)
	    MHrate += dinvgamma(Klambdas[ih], candi, gam,1) - dinvgamma(Klambdas[ih], dum, gam,1);
	  if ( runif(0,1.0) < exp(MHrate) )
	    {
	      dum = candi;
	      ct_acc[1]++;
	    }
	}
      // Update D
      // D ~ gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  

      MHrate =0.0;
      candi = rnorm(D,can[3]);
      ct_tot[3]++;
      if (a_D <= candi && candi <= ib_D )
	{
	  for (ih = 0; ih < M-1 ; ih++)
	    MHrate += dbeta(vs[ih], 1.0, candi,1) - dbeta(vs[ih], 1.0, D,1);
	  if ( runif(0,1.0) < exp(MHrate) )
	    {
	      D = candi;
	      ct_acc[3]++;
	    }
	}

      if (iB < burnin)
	{
	  for (ih = 0; ih < Nacc ; ih++)
	    {
	      w = ct_acc[ih]/ct_tot[ih];
	      if (can[ih] > 0.001 && can[ih] < maxCan[ih])
		{
		  if (w < 0.2) can[ih] *= 0.7;
		  else if (w > 0.5) can[ih] *=1.3;
		}
	      else if (can[ih] < 0.001) can[ih] = 0.001;
	      else if (can[ih] > maxCan[ih]) can[ih] = maxCan[ih];
	    }
	}
      else if (iB==burnin)
	{
	  for (ih = 0; ih < Nacc ; ih++)
	    {
	      ct_tot[ih] = 0.0;
	      ct_acc[ih] = 0.0;
	    }
	}

      // Store the result from this iteration
      for (iN = 0 ; iN < N ; iN++ )
	{
	  idh = iN + iB * N;
	  INTEGER(post_S)[idh] = S[iN];
	  REAL(post_betas)[idh] = betas[S[iN]];
	  REAL(post_Klambdas)[idh] = Klambdas[S[iN]]; 
	}
      for (ih = 0 ; ih < M; ih++ )
	{
	  idh = ih + iB * M;
	  REAL(post_betas_uni)[idh] = betas[ih];
	  REAL(post_Klambdas_uni)[idh] = Klambdas[ih]; 
	  REAL(post_weightS)[idh] = weightS[ih];
	  REAL(post_vs)[idh] = vs[ih]; 
	}
      REAL(post_D)[iB] = D;
      REAL(post_phi)[iB] = phi;
      REAL(post_gam)[iB] = gam;
      REAL(post_dum)[iB] = dum;
	
      
      if(iB % printFreq == 0) Rprintf(" %d iterations are done...   ", iB);
      
    }
  for (ih = 0 ; ih < Nacc ; ih ++)
    {
      REAL(post_AR)[ih] = ct_acc[ih]/ct_tot[ih];
      REAL(post_can)[ih] = can[ih];
    }
  PutRNGstate();
  UNPROTECT(1);
  return res;

}




// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP gibbsAllUnif( SEXP inpts, SEXP inpB, SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inploc_phi, 
// 		     SEXP inpaGam, SEXP inploc_bGam, 
// 		     //SEXP inpmax_dum, 
// 		     SEXP inpburnin,  SEXP inpprintFreq) ;
// }



SEXP gibbsAllUnif( SEXP inpts, SEXP inpB, SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inploc_phi, 
		   SEXP inpaGam, SEXP inploc_Gam, 
		   // SEXP inpmax_dum, 
		   SEXP inpburnin,  SEXP inpprintFreq){

  /* ==== :Model description: =====
     ti | beta_i Klambda_i ~ weibull(lambdai,beta_i)
     beta_i, Klambda_i | G ~ G
     G ~ DP(D,G0)
     where G0 =  unif(beta; 0.0001,phi) * unif(Klambda;0.001,scale=gam) 
     phi     ~  pareto(shape=aphi*,location=loc_phi*) 
     gam     ~  pareto(shape=aGam*,location=loc_Gam*)
     D       ~  unif(a_D*,ib_D*)
     
     * are hyperparameters 
     Kottas recommended dum* = 2 (this code treat dum as random quantity), aGam* = 1, aphi* = 2, a_D = 2, and ib_D = 0.9  
     
     The weibull pdf is parametrized as:
     k_W(t|Klambda,beta) = beta*t^{beta-1}/lambda * exp(- t^{beta}/Klambda)
     　　　
     　　The relationship with the standard weibull Weibul(shape=beta,scale=lambda) is: lambda = Klambda^{1/beta}
  */
  
  // === Importing values === // 
  const double *ts=REAL(inpts); // array of length N
  const int N = length(inpts);
  const int B = INTEGER(inpB)[0]; 
  const int M = INTEGER(M_)[0]; // step 1
  const double a_D = REAL(inpa_D)[0];                   // step 2
  const double ib_D = REAL(inpib_D)[0]; 
  const double aphi = REAL(inpaphi)[0];               // step 3
  const double loc_phi = REAL(inploc_phi)[0];               
  const double aGam = REAL(inpaGam)[0];               // step 4  
  const double loc_Gam = REAL(inploc_Gam)[0]; 
  //const double max_dum = REAL(inpmax_dum)[0];   // default value should be 2 by Kottas
  const int burnin = INTEGER(inpburnin)[0];
  const int printFreq = INTEGER(inpprintFreq)[0];

 
  double phi=0, gam=0, D=0, beta=0, candi=0,tempDbl1=0, w=0,i1=0, i2=0,MHrate=0.0; 
  const int Nacc = 3;
  // beta, Lambda, D 
  double can[Nacc], ct_tot[Nacc],ct_acc[Nacc];
  can[0]=0.5;can[1]=0.5;can[2]=0.5;
  ct_tot[0]=1.0;ct_tot[1]=1.0;ct_tot[2]=1.0;
  ct_acc[0]=1.0;ct_acc[1]=1.0;ct_acc[2]=1.0;
  int iN=0,ih=0,idh=0, S[N];
  double Lambdas[M],betas[M],vs[M],postprob[M],weightS[M]; // change the size 
  const double smallValue = 0.1;  
  double maxCan[Nacc]; // beta, Lambda, D 
  maxCan[0]=10;maxCan[1]=10;maxCan[2]=1.2;
  double minCan[Nacc]; // beta, Lambda, D 
  minCan[0]=0.001;minCan[1]=0.001;minCan[2]=0.1;
  // beta is supposed to be from uniform (0,max_gum), however if beta<0.01, dwei(shape=beta,scale=Lambda) becomes undefined
  // for Lambda > 20,000. As our model should favor the set (beta, Lambda) = about (3,20,000), we must restrict any possible values of beta
  // to accept Lambda = 20,000. For this reason, we turncate the left hand extreme value of beta
  // Create initial values of Lambdas, betas:

  // Unlike gibbs code 
  // - D ~ unif(a_D,ib_D)
  // - gam ~ unif(aGam,loc_Gam)
  D = runif(
	    a_D, // min 
	    ib_D  // max
	    );
    
  for (ih = 0 ; ih < (M-1) ; ih ++) vs[ih] = rbeta(1.0,D);
  vs[M-1]=1.0;
  
  weightS[0] = vs[0];
  w = 1.0 ;
  for (ih=1; ih < M; ih++)
    {
      w *= (1.0-vs[ih-1]);
      weightS[ih] = vs[ih]*w;
    }


  gam = rpareto(loc_Gam,
		aGam);

  phi = rpareto(loc_phi, // location > 0
		aphi  // shape >  0
		); // returns double


  for (iN = 0 ; iN < N; iN++ ) S[iN] = 0;

  for (ih = 0 ; ih < M ; ih ++)
    {
      Lambdas[ih] =  runif(smallValue,
			    gam); // return double
      betas[ih] = runif(
			smallValue, // min
			phi // max
			);
    }

  // Storage



  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP,12)); // result is stored in res
  SEXP post_S = allocVector(INTSXP, B*N); 
  SET_VECTOR_ELT(res, 0, post_S); 
  SEXP post_betas = allocVector(REALSXP, B*N); 
  SET_VECTOR_ELT(res, 1, post_betas); 
  SEXP post_Lambdas = allocVector(REALSXP, B*N);
  SET_VECTOR_ELT(res, 2, post_Lambdas);
 
  SEXP post_betas_uni = allocVector(REALSXP, B*M); 
  SET_VECTOR_ELT(res, 3, post_betas_uni); 
  SEXP post_Lambdas_uni = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 4, post_Lambdas_uni); 
  SEXP post_vs = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 5, post_vs); 
  SEXP post_weightS = allocVector(REALSXP, B*M);
  SET_VECTOR_ELT(res, 6, post_weightS); 

  SEXP post_D = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 7, post_D); 
  SEXP post_phi = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 8, post_phi); 
  SEXP post_gam = allocVector(REALSXP, B);
  SET_VECTOR_ELT(res, 9, post_gam); 
  //SEXP post_dum = allocVector(REALSXP, B);
  //SET_VECTOR_ELT(res, 10, post_dum); 

  SEXP post_AR = allocVector(REALSXP,Nacc); // ct_integ_tot:
  SET_VECTOR_ELT(res, 10, post_AR);
  
  SEXP post_can = allocVector(REALSXP,Nacc); // ct_integ_tot:
  SET_VECTOR_ELT(res, 11, post_can);
  
  GetRNGstate();
  
  // Let's Roll! MCMC! 
  for (int iB = 0 ; iB < B ; iB++)
    {
      
      //  Update S, a vector of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of S has a categorical distribution with M categories. 
      //  if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
      for (iN = 0 ; iN < N ; iN++)
	{ 
	  // Given patient iN, compute P(H_1i=ih;-) for each ih 
	  // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	  for (ih=0; ih < M; ih++ ) postprob[ih] = dweibull(ts[iN],betas[ih],Lambdas[ih],0)*weightS[ih];
	  S[iN] = sample(postprob,M); 
	  
	}
      // Update vs, a vector of length M, containing the latent variable v, used to construct pis. 
      // Once vs[ih] is updated, update the probabilities of the categorical distribution of H1 based on the formula:
      // weightS[0] = vs[0]
      // weightS[ih]= vs[ih]*(1-vh[ih-1])**(1-vh[ih-2])*...**(1-vh[0]) for ih =1,2,...M-1
      w = 1.0;
      for (ih = 0 ; ih < (M-1) ; ih ++) // vs[M-1] = 1 always!
	{
	  i1 = 0.0;    // # patients with S == ih
	  i2 = 0.0; // # patients with S > ih
	  for (iN = 0 ; iN < N; iN++ ) 
	    {
	      if ( S[iN]==ih ) i1++ ;
	      if ( S[iN] > ih ) i2++ ;
	    }
	  vs[ih] = rbeta(
			 1.0 + i1, 
			 D + i2 
			 );
	  //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+i1+idh ,  D + (double)i2 + (double)ivec);
	  if (ih > 0) w *= (1.0-vs[ih-1]);
	  weightS[ih] = vs[ih]*w;
	}
      weightS[M-1] = w*(1-vs[M-2]);


      // Updating Lambdas[ih] and betas[ih]
      for (ih=0;ih < M ; ih ++)
	{
	  // Update Lambdas M-H
	  candi = rnorm(Lambdas[ih],can[1]);
	  if (weightS[ih] > 0.1) ct_tot[1]++;
	  // the algorithm requires to compute Lambda hence Lambda and beta must not be too small
	  if (candi > smallValue && candi < gam)  
	    {
	      MHrate = 0.0;  // If there is no obs to h^th component, candi is always accepted if it lies within [0.0001, phi]
	      for (iN = 0 ; iN < N ; iN++)
		{
		  if (S[iN]==ih)
		    MHrate += dweibull(ts[iN],betas[ih],candi,1) - dweibull(ts[iN],betas[ih],Lambdas[ih],1);
		}
	      // If the proposed beta does not return finite value of density of Weibull, such proposal is not accepted
	      if (runif(0.0,1.0) < exp(MHrate) && R_FINITE(dweibull(6.0,candi,Lambdas[ih],0)) )
		{
		  if (weightS[ih] > 0.1) ct_acc[1]++;
		  Lambdas[ih] = candi;
		}
	    }

	  // Update betas M-H
	  candi = rnorm(betas[ih],can[0]); 
	  if (weightS[ih] > 0.1) ct_tot[0]++;
	  // the algorithm requires to compute Lambda^beta hence Lambda and beta must not be too small
	  if (candi > smallValue && candi < phi)  
	    {
	      MHrate = 0.0;  // If there is no obs to h^th component, candi is always accepted if it lies within [0.0001, phi]
	      for (iN = 0 ; iN < N ; iN++)
		{
		  if (S[iN]==ih)
		    MHrate += dweibull(ts[iN],candi,Lambdas[ih],1) - dweibull(ts[iN],betas[ih],Lambdas[ih],1);
		}
	      // If the proposed beta does not return finite value of density of Weibull, such proposal is not accepted
	      if (runif(0.0,1.0) < exp(MHrate) && R_FINITE(dWEI2(6.0,candi,Lambdas[ih],0)) )
		{
		  if (weightS[ih] > 0.1) ct_acc[0]++;
		  betas[ih] = candi;
		}
	    }
	}
      // Update phi 
      /* Step 3: update phi
	 phi  ~ pareto(shape=aphi*,location=loc_phi*) 
	 => phi |  betas ~ pareto(shape=aphi+Nuniq,location=max(betas,loc_phi)) */
      phi = rpareto(max_VecScal(betas,loc_phi,M), //  const double location, // location > 0
		    M + aphi             //  const double shape // shape >  0
		    );

      /* Update gamma:
	 gamma ~pareto(shape=aGam,location=loc_Gam) via M-H
      */
      gam = rpareto(max_VecScal(Lambdas,loc_Gam,M), //  const double location, // location > 0
		    M + aGam             //  const double shape // shape >  0
		    );

      // Update D
      // D ~ gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
 
      MHrate =0.0;
      candi = rnorm(D,can[2]);
      ct_tot[2]++;
      if (a_D <= candi && candi <= ib_D )
	{
	  for (ih = 0; ih < M-1 ; ih++)
	    MHrate += dbeta(vs[ih], 1.0, candi,1) - dbeta(vs[ih], 1.0, D,1);
	  if ( runif(0,1.0) < exp(MHrate) )
	    {
	      D = candi;
	      ct_acc[2]++;
	    }
	}

      if (iB < burnin)
	{
	  for (ih = 0; ih < Nacc ; ih++)
	    {
	      w = ct_acc[ih]/ct_tot[ih];
	      if (can[ih] > minCan[ih] && can[ih] < maxCan[ih])
		{
		  if (w < 0.3) can[ih] *= 0.97;
		  else if (w > 0.6) can[ih] *=1.01;
		}
	      else if (can[ih] < minCan[ih]) can[ih] = minCan[ih];
	      else if (can[ih] > maxCan[ih]) can[ih] = maxCan[ih];
	    }
	}
      else if (iB==burnin)
	{
	  for (ih = 0; ih < Nacc ; ih++)
	    {
	      ct_tot[ih] = 0.0;
	      ct_acc[ih] = 0.0;
	    }
	}

      // Store the result from this iteration
      for (iN = 0 ; iN < N ; iN++ )
	{
	  idh = iN + iB * N;
	  INTEGER(post_S)[idh] = S[iN];
	  REAL(post_betas)[idh] = betas[S[iN]];
	  REAL(post_Lambdas)[idh] = Lambdas[S[iN]]; 
	}
      for (ih = 0 ; ih < M; ih++ )
	{
	  idh = ih + iB * M;
	  REAL(post_betas_uni)[idh] = betas[ih];
	  REAL(post_Lambdas_uni)[idh] = Lambdas[ih]; 
	  REAL(post_weightS)[idh] = weightS[ih];
	  REAL(post_vs)[idh] = vs[ih]; 
	}
      REAL(post_D)[iB] = D;
      REAL(post_phi)[iB] = phi;
      REAL(post_gam)[iB] = gam;
      //REAL(post_dum)[iB] = dum;
	
      
      if(iB % printFreq == 0) Rprintf(" %d iterations are done...   ", iB);
      
    }
  for (ih = 0 ; ih < Nacc ; ih ++)
    {
      REAL(post_AR)[ih] = ct_acc[ih]/ct_tot[ih];
      REAL(post_can)[ih] = can[ih];
    }
  PutRNGstate();
  UNPROTECT(1);
  return res;

}







double mixtureP(double pr[],const double ti, const double betas[], const double Klambdas[],const int iB, const int M,const double alpha)
{
  // cdf of mixture 
  double temp = 0;
  int idh=0;
  for (int ih =0; ih < M ; ih++) 
    //rN[0]/L*dWEI2(ts[it],Klambdas[M*iB+0],betas[M*iB+0],0))+...+rN[M-1]/L*dWEI2(ts[it],Klambdas[M*iB+M-1],betas[M*iB+M-1],0))
    {
      idh = ih+M*iB;
      if (pr[ih]>0)
	temp += (double) pr[ih]*pWEI2(ti,betas[idh],Klambdas[idh],0);
    }
  return temp-alpha;
}



// extern "C"{
//   SEXP getDens(SEXP ts_, SEXP prob_, SEXP Klambdas_, SEXP betas_, SEXP M_, SEXP B_,SEXP alpha_,SEXP densi_);
// }

SEXP getDens(SEXP ts_, SEXP weightS_, SEXP Klambdas_, SEXP betas_, SEXP M_, SEXP B_,SEXP alpha_,SEXP densi_)
{
  const double *ts= REAL(ts_);
  const int lt = length(ts_);
  const int M = INTEGER(M_)[0]; // truncation value 
  const int B = INTEGER(B_)[0]; // truncation value 
  const double *Klambdas = REAL(Klambdas_); // vector of length B*M
  const double *betas = REAL(betas_); // vector of length B*M
  const double *weightS = REAL(weightS_);
  const double alpha = REAL(alpha_)[0];
  const int densi = INTEGER(densi_)[0];
  //int rN[M], 
  int idh, OK;
  double pr[M],temp=0, X1=0.0001, X2=10.0, Y1,Y2,X,Y;

  SEXP res = PROTECT(allocVector(VECSXP,3)); // result is stored in res
  SEXP mat = allocVector(REALSXP, B*lt); 
  SET_VECTOR_ELT(res, 0, mat); 
  SEXP quant = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, quant); 
  SEXP NconvItr = allocVector(INTSXP, B); 
  SET_VECTOR_ELT(res, 2, NconvItr); 

  GetRNGstate();

  for (int iB=0; iB < B; iB++)
    {
      for (int ih=0;ih < M ; ih++) pr[ih]=weightS[ih+M*iB];
      // rN is a vector of length M containing the number of draws from each of the M components
      // sum(rN) = L;

      //rmultinom(L,pr,M,rN);

      if (densi||alpha==0) // approximate the density curve
	{
	  for (int it = 0 ; it < lt ; it++)
	    {
	      temp = 0;
	      for (int ih =0; ih < M ; ih++) 
		//rN[0]/L*dWEI2(ts[it],Klambdas[M*iB+0],betas[M*iB+0],0))+...+rN[M-1]/L*dWEI2(ts[it],Klambdas[M*iB+M-1],betas[M*iB+M-1],0))
		{
		  idh = ih+M*iB;
		  if (pr[ih]>0)
		    temp += (double) pr[ih]*dWEI2(ts[it],betas[idh],Klambdas[idh],0);
		}
	      REAL(mat)[it+iB*lt]=temp; ///((double)L);
	    }
	}
      if (alpha!=0)
	{
	  // bisection algorithm
	  X1=0.0001; // Left hand side initial value
	  X2=10.0;   // Right hand side initial value 
	  Y1=mixtureP(pr,X1,betas,Klambdas,iB, M,alpha);
	  Y2=mixtureP(pr,X2,betas,Klambdas,iB, M,alpha);
	  OK = (Y1 > 0) != (Y2 > 0);
	  if (OK) 
	    {
	      idh = 0;
	      while ((fabs (X1 - X2) > 1.0e-10) && (idh < 1000))
		{
		  //Rprintf("         Y1 %f Y2 %f X1 %f X2 %f ",Y1,Y2,X1,X2);
		  X = (X1 + X2) / 2.0;
		  Y = mixtureP(pr,X,betas,Klambdas,iB,M,alpha);
		  if ((Y1 < 0) && (Y < 0)) {
		    X1 = X;
		    Y1 = Y;
		  } else {
		    X2 = X;
		    Y2 = Y;
		  }
		  ++idh;
		}
	    }
	  REAL(quant)[iB]=X1;
	  INTEGER(NconvItr)[iB]=idh;
	}
    }
  

  PutRNGstate();
  UNPROTECT(1);
  return res;
}






double mixtureP_weib(double pr[],const double ti, const double betas[], const double Lambdas[],const int iB, const int M,const double alpha)
{
  // cdf of mixture 
  double temp = 0;
  int idh=0;
  for (int ih =0; ih < M ; ih++) 
    //rN[0]/L*dweibull(ts[it],Lambdas[M*iB+0],betas[M*iB+0],0))+...+rN[M-1]/L*dweibull(ts[it],Lambdas[M*iB+M-1],betas[M*iB+M-1],0))
    {
      idh = ih+M*iB;
      if (pr[ih]>0)
	temp += (double) pr[ih]*pweibull(ti,betas[idh],Lambdas[idh],
					 1, // Lowertail
					 0); // log
    }
  return temp-alpha;
}



// extern "C"{
//   SEXP getDens_weib(SEXP ts_, SEXP prob_, SEXP Lambdas_, SEXP betas_, SEXP M_, SEXP B_,SEXP alpha_,SEXP densi_);
// }

SEXP getDens_weib(SEXP ts_, SEXP weightS_, SEXP Lambdas_, SEXP betas_, SEXP M_, SEXP B_,SEXP alpha_,SEXP densi_)
{
  const double *ts= REAL(ts_);
  const int lt = length(ts_);
  const int M = INTEGER(M_)[0]; // truncation value 
  const int B = INTEGER(B_)[0]; // truncation value 
  const double *Lambdas = REAL(Lambdas_); // vector of length B*M
  const double *betas = REAL(betas_); // vector of length B*M
  const double *weightS = REAL(weightS_);
  const double alpha = REAL(alpha_)[0];
  const int densi = INTEGER(densi_)[0];
  //int rN[M], 
  int idh, OK;
  double pr[M],temp=0, X1=0.0001, X2=10.0, Y1,Y2,X,Y;

  SEXP res = PROTECT(allocVector(VECSXP,3)); // result is stored in res
  SEXP mat = allocVector(REALSXP, B*lt); 
  SET_VECTOR_ELT(res, 0, mat); 
  SEXP quant = allocVector(REALSXP, B); 
  SET_VECTOR_ELT(res, 1, quant); 
  SEXP NconvItr = allocVector(INTSXP, B); 
  SET_VECTOR_ELT(res, 2, NconvItr); 

  GetRNGstate();

  for (int iB=0; iB < B; iB++)
    {
      for (int ih=0;ih < M ; ih++) pr[ih]=weightS[ih+M*iB];
      // rN is a vector of length M containing the number of draws from each of the M components
      // sum(rN) = L;

      //rmultinom(L,pr,M,rN);

      if (densi||alpha==0) // approximate the density curve
	{
	  for (int it = 0 ; it < lt ; it++)
	    {
	      temp = 0;
	      for (int ih =0; ih < M ; ih++) 
		//rN[0]/L*dweibull(ts[it],Lambdas[M*iB+0],betas[M*iB+0],0))+...+rN[M-1]/L*dweibull(ts[it],Lambdas[M*iB+M-1],betas[M*iB+M-1],0))
		{
		  idh = ih+M*iB;
		  if (pr[ih]>0)
		    temp += (double) pr[ih]*dweibull(ts[it],betas[idh],Lambdas[idh],0);
		}
	      REAL(mat)[it+iB*lt]=temp; ///((double)L);
	    }
	}
      if (alpha!=0)
	{
	  // bisection algorithm
	  X1=0.0001; // Left hand side initial value
	  X2=10.0;   // Right hand side initial value 
	  Y1=mixtureP_weib(pr,X1,betas,Lambdas,iB, M,alpha);
	  Y2=mixtureP_weib(pr,X2,betas,Lambdas,iB, M,alpha);
	  OK = (Y1 > 0) != (Y2 > 0);
	  if (OK) 
	    {
	      idh = 0;
	      while ((fabs (X1 - X2) > 1.0e-10) && (idh < 1000))
		{
		  //Rprintf("         Y1 %f Y2 %f X1 %f X2 %f ",Y1,Y2,X1,X2);
		  X = (X1 + X2) / 2.0;
		  Y = mixtureP_weib(pr,X,betas,Lambdas,iB,M,alpha);
		  if ((Y1 < 0) && (Y < 0)) {
		    X1 = X;
		    Y1 = Y;
		  } else {
		    X2 = X;
		    Y2 = Y;
		  }
		  ++idh;
		}
	    }
	  REAL(quant)[iB]=X1;
	  INTEGER(NconvItr)[iB]=idh;
	}
    }
  

  PutRNGstate();
  UNPROTECT(1);
  return res;
}






// // OLD C++ codes
// //  ========== Old codes ======= =======
// extern "C"{
//   SEXP multin(SEXP L);     
// }
// SEXP multin(SEXP L)
// {
//   int M = 10;
//   double prob[M];
//   int rN[M];
//   SEXP res = PROTECT(allocVector(VECSXP,1)); // result is stored in res
//   SEXP post_S = allocVector(INTSXP, M); 
//   SET_VECTOR_ELT(res, 0, post_S); 
//   GetRNGstate();
  
//   for (int i = 0 ; i< M; i++) prob[i]=(double) i/45;
//   for (int i = 0 ; i< M; i++) Rprintf("  %f",prob[i]);
//   rmultinom(INTEGER(L)[0],prob,M,rN);
  
//   for (int i = 0 ; i< M; i++) INTEGER(post_S)[i] = rN[i];

//   PutRNGstate();
//   UNPROTECT(1);
//   return res;
// }

// // OLD C++ codes
// // ====== OLD codes for collapsed Gibbs sampling =======
// // 2
// void getQis(
// 	    vector<double> &qis,
// 	    const int iN, 
// 	    const int S[], 
// 	    const double ti, 
// 	    const vector<double> betas, 
// 	    const vector<double> Klambdas,
// 	    const int N
// 	    )
// {
//   // compute the prob of drawing from the existing set
//   // get qi's: vector of length n.uniq
//   // 
//   // qis:          A vector of size betas.size(), containing the probabilities that the i^th obs
//   //               to belong to the existing clusters.
//   // tempNV:       A NumericVector of size 1, containing temporary values
//   // S:            A IntegerVector of size N, containing class labels of observations.
//   // ti:           A NumericVector of size 1, containing the observed value of the i^th individual
//   //               It has to be NumericVector because it will be used as the first arg of dweibull
//   // betas:        A vector of size # uniq clusters from the previous Gibbs steps, containing the uniq values of betas
//   // Klambdas:     A vector of size # uniq clusters from the previous Gibbs steps, containing the uniq values of Klambdas
//   //               { betas[i],Klambda[i]} create a cluster

//   // resize the qis 
//   // for (tempInt=0;tempInt < 5; tempInt++) Rprintf("S: %d",S[tempInt]);

//   qis.resize(betas.size());
//   // Then Klambdas.size()==betas.size()==qis.size()
//   // 
//   int iS = 0, n = 0;
//   for (vector<double>::size_type ibeta = 0; ibeta < betas.size(); ibeta ++)
//     {
//       // if the i^th obs was creating a parameter set by itself
//       // then prob that the new drawn parameter set to be the
//       // same set is 0.
//       // This is because n = 0 if iN belongs to the uniq cluster S[iN] 
//       // due to (!)
	
//       // n counts the number of observations which belongs to ibeta.
//       // excluding the iN
//       n = 0;
//       // iS run through 1 to N (sample size) 
//       for (iS=0; iS < N; iS++) 
// 	{
// 	  if ((S[iS] == ibeta)&&(iS!=iN)) n++;// (!)
// 	}
//       // Rprintf("\v  the number of obs other than %d (=iN) which is in cluster %d is %d",iN, ibeta, n);
//       qis[ibeta] =  dWEI2(             
// 			  ti,              // A value at which the density is evaluated at:
// 			  betas[ibeta],    // shape 
// 			  Klambdas[ibeta], // Klambda
// 			  0                // log = FALSE
// 				       )*n;
//     }
//   // function returns qis: a vector of length # unique clusters from the previous iterations, containing the 
//   //                       probability that qis takes its value.
// }

// // 3
// int sample2(vector<double> &prob, // A vector of any length, the elements do not have to add up to 1              
// 	    const int tempInt2 // If == 0 then probability will be adjusted to be summed up to 10000
// 	    //,const bool uniqVanish
// 	    )
// {

//   /* generate a random number from 0 to prob.size()-1.
//      the probability of i^th obj to be selected is proportional to prob[i-1] */
//   if (tempInt2==0)
//     {
//       /* Step1: calculate the sum of all the weights */
//       double tempDbl = 0.0;
//       for( vector<double>::size_type i = 0; i < prob.size(); i++ )
// 	tempDbl += prob[i];
//       if (tempDbl == 0)
// 	{
// 	  Rprintf("bug! subroutine:sample, prob: ");
// 	  for( vector<double>::size_type i = 0; i < prob.size(); i++ )
// 	    Rprintf(" %f", prob[i]);
// 	  error("error");
// 	}
//       /* Step2: Each entry of prob is divided by the sum then multiplied by 10000
// 	 This will make the entries of prob to add up to 1000 and all greater or equal to 0. */
//       for( vector<double>::size_type i = 0; i < prob.size(); i++ )
// 	prob[i] *=10000/tempDbl;
//     }
 
 
//   /* Step3: generate single random number from 1 to 10,000 */
//   int intrandomN = rand()%10000 + 1;
//   double randomN = (double) intrandomN ;  // type casting
//   randomN = randomN - 0.00001; 
//   // there is a chance that randomN == 10000, then the equality randomN < prob[i]
//   // will never be hold. To avoid this issue, we make randomN = randomN - 0.00001;  
//   // Rprintf("\v random number: %d",randomN);
//   // must be double because of the subtraction operation below
  
//   // if (uniqVanish) 
//   //  {
//   //    Rprintf("sample, randomN %d prob: ",randomN);
//   //    for( vector<double>::size_type i = 0; i < prob.size(); i++ )
//   // 	Rprintf(" %f", prob[i]);
//   //  }
//   // to find bugs

//   /* Step4:go through the items one at a time, 
//      subtracting their weight from your random number, 
//      until you get the item where the random number is less than that item's weight*/
//   for( vector<double>::size_type i = 0; i < prob.size(); i++ ) // THE i MUST SHOULD HAVE TYPE int!
//     {
//       if(randomN < prob[i])
// 	{
// 	  //Rprintf("\v selected cluster %d",i);
// 	  return i;
// 	}
//       randomN -= prob[i]; 
//     }
//   Rprintf("\v error in sample subroutine \v probabilities:");
//   for ( vector<double>::size_type i = 0; i < prob.size();i ++ ) Rprintf(" %f",prob[i]);
//   Rprintf("\v selected random number %f",randomN);
//   error("error");
// }


// // 8
// double f(const double beta, // want to integrate over beta
// 	 const double ti, const double gam, const double dum)
// {
//   // take an integral of this function over 0 and phi
//   return beta*R_pow(ti,beta)/R_pow(gam+R_pow(ti,beta),dum+1);
// }






// //
// // Recursive auxiliary function for adaptiveSimpsons() function below
// // 9
// double adaptiveSimpsonsAux(double a, double b, double epsilon,                 
// 			   double S, double fa, double fb, double fc, int bottom,
// 			   const double ti, 
// 			   const double gam, 
// 			   const double dum
// 			   ,unsigned long &count_integ, 
// 			   unsigned long &count_integ_tot
// 			   )
// {          
//   count_integ_tot++;
//   double c = (a + b)/2, h = b - a;                                                                  
//   double d = (a + c)/2, e = (c + b)/2;                                                              
//   double fd = f(d,ti, gam, dum), fe = f(e,ti, gam, dum);                                                                   
//   double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
//   double Sright = (h/12)*(fc + 4*fe + fb);                                                          
//   double S2 = Sleft + Sright;    
//   if (bottom <= 0 )
//     {
//       count_integ++;
//       //Rprintf("!! The numerical approx. of q0 might be inaccurate! %d out of %d",count_integ,count_integ_tot);
//     }
//   if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)                                                    
//     return S2 + (S2 - S)/15;                                                                        
//   return adaptiveSimpsonsAux(a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1,ti, gam, dum
// 			     ,count_integ,count_integ_tot
// 			     ) +                    
//     adaptiveSimpsonsAux(c, b, epsilon/2, Sright, fc, fb, fe, bottom-1,ti, gam, dum
// 			,count_integ,count_integ_tot
// 			) ;                     
// }         
// // 10
// double adaptiveSimpsons( double a, 
// 			 double b,  // interval [a,b]
// 			 double epsilon,  // error tolerance
// 			 int maxRecursionDepth,
// 			 const double ti, 
// 			 const double gam, 
// 			 const double dum
// 			 ,unsigned long &count_integ,
// 			 unsigned long &count_integ_tot
// 			 ){   // recursion cap        
//   double c = (a + b)/2, h = b - a;                                                                  
//   double fa = f(a,ti, gam, dum), fb = f(b,ti, gam, dum), fc = f(c,ti, gam, dum);                                                           
//   double S = (h/6)*(fa + 4*fc + fb);                                                                
//   return adaptiveSimpsonsAux(a, b, epsilon, S, fa, fb, fc, maxRecursionDepth,ti, gam, dum
// 			     ,count_integ,count_integ_tot
// 			     );                   
// }                                                                                                   
 


// //11
// void getQi0(
// 	    double &qi0,
// 	    const int maxRecursion,
// 	    const double epsilon,
// 	    const double phi,
// 	    const double ti,
// 	    const double gam,
// 	    const double dum,
// 	    const double V
// 	    ,unsigned long &count_integ,
// 	    unsigned long &count_integ_tot
// 	    )
// {
 
//   // 300 is the number of intervals. It must be an even number 
//   // This number should vary based on other parameters to get a reliable integral estimate
//   // Taking integration over 0 to phi
//   qi0 = adaptiveSimpsons(  
// 			 0,            // double a, 
// 			 phi,          // double &b,  // interval [a,b]
// 			 epsilon,      // double &epsilon,  // error tolerance
// 			 maxRecursion, // const int &maxRecursionDepth,
// 			 // variables for adaptive sympson's rule
// 			 ti,        // double &ti, 
// 			 gam,          // double &gam, 
// 			 dum          // double &dum
// 			 ,count_integ, 
// 			 count_integ_tot
// 			   );
//   qi0 *= (V*dum*R_pow(gam,dum))/(phi*ti);
// }




// // 12
// // 5
// double max_VecScal2(const vector<double> &vec,
// 		    const double scal
// 		    )  
// {
 
//   double tempDbl1 = scal; // This may need to be initalised different base on your data
//   for (vector<double>::size_type i = 0; i < vec.size(); i++)
//     {
//       if (tempDbl1 < vec[i]) tempDbl1 = vec[i]; // vec.at(i) = vec[i] 
//     }
//   return tempDbl1;
// }



// double rWEI2 (const double beta, const double Klambda) 
// {

//   if (Klambda <= 0) 
//     {
//       error("subroutine:rWEI2, Klambda must be positive");
//     }
//   if (beta <= 0)
//     {
//       error("subroutine:rWEI2, beta must be positive");
//     }
 
//   return(rweibull(
// 		  beta, // shape
// 		  R_pow(Klambda,(1/beta)) //scale
// 		  )
// 	 );
    
// }




// double sum_inv2(const vector<double> &Klambdas)
// {
//   double tempDbl1 = 0;
//   for (vector<double>::size_type i=0; i < Klambdas.size();i++)
//     {
//       tempDbl1 += 1/Klambdas[i];
//     }
//   return tempDbl1;
// }





// // declarations (necessary for C++ codes)
// extern "C" {
//   SEXP gibbsPU(SEXP inpts, SEXP inpB, SEXP inpprp_sd_beta, SEXP inpaV, SEXP inpbV, SEXP inpaphi, SEXP inpbphi, SEXP inpaGam, SEXP inpbGam, SEXP inpdum, SEXP inpmaxRecursion, SEXP inpepsilon, SEXP inpburnin, SEXP inpL, SEXP inppred, SEXP inpprintFreq,SEXP alpha_) ;
// }



// SEXP gibbsPU( SEXP inpts, SEXP inpB, SEXP inpprp_sd_beta, SEXP inpaV, SEXP inpbV, SEXP inpaphi, SEXP inpbphi, SEXP inpaGam, SEXP inpbGam, SEXP inpdum, SEXP inpmaxRecursion, SEXP inpepsilon, SEXP inpburnin, SEXP inpL, SEXP inppred, SEXP inpprintFreq, SEXP alpha_){

//   /* ==== :Model description: =====
//      ti | beta_i Klambda_i ~ kW(lambdai,beta_i)
//      beta_i, Klambda_i | G ~ G
//      G ~ DP(V,D)
//      where D =  unif(beta; 0.0001,phi) * igamma(Klambda;shape=dum*,scale=gam) 
//      gam     ~  gamma(shape=aGam*,rate=bGam*)
//      phi     ~  pareto(shasssspe=aphi*,location=bphi*) 
//      V       ~  gamma(shape=aV*,scale=1/bV*)
     
//      * are hyperparameters 
//      Kottas recommended dum* = 2, aGam* = 1, aphi* = 2, aV = 2, and bV = 0.9 
     
//      The weibull pdf is parametrized as:
//      k_W(t|Klambda,beta) = beta*t^{beta-1}/lambda * exp(- t^{beta}/Klambda)
//      　　　
//      　　The relationship with the standard weibull Weibul(shape=beta,scale=lambda) is: lambda = Klambda^{1/beta}
//   */
  
//   // === Importing values === // 
//   const double *ts=REAL(inpts); // array of length N
//   const int B = INTEGER(inpB)[0]; 
//   const double prp_sd_beta = REAL(inpprp_sd_beta)[0]; // step 1
//   const double aV = REAL(inpaV)[0];                   // step 2
//   const double bV = REAL(inpbV)[0]; 
//   const double aphi = REAL(inpaphi)[0];               // step 3
//   const double bphi = REAL(inpbphi)[0];               
//   const double aGam = REAL(inpaGam)[0];               // step 4  
//   const double bGam = REAL(inpbGam)[0]; 
//   const double dum = REAL(inpdum)[0];   // default value should be 2 by Kottas
//   const int maxRecursion = INTEGER(inpmaxRecursion)[0];
//   const double epsilon = REAL(inpepsilon)[0];
//   const int burnin = INTEGER(inpburnin)[0];
//   const int L = INTEGER(inpL)[0];
//   const int pred = INTEGER(inppred)[0];
//   const int printFreq = INTEGER(inpprintFreq)[0];
//   const int N = length(inpts);
 
//   double phi=0, gam=0, V=0, beta=0, nalda = 0.0, qi0=0, candi=0,tempDbl1=0, ti = 0, Newb=0, NewK=0;
//   int iold = 0, Nuniq = 0, inew=0, ct_new=0,ct_accept=0, tempInt1 = 0, iN=0,iL=0,idh, S[N];
//   unsigned long ct_integ=0,ct_integ_tot=0;
//   bool uniqVanish = 0,bool_checknan = 0;
//   vector<double> Klambdas(1),betas(1),qis(1),probs(2),Nset(1); // change the size 
//   vector<int> N_Newbk;
  
//   // Create initial values of Klambdas, betas:
//   V = rgamma(
// 	     aV, // shape 
// 	     1/bV  // scale
// 	     );

//   gam = rgamma(aGam,1/bGam);

//   phi = rpareto(bphi, // location > 0
// 		aphi  // shape >  0
// 		); // returns double


//   for (iN = 0 ; iN < N; iN++ ) S[iN] = 0;

//   //Klambdas.resize (1);
//   //betas.resize (1);
//   Klambdas[0] =  rigamma(dum,gam); // return double
//   betas[0] = runif(
// 		   0.1, // min
// 		   phi // max
// 		   );
//   Nuniq = betas.size();

//   // Storage



//   // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
//   SEXP res = PROTECT(allocVector(VECSXP,10)); // result is stored in res
//   SEXP post_S = allocVector(INTSXP, B*N); 
//   SET_VECTOR_ELT(res, 0, post_S); 
//   SEXP post_betas = allocVector(REALSXP, B*N); 
//   SET_VECTOR_ELT(res, 1, post_betas); 
//   SEXP post_Klambdas = allocVector(REALSXP, B*N);
//   SET_VECTOR_ELT(res, 2, post_Klambdas); 
//   SEXP post_Nuniq = allocVector(INTSXP, B);
//   SET_VECTOR_ELT(res, 3, post_Nuniq); 
//   SEXP post_V = allocVector(REALSXP, B);
//   SET_VECTOR_ELT(res, 4, post_V); 
//   SEXP post_phi = allocVector(REALSXP, B);
//   SET_VECTOR_ELT(res, 5, post_phi); 
//   SEXP post_gam = allocVector(REALSXP, B);
//   SET_VECTOR_ELT(res, 6, post_gam); 
//   SEXP post_Eta=allocVector(REALSXP,(B-burnin));
//   SET_VECTOR_ELT(res, 7, post_Eta);
//   SEXP post_ct_integ_tot = allocVector(REALSXP,2); // ct_integ_tot:
//   SET_VECTOR_ELT(res,8,post_ct_integ_tot);
//   SEXP post_ct_accept = allocVector(REALSXP,2); // ct_integ_tot:
//   SET_VECTOR_ELT(res,9,post_ct_accept);

//   GetRNGstate();
  
//   // MCMC! 
//   for (int iB = 0 ; iB < B ; iB++)
//     {

//       // STEP 1: Update betas, Klambdas, S and Nuniq
//       for ( iN=0; iN < N ; iN++)
// 	{

// 	  if (Klambdas.size() != betas.size() || Klambdas.size() != Nuniq) 
// 	    {
// 	      Rprintf("bug! b.s(): %d",betas.size());
// 	      Rprintf("bug! Nuniq: %d",Nuniq);
// 	      error("error");
// 	    }

// 	  iold = S[iN];     // class belonging of the i^th obs from the previous step
// 	  ti = ts[iN];      // the i^th observation
// 	  /* Check if the i^th obs is creating a set by itself
// 	     if so, the unique parameter set corresponding the i^th obs
// 	     will disappear after this iteration */
// 	  uniqVanish = 1;
// 	  for (tempInt1=0; tempInt1 < N; tempInt1++ ) 
// 	    {
// 	      if (S[tempInt1]==iold && tempInt1!=iN)
// 		{
// 		  /* tempInt1^th patient also belongs to iold^th cluster
// 		     => iN^th patient is NOT creating a uniqe cluster   */
// 		  uniqVanish = 0; 
// 		  break; 
// 		  // if we find at least one observation which  belongs to the same cluster as iN, uniqVanish = 0 and we are done.
// 		}
// 	    }
	  
// 	  // prob( the i^th observation belongs to the existing cluster k  ); k = 1,..., Nuniq 
// 	  getQis(
// 		 qis,       // vector<double>  // qis changes the value!
// 		 iN,        // const int 
// 		 S,         // const IntegerVector
// 		 ti,        // const NumericVector &ti,
// 		 betas,     // const vector<double> 
// 		 Klambdas,   // const vector<double> 
// 		 N
// 		 );
// 	  // prob(  the i^th is from the new cluster  ) 
// 	  getQi0(
// 		 qi0, // double &qi0,
// 		 maxRecursion, // const int &maxRecursion,
// 		 epsilon,  // const double &epsilon,
// 		 phi,      // const double &phi,
// 		 ti,       // const NumericVector &ti,
// 		 gam,      // const double &gam,
// 		 dum,      // const double &dum,
// 		 V       //   const double &V,
// 		 ,ct_integ, // unsigned long &ct_integ,
// 		 ct_integ_tot //unsigned long &ct_integ_tot
// 		 );

// 	  /* qis = c(qis, q0) */
// 	  qis.push_back(qi0);
	  
	  
// 	  if (qis.size()!= Nuniq+1) Rprintf("bug!   qis.s(): %d", qis.size()); // bug trap
	  
// 	  if (uniqVanish == 1 && qis[iold] != 0 )  // bug trap
// 	    {
// 	      Rprintf(":::bug::: iN: %d uV: %d iold: %d qis[iold]: %f",iN,uniqVanish, iold, qis[iold]);
// 	      Rprintf(" qis:");
// 	      for (int iqis = 0; iqis< qis.size();iqis++) Rprintf(" %f",qis[iqis]);
// 	      Rprintf(" S:");
// 	      for (int iSS = 0; iSS< N ;iSS++) Rprintf(" %f",S[iSS]);
// 	      Rprintf(" Klambdas:");
// 	      for (int iSS = 0; iSS< Klambdas.size();iSS++) Rprintf(" %f",Klambdas[iSS]);
// 	      Rprintf(" betas:");
// 	      for (int iSS = 0; iSS< betas.size();iSS++) Rprintf(" %f",betas[iSS]);
// 	      error("error");
// 	    }
// 	  /* S[iN]: A integer between 0 and Nuniq. If it is between 0 and beta.size()-1 then 
// 	     The iN^th observation belongs to the existing cluster.
// 	     If it is beta.size(), then new cluster will be created.  */

// 	  S[iN]= sample2(qis,        // vector<double> &prob, // A vector of any length, the elements do not have to add up to 1             
// 			 0         // int &tempInt2,  // tempInt2 == 0 then probability will be adjusted to be summed up to 10000
// 			 );


// 	  if (S[iN]==Nuniq) 
// 	    {
// 	      //Rprintf(" IN! generate samples from the base distribution ");
// 	      ct_new++;
	    
// 	      /* The new parameter set (lambda,beta) is drawn
// 		 from the base distribution.
// 		 Input
// 		 STEP 1-1: generate a new lambda value from its full conditional: ivgamma 
// 		 lambda ~ igamma(shape=dum,scale=gam) => lambda|t1 ~ igamma(shape=dum+1,scale=gam+t1^beta) 
// 		 = 1/lambda ~ gamma(shape=dum+1,scale=1/(gam+t1^beta))  */
// 	      if (iold < 0 || iold >= betas.size() ){
// 		Rprintf("\v bug trap 9! iold=",iold);
// 		for (vector<double>::size_type i_st=0; i_st <betas.size(); i_st++)
// 		  Rprintf("\v betas[ii]: %f",betas[i_st]);
// 		error("error");
// 	      }

// 	      //Rprintf("igamma %f2.2", tempDbl1);
// 	      /* At this point Klambda.size() = Nuniq - 1. and Nuniq is defined from the previous step
// 		 The next line adds the new value to Klambda[Nuniq] */
// 	      Klambdas.push_back(rigamma(dum+1, // const double shape,
// 					 gam+R_pow(ti,betas[iold]) //  const double scale,
// 					 ) 
// 				 ); 
	     
	  
// 	      /* STEP 1-2: generate a beta via M-H
// 		 generate a canbeta from the normal proposal distribution */
// 	      beta = betas[iold];

// 	      candi = rnorm(beta,prp_sd_beta);
// 	      /* Left and Right truncation to avoid prior density = 0
// 		 and for Weibull parameters to be well defined (i.e. beta > 0) */

// 	      if (Nuniq != Klambdas.size()-1) // bug trap
// 		{
// 		  Rprintf("\v bug trap8!: Nuniq %d",Nuniq );
// 		  for (vector<double>::size_type i_st=0; i_st <Klambdas.size(); i_st++)
// 		    Rprintf("\v Klambdas[i]: %f",Klambdas[i_st]);
// 		  error("error");
// 		}

// 	      if (candi < phi && 0.0001 < candi)
// 		{
// 		  tempDbl1 = dWEI2(ti,candi,Klambdas[Nuniq],1)-dWEI2(ti,beta,Klambdas[Nuniq],1);
// 		  // const NumericVector &x,// const double &beta,
// 		  //const double &Klambda,//const bool &logTF 
			 
// 		  if (R_FINITE(tempDbl1)) // !tempLV1[0] = TRUE means that tempDbl1 is not NA
// 		    {
// 		      if (runif(0.0,1.0) < exp(tempDbl1))
// 			{
// 			  beta=candi;
// 			  ct_accept++;
// 			}
// 		    }else{
// 		    Rprintf("\v checknan is TRUE! ti %f candi %f Klambdas[Nuniq] %f beta %f",ti,candi,Klambdas[Nuniq],beta);
// 		  } // bug trap
// 		}
// 	      betas.push_back(beta);

// 	    } // if (S[iN]==Nuniq)

// 	  if(uniqVanish)
// 	    {
// 	      /* Since the number of unique sets of (lambda,beta) decreases,
// 		 We must trim the vanished (lambda,beta) from setlb         */
	      
// 	      //Rprintf("!uV! iold %d Kl.s %d",iold,Klambdas.size());

// 	      Klambdas.erase(Klambdas.begin()+iold); // Remove element at position iold
// 	      betas.erase(betas.begin()+iold);       // Remove element at position iold
// 	      for (tempInt1=0;tempInt1 < N;tempInt1++)
// 		{
// 		  if (S[tempInt1] > iold) S[tempInt1] -=1;
// 		  if (S[tempInt1] < 0) { 
// 		    Rprintf("bag trap10! tempInt1 %d and S[tempInt1] %d",tempInt1,S[tempInt1]);
// 		    error("error");}
// 		}
// 	    } 
	
// 	  //for (tempInt1=0; tempInt1 < Klambdas.size(); tempInt1++) 
// 	  //Rprintf(" new Klambdas %f ",Klambdas[tempInt1]);
// 	  Nuniq = betas.size();  
// 	} // Loop iN=0; iN < N ; iN++

//       if (Klambdas.size()!=betas.size()){Rprintf("  bug trap 1! ");error("error");}



//       // Step 2:  update V, the precision parameter
    
//       nalda = rbeta(V+1, N);
//       probs[0]=(aV+Nuniq-1)/(N*(bV-log(nalda))+aV+Nuniq-1);
//       probs[1] = 1-probs[0];
	
//       tempInt1 = sample2(probs,    // vector<double> &prob, // A vector of any length, the elements do not have to add up to 1
// 			 0         // int &tempInt1, tempInt2 == 0 then probability will be adjusted to be summed up to 10000        
// 			 );
//       // double &tempDbl)
//       if (tempInt1==0)
// 	{
// 	  V = rgamma(
// 		     aV+Nuniq, // shape
// 		     1/(bV-log(nalda)) // scale
// 		     ); 
// 	}else if(tempInt1==1){
// 	V = rgamma( 
// 		   aV+Nuniq-1,  // shape
// 		   1/(bV-log(nalda))// scale
// 		    ); 
//       }
	
//       /* Step 3: update phi
// 	 phi  ~ pareto(shape=aphi*,location=bphi*) 
// 	 => phi |  betas ~ pareto(shape=aphi+Nuniq,location=max(betas,bphi)) */
//       phi = rpareto(max_VecScal2(betas,bphi), //  const double location, // location > 0
// 		    Nuniq + aphi             //  const double shape // shape >  0
// 		    );

//       /* Step 4: update gamma:
// 	 gam     ~ gamma(shape=aGam*,rate=bGam*)
// 	 => gam| lambda ~ gamma(shape=dum+) */

//       gam = rgamma(
// 		   aGam+dum*Nuniq, // shape 
// 		   1/(bGam+sum_inv2(Klambdas))  // scale  // sum_inv(Klambda) = sum(1/Klambda)
// 		   );

//       if (Klambdas.size()!=betas.size()){Rprintf("  bug trap 2! ");error("error");}
      
	

//       // Store the result from this iteration
//       for (iN = 0 ; iN < N ; iN++ )
// 	{
// 	  idh = iN + iB * N;
// 	  INTEGER(post_S)[idh] = S[iN];
// 	  REAL(post_betas)[idh] = betas[S[iN]];
// 	  REAL(post_Klambdas)[idh] = Klambdas[S[iN]]; 
// 	}
//       INTEGER(post_Nuniq)[iB] = Nuniq;
//       REAL(post_V)[iB] = V;
//       REAL(post_phi)[iB] = phi;
//       REAL(post_gam)[iB] = gam;

//       if (iB >= burnin)
// 	{
// 	  // Nset is a vector<double> of size Nuniq
// 	  // The first 0 to Nuniq-1 entries will contain the number of observations in each clusters.
// 	  // and the Nuniq^th entry contains V. The rest will contain zero.
// 	  // Nuniq does not exceed N in this model
// 	  Nset.resize(Nuniq+1);
// 	  for (tempInt1 = 0; tempInt1 < Nuniq; tempInt1++) 
// 	    {
// 	      Nset[tempInt1] = 0.0;
// 	      for (iN = 0 ; iN < N ; iN++ )
// 		{
// 		  if (S[iN]==tempInt1) ++Nset[tempInt1];
// 		}
// 	    }
// 	  Nset[Nuniq] = V;
	  
// 	  tempInt1 = sample2(Nset,    // vector<double> &prob, 
// 			     // A vector of any length, the elements do not have to add up to 1      
// 			     0   // int &tempInt2, tempInt2 == 0 then probability will be adjusted to be summed up to 10000
// 			     );
// 	  if (tempInt1 < 0 ){Rprintf("  bug trap 3! ");error("error");}
// 	  else if (tempInt1 > Nuniq ){Rprintf("  bug trap 11! Nuniq %d tempInt1 %d ", Nuniq,tempInt1);error("error");}
	  
// 	  // But probability will NOT be again computed to adjust to be summed up to 10000 this time because it is already.
// 	  if (tempInt1==Nuniq)
// 	    {
// 	      Newb = runif(0.0001,phi);       // const double &beta shape
// 	      NewK = rigamma(dum,gam);     // const double &Klambda
// 	    }
// 	  else if(tempInt1 < Nuniq || tempInt1 >= 0 )
// 	    {
// 	      Newb = betas[tempInt1];       // const double &beta shape
// 	      NewK = Klambdas[tempInt1];     // const double &Klambda
// 	    }else{error("bug at predictive dist!");}  
// 	  REAL(post_Eta)[iB-burnin]=pow(-NewK*log(1-REAL(alpha_)[0]),1/Newb); 
// 	}

      
//       if(iB % printFreq == 0) Rprintf(" %d iterations are done...   ", iB);
      
//     }
//   // MCMC is done: Based on the predictive sample2s of (beta^{b},Klambda^{b}),b=1,...B, 
//   // we approximate the predictive distribution of the alpha^{th} quantile as:
//   // Pr(\tilde{eta}_{\alpha} <= t) = Pr( {- Klambda*log(1-alpha)}^{1/beta} <= t)
//   //                               = E( Pr( {- Klambda*log(1-alpha)}^{1/beta} <= t| Klambda, beta) )
//   //                         \approx sum_{b=1}^B Pr( {- Klambda^{b}*log(1-alpha)}^{1/beta^{b}} <= t| Klambda^{b}, beta^{b})
//   //
//   REAL(post_ct_integ_tot)[0] =ct_integ;
//   REAL(post_ct_integ_tot)[1] =ct_integ_tot;
//   REAL(post_ct_accept)[0] =ct_accept;
//   REAL(post_ct_accept)[1] =ct_new;

//   PutRNGstate();
//   UNPROTECT(1);
//   return res;

// }

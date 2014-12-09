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
  if (p_tot <0) 
    {
      error("At least one prob must be positive!");
    }
  for (ih = 0; ih < size; ih ++) prob[ih] /= p_tot;

  rmultinom(1, prob, size, rN);

  for ( ih=0; ih < size; ih++)
    {
      if(rN[ih]==1) return(ih);
    }
  return 1;
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



double max_MatScal(const int K, const int m, const double mat[K][m],const double scal)
{
  double tempDbl1 = scal;
  for (int ik = 0; ik < K; ik++)
    for (int im = 0; im < m; im++)
      if( tempDbl1 < mat[K][m]) tempDbl1 = mat[K][m];
  return tempDbl1;
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

 
  double phi=0, gam=0, D=0, candi=0, w=0,i1=0, i2=0,MHrate=0.0,dum=(max_dum+1.0)/2.0; 
  double can[2]={2.0,2.0}, ct_tot[2]={1.0,1.0},ct_acc[2]={1.0,1.0};
  int iN=0,ih=0,idh=0, S[N];
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

 
  double phi=0, gam=0, D=0,  candi=0, w=0,i1=0, i2=0,MHrate=0.0,dum=(max_dum+1.0)/2.0; 
  const int Nacc = 4;
  double can[Nacc], ct_tot[Nacc],ct_acc[Nacc];
  can[0]=2.0;can[1]=2.0;can[2]=15.0;can[3]=2.0;
  ct_tot[0]=1.0;ct_tot[1]=1.0;ct_tot[2]=1.0;ct_tot[3]=1.0;
  ct_acc[0]=1.0;ct_acc[1]=1.0;ct_acc[2]=1.0;ct_acc[3]=1.0;

  int iN=0,ih=0,idh=0, S[N];
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






SEXP gibbsAllUnif( SEXP inpts, SEXP inpB, SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inploc_phi, 
		   SEXP inpaGam, SEXP inploc_Gam, 
		   // SEXP inpmax_dum, 
		   SEXP inpburnin,  SEXP inpprintFreq){
  R_CheckUserInterrupt(); 
  R_ProcessEvents();
  R_FlushConsole(); 
  GetRNGstate();
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

 
  double phi=0, gam=0, D=0, candi=0, w=0,i1=0, i2=0,MHrate=0.0; 
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
  R_CheckUserInterrupt(); 
  R_ProcessEvents();
  R_FlushConsole(); 


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
      R_CheckUserInterrupt(); 
      R_ProcessEvents();
      R_FlushConsole(); 
      if(iB % printFreq == 0) 
	{
	  Rprintf("\v %d iterations are done...   ", iB);
	  //R_ProcessEvents();
	}
      
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




SEXP gibbsAllUnifAllSample( SEXP inpts, 
			    SEXP inpB, 
			    
			    SEXP M_, SEXP inpa_D, SEXP inpib_D, SEXP inpaphi, SEXP inploc_phi, 
			    SEXP inpaGam, SEXP inploc_Gam, 

			    SEXP inpburnin,  SEXP inpprintFreq,
			    SEXP K_, 
			    SEXP mks_,
			    SEXP max_mk_)
{
  // K = length(mks) the number of species
  // mks = the number of samples from each species
  // inpts = a vector of length sum(mks) containing all the species samples 
  R_CheckUserInterrupt(); 
  R_ProcessEvents();
  R_FlushConsole(); 
  GetRNGstate();
  /* ==== :Model description: =====
     tki | beta_ki lambda_ki ~ weibull(lambda_ki,beta_ki)
     beta_ki, lambda_ki | Gk ~ Gk
     Gk ~ DP(D,G0)
     where G0 =  unif(beta; 0.0001,phi) * unif(lambda;0.001,scale=gam) 
     phi     ~  pareto(shape=aphi*,location=loc_phi*) 
     gam     ~  pareto(shape=aGam*,location=loc_Gam*)
     D       ~  unif(a_D*,ib_D*)
  */
  
  // === Importing values === // 

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
  const int * mks = INTEGER(mks_);

  const int max_mk = INTEGER(max_mk_)[0];
  const int K = INTEGER(K_)[0];
  double ts[K][max_mk];
  int m_temp = 0,cumsum_m[K];
  for (int ik = 0 ;ik < K; ik ++)
    {
      cumsum_m[ik] = m_temp; // contains 0, mks[0], mks[0] + mks[1], mks[0] + mks[1] + mks[2],...
      for (int ipat = 0 ; ipat < mks[ik] ; ipat++)
	{
	  ts[ik][ipat] = REAL(inpts)[ipat + m_temp]; // array of length N
	}
      m_temp += mks[ik];
    }
  const double tot_m = m_temp;
  // m_sum contains the total number of samples
  
  double phi=0, gam=0, D=0, candi=0, w=0,i1=0, i2=0,MHrate=0.0; 
  const int Nacc = 3;
  // beta, Lambda, D 
  double can[Nacc], ct_tot[Nacc],ct_acc[Nacc];
  can[0]=0.5;can[1]=0.5;can[2]=0.5;
  ct_tot[0]=1.0;ct_tot[1]=1.0;ct_tot[2]=1.0;
  ct_acc[0]=1.0;ct_acc[1]=1.0;ct_acc[2]=1.0;
  int iN=0,ih=0,ik=0, idh=0, S[K][max_mk],reasonable;
  double Lambdas[K][M],betas[K][M],vs[K][M],postprob[K][M],weightS[K][M]; // change the size 
  const double smallValue = 0.1;  
  double maxCan[Nacc]; // beta, Lambda, D 
  maxCan[0]=10;maxCan[1]=10;maxCan[2]=1.2;
  double minCan[Nacc]; // beta, Lambda, D 
  minCan[0]=0.001;minCan[1]=0.001;minCan[2]=0.1;
  // beta is supposed to be from uniform (0,max_gum), however if beta<0.01, dwei(shape=beta,scale=Lambda) becomes undefined
  // for Lambda > 20,000. As our model should favor the set (beta, Lambda) = about (3,20,000), we must restrict any possible values of beta
  // to accept Lambda = 20,000. For this reason, we turncate the left hand extreme value of beta
  // Create initial values of Lambdas, betas:
  // ---------------------------------- //
  // Initializations of D, gam and phi  //
  // ---------------------------------- //
  D = runif(
	    a_D, // min 
	    ib_D  // max
	    );
  gam = rpareto(loc_Gam,
		aGam);

  phi = rpareto(loc_phi, // location > 0
		aphi  // shape >  0
		); // returns double
  // ------------------------------------------------- //
  // initializations of vs, weightS, Lambdas, betas, S //
  // ------------------------------------------------- //
  for (ik = 0 ; ik < K ; ik ++)
    {
      for (ih = 0 ; ih < (M-1) ; ih ++) vs[ik][ih] = rbeta(1.0,D);
      vs[ik][M-1]=1.0;
  
      weightS[ik][0] = vs[ik][0];
      w = 1.0 ;
      for (ih=0; ih < M; ih++)
	{
	  if (ih > 0){
	    w *= (1.0-vs[ik][ih-1]);
	    weightS[ik][ih] = vs[ik][ih]*w;
	  }
	  Lambdas[ik][ih] =  runif(smallValue,
				   gam); // return double
	  betas[ik][ih] = runif(
				smallValue, // min
				phi // max
				);

	  //Rprintf("\n ik %d ih %d beta %f Lambda %f",ik,ih,betas[ik][ih],Lambdas[ik][ih]);
	}
      for (iN = 0 ; iN < mks[ik]; iN++ ) S[ik][iN] = 0;
    }






  // Storage
  R_CheckUserInterrupt(); 
  R_ProcessEvents();
  R_FlushConsole(); 


  // http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
  SEXP res = PROTECT(allocVector(VECSXP,12)); // result is stored in res
  SEXP post_S = allocVector(INTSXP, B*tot_m); 
  SET_VECTOR_ELT(res, 0, post_S); 
  SEXP post_betas = allocVector(REALSXP, B*tot_m); 
  SET_VECTOR_ELT(res, 1, post_betas); 
  SEXP post_Lambdas = allocVector(REALSXP, B*tot_m);
  SET_VECTOR_ELT(res, 2, post_Lambdas);
 
  SEXP post_betas_uni = allocVector(REALSXP, B*M*K); 
  SET_VECTOR_ELT(res, 3, post_betas_uni); 
  SEXP post_Lambdas_uni = allocVector(REALSXP, B*M*K);
  SET_VECTOR_ELT(res, 4, post_Lambdas_uni); 
  SEXP post_vs = allocVector(REALSXP, B*M*K);
  SET_VECTOR_ELT(res, 5, post_vs); 
  SEXP post_weightS = allocVector(REALSXP, B*M*K);
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
  

  
  // Let's Roll! MCMC! 
  for (int iB = 0 ; iB < B ; iB++)
    {
      //  Update S, a vector of length N, containing the cluster index of random effect 1 of each patient.
      //  the cluster index is between 1 and M.
      //  conjugate: posterior distribution of S has a categorical distribution with M categories. 
      //  if (printing && iB%printEvery==0)  Rprintf("\v\v Step 4");
      for (ik = 0 ; ik < K ; ik++)
	{
	  for (iN = 0 ; iN < mks[ik] ; iN++)
	    { 
	      // Given patient iN, compute P(H_1i=ih;-) for each ih 
	      // P(H_1i=ih;-) = prod_{i=1}^{# repeated measure} dnbinom(yij; size=exp(b0+b1*x1i+..), prob=1/(gij+1))
	      for (ih=0; ih < M; ih++ ) postprob[ik][ih] = dweibull(ts[ik][iN],betas[ik][ih],Lambdas[ik][ih],0)*weightS[ik][ih];
	      S[ik][iN] = sample(postprob[ik],M); 
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
	      for (iN = 0 ; iN < mks[ik]; iN++ ) 
		{
		  if ( S[ik][iN] == ih ) i1++ ;
		  if ( S[ik][iN] > ih ) i2++ ;
		}
	      vs[ik][ih] = rbeta(
				 1.0 + i1, 
				 D + i2 
				 );
	      //Rprintf("ih %d sp1 %f sp2 %f",ih, (double) 1+i1+idh ,  D + (double)i2 + (double)ivec);
	      if (ih > 0) w *= (1.0-vs[ik][ih-1]);
	      weightS[ik][ih] = vs[ik][ih]*w;
	    }
	  weightS[ik][M-1] = w*(1-vs[ik][M-2]);
	  
	  // Updating Lambdas[ih] and betas[ih]
	  for (ih=0;ih < M ; ih ++)
	    {
	      // Update Lambdas M-H
	      candi = rnorm(Lambdas[ik][ih],can[1]);
	      if (weightS[ik][ih] > 0.1) ct_tot[1]++;
	      // the algorithm requires to compute Lambda hence Lambda and beta must not be too small
	      if (candi > smallValue && candi < gam)  
		{
		  MHrate = 0.0;  // If there is no obs to h^th component, candi is always accepted if it lies within [0.0001, phi]
		  for (iN = 0 ; iN < mks[ik] ; iN++)
		    {
		      if (S[ik][iN]==ih)
			MHrate += dweibull(ts[ik][iN],betas[ik][ih],candi,1) - dweibull(ts[ik][iN],betas[ik][ih],Lambdas[ik][ih],1);
		    }
		  
		  reasonable = R_FINITE(dweibull(6.0,betas[ik][ih],candi,0));
		  /* if(!reasonable)  
		     Rprintf("\n ik %d ih %d betas %f candi %f smallValue %f 
		     gam %f \n candi > smallValue %d candi > smallValue && candi < gam %d \n", 
		     ik,ih,betas[ik][ih],candi,smallValue, gam,candi > smallValue, 
		     candi > smallValue && candi < gam); */
		  
		  // If the proposed beta does not return finite value of density of Weibull, such proposal is not accepted
		  if (runif(0.0,1.0) < exp(MHrate) && reasonable )
		    {
		      if (weightS[ik][ih] > 0.1) ct_acc[1]++;
		      Lambdas[ik][ih] = candi;
		    }
		}
	      
	      // Update betas M-H
	      candi = rnorm(betas[ik][ih],can[0]); 
	      if (weightS[ik][ih] > 0.1) ct_tot[0]++;
	      // the algorithm requires to compute Lambda^beta hence Lambda and beta must not be too small
	      if (candi > smallValue && candi < phi)  
		{
		  MHrate = 0.0;  // If there is no obs to h^th component, candi is always accepted if it lies within [0.0001, phi]
		  for (iN = 0 ; iN < mks[ik] ; iN++)
		    {
		      if (S[ik][iN]==ih)
			MHrate += dweibull(ts[ik][iN],candi,Lambdas[ik][ih],1) -
			  dweibull(ts[ik][iN],betas[ik][ih],Lambdas[ik][ih],1);
		    }
		  // candi could be extreme values when no one belong to that class
		  // necessary to check that (candi,Lambdas[ik][ih]) returns finite value
		  reasonable = R_FINITE(dweibull(6.0,candi,Lambdas[ik][ih],0));
		  // if(!reasonable) Rprintf("\n ik%d ih%d candi %f Lambda %f",ik,ih,candi,Lambdas[ik][ih]);

		  // If the proposed beta does not return finite value of density of Weibull, such proposal is not accepted
		  if (runif(0.0,1.0) < exp(MHrate) && reasonable)
		    {
		      if (weightS[ik][ih] > 0.1) ct_acc[0]++;
		      betas[ik][ih] = candi;
		    }
		}
	    }
	}


      
      // Update phi 
      /* Step 3: update phi
	 phi  ~ pareto(shape=aphi*,location=loc_phi*) 
	 => phi |  betas ~ pareto(shape=aphi+Nuniq,location=max(betas,loc_phi)) */
      phi = rpareto(max_MatScal(K, M, betas,loc_phi), //  const double location, // location > 0
		    M*K + aphi             //  const double shape // shape >  0
		    );

      /* Update gamma:
	 gamma ~pareto(shape=aGam,location=loc_Gam) via M-H
      */
      gam = rpareto(max_MatScal(K, M, Lambdas,loc_Gam), //  const double location, // location > 0
		    M*K + aGam             //  const double shape // shape >  0
		    );

      // Update D
      // D ~ gamma(scale=a.D,shape=b.D) = gamma(scale=a.D,shape=1/ib.D)  
 
      MHrate = 0.0;
      candi = rnorm(D,can[2]);
      ct_tot[2]++;
      if (a_D <= candi && candi <= ib_D )
	{
	  for (ik = 0 ; ik < K; ik++)
	    for (ih = 0; ih < M-1 ; ih++)
	      MHrate += dbeta(vs[ik][ih], 1.0, candi,1) - dbeta(vs[ik][ih], 1.0, D,1);
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
      for (ik = 0 ; ik < K ; ik++ )
	{
	  for (iN = 0 ; iN < mks[ik] ; iN++ )
	    {
	      idh = iN + cumsum_m[ik] + iB * tot_m; // need to think here
	      INTEGER(post_S)[idh] = S[ik][iN];
	      REAL(post_betas)[idh] = betas[ik][S[ik][iN]];
	      REAL(post_Lambdas)[idh] = Lambdas[ik][S[ik][iN]]; 
	    }
	  for (ih = 0 ; ih < M; ih++ )
	    {
	      idh = ih + ik*M + iB * K * M;
	      REAL(post_betas_uni)[idh] = betas[ik][ih];
	      REAL(post_Lambdas_uni)[idh] = Lambdas[ik][ih]; 
	      REAL(post_weightS)[idh] = weightS[ik][ih];
	      REAL(post_vs)[idh] = vs[ik][ih]; 
	    }
	}

      REAL(post_D)[iB] = D;
      REAL(post_phi)[iB] = phi;
      REAL(post_gam)[iB] = gam;
      //REAL(post_dum)[iB] = dum;
      R_CheckUserInterrupt(); 
      R_ProcessEvents();
      R_FlushConsole(); 
      if(iB % printFreq == 0) 
	{
	  Rprintf("\v %d iterations are done...   ", iB);
	  //R_ProcessEvents();
	}
      
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
  GetRNGstate();
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

/* // sprodC2 compute prod_{j \ne k } Pr(eta_j > t | data) */
/* SEXP sprobC2(SEXP sortUs_, SEXP tsGrid_, SEXP CDF_, SEXP TsOthers_,  */
/* 	     SEXP ms_, SEXP par_, SEXP c0s_, SEXP alpha_,SEXP max_m_) */
/* { */
/*   GetRNGstate(); */

/*   const double *sortUs = REAL(sortUs_); */
/*   const double *tsGrid = REAL(tsGrid_); */
/*   const double *CDF = REAL(CDF_); // the values of CDF of this specie (corresponding to ts) */
/*   const double *TsOthers = REAL(TsOthers_); // the values of CDF of other species (corresponding to ts) */
/*   const int *ms = INTEGER(ms_); */
/*   const double alpha = REAL(alpha_)[0]; */
/*   const double *c0s = REAL(c0s_); */
/*   const int NMC = length(sortUs_); */
/*   const int Ngrid = length(tsGrid_); */
/*   const int K1 = length(ms_); */
/*   const int max_m = INTEGER(max_m_)[0]; */

/*   int iUpper, iLower=0,iMiddle; */
/*   double sampleQ, u, prod, zeta, temp = 0; */

/*   SEXP res = PROTECT(allocVector(VECSXP,1)); // result is stored in res */
/*   SEXP mat = allocVector(REALSXP, 1 );  */
/*   SET_VECTOR_ELT(res, 0, mat);  */


/*   int i2,ik; */
/*   double par[K1][2]; */
/*   Rprintf("\n max_n %d",max_m); */
/*   for (ik = 0; ik < K1 ; ik++){ */
/*     par[ik][0] = REAL(par_)[ik*2]; */
/*     par[ik][1] = REAL(par_)[ik*2+1]; */
/*     Rprintf("\n par[%d][0] %f",ik,par[ik][0]); */
/*     Rprintf("\n par[%d][1] %f",ik,par[ik][1]); */
/*     Rprintf("\n c0s[%d] %f",ik,c0s[ik]); */
/*   } */

/*   Rprintf("NMC %d Ngrid %d K1 %d",NMC, Ngrid, K1); */
/*   Rprintf("\n length(TsOthers_) %d ",length(TsOthers_)); */

/*   int empirical[K1]; */
/*   for (i2 = 0 ; i2 < K1 ; i2++) empirical[i2] = 0; */
/*   for (int i=0 ; i < NMC; i++) // compute the  */
/*     { */
/*       R_CheckUserInterrupt();  */
/*       R_ProcessEvents(); */
/*       R_FlushConsole();  */
/*       // STEP 1: generate sample of size N.MC from uniform dist'n */
/*       //Rprintf("\n\n i %d",i); */
/*       // t[w > p[i]] is the values of F; F(t) > p[i] the smallest of such is the p[i]^th quantile  */
/*       u = sortUs[i]; */
/*       iUpper = Ngrid-1; */
      
/*       // iLower is updated */
/*       while (iLower + 1 < iUpper) */
/*       	{ */
/*       	  iMiddle = (int) (iLower + iUpper)/2; */
/* 	  //Rprintf("\n i %d iLower %d iUpper %d u %f",i, iLower,iUpper, u); */
/* 	  //Rprintf("\n iMiddle %d CDF[ iMiddle ] %f",iMiddle,CDF[ iMiddle ]); */
/*        	  if (CDF[ iMiddle ] < u){ */
/*       	    iLower = iMiddle; */
/*       	  }else{ // CDF[ iMiddle ] >= u */
/*       	    iUpper = iMiddle; */
/*       	  } */
/*       	} */
/*       sampleQ = tsGrid[iUpper];// sampled quantile from Pr(eta_{k,alpha} <= t | data) */
/*       // Get H */
/*       prod = 1; */
/*       for (ik = 0; ik < K1; ik++) */
/* 	{ */
/* 	  i2 = empirical[ik]; */
/* 	  while(TsOthers[ik*max_m + i2] <= sampleQ) i2++; */
	    
/* 	  // The number of Ts from the species i2 such that Ts <= sampleQ minimum is zero */
/* 	  empirical[ik] = i2; */
	  
/* 	  zeta = empirical[ik] + c0s[ik]*pweibull(sampleQ,par[ik][0],par[ik][1],1,0); // shape scale zetaive_log */
/* 	  prod *= pbeta(alpha, zeta, c0s[ik] + ms[ik] - zeta, 1, 0); */
	  
/* 	  //Rprintf("\n alpha %f sampleQ %f prod %f zeta %f c+m-zeta %f empirical[%d] %d      ms[%d] %d    c0s %f", */
/* 	  //alpha, sampleQ,   prod,   zeta,   c0s[ik] + ms[ik] - zeta,  ik, empirical[ik],    ik, ms[ik] ,    c0s[ik]); */
/* 	} */
/*       temp = temp + prod; */
/*     } */
/*   REAL(mat)[0] = temp/NMC; */
/*   // returns sum_j prod_{k !=s } P( Fs^{-1}(alpha)^(j) <= Fk^{-1}(alpha) | data) / N.MC */
/*   PutRNGstate(); */
/*   UNPROTECT(1); */
/*   return res; */
/* } */




double nu(const double t,const double data[] ,const double c0, const double shape, const double scale, int *icountprev, const int n)
{
  double F0 = pweibull(t,shape,scale,1,0);
  int icount = *icountprev;

  while( data[icount] <= t & icount <= n)
    {
      //Rprintf("\n t %f data[%d] %f",t,icount,data[icount]);
      icount++;
    }
  *icountprev = icount;
  return (icount + F0*c0);
}

double petapost(const double t,const double shape,const double scale,const double alpha,const double data[], 
		const double c0,int *icountprev, const int n)
    {
      double shape1 = nu(t,data,c0,shape,scale,icountprev,n);
      double shape2 = c0 + n - shape1;
      if (shape1 < 1E-20) shape1 = 1E-20;
      if (shape2 < 1E-20) shape2 = 1E-20;
      return (1 - pbeta(alpha,shape1,shape2,1,0)); 
    }

SEXP cdf_sub(SEXP sortTs,SEXP shape, SEXP scale, SEXP alpha, SEXP dat,SEXP c0)
{
  
  GetRNGstate();

  const int n = length(dat);
  const int M = length(sortTs);
  SEXP res = PROTECT(allocVector(VECSXP,1)); // result is stored in res
  SEXP mat = allocVector(REALSXP, M); 
  SET_VECTOR_ELT(res, 0, mat); 
  
  int icountprev = 0;
  for (int it = 0 ; it < M; it++ ){
    R_CheckUserInterrupt(); 
    R_ProcessEvents();
    R_FlushConsole(); 
    REAL(mat)[it] = petapost(REAL(sortTs)[it] , 
			     REAL(shape)[0], REAL(scale)[0],
			     REAL(alpha)[0],REAL(dat), 
			     REAL(c0)[0],&icountprev, n);
    //Rprintf("mF %d ",icountprev);

  }
  PutRNGstate();
  UNPROTECT(1);
  return res;

}

/* SEXP sprobC_unsort(SEXP NMC_, SEXP ts_, SEXP CDF_, SEXP CDFs_) */
/* { */
/*   GetRNGstate(); */

/*   const long int NMC = INTEGER(NMC_)[0]; */
/*   const double *ts = REAL(ts_); */
/*   const double *CDF = REAL(CDF_); // the values of CDF of this specie (corresponding to ts) */
/*   const double *CDFs = REAL(CDFs_); // the values of CDF of other species (corresponding to ts) */

/*   const int Ngrid = length(ts_); */
/*   const int K1 = length(CDFs_)/length(CDF_); */

/*   double * sortUs = malloc(NMC * sizeof(double)); */

  
/*   for (int i = 0 ;i < NMC ; i++) */
/*     { */
/*       sortUs[i]=unif_rand(); */
/*     } */
/*   R_rsort(sortUs, NMC);  */
/*   //for (int i = 0 ;i < NMC; i++)  */
/*   //    Rprintf("\n %1.3f ",sortUs[i]);  */


/*   int iUpper, iLower=0,iMiddle; */
/*   double sampleQ, u, prod, temp = 0; */

/*   SEXP res = PROTECT(allocVector(VECSXP,1)); // result is stored in res */
/*   SEXP mat = allocVector(REALSXP, 1 );  */
/*   SET_VECTOR_ELT(res, 0, mat);  */

/*   Rprintf("NMC %d Ngrid %d K1 %d",NMC, Ngrid, K1); */

/*   for (int i=0 ; i < NMC; i++) // compute the */
/*     { */
/*       R_CheckUserInterrupt(); */
/*       R_ProcessEvents(); */
/*       R_FlushConsole(); */
/*       // STEP 1: generate sample of size N.MC from uniform dist'n */
/*       // Rprintf("i %d",i); */
/*       // t[w > p[i]] is the values of F; F(t) > p[i] the smallest of such is the p[i]^th quantile */
/*       u = sortUs[i]; */
/*       iUpper = Ngrid-1; */

/*       // iLower is updated */
/*       while (iLower + 1 < iUpper) */
/*       	{ */
/*       	  iMiddle = (int) (iLower + iUpper)/2; */
/*   	  //Rprintf("\n i %d iLower %d iUpper %d u %f",i, iLower,iUpper, u); */
/*   	  //Rprintf("\n iMiddle %d CDF[ iMiddle ] %f",iMiddle,CDF[ iMiddle ]); */
/*        	  if (CDF[ iMiddle ] < u){ */
/*       	    iLower = iMiddle; */
/*       	  }else{ // CDF[ iMiddle ] >= u */
/*       	    iUpper = iMiddle; */
/*       	  } */
/*       	} */
/*       sampleQ = ts[iUpper];// sampled quantile from Pr(eta_{k,alpha} <= t | data) */
/*       // Get H */
/*       prod = 1; */
/*       for (int ik=0; ik < K1; ik++){ */
/*       	prod *= (1 - CDFs[iUpper + Ngrid*ik]); */
/*       } */
/*       temp = temp + prod; */
/*     } */
/*   REAL(mat)[0] = (temp/NMC); */
/*   // returns sum_j prod_{k !=s } P( Fs^{-1}(alpha)^(j) <= Fk^{-1}(alpha) | data) / N.MC */
/*   PutRNGstate(); */
/*   UNPROTECT(1); */
/*   free(sortUs); /\* releases the memory for other applications *\/ */
/*   return res; */
/* } */


SEXP sprobC(SEXP sortUs_, SEXP ts_, SEXP CDF_, SEXP CDFs_)
{
  GetRNGstate();

  //const long int NMC = INTEGER(NMC_)[0];
  const double *ts = REAL(ts_);
  const double *CDF = REAL(CDF_); // the values of CDF of this specie (corresponding to ts)
  const double *CDFs = REAL(CDFs_); // the values of CDF of other species (corresponding to ts)

  const int Ngrid = length(ts_);
  const int K1 = length(CDFs_)/length(CDF_);

  const double *sortUs = REAL(sortUs_);
  const int NMC = length(sortUs_);
  /* double * sortUs = malloc(NMC * sizeof(double)); */  
  /* for (int i = 0 ;i < NMC ; i++) */
  /*   { */
  /*     sortUs[i]=unif_rand(); */
  /*   } */
  /* R_rsort(sortUs, NMC);  */
  //for (int i = 0 ;i < NMC; i++) 
  //    Rprintf("\n %1.3f ",sortUs[i]); 


  int iUpper, iLower=0,iMiddle;
  double sampleQ, u, prod, temp = 0;

  SEXP res = PROTECT(allocVector(VECSXP,1)); // result is stored in res
  SEXP mat = allocVector(REALSXP, 1 ); 
  SET_VECTOR_ELT(res, 0, mat); 
  
  Rprintf("NMC %d Ngrid %d K1 %d",NMC, Ngrid, K1);

  for (int i=0 ; i < NMC; i++) // compute the
    {
      R_CheckUserInterrupt();
      R_ProcessEvents();
      R_FlushConsole();
      // STEP 1: generate sample of size N.MC from uniform dist'n
      // Rprintf("i %d",i);
      // t[w > p[i]] is the values of F; F(t) > p[i] the smallest of such is the p[i]^th quantile
      u = sortUs[i];
      iUpper = Ngrid-1;

      // iLower is updated
      while (iLower + 1 < iUpper)
      	{
      	  iMiddle = (int) (iLower + iUpper)/2;
  	  //Rprintf("\n i %d iLower %d iUpper %d u %f",i, iLower,iUpper, u);
  	  //Rprintf("\n iMiddle %d CDF[ iMiddle ] %f",iMiddle,CDF[ iMiddle ]);
       	  if (CDF[ iMiddle ] < u){
      	    iLower = iMiddle;
      	  }else{ // CDF[ iMiddle ] >= u
      	    iUpper = iMiddle;
      	  }
      	}
      // iUpper is the smallest index such that CDF[index] >= u i.e. ts[index] is the smallest t such that CDF(t) >= u
      sampleQ = ts[iUpper];// sampled quantile from Pr(eta_{k,alpha} <= t | data)
      // Get H
      prod = 1;
      for (int ik=0; ik < K1; ik++){
      	prod *= (1 - CDFs[iUpper + Ngrid*ik]);
      }
      temp = temp + prod;
    }
  REAL(mat)[0] = (temp/NMC);
  // returns sum_j prod_{k !=s } P( Fs^{-1}(alpha)^(j) <= Fk^{-1}(alpha) | data) / N.MC
  PutRNGstate();
  UNPROTECT(1);
  //free(sortUs); /* releases the memory for other applications */
  return res;
}




SEXP sprobC2(SEXP sortUs_, SEXP ts_, SEXP CDF_, SEXP CDFs_, SEXP EXCLUDE_k)
{
  GetRNGstate();

  //const long int NMC = INTEGER(NMC_)[0];
  const double *ts = REAL(ts_);
  const double *CDF = REAL(CDF_); // the values of CDF of this specie (corresponding to ts)
  const double *CDFs = REAL(CDFs_); // the values of CDF of other species (corresponding to ts)

  const int Ngrid = length(ts_);
  const int K1 = length(CDFs_)/length(CDF_);

  const double *sortUs = REAL(sortUs_);
  const int NMC = length(sortUs_);
  /* double * sortUs = malloc(NMC * sizeof(double)); */  
  /* for (int i = 0 ;i < NMC ; i++) */
  /*   { */
  /*     sortUs[i]=unif_rand(); */
  /*   } */
  /* R_rsort(sortUs, NMC);  */
  //for (int i = 0 ;i < NMC; i++) 
  //    Rprintf("\n %1.3f ",sortUs[i]); 


  int iUpper, iLower=0,iMiddle;
  double sampleQ, u, prod1, prod2, temp1 = 0, temp2 = 0, temp3 = 0;

  SEXP res = PROTECT(allocVector(VECSXP,3)); // result is stored in res
  SEXP mat1 = allocVector(REALSXP, 1 ); 
  SET_VECTOR_ELT(res, 0, mat1); 
  SEXP mat2 = allocVector(REALSXP, 1 ); 
  SET_VECTOR_ELT(res, 1, mat2);
  SEXP eta = allocVector(REALSXP, 1);
  SET_VECTOR_ELT(res, 2, eta);
  //Rprintf("NMC %d Ngrid %d K1 %d",NMC, Ngrid, K1);

  for (int i=0 ; i < NMC; i++) // compute the
    {
      R_CheckUserInterrupt();
      R_ProcessEvents();
      R_FlushConsole();
      // STEP 1: generate sample of size N.MC from uniform dist'n
      // Rprintf("i %d",i);
      // t[w > p[i]] is the values of F; F(t) > p[i] the smallest of such is the p[i]^th quantile
      u = sortUs[i];
      iUpper = Ngrid-1;

      // iLower is updated
      while (iLower + 1 < iUpper)
      	{
      	  iMiddle = (int) (iLower + iUpper)/2;
  	  //Rprintf("\n i %d iLower %d iUpper %d u %f",i, iLower,iUpper, u);
  	  //Rprintf("\n iMiddle %d CDF[ iMiddle ] %f",iMiddle,CDF[ iMiddle ]);
       	  if (CDF[ iMiddle ] < u){
      	    iLower = iMiddle;
      	  }else{ // CDF[ iMiddle ] >= u
      	    iUpper = iMiddle;
      	  }
      	}
      // iUpper is the smallest index such that CDF[index] >= u i.e. ts[index] is the smallest t such that CDF(t) >= u
      sampleQ = ts[iUpper];// sampled quantile from Pr(eta_{k,alpha} <= t | data)
      // Get H
      prod1 = 1;
      prod2 = 1;
      for (int ik=0; ik < K1; ik++){
      	prod1 *= (1 - CDFs[iUpper + Ngrid*ik]);
	if (ik != INTEGER(EXCLUDE_k)[0]) prod2 *= (1 - CDFs[iUpper + Ngrid*ik]);
      }
      temp1 += prod1;
      temp2 += prod2;
      temp3 += sampleQ;
    }
  REAL(mat1)[0] = temp1/NMC;
  if (INTEGER(EXCLUDE_k)[0] < 0 ){
    REAL(mat2)[0] = R_NaReal;
  }else{
    REAL(mat2)[0] = temp2/NMC;
  }
  REAL(eta)[0] = temp3/NMC;
  // returns sum_j prod_{k !=s } P( Fs^{-1}(alpha)^(j) <= Fk^{-1}(alpha) | data) / N.MC
  PutRNGstate();
  UNPROTECT(1);
  //free(sortUs); /* releases the memory for other applications */
  return res;
}





SEXP getDens_weib(SEXP ts_, SEXP weightS_, SEXP Lambdas_, SEXP betas_, SEXP M_, SEXP B_,SEXP alpha_,SEXP densi_)
{
  GetRNGstate();
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



  for (int iB=0; iB < B; iB++)
    {
      R_CheckUserInterrupt(); 
      R_ProcessEvents();
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
		  idh = ih + M*iB;
		  if (pr[ih]>0)
		    temp += (double) pr[ih]*dweibull(ts[it],betas[idh],Lambdas[idh],0);
		  
		  if (!R_FINITE(temp)) Rprintf("\ndens %f ts %f betas %f Lambdas %f",
					       temp,ts[it],betas[idh],Lambdas[idh]);
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


functions {
  real pert_dist_lpdf(real theta, real lb, real md, real ub, real gama){
    real x1;
    real x2;
    real x3;
    real x4;
    real alpha; 
    real beta;
    
    alpha = 1 + gama * (md - lb)/(ub - lb);
    beta = 1 + gama * (ub - md)/(ub - lb);
  
    x1 = (alpha-1) * log((theta - lb));
    x2 = (beta-1) * log((ub - theta));
    x3 = (alpha+beta-1) * log((ub - lb));
    x4 = lbeta(alpha, beta);
    return( x1 + x2 - x3 - x4);
  }
}
data{
  int N;  // the total number of distinct dose group
  vector[N] n;  // the sample size for each dose group
  vector[N] x;  // the dose level of each dose group
  vector[N] y;  // the number of adverse events for each dose group
  real q;       // the BMR
  vector[4] priormuQ; 
  real priorlb[2]; //lower bound
  real priorub[2]; //upper bound
  real priorgama[2];
  real eps;
  cov_matrix[3] priorSigmaQ;
  real truncdQ;
  int<lower=0, upper=1> is_bin;       //model type 1 = Binomial 0 = otherwise
  int<lower=0, upper=1> is_betabin;  //model type 1 = Beta-Binomial 0 = otherwise
}
parameters{
  real<lower=0, upper=1> par1; //a 
  real<lower=0> par2; //BMD
  real par3; // d on a log scale
  real rho[is_betabin];
}
transformed parameters{
  real a;
  real b;
  real d;
  real k;
  real m[N];
  real abet[N];
  real bbet[N];
  real<lower=0> BMD;
  BMD = par2;
  a = par1;
  d = exp(par3);
  k = log(par2);
  b = (-log(1-q)) / (exp(k)+(((exp(k)*(exp(k)-1))/exp(d))));
  
  for(i in 1:N){
    if(x[i] == 0){
      m[i] = a;
     } else if(x[i] > 0) {
      m[i] = a + (1 - a)*(1 - exp(-b*x[i] - ( (b/exp(d)) * x[i] * (x[i] - 1) )));
     }
  }
    
    
  if(is_bin == 0) {
    for(i in 1:N){
      abet[i] = m[i]*((1/rho[is_betabin])-1);
      bbet[i] = (1.0 - m[i])*((1/rho[is_betabin])-1);
    }
  } else {
    for(i in 1:N){
      abet[i] = 0.0;
      bbet[i] = 0.0;
    }
  }
      
 
}
model{
    par1 ~ pert_dist(priorlb[1], priormuQ[1], priorub[1], priorgama[1]); //prior for a
    par2 ~ pert_dist(priorlb[2], priormuQ[2], priorub[2], priorgama[2]); //prior for BMD
    par3 ~ normal(priormuQ[3], priorSigmaQ[3,3])T[, truncdQ]; //prior for d

    //1/exp(d) ~ exponential(1);
    //1/exp(d) ~ uniform(0,1);
    
   if(is_bin==1) {
    
      for(i in 1:N){
        target += lchoose(n[i], y[i]) + y[i]*log(m[i]+eps) + (n[i] - y[i])*log(1 - m[i]+eps);
       // target += par3 - exp(par3);
      }
    
    } else {
      rho[is_betabin] ~ pert_dist(0.0, priormuQ[4], 1.0, 4.0);
      for(i in 1:N){
        target += lchoose(n[i], y[i]) + lgamma(abet[i]+y[i]+eps) + lgamma(bbet[i]+n[i]-y[i]+eps) - 
                  lgamma(abet[i]+bbet[i]+n[i]+eps) - lgamma(abet[i]+eps) - lgamma(bbet[i]+eps) + 
                  lgamma(abet[i]+bbet[i]+eps);
       // target += par3 - exp(par3);
      }
    }
}

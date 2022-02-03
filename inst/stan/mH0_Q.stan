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
  vector[N] y;  // the arithmetic mean of the response values for each dose group
  real priormu;
  real<lower=0> priorlb; //lower bound
  real<upper=1> priorub; //upper bound
  real priorgama;
  real eps;
  int<lower=0, upper=1> is_bin;  //model type 1 = Binomial 0 = otherwise
  int<lower=0, upper=1> is_betabin;  //model type 1 = Beta-Binomial 0 = otherwise
}
parameters{
  real<lower=0,upper=1> par; // par[1]=overall log(mean), par[2]=log(invsigma2)
  real etarho[is_betabin]; //will be defined if beta-binomial is to be fitted
}
transformed parameters{
  real a;
  real rho[is_betabin];
  real abet[N];
  real bbet[N];

  a=par;

  if(is_betabin == 1) {
    rho[is_betabin] = (exp(etarho[is_betabin])-1.0)/(exp(etarho[is_betabin])+1.0);
    for(i in 1:N){
      abet[i] = a*((1/rho[is_betabin])-1.0);
      bbet[i] = (1.0 - a)*((1.0/rho[is_betabin])-1);
    }
  } else if(is_bin == 1) {
    for(i in 1:N){
      abet[i] = 0.0;
      bbet[i] = 0.0;
    }
  }


}
model{
   par ~ pert_dist(priorlb, priormu, priorub, priorgama);

   if(is_bin==1) {

      for (i in 1:N){
        target += lchoose(n[i], y[i]) + y[i]*log(a+eps) + (n[i] - y[i])*log(1 - a+eps);
      }

    } else if(is_betabin==1){
      etarho ~ normal(0,1);
      for(i in 1:N){
        target += lchoose(n[i], y[i]) + lgamma(abet[i]+y[i]) + lgamma(bbet[i]+n[i]-y[i]) -
                  lgamma(abet[i]+bbet[i]+n[i]) - lgamma(abet[i]) - lgamma(bbet[i]) +
                  lgamma(abet[i]+bbet[i]);
      }
    }


}

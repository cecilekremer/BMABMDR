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
  vector[N] priormu;
  real priorlb; //lower bound
  vector[2] priorub; //upper bound
  real priorgama;
  real eps;
  int<lower=0, upper=1> is_bin;  //model type 1 = Binomial 0 = otherwise
  int<lower=0, upper=1> is_betabin;  //model type 1 = Beta-Binomial 0 = otherwise
}
parameters{
  vector[N] par; // par[1]=background, par[2:N]=increment per dose group
  real etarho[is_betabin]; //will be defined if beta-binomial is to be fitted
}
transformed parameters{
  vector[N] a;
  real rho[is_betabin];
  real abet[N];
  real bbet[N];


  a[1] = par[1];
  for(k in 2:N){
      a[k] = a[k-1] + par[k];
      if(a[k] <= 0){
        a[k] = 0.01;
      }
      if(a[k] >= 1){
        a[k] = 0.99;
      }
  }


  if(is_betabin == 1) {
    rho[is_betabin] = (exp(etarho[is_betabin])-1.0)/(exp(etarho[is_betabin])+1.0);
    for(i in 1:N){
      abet[i] = a[i]*((1/rho[is_betabin])-1.0);
      bbet[i] = (1.0 - a[i])*((1.0/rho[is_betabin])-1);
    }
  } else if(is_bin == 1) {
    for(i in 1:N){
      abet[i] = 0.0;
      bbet[i] = 0.0;
    }
  }
}
model{


    if(is_bin==1) {

      for (i in 1:N){
       par[1] ~ pert_dist(priorlb, priormu[1], priorub[1], priorgama);

         for(k in 2:N){

            par[k] ~ uniform(-priorub[2], priorub[2]);

            }

       target += lchoose(n[i], y[i]) + y[i]*log(a[i]+eps) + (n[i] - y[i])*log(1 - a[i]+eps);

      }

    } else if(is_betabin==1){
    etarho ~ normal(0,1);
      for (i in 1:N){
       par[1] ~ pert_dist(priorlb, priormu[1], priorub[1], priorgama);

         for(k in 2:N){

            par[k] ~ uniform(-priorub[2], priorub[2]);

            }
        target += lchoose(n[i], y[i]) + lgamma(abet[i]+y[i]) + lgamma(bbet[i]+n[i]-y[i]) -
                  lgamma(abet[i]+bbet[i]+n[i]) - lgamma(abet[i]) - lgamma(bbet[i]) +
                  lgamma(abet[i]+bbet[i]);
      }
    }


}

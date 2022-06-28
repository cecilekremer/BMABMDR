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
  vector[4] priormu;
  real priorlb[2]; //lower bound
  real priorub[2]; //upper bound
  real priorgama[2];
  real eps;
  cov_matrix[3] priorSigma;
  real truncd;
  int<lower=0, upper=1> is_bin;  //model type 1 = Binomial 0 = otherwise
  int<lower=0, upper=1> is_betabin;  //model type 1 = Beta-Binomial 0 = otherwise
}
parameters{
  //vector[3] par; // at=log(a), BMD on log-scale, ct=log(c-1), dt=log(d), invsigma2t=log(invsigma2)
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
  k = log(par2);
  d = exp(par3);
  b = exp(inv_Phi(q)-(k*d));

  for(i in 1:N){
    if(x[i] == 0){
      m[i] = a;
    } else if(x[i] > 0) {
      m[i] = a + (1 - a)*Phi(log(b) + d*log(x[i]));
    }
  }


  if(is_bin == 0) {
    for(i in 1:N){
      abet[i] = m[i]*((1.0/rho[1])-1.0);
      bbet[i] = (1.0 - m[i])*((1.0/rho[1])-1.0);
    }
  } else {
    for(i in 1:N){
      abet[i] = 0.0;
      bbet[i] = 0.0;
    }
  }

}
model{
    par1 ~ pert_dist(priorlb[1], priormu[1], priorub[1], priorgama[1]); //prior for a
    par2 ~ pert_dist(priorlb[2], priormu[2], priorub[2], priorgama[2]); //prior for BMD
    par3 ~ normal(priormu[3], priorSigma[3,3])T[,truncd]; //prior for d

  if(is_bin==1) {

      for(i in 1:N){
        target += lchoose(n[i], y[i]) + y[i]*log(m[i]+eps) + (n[i] - y[i])*log(1 - m[i]+eps);
      }

    } else {

      rho[is_betabin] ~ pert_dist(0.0, priormu[4], 1.0, 4.0);
      for(i in 1:N){
        target += lchoose(n[i], y[i]) + lgamma(abet[i]+y[i]+eps) + lgamma(bbet[i]+n[i]-y[i]+eps) -
                  lgamma(abet[i]+bbet[i]+n[i]+eps) - lgamma(abet[i]+eps) - lgamma(bbet[i]+eps) +
                  lgamma(abet[i]+bbet[i]+eps);
      }
    }
}

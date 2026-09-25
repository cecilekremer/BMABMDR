functions{
  vector algebra_system(vector yG,        // unknowns
  vector theta,    // parameters
        real q){  // data (integer)
  vector[1] x;
  if(yG[1]>0) x[1] = gamma_p(theta[2],yG[1]*theta[1]) - q;
  else if(yG[1]<=0) x[1]=1;
  return x;
  }

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
  real init_b;
  vector[4] priormu;
  array[2] real priorlb; // lower bound
  array[2] real priorub; // upper bound
  array[2] real priorgama;
  real eps;
  cov_matrix[3] priorSigma;
  real truncd;
  int<lower=0, upper=1> is_bin;       //model type 1 = Binomial 0 = otherwise
  int<lower=0, upper=1> is_betabin;  //model type 1 = Beta-Binomial 0 = otherwise
}
parameters{
  real<lower=0, upper=1> par1; //a
  real<lower=0> par2; //BMD
  real par3; // d on a log scale
  array[is_betabin] real rho; //will be defined if beta-binomial is to be fitted
}
transformed parameters{
  real a;
  real b;
  real d;
  real k;
  array[N] real m;
  array[N] real abet;
  array[N] real bbet;
  vector[2] theta;
  vector[1] y_guess;
  vector[1] yG;
  real<lower=0> BMD;
  BMD = par2;
  a = par1;
  d = exp(par3);
  k = log(par2);
  theta[1] = BMD;
  theta[2] = d;
  y_guess[1] = init_b;

  yG = solve_powell_tol(algebra_system, y_guess, 1e-10, positive_infinity(), 1000, theta, q);
  b = yG[1];

  for(i in 1:N){
    if(x[i] == 0){
      m[i] = a;
    } else if(x[i] > 0) {
      m[i] = a + (1 - a)*gamma_cdf(x[i] | d, b);
    }
  }


  if(is_bin == 0) {

    for(i in 1:N){
      abet[i] = m[i]*((1.0/rho[is_betabin])-1.0);
      bbet[i] = (1.0 - m[i])*((1/rho[is_betabin])-1.0);
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


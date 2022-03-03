functions{
   vector algebra_system(vector yG,        // unknowns
               vector theta,    // parameters
               real[] x_r,      // data (real)
               int[] x_i) {     // data (integer)
   vector[1] x;
   real q = x_r[1];
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
  real priorlb[2]; //lower bound
  real priorub[2]; //upper bound
  real priorgama[2];
  real eps;
  cov_matrix[3] priorSigma;
  int<lower=0, upper=1> is_bin;       //model type 1 = Binomial 0 = otherwise
  int<lower=0, upper=1> is_betabin;  //model type 1 = Beta-Binomial 0 = otherwise
 }
 transformed data{
   real x_r[1] = {q};
   int x_i[0];
 }
 parameters{
  real<lower=0, upper=1> par1; //a
  real<lower=0> par2; //BMD
  real par3; // d on a log scale
  real rho[is_betabin];
}
 transformed parameters{
   real a;
   real d;
   real b;
   real k;
   vector[2] theta;
   vector[1] y_guess;
   vector[1] yG;
   real m[N];
   real abet[N];
   real bbet[N];
   real<lower=0> BMD;
   BMD = par2;
   a = par1;
   d = exp(par3);
   k = log(par2);
   theta[1] = BMD;
   theta[2] = d;
   y_guess[1] = init_b;

   yG = algebra_solver(algebra_system, y_guess, theta, x_r, x_i, 1e-10, positive_infinity(), 1e3);
   b = yG[1];

   for(i in 1:N){
    if(x[i] == 0){
      m[i] = a;
    } else if(x[i] > 0) {
      m[i] = a + (1 - a)*gamma_cdf(x[i], d, b);
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
    par3 ~ normal(priormu[3], priorSigma[3,3]); //prior for d

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

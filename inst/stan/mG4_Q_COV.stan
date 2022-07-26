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
  int nlevels;
  int nlevels_a;
  int nlevels_BMD;
  int nlevels_d;
  int nlevels_b;
  vector[N] n;  // the sample size for each dose group
  vector[N] x;  // the dose level of each dose group
  vector[N] y;  // the number of adverse events for each dose group
  matrix[N, nlevels] trt_ind;
  real q;       // the BMR
  real init_b;
  matrix[4, nlevels] priormu;
  matrix[4, nlevels] priorlb; //lower bound
  matrix[4, nlevels] priorub; //upper bound
  matrix[2, nlevels] priorgama;
  real eps;
  cov_matrix[3] priorSigma;
  real truncd;
 }
transformed data{
   real x_r[1] = {q};
   int x_i[0];
 }
parameters{
  real<lower=0, upper=1> par1[nlevels_a]; //a
  real<lower=0> par2[nlevels_BMD]; //BMD
  real par3[nlevels_d]; // d on a log scale
}
 transformed parameters{
  real a[nlevels_a];
  real b[nlevels_b];
  real d[nlevels_d];
  real k[nlevels_BMD];
//  real m[N];
  real<lower=0> BMD[nlevels_BMD];
  vector[2] theta;
  vector[1] y_guess;
  vector[1] yG;

  for(mn in 1:nlevels_BMD){
    BMD[mn] = par2[mn];
    k[mn] = log(par2[mn]);
  }

  for(mn in 1:nlevels_a){
    a[mn] = par1[mn];
  }

  for(mn in 1:nlevels_d){
    d[mn] = exp(par3[mn]);
  }

  if(nlevels_a == 1 && nlevels_d == 1 && nlevels_BMD == 1){
    for(mn in 1:nlevels_b){
      theta[1] = BMD[1];
      theta[2] = d[1];
      y_guess[1] = init_b;
      yG = algebra_solver(algebra_system, y_guess, theta, x_r, x_i, 1e-10, positive_infinity(), 1e3);
      b[mn] = yG[1];
    }
  }else if(nlevels_a > 1 && nlevels_d == 1 && nlevels_BMD == 1){
    for(mn in 1:nlevels_b){
      theta[1] = BMD[1];
      theta[2] = d[1];
      y_guess[1] = init_b;
      yG = algebra_solver(algebra_system, y_guess, theta, x_r, x_i, 1e-10, positive_infinity(), 1e3);
      b[mn] = yG[1];
    }
  }else if(nlevels_a == 1 && nlevels_d > 1 && nlevels_BMD > 1){
    for(mn in 1:nlevels_b){
      theta[1] = BMD[mn];
      theta[2] = d[mn];
      y_guess[1] = init_b;
      yG = algebra_solver(algebra_system, y_guess, theta, x_r, x_i, 1e-10, positive_infinity(), 1e3);
      b[mn] = yG[1];
    }
  }else if(nlevels_a > 1 && nlevels_d > 1 && nlevels_BMD > 1){
    for(mn in 1:nlevels_b){
      theta[1] = BMD[mn];
      theta[2] = d[mn];
      y_guess[1] = init_b;
      yG = algebra_solver(algebra_system, y_guess, theta, x_r, x_i, 1e-10, positive_infinity(), 1e3);
      b[mn] = yG[1];
    }
  }

}
 model{

  // a prior
  for(i in 1:nlevels_a){
    par1[i] ~ pert_dist(priorlb[1,i], priormu[1, i], priorub[1, i], priorgama[1, i]);
  }

  // BMD prior
  for(i in 1:nlevels_BMD){
    par2[i] ~ pert_dist(priorlb[2, i], priormu[2, i], priorub[2, i], priorgama[2, i]);
  }

  // d prior
  for(i in 1:nlevels_d){
    par3[i] ~ normal(priormu[3, i], priorSigma[3, 3])T[,truncd];
  }

  if(nlevels_d > 1 && nlevels_BMD > 1 && nlevels_a == 1){
    for(i in 1:N){
      for(mn in 1:nlevels){
        if(x[i] == 0) { target += ( lchoose(n[i], y[i]) +
                  y[i] * log(a[1]+eps) +
                  (n[i] - y[i]) * log(1 - a[1]+eps)) * trt_ind[i, mn];
        }else if(x[i] > 0) { target += ( lchoose(n[i], y[i]) +
                                         y[i] * log((a[1] + (1-a[1])*gamma_cdf(x[i], d[mn], b[mn])) +eps) +
                                        (n[i] - y[i]) * log(1 - (a[1] + (1-a[1])*gamma_cdf(x[i], d[mn], b[mn]))+eps)) * trt_ind[i, mn];
        }
      }
    }
  }else if(nlevels_a > 1 && nlevels_d == 1 && nlevels_BMD == 1){
    for(i in 1:N){
      for(mn in 1:nlevels){
        if(x[i] == 0) { target += ( lchoose(n[i], y[i]) +
                  y[i] * log(a[mn]+eps) +
                  (n[i] - y[i]) * log(1 - a[mn]+eps)) * trt_ind[i, mn];
        }else if(x[i] > 0) { target += ( lchoose(n[i], y[i]) +
                                         y[i] * log((a[mn] + (1-a[mn])*gamma_cdf(x[i], d[1], b[1])) +eps) +
                                        (n[i] - y[i]) * log(1 - (a[mn] + (1-a[mn])*gamma_cdf(x[i], d[1], b[1]))+eps)) * trt_ind[i, mn];
        }
      }
    }
  }else if(nlevels_a > 1 && nlevels_d > 1 && nlevels_BMD > 1){
    for(i in 1:N){
      for(mn in 1:nlevels){
        if(x[i] == 0) { target += ( lchoose(n[i], y[i]) +
                  y[i] * log(a[mn]+eps) +
                  (n[i] - y[i]) * log(1 - a[mn]+eps)) * trt_ind[i, mn];
        }else if(x[i] > 0) { target += ( lchoose(n[i], y[i]) +
                                         y[i] * log((a[mn] + (1-a[mn])*gamma_cdf(x[i], d[mn], b[mn])) +eps) +
                                        (n[i] - y[i]) * log(1 - (a[mn] + (1-a[mn])*gamma_cdf(x[i], d[mn], b[mn]))+eps)) * trt_ind[i, mn];
        }
      }
    }
  }else if(nlevels_a == 1 && nlevels_d == 1 && nlevels_BMD == 1){
    for(i in 1:N){
      for(mn in 1:nlevels){
        if(x[i] == 0) { target += ( lchoose(n[i], y[i]) +
                  y[i] * log(a[1]+eps) +
                  (n[i] - y[i]) * log(1 - a[1]+eps)) * trt_ind[i, mn];
        }else if(x[i] > 0) { target += ( lchoose(n[i], y[i]) +
                                         y[i] * log((a[1] + (1-a[1])*gamma_cdf(x[i], d[1], b[1])) +eps) +
                                        (n[i] - y[i]) * log(1 - (a[1] + (1-a[1])*gamma_cdf(x[i], d[1], b[1]))+eps)) * trt_ind[i, mn];
        }
      }
    }
  }

 }

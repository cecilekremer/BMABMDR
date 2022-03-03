functions{
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
  vector[N] m;  // the arithmetic mean of the response values for each dose group
  vector[N] s2;  // the arithmetic variance of the response values for each dose group
  real shift; // data are shifted if different from 0
  real priorlb;
  real priorg;
  vector[2] priorub;
  vector[N+1] priormu;
  cov_matrix[N+1] priorSigma;
  int data_type; // data_type; 1 = increasing N, 2 = increasing LN, 3 = decreasing N, 4 = decreasing LN
}
parameters{
  vector[N+1] par; // par[1]=background, par[2:N]=increment per dose group, par[N+1]=log(invsigma2)
}
transformed parameters{
  vector[N] a;
  real invsigma2;
  vector[N] mu;

  // means per dose group
  mu[1] = par[1];
  for(k in 2:N){
    mu[k] = mu[k-1] + par[k];
  }

  // parameter a
  if(data_type == 1 || data_type == 3){ // normal dist
    a = mu;
  }else if(data_type == 2 || data_type == 4){ // lognormal dist
    a = log(mu) - shift;
  }

  invsigma2=exp(par[N+1]);
}
model{
   par[1] ~ pert_dist(priorlb, priormu[1], priorub[1], priorg);

   for(k in 2:N){

    par[k] ~ uniform(-priorub[2], priorub[2]);

   }

   par[N+1] ~ normal(priormu[N+1],priorSigma[(N+1),(N+1)]);

  if(data_type == 1 || data_type == 3){ // normal dist
   for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a[i])*invsigma2;
   }
  }else if(data_type == 2 || data_type == 4){ // lognormal dist
   for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a[i])*invsigma2 - m[i]*n[i];
   }
  }

}


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
  vector[2] priormu;
  cov_matrix[2] priorSigma;
  real priorlb; //lower bound
  real priorub; //upper bound
  real priorg;
  int data_type; // data_type; 1 = increasing N, 2 = increasing LN, 3 = decreasing N, 4 = decreasing LN
}
parameters{
  vector[2] par; // par[1]=overall mean, par2=log(invsigma2)
}
transformed parameters{
  real a;
  real invsigma2;

  if(data_type == 1 || data_type == 3){
    a = par[1];
  }else if(data_type == 2 || data_type == 4){
    a = log(par[1]) - shift;
  }

  invsigma2=exp(par[2]);
}
model{
   par[1] ~ pert_dist(priorlb, priormu[1], priorub, priorg);
   par[2] ~ normal(priormu[2],priorSigma[2,2]);
   if(data_type == 1 || data_type == 3){
    for (i in 1:N){
       target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a)*invsigma2;
    }
   }else if(data_type == 2 || data_type == 4){
    for (i in 1:N){
      target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a)*invsigma2 - m[i]*n[i];
    }
   }

}

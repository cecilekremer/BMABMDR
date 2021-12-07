data{
  int N;  // the total number of distinct dose group
  vector[N] n;  // the sample size for each dose group
  vector[N] m;  // the arithmetic mean of the response values for each dose group
  vector[N] s2;  // the arithmetic variance of the response values for each dose group
  vector[2] priormu;
  cov_matrix[2] priorSigma;
}
parameters{
  vector[2] par; // par[1]=overall mean, par[2]=log(invsigma2)
}
transformed parameters{
  real a;
  real invsigma2;
  a=par[1];
  invsigma2=exp(par[2]);
}
model{
   par ~ multi_normal(priormu,priorSigma);
   for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a)*invsigma2;
   }
}

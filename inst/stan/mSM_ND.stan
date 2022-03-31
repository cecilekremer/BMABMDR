data{
  int N;  // the total number of distinct dose group
  vector[N] n;  // the sample size for each dose group
  vector[N] m;  // the arithmetic mean of the response values for each dose group
  vector[N] s2;  // the arithmetic variance of the response values for each dose group
  vector[N+1] priormu;
  cov_matrix[N+1] priorSigma;
}
parameters{
  vector[N+1] par; // par[1..N]=local log means, par[N+1]=log(invsigma2)
}
transformed parameters{
  vector[N] a;
  real invsigma2;
  a=exp(par[1:N]);
  invsigma2=exp(par[N+1]);
}
model{
   par ~ multi_normal(priormu, priorSigma);
   for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a[i])*invsigma2;
   }
}

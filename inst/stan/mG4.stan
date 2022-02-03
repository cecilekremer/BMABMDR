 functions{
   vector algebra_system(vector y,        // unknowns
               vector theta,    // parameters
               real[] x_r,      // data (real)
               int[] x_i) {     // data (integer)
   vector[1] x;
   real q = x_r[1];
   real data_type = x_i[1];
   if(data_type == 1){
       if(y[1]>0) x[1] = gamma_p(theta[3],y[1]*theta[1]) - (q/(theta[2]-1));
       else if(y[1]<=0) x[1]=1;
   }else if(data_type == 2){
       if(y[1]>0) x[1] = gamma_p(theta[3],y[1]*theta[1]) - (log(1+q)/(theta[4]*(theta[2]-1)));
       else if(y[1]<=0) x[1]=1;
   }else if(data_type == 3){
       if(y[1]>0) x[1] = gamma_p(theta[3],y[1]*theta[1]) - ((-q)/(theta[2]-1));
       else if(y[1]<=0) x[1]=1;
   }else if(data_type == 4){
      if(y[1]>0) x[1] = gamma_p(theta[3],y[1]*theta[1]) - (log(1-q)/(theta[4]*(theta[2]-1)));
      else if(y[1]<=0) x[1]=1;
   }
   return x;
   }

  real pert_dist_lpdf(real theta, real alpha, real beta, real lb, real ub){
    real x1;
    real x2;
    real x3;
    real x4;

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
   vector[N] m;  // the arithmetic mean of the response values for each dose group
   vector[N] s2;  // the arithmetic variance of the response values for each dose group
   real q;       // the BMR
   real init_b; // initial value for b
   vector[5] priormu;
   vector[5] priorlb;
   vector[5] priorub;
   vector[5] shape1;
   vector[5] shape2;
   cov_matrix[5] priorSigma;
  int data_type; // data data_type; 1 = increasing N, 2 = increasing LN, 3 = decreasing N, 4 = decreasing LN
  int<lower=0, upper=1> is_increasing; // indicator for increasing data
  real L; //lower bound for increasing data
  int<lower=0, upper=1> is_decreasing; // indicator for decreasing data
  real U; // upper bound for decreasing data
 }
 transformed data{
   real x_r[1] = {q};
   int x_i[1] = {data_type};
 }
 parameters{
 real<lower=0> par1;
 real<lower=0, upper=1> par2;
 real<lower=0> pars3i[is_increasing]; // will be size one if is_increasing
 real<lower=0, upper=1> pars3d[is_decreasing]; // will be size one if is_decreasing
 real par4;
 real par5;
}
 transformed parameters{
  real b;
  real a;
  real c;
  real par3;
  real d;
  real k;
  real mu_inf;
  real invsigma2;
  real mu_0;
  vector[4] theta;
  vector[1] y_guess;
  vector[1] y;

  mu_0 = par1;

  if(data_type == 1 || data_type == 3){
    a = par1;
  }else if(data_type == 2 || data_type == 4){
    a = log(par1);
  }

  if(is_increasing){
    par3 = L + pars3i[1];
  }else if(is_decreasing){
    par3 = L + (U - L) .* pars3d[1];
  }

  mu_inf = par1*par3;
  if(data_type == 1 || data_type == 3){
    c = mu_inf/mu_0;
  }else if(data_type == 2 || data_type == 4){
    c = log(mu_inf)/log(mu_0);
  }
  d = exp(par4);
  k = log(par2);
  invsigma2=exp(par5);

   theta[1] = exp(k);
   theta[2] = c;
   theta[3] = d;
   theta[4] = a;
   y_guess[1] = init_b;

   y = algebra_solver(algebra_system, y_guess, theta, x_r, x_i, 1e-10, positive_infinity(), 1e3);
   b = y[1];

 }
 model{
    par1 ~ pert_dist(shape1[1], shape2[1], priorlb[1], priorub[1]);
    par2 ~ pert_dist(shape1[2], shape2[2], priorlb[2], priorub[2]);
    par3 ~ pert_dist(shape1[3], shape2[3], priorlb[3], priorub[3]);
    par4 ~ normal(priormu[4],priorSigma[4,4]);
    par5 ~ normal(priormu[5],priorSigma[5,5]);

  if(data_type == 1 || data_type == 3){
    for (i in 1:N){
      target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a-a*(c-1)*
      gamma_cdf(x[i],d,b))*invsigma2;
     }
  }else if(data_type == 2 || data_type == 4){
      for (i in 1:N){
      target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*square(m[i]-a-a*(c-1)*
      gamma_cdf(x[i],d,b))*invsigma2 - m[i]*n[i];
     }
  }

 }

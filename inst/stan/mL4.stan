functions {
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
  real shift; // data are shifted if different from 0
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

  mu_0 = par1;

  if(is_increasing == 1){
    par3 = L + pars3i[1];
  }else if(is_decreasing == 1){
    par3 = L + (U - L) .* pars3d[1];
  }

  mu_inf = par1*par3;

  if(data_type == 1 || data_type == 2){
    if(data_type == 1){
      a = mu_inf;
    }else if(data_type == 2){
      a = log(mu_inf) - shift;
    }
  }else if(data_type == 3 || data_type == 4){
    if(data_type == 3){
      a = mu_0;
    }else if(data_type == 4){
      a = log(mu_0) - shift;
    }
  }

  if(data_type == 1 || data_type == 2){
    if(data_type == 1){
      c = pow(3,0.5)/pi()*logit(mu_0/mu_inf);
    }else if(data_type == 2){
      c = pow(3,0.5)/pi()*logit((log(mu_0)-shift)/(log(mu_inf)-shift));
    }
  }else if(data_type == 3 || data_type == 4){
    if(data_type == 3){
      c = pow(3,0.5)/pi()*logit(mu_inf/mu_0);
    }else if(data_type == 4){
      c = pow(3,0.5)/pi()*logit((log(mu_inf)-shift)/(log(mu_0)-shift));
    }
  }

  d = exp(par4);
  k = log(par2);

  if(data_type == 1){
     b=exp(-k*d)*(pow(3,0.5)/pi()*logit(inv_logit(pi()/pow(3,0.5)*c)*(1+q))-c);
  }else if(data_type == 2){
    b=exp(-k*d)*(pow(3,0.5)/pi()*logit(inv_logit(pi()/pow(3,0.5)*c)+log(1+q)/a)-c);
  }else if(data_type == 3){
    b=exp(-k*d)*(pow(3,0.5)/pi()*logit(inv_logit(pi()/pow(3,0.5)*c)+q)-c);
  }else if(data_type == 4){
    b=exp(-k*d)*(pow(3,0.5)/pi()*logit(inv_logit(pi()/pow(3,0.5)*c)-log(1-q)/a)-c);
  }

  invsigma2=exp(par5);
}
model{
    par1 ~ pert_dist(shape1[1], shape2[1], priorlb[1], priorub[1]);
    par2 ~ pert_dist(shape1[2], shape2[2], priorlb[2], priorub[2]);
    par3 ~ pert_dist(shape1[3], shape2[3], priorlb[3], priorub[3]);
    par4 ~ normal(priormu[4],priorSigma[4,4]);
    par5 ~ normal(priormu[5],priorSigma[5,5]);

  if(is_increasing == 1){
    if(data_type == 1){
       for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*
     square(m[i]-(a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d)))))*invsigma2;
      }
    }else if(data_type == 2){
      for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*
     square(m[i]-(a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d)))))*invsigma2 - m[i]*n[i];
      }
    }

  }else if(is_decreasing == 1){
    if(data_type == 3){
      for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*
     square(m[i]-((a*(1+inv_logit(pi()/pow(3,0.5)*c)))-(a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d))))))*invsigma2;
      }
    }else if(data_type == 4){
      for (i in 1:N){
     target += -0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2)-0.5*(n[i]-1)*s2[i]*invsigma2-0.5*n[i]*
     square(m[i]-((a*(1+inv_logit(pi()/pow(3,0.5)*c)))-(a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d))))))*invsigma2 - m[i]*n[i];
      }
    }

  }

}

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
  int nlevels; //number of levels for the covariate
  int nlevels_a;
  int nlevels_BMD;
  int nlevels_d;
  int nlevels_sigma;
  int nlevels_b;
  vector[N] n;  // the sample size for each dose group
  vector[N] x;  // the dose level of each dose group
  vector[N] m;  // the arithmetic mean of the response values for each dose group
  vector[N] s2;  // the arithmetic variance of the response values for each dose group
  matrix[N, nlevels] trt_ind;
  real q;       // the BMR
  real shift; // data are shifted if different from 0
  matrix[5,nlevels] priormu;
  matrix[5,nlevels] priorlb; //lower bound
  matrix[5,nlevels] priorub; //upper bound
  matrix[5,nlevels] shape1; //alpha
  matrix[5,nlevels] shape2; //beta
  cov_matrix[5] priorSigma;
  int data_type; // data_type; 1 = increasing N, 2 = increasing LN, 3 = decreasing N, 4 = decreasing LN
  int<lower=0, upper=1> is_increasing; // indicator for increasing data
  real L; //lower bound for increasing data
  int<lower=0, upper=1> is_decreasing; // indicator for decreasing data
  real U; // upper bound for decreasing data
  real truncd;
}
parameters{
  real<lower=0> par1[nlevels_a];
  real<lower=0> par2[nlevels_BMD]; // BMD
  real<lower=0> pars3i[is_increasing]; // will be size one if is_increasing
  real<lower=0, upper=1> pars3d[is_decreasing]; // will be size one if is_decreasing
  real par4[nlevels_d];
  real par5[nlevels_sigma];
}
transformed parameters{
  real b[nlevels_b];
  real a[nlevels_a];
  real c[nlevels_a];
  real par3;
  real d[nlevels_d];
  real k[nlevels_BMD];
  real mu_inf[nlevels_a];
  real invsigma2[nlevels_sigma];
  real mu_0[nlevels_a];

  //d
  for(mn in 1:nlevels_d){
    d[mn] = exp(par4[mn]);
  }

  //a and c
  if(is_increasing){
    par3 = L + pars3i[1];
  }else if(is_decreasing){
    par3 = L + (U - L) .* pars3d[1];
  }

  for(mn in 1:nlevels_a){
    mu_0[mn] = par1[mn];

    if(data_type == 1 || data_type == 3){
      a[mn] = par1[mn];
    }else if(data_type == 2 || data_type == 4){
      a[mn] = log(par1[mn]) - shift;
    }
    mu_inf[mn] = par1[mn]*par3;

  }

  // sigma2
  for(mn in 1:nlevels_sigma)
  invsigma2[mn]=exp(par5[mn]);

  //BMD
  for(mn in 1:nlevels_BMD)
  k[mn] = log(par2[mn]);

  //c
  for(mn in 1:nlevels_a){
    if(data_type == 1 || data_type == 3){
      c[mn] = mu_inf[mn]/mu_0[mn];
    }else if(data_type == 2 || data_type == 4){
      c[mn] = (log(mu_inf[mn])-shift)/(log(mu_0[mn])-shift);
    }
  }

  //b
  if(nlevels_a == 1 && nlevels_d == 1 && nlevels_BMD == 1) {

    for(mn in 1:nlevels_b){

      if(data_type == 1){
        b[1]=exp(inv_Phi(q/(c[1]-1))-(k[1]*d[1]));
      }else if(data_type == 2){
        b[1]=exp(inv_Phi((log(1+q))/(a[1]*(c[1]-1)))-(k[1]*d[1]));
      }else if(data_type == 3){
        b[1]=exp(inv_Phi(-q/(c[1]-1))-(k[1]*d[1]));
      }else if(data_type == 4){
        b[1]=exp(inv_Phi((log(1-q))/(a[1]*(c[1]-1)))-(k[1]*d[1]));
      }

    }

  } else if(nlevels_a > 1 && nlevels_d == 1 && nlevels_BMD == 1) {

    for(mn in 1:nlevels_b){

      if(data_type == 1){
        b[mn]=exp(inv_Phi(q/(c[mn]-1))-(k[1]*d[1]));
      }else if(data_type == 2){
        b[mn]=exp(inv_Phi((log(1+q))/(a[mn]*(c[mn]-1)))-(k[1]*d[1]));
      }else if(data_type == 3){
        b[mn]=exp(inv_Phi(-q/(c[mn]-1))-(k[1]*d[1]));
      }else if(data_type == 4){
        b[mn]=exp(inv_Phi((log(1-q))/(a[mn]*(c[mn]-1)))-(k[1]*d[1]));
      }

    }

  } else if(nlevels_a == 1 && nlevels_d > 1 && nlevels_BMD > 1) {

    for(mn in 1:nlevels_b){

      if(data_type == 1){
        b[mn]=exp(inv_Phi(q/(c[1]-1))-(k[mn]*d[mn]));
      }else if(data_type == 2){
        b[mn]=exp(inv_Phi((log(1+q))/(a[1]*(c[1]-1)))-(k[mn]*d[mn]));
      }else if(data_type == 3){
        b[mn]=exp(inv_Phi(-q/(c[1]-1))-(k[mn]*d[mn]));
      }else if(data_type == 4){
        b[mn]=exp(inv_Phi((log(1-q))/(a[1]*(c[1]-1)))-(k[mn]*d[mn]));
      }

    }


  } else if(nlevels_a > 1 && nlevels_d > 1 && nlevels_BMD > 1) {

    for(mn in 1:nlevels_b){

      if(data_type == 1){
        b[mn]=exp(inv_Phi(q/(c[mn]-1))-(k[mn]*d[mn]));
      }else if(data_type == 2){
        b[mn]=exp(inv_Phi((log(1+q))/(a[mn]*(c[mn]-1)))-(k[mn]*d[mn]));
      }else if(data_type == 3){
        b[mn]=exp(inv_Phi(-q/(c[mn]-1))-(k[mn]*d[mn]));
      }else if(data_type == 4){
        b[mn]=exp(inv_Phi((log(1-q))/(a[mn]*(c[mn]-1)))-(k[mn]*d[mn]));
      }

    }

  }

}
model{

  //a prior
  for(i in 1:nlevels_a)
  par1[i] ~ pert_dist(shape1[1,i], shape2[1,i], priorlb[1,i], priorub[1,i]);

  //BMD prior
  for(i in 1:nlevels_BMD)
  par2[i] ~ pert_dist(shape1[2,i], shape2[2,i], priorlb[2,i], priorub[2,i]);

  // sigma prior
  for(i in 1:nlevels_sigma)
  par5[i] ~ normal(priormu[5,i],priorSigma[5,5]);

  //d prior
  for(i in 1:nlevels_d)
  par4[i] ~ normal(priormu[4,i],priorSigma[4,4])T[,truncd];

  // c prior (always constant)
  par3 ~ pert_dist(shape1[3,1], shape2[3,1], priorlb[3,1], priorub[3,1]);

  if(nlevels_d > 1 && nlevels_BMD > 1 && nlevels_a == 1 && nlevels_sigma == 1){

    //BMD & d on covariate
    if(data_type == 1 || data_type == 3){
      for (i in 1:N){
        for(mn in 1:nlevels){

          if (x[i]==0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*
          square(m[i]-a[1])*
          invsigma2[1])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*square(m[i]-a[1]-a[1]*(c[1]-1)*Phi(log(b[mn])+(d[mn]*log(x[i]))))*
          invsigma2[1])*trt_ind[i,mn];}

        }
      }
    }else if(data_type == 2 || data_type == 4){
      for (i in 1:N){
        for(mn in 1:nlevels){
          if (x[i]==0){target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*square(m[i]-a[1])*
          invsigma2[1])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*square(m[i]-a[1]-a[1]*(c[1]-1)*Phi(log(b[mn])+(d[mn]*log(x[i]))))*
          invsigma2[1] - m[i]*n[i])*trt_ind[i,mn];}
        }
      }
    }

  } else if(nlevels_a > 1 && nlevels_sigma > 1 && nlevels_BMD == 1 && nlevels_d == 1) {

    //a & sigma on covariate
    if(data_type == 1 || data_type == 3){
      for (i in 1:N){
        for(mn in 1:nlevels){

          if (x[i]==0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*
          square(m[i]-a[mn])*invsigma2[mn])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*square(m[i]-a[mn]-a[mn]*(c[mn]-1)*Phi(log(b[mn])+(d[1]*log(x[i]))))*
          invsigma2[mn])*trt_ind[i,mn];}

        }
      }
    } else if(data_type == 2 || data_type == 4){
      for (i in 1:N){
        for(mn in 1:nlevels){

          if (x[i]==0){target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*
          square(m[i]-a[mn])*invsigma2[mn] - m[i]*n[i])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*square(m[i]-a[mn]-a[mn]*(c[mn]-1)*Phi(log(b[mn])+(d[1]*log(x[i]))))*
          invsigma2[mn] - m[i]*n[i])*trt_ind[i,mn];}
        }
      }
    }


  } else if( nlevels_d > 1 && nlevels_BMD > 1 && nlevels_a > 1 && nlevels_sigma > 1) {

    //a, BMD, d & sigma on covariate
    if(data_type == 1 || data_type == 3){
      for (i in 1:N){
        for(mn in 1:nlevels){

          if (x[i]==0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*
          square(m[i]-a[mn])*invsigma2[mn])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*square(m[i]-a[mn]-a[mn]*(c[mn]-1)*Phi(log(b[mn])+(d[mn]*log(x[i]))))*
          invsigma2[mn])*trt_ind[i,mn];}

        }
      }
    } else if(data_type == 2 || data_type == 4){
      for (i in 1:N){
        for(mn in 1:nlevels){

          if (x[i]==0){target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*
          square(m[i]-a[mn])*invsigma2[mn] - m[i]*n[i])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[mn])-
          0.5*(n[i]-1)*s2[i]*invsigma2[mn]-0.5*n[i]*square(m[i]-a[mn]-a[mn]*(c[mn]-1)*Phi(log(b[mn])+(d[mn]*log(x[i]))))*
          invsigma2[mn] - m[i]*n[i])*trt_ind[i,mn];}
        }
      }
    }

  } else {
    //no dependence on covariates
    if(data_type == 1 || data_type == 3){
      for (i in 1:N){
        for(mn in 1:nlevels){

          if (x[i]==0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*
          square(m[i]-a[1])*invsigma2[1])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*square(m[i]-a[1]-a[1]*(c[1]-1)*Phi(log(b[1])+(d[1]*log(x[i]))))*
          invsigma2[1])*trt_ind[i,mn];}

        }
      }
    } else if(data_type == 2 || data_type == 4){
      for (i in 1:N){
        for(mn in 1:nlevels){

          if (x[i]==0){target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*
          square(m[i]-a[1])*invsigma2[1] - m[i]*n[i])*trt_ind[i,mn];}

          if (x[i]>0) {target += (-0.5*n[i]*log(2*pi())+0.5*n[i]*log(invsigma2[1])-
          0.5*(n[i]-1)*s2[i]*invsigma2[1]-0.5*n[i]*square(m[i]-a[1]-a[1]*(c[1]-1)*Phi(log(b[1])+(d[1]*log(x[i]))))*
          invsigma2[1] - m[i]*n[i])*trt_ind[i,mn];}
        }
      }
    }
  }
}

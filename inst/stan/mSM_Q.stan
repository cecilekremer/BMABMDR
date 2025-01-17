functions {
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
  int N;  // the total number of observations
  int Ndose; // the total number of distinct dose group
  int n_litter[Ndose]; // number of litters per dose group
  int maxl; // max number of litters per dose

  // matrix[Ndose, maxl] n;
  // matrix[Ndose, maxl] y;
  int n[Ndose, maxl];
  int y[Ndose, maxl];
  int use_data[Ndose, maxl];

  // int n[N];  // the sample size for each dose group
  // int y[N];  // the arithmetic mean of the response values for each dose group

  real priormu[2];
  real priorlb; //lower bound
  real priorub[2]; //upper bound
  real priorgama;
  real eps;
  int<lower=0, upper=1> is_bin;  //model type 1 = Binomial 0 = otherwise
  int<lower=0, upper=1> is_betabin;  //model type 1 = Beta-Binomial 0 = otherwise
  int<lower=0, upper=1> force_monotone;
}
parameters{
  row_vector[Ndose] par; // par[1]=background, par[2:N]=increment per dose group
  real<lower=0, upper=1> rho[is_betabin]; //will be defined if beta-binomial is to be fitted
}
transformed parameters{
  // vector[Ndose] a; // at the dose level
  matrix[Ndose, maxl] a;

  // real abet[N];
  // real bbet[N];
  row_vector[N] abet;
  row_vector[N] bbet;

  // mean in lowest dose
  // a[1] = par[1];
  //
  // for(k in 2:Ndose){
    //   a[k] = a[k-1] + par[k];
    //   if(a[k] <= 0){
      //     a[k] = 0.0001;
      //   }
      //   if(a[k] >= 1){
        //     a[k] = 0.9999;
        //   }
        // }
        a[1, ] = rep_row_vector(par[1], n_litter[1]);

        for(k in 2:Ndose){
          a[k, ] = a[k-1, ] + par[k];
          // if(a[k, 1] <= 0){
          //   a[k, ] = rep_row_vector(0.0001, n_litter[k]);
          // }
          // if(a[k, 1] >= 1){
          //   a[k, ] = rep_row_vector(0.9999, n_litter[k]);
          // }
        }

        if(is_bin == 0) {

          int j = 1;
          for(i in 1:Ndose){
            for(k in 1:n_litter[i]){
              abet[j] = (a[i,k])*((1/rho[is_betabin])-1.0);
              bbet[j] = (1.0 - (a[i,k]))*((1.0/rho[is_betabin])-1);
              j = j + 1;
            }
          }

        }else {
          for(i in 1:N){
            abet[i] = 0.0;
            bbet[i] = 0.0;
          }
        }


    // Adjust parameter matrix for data not used
    for(i in 1:Ndose){
      for(j in 1:maxl){
        if(use_data[i, j] == 0){
          a[i, j] = -1;
        }else{
          if(a[i, j] <= 0){
            a[i, j] = 0.0001;
          }
          if(a[i, j] >= 1){
            a[i, j] = 0.9999;
          }
        }
      }
    }

}
model{

  par[1] ~ pert_dist(priorlb, priormu[1], priorub[1], priorgama);

  for(k in 2:Ndose){
    if(force_monotone == 1){
      par[k] ~ uniform(0, priorub[2]);
      // par[k] ~ uniform(0, 1);
    }else{
      par[k] ~ uniform(-priorub[2], priorub[2]);//pert_dist(priorlb, priormu[1], priorub, 4);
      // par[k] ~ uniform(-1, 1);
    }
  }

  if(is_bin == 1) {

    // for (i in 1:N){
      //   // target += lchoose(n[i], y[i]) + y[i]*log(a[i]+eps) + (n[i] - y[i])*log(1 - a[i]+eps);
      //   target += binomial_lpmf(y[i] | n[i], a[i]);
      //
      // }
      for(k in 1:Ndose){
        for(i in 1:n_litter[k]){
          if(use_data[k, i] == 1){
            target += binomial_lpmf(y[k, i] | n[k, i], a[k, i]);
          }
        }
      }

  } else if(is_betabin == 1) {

    rho[is_betabin] ~ pert_dist(0.0, priormu[2], 1.0, 4.0);
    // for (i in 1:N){
      //   // target += lchoose(n[i], y[i]) + lgamma(abet[i]+y[i]+eps) + lgamma(bbet[i]+n[i]-y[i]+eps) -
      //   // lgamma(abet[i]+bbet[i]+n[i]+eps) - lgamma(abet[i]+eps) - lgamma(bbet[i]+eps) +
      //   // lgamma(abet[i]+bbet[i]+eps);
      //
      //   target += beta_binomial_lpmf(y[i] | n[i], abet[i], bbet[i]);
      // }

      int j = 1;
      for(k in 1:Ndose){
        for(i in 1:n_litter[k]){
          if(use_data[k, i] == 1){
            target += beta_binomial_lpmf(y[k, i] | n[k, i], abet[j], bbet[j]);
            j = j + 1;
          }
        }
      }
  }


}

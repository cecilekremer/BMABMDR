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
  int N;  // the total number of distinct dose groups
  int n[N]; // number of litters per dose group
  int nc; // number of unique dose x litter combinations (i.e. clusters)
  int maxN; // max number of obs per cluster
  int maxNc; // max number of litters per dose group
  int nij[N, maxNc]; // dose x litter matrix with the number of fetuses for each combination
  matrix[nc, maxN] y; // responses
  vector[N] x;  // the dose level of each dose group
  real q;       // the BMR
  real shift; // data are shifted if different from 0
  vector[6] priormu;
  vector[6] priorlb; //lower bound
  vector[6] priorub; //upper bound
  vector[6] shape1; //alpha
  vector[6] shape2; //beta
  cov_matrix[5] priorSigma;
  real truncd;
  int data_type; // data_type; 1 = increasing N, 2 = increasing LN, 3 = decreasing N, 4 = decreasing LN
  int<lower=0, upper=1> is_increasing; // indicator for increasing data
  real L; //lower bound for increasing data
  int<lower=0, upper=1> is_decreasing; // indicator for decreasing data
  real U; // upper bound for decreasing data
}
parameters{
  real<lower=0> par1;
  real<lower=0> par2; // BMD
  real<lower=0> pars3i[is_increasing]; // will be size one if is_increasing
  real<lower=0, upper=1> pars3d[is_decreasing]; // will be size one if is_decreasing
  real par4;
  real par5; // variance constant across cluster and dose groups
  real par6; // correlation parameter rho
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
  real rho_cluster;

  rho_cluster = par6;

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
  par4 ~ normal(priormu[4],priorSigma[4,4])T[,truncd];
  par5 ~ normal(priormu[5],priorSigma[5,5]);
  par6 ~ pert_dist(shape1[6], shape2[6], priorlb[6], priorub[6]);


  if(data_type == 1 || data_type == 3){

    int cnt;

    cnt = 1; // first cluster (ID)


    for(i in 1:N){ // for each dose group

    int nl;
    real mx;

    if(data_type == 1){
      mx = a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d)));
    }else if(data_type == 3){
      mx = (a*(1+inv_logit(pi()/pow(3,0.5)*c)))-(a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d))));
    }

    nl = n[i]; // number of litters in dose group i

    for(j in 1:nl){ // for each litter in dose group i

    int lt = nij[i, j]; // litter size

      row_vector[lt] resp;
      row_vector[lt] m;
      matrix[lt, lt] P;
      matrix[lt, lt] Sigma;

      resp = y[cnt, 1:lt];

      m = rep_row_vector(mx, lt);

      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1;
          }else{
            P[id1, id2] = rho_cluster;
          }
        }
      }

      Sigma = (1/invsigma2)*P;

      target += multi_normal_lpdf(resp | m, Sigma);


    cnt = cnt + 1; // next cluster

    }


    }

  }else if(data_type == 2 || data_type == 4){

    int cnt;

    cnt = 1; // first cluster (ID)


    for(i in 1:N){ // for each dose group

    int nl;
    real mx;

    if(data_type == 2){
      mx = a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d)));
    }else if(data_type == 4){
      mx = (a*(1+inv_logit(pi()/pow(3,0.5)*c)))-(a*inv_logit(pi()/pow(3,0.5)*(c+b*pow(x[i],d))));
    }

    nl = n[i]; // number of litters in dose group i

    for(j in 1:nl){ // for each litter in dose group i

    int lt = nij[i, j]; // litter size

      row_vector[lt] resp;
      row_vector[lt] m;
      matrix[lt, lt] P;
      matrix[lt, lt] Sigma;

      resp = y[cnt, 1:lt];

      m = rep_row_vector(mx, lt);

      for(id1 in 1:lt){
        for(id2 in 1:lt){
          if(id1 == id2){
            P[id1, id2] = 1;
          }else{
            P[id1, id2] = rho_cluster;
          }
        }
      }

      Sigma = (1/invsigma2)*P;

      target += multi_normal_lpdf(resp | m, Sigma);

      target += -(sum(resp));


    cnt = cnt + 1; // next cluster

    }


    }
  }
}

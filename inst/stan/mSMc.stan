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
  int N;  // the total number of distinct dose groups
  int n[N]; // number of litters per dose group
  int nc; // number of unique dose x litter combinations (i.e. clusters)
  int maxN; // max number of obs per cluster
  int maxNc; // max number of litters per dose group
  int nij[N, maxNc]; // dose x litter matrix with the number of fetuses for each combination
  matrix[nc, maxN] y; // responses
  // vector[N] x;  // the dose level of each dose group
  real q;       // the BMR
  real shift; // data are shifted if different from 0
  vector[N+2] priormu;
  cov_matrix[N+1] priorSigma;
  real priorlb;
  vector[2] priorub;
  real priorg;
  int data_type; // data_type; 1 = increasing N, 2 = increasing LN, 3 = decreasing N, 4 = decreasing LN
}
parameters{
  vector[N+2] par; // par[1] = background, par[2:N]=increment per dose group, par[N+1]=log(invsigma2), par[N+2]=rho_cluster
}
transformed parameters{
  vector[N] a;
  vector[N] mu;
  real invsigma2;
  real rho_cluster;

  rho_cluster = par[N+2];

  // means per dose group
  mu[1] = par[1];
  for(k in 2:N){
    mu[k] = mu[k-1] + par[k];
  }

  // parameter a
  if(data_type == 1 || data_type == 3){
    a = mu;
  }else if(data_type == 2 || data_type == 4){
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
  par[N+2] ~ pert_dist(0, 0.5, 1, 0.0001);

  if(data_type == 1 || data_type == 3){

    int cnt;

    cnt = 1; // first cluster (ID)

    for(i in 1:N){ // for each dose group

    int nl;
    // real mx;

    // mx = a;

    nl = n[i]; // number of litters in dose group i

    for(j in 1:nl){ // for each litter in dose group i

    int lt = nij[i, j]; // litter size

    row_vector[lt] resp;
    row_vector[lt] m;
    matrix[lt, lt] P;
    matrix[lt, lt] Sigma;

    resp = y[cnt, 1:lt];

    // m = rep_row_vector(mx, lt);
    m = rep_row_vector(a[i], lt);

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

    // mx = a;

    nl = n[i]; // number of litters in dose group i

    for(j in 1:nl){ // for each litter in dose group i

    int lt = nij[i, j]; // litter size

    row_vector[lt] resp;
    row_vector[lt] m;
    matrix[lt, lt] P;
    matrix[lt, lt] Sigma;

    resp = y[cnt, 1:lt];

    m = rep_row_vector(a[i], lt);

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

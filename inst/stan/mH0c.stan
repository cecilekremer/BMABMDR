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
  vector[2] priormu;
  cov_matrix[2] priorSigma;
  real priorlb;
  real priorub;
  real priorg;
  int data_type; // data_type; 1 = increasing N, 2 = increasing LN, 3 = decreasing N, 4 = decreasing LN
}
parameters{
  vector[3] par; // par[1] = overall mean, par[2] = log(invsigma2), par[3] = rho
}
transformed parameters{
  real a;
  real invsigma2;
  real rho_cluster;

  rho_cluster = par[3];

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
  par[3] ~ pert_dist(0, 0.5, 1, 0.0001);

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
    m = rep_row_vector(a, lt);

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
  //  real mx;

  //  mx = a;

    nl = n[i]; // number of litters in dose group i

    for(j in 1:nl){ // for each litter in dose group i

    int lt = nij[i, j]; // litter size

    row_vector[lt] resp;
    row_vector[lt] m;
    matrix[lt, lt] P;
    matrix[lt, lt] Sigma;

    resp = y[cnt, 1:lt];

//    m = rep_row_vector(mx, lt);
    m = rep_row_vector(a, lt);
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

data {
  int<lower=1> N;  // total number of observations
  // array[N] int Y;  // response variable
  // array[N] int trials;  // number of trials
  int Y[N];
  int trials[N];
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
}
model {
  // likelihood including constants
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept;
    target += binomial_logit_lpmf(Y | trials, mu);
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}

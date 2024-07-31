data {
  int<lower=1> N;
  int<lower=1> P;
  vector[N] ybar;
  vector[N] se;
  matrix[N, P] X;  // design matrix for fixed + random effects
  int<lower=1> N_pred;
  matrix[N_pred, P] X_pred;
  // vector[N_pred] se_pred;
}

transformed data {
real L = 0.;
real D = 1.;
}

parameters {
  real alpha;
  vector[P] beta;
  // real<lower=0> D;
  // real L;
  real lognu;
}

transformed parameters {
  vector[N] eta = alpha + X*beta;
  real<lower=0> nu = exp(lognu);
  vector[N] yhat;
  for (i in 1:N) yhat[i] = L + D*student_t_cdf(eta[i], nu, 0, 1);
}

model {
// {
//   vector[N] std;
//   for (i in 1:N) std[i] = sqrt(pow(se[i],2) + pow(sigma,2));
// }
  yhat ~ normal(ybar, se);

//   alpha ~ normal(0, 0.5);
//   beta ~ normal(0, 0.25);
//   lognu ~ normal(0, 2.5);
//   L ~ normal(0, 0.2);
//   D ~ normal(1, 0.2);
}

generated quantities {
  vector[N_pred] ypred;
  {
    vector[N_pred] eta_pred = alpha + X_pred*beta;
    for (i in 1:N_pred) ypred[i] = L + D*student_t_cdf(eta_pred[i], nu, 0, 1);
  }
}
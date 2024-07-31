data {
  int<lower=1> N;
  int<lower=1> P;
  vector[N] ybar;
  vector[N] M;
  matrix[N, P] X;  // design matrix for fixed + random effects
  int<lower=1> N_pred;
  matrix[N_pred, P] X_pred;
  real<lower=0> m0;
}

transformed data {
  int L = 0;
  int D = 1;
}

parameters {
  real alpha;
  vector[P] beta;
  // real<lower=0> D;
  // real L;
  real lognu;
  // real<lower=0> sigma;
}

transformed parameters {
  vector[N] eta = alpha + X*beta;
  real<lower=0> nu = exp(lognu);
  vector[N] yhat;
  for (i in 1:N) yhat[i] = L + D*student_t_cdf(eta[i], nu, 0, 1);
}

model {
  yhat ~ beta(ybar.*M + m0, (1-ybar).*M + m0);

  // alpha ~ normal(0, 1);
  // beta ~ normal(0, 1);
  // lognu ~ normal(0, 2.5);
  // L ~ normal(0, 0.2);
  // D ~ normal(1, 0.2);
}

generated quantities {
  vector[N_pred] ypred;
  // real pp_lr = 0;
  {
    vector[N_pred] eta_pred = alpha + X_pred*beta;
    for (i in 1:N_pred) ypred[i] = L + D*student_t_cdf(eta_pred[i], nu, 0, 1);
    // pp_lr = normal_lpdf(ybar | yhat, se_pred) - normal_lpdf(ypred | yhat, se_pred);
  }
}

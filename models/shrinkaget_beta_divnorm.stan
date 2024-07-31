data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> Q;
  vector[N] ybar;
  vector[N] M;
  matrix[N, P] X;  // design matrix for fixed + random effects
  matrix[N, Q] Z;  // design matrix for divnorm
  real<lower=0> m0;
}

transformed data {
  real L = 0.;
  real D = 1.;
  // real nu = 100;
}

parameters {
  real alpha;
  vector[P] beta;
  // real<lower=0> D;
  // real L;
  real lognu;
  vector[Q] gamma;
  // real kappa;
}

transformed parameters {
  vector[N] yhat;
  {
    vector[N] sigma = exp(Z*gamma);
    vector[N] eta = alpha + X*beta;
    real nu = exp(lognu);
    for (i in 1:N) yhat[i] = L + D*student_t_cdf(eta[i], nu, 0, sigma[i]);
  }
}

model {
// {
//   vector[N] std;
//   for (i in 1:N) std[i] = sqrt(pow(se[i],2) + pow(sigma,2));
// }
    yhat ~ beta(ybar.*M + m0, (1-ybar).*M + m0);

//   alpha ~ normal(0, 0.5);
//   beta ~ normal(0, 0.25);
//   lognu ~ normal(0, 2.5);
//   L ~ normal(0, 0.2);
//   D ~ normal(1, 0.2);
}

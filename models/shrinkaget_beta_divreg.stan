data {
  int<lower=1> N;
  int<lower=1> P;
  vector[N] ybar;
  vector[N] M;
  matrix[N, P] X;  // design matrix for fixed + random effects
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
  real l0;
  real l1;
  real lognu;
}

transformed parameters {
  vector[N] yhat;
  { 
    vector[P] gamma = l0 + l1.*abs(beta);
    vector[N] sigma = exp(X*gamma);
    vector[N] eta = alpha + X*beta;
    real nu = exp(lognu);
    for (i in 1:N) yhat[i] = L + D*student_t_cdf(eta[i], nu, 0, sigma[i]);
  }
}

model {
  
  yhat ~ beta(ybar.*M + m0, (1-ybar).*M + m0);
}

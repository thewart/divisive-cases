data {
  int<lower=1> N;
  int<lower=1> P;
  vector[N] ybar;
  vector[N] se;
  matrix[N, P] X;  // design matrix for fixed + random effects
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
  real l0;
  real l1;
  real lognu;
  // real kappa;
}

transformed parameters {
  vector[N] yhat;
  { 
    vector[P] gamma = l0 + l1*abs(beta);
    vector[N] sigma = exp(X*gamma);
    vector[N] eta = alpha + X*beta;
    real nu = exp(lognu);
    for (i in 1:N) yhat[i] = L + D*student_t_cdf(eta[i], nu, 0, sigma[i]);
  }
}

model {
  
  yhat ~ normal(ybar, se);
}

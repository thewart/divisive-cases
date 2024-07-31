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
array[11] int backmap = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 3};
// real nu = 100;
}

parameters {
  real alpha;
  vector[P] beta;
  // real<lower=0> D;
  // real L;
  vector[3] l0;
  vector[2] l1;
  real lognu;
  // real kappa;
}

transformed parameters {
  vector[N] yhat;
  { 
    vector[3] lambda_0 = [l0[1]+l0[2], l0[2], l0[3]+l0[2]]';
    vector[3] lambda_1 = [l1[1], 0, l1[2]+l1[1]]';
    vector[P] gamma = lambda_0[backmap] + lambda_1[backmap].*abs(beta);
    vector[N] sigma = exp(X*gamma);
    vector[N] eta = alpha + X*beta;
    real nu = exp(lognu);
    for (i in 1:N) yhat[i] = L + D*student_t_cdf(eta[i], nu, 0, sigma[i]);
  }
}

model {
  
  yhat ~ normal(ybar, se);
}

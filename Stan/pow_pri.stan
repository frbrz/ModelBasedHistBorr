// This Stan file defines a simple power prior

data {
  int<lower = 0> N; // Number of OBS
  int<lower = 1> M; // Number of covariates 
  real<lower = 0, upper = 1> alpha; //power prior
  vector[N] y;  // observation
  matrix[N, M] X; // matrix of covariates, first col is intercept
}

parameters {
  vector[M] beta;   // covariate coeff, beta[1] intercept, i.e. mean trt eff
  real<lower=0> sigma; // btwn trial variability
}

model {
  // expected value for each observation
  vector[N] mu = X * beta;

  // Log-ikelihood
  // alpha is the power (that multiplies cos on log)
  target += alpha * normal_lpdf(y | mu, sigma);
  
  // Weakly info priors
  beta[1] ~ normal(10, 5); // intercept
  beta[2] ~ normal(0, 5);  // covariate effect
  sigma ~ normal(0, 10);
}


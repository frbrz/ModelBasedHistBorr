// This Stan file defines a simple individual level
// meta-analytic prior

data {
  int<lower=0> N; // Number of OBS
  int<lower=1> M; // Number of covariates 
  int<lower=1> S; // Number of studies
  vector[N] y;  // observation
  int<lower = 1, upper = S> std[N]; // study index for each observation
  matrix[N, M] X; // matrix of covariates, first col is intercept
}

parameters {
  vector[M] beta;   // covariate coeff, beta[1] intercept, i.e. mean trt eff
  real<lower=0> tau; //btwn trial variability
  real<lower=0> sigma; // btwn trial variability
  vector[S] eta_raw; // vector of random effects
}

transformed parameters{
  // Non-centered parameterization for eta
  vector[S] eta = eta_raw * tau;
}

model {
  // expected value for each observation
  vector[N] mu = X * beta + eta[std];

  // Likelihood
  y ~ normal(mu, sigma);
  
  // Weakly info priors
  eta_raw ~ std_normal();
  beta[1] ~ normal(10, 5); // intercept
  beta[2] ~ normal(0, 5);  // covariate effect
  sigma ~ normal(0, 10);
  tau ~ normal(0, 5);
}

generated quantities {
  real theta_star = normal_rng(beta[1], tau);
}
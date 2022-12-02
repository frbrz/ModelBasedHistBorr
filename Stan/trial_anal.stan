// Stan model for trial analysis of a Normal endpoint

data {
  int<lower=1> N;  // number of observations
  vector[N] Y;     // response variable (CVA)
  matrix[N, 2] X;  // design matrix
                   // col1 = TRT allocation
                   // col2 = X_VA0 (centered)
  
  // // Prior hyperparameters
  
  // Mixture priors parameters for SOC term
  // 2 terms = MAP prior
  real<lower = 0, upper = 1> mix_w;
  
  // Means
  vector[2] mu_soc;         // SOC effect 
  real mu_trt;              // Treatment effect
  real mu_va0_cov;          // VA0 covariate effect
  real<lower = 0> mu_sig;   // Observation error

  // Standard deviations
  vector<lower = 0>[2] sd_soc;  // SOC effect
  real<lower = 0> sd_trt;       // Treatment effect
  real<lower = 0> sd_va0_cov;   // VA0 covariate effect
  real<lower = 0> sd_sig;       // Observation error
}

parameters {
  real soc_eff;  // intercept (i.e. SOC treatment effect)
  real trt_eff;
  real va0_eff;
  real<lower=0> sigma_raw;  // residual SD
}

transformed parameters {
  vector[2] beta;  // vector of model params (to use efficient normal_id_glm)
                   // beta[1] treatment effect
                   // beta[2] covariate (X_VA0) effect
  real<lower = 0> sigma; // observation noise sd  
  
  beta[1] = trt_eff;
  beta[2] = va0_eff;
  
  // Non centered parameerization of sigma
  sigma = mu_sig + sigma_raw * sd_sig;  
}

model {
  // model is:
  // alpha + beta_1 * TRT + beta_2 * X_VA0
  // for computational efficacy rewritten with beta as vectors
  // alpha + beta * X 
  // X[, 1]: TRT and X[, 2] X_VA0
  
  // Likelihood 
  Y ~ normal_id_glm(X, soc_eff, beta, sigma);
  
  // Priors
  sigma_raw ~ std_normal();
  trt_eff ~ normal(mu_trt, sd_trt);
  va0_eff ~ normal(mu_va0_cov, sd_va0_cov);
  
  // Prior on SOC effect can be a mixture (MAP prior)
  // code to be efficient needs if-else
  // for 1 (i.e. no mixture) vs 2 components mixtures
  if (mix_w == 1) {
    // If no-mixture efficient coding
    soc_eff ~ normal(mu_soc[1], sd_soc[1]);
  } else {
    // If mixture less efficient coding
    target += log_mix(
      mix_w,
      normal_lpdf(soc_eff | mu_soc[1], sd_soc[1]),
      normal_lpdf(soc_eff | mu_soc[2], sd_soc[2])
    );
  }
}

// generated quantities needed for power-scale sensitivity
generated quantities {
  vector [N] log_lik ; // log likelihood
  real lprior_soc ; // temporary log prior for SOC parameter
  real lprior ; // joint log prior
  
  // log likelihood (individual)
  for (n in 1:N) { 
    log_lik[n] = normal_lpdf(Y[n]| X[n] * beta + soc_eff, sigma);
  }
  
  // log prior
  if (mix_w == 1) {
    // If no-mixture efficient coding
    lprior_soc = normal_lpdf(soc_eff | mu_soc[1], sd_soc[1]);
  } else {
    // If mixture less efficient coding
    lprior_soc = log_mix(
      mix_w,
      normal_lpdf(soc_eff | mu_soc[1], sd_soc[1]),
      normal_lpdf(soc_eff | mu_soc[2], sd_soc[2])
    );
  }
//  lprior = lprior_soc 
//         + normal_lpdf(trt_eff|mu_trt, sd_trt) 
//         + normal_lpdf(va0_eff | mu_va0_cov, sd_va0_cov)
//         + normal_lpdf(sigma|mu_sig, sd_sig);
  
  // Only care about power scaling soc_eff prior
  lprior = lprior_soc;
}


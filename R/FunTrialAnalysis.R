## AUTHOR: Francesco Brizzi
## AIM: Helper functions to run trial simulations

# Libraries needed -------------------------------------------------------------

on_cluster <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

if (on_cluster) {
  # Reduce amount of output printed on cluster
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(checkmate))
  suppressPackageStartupMessages(library(cmdstanr))
  suppressPackageStartupMessages(library(posterior))
  suppressPackageStartupMessages(library(priorsense))
} else {
  library(tidyverse)
  library(checkmate)
  library(cmdstanr)
  library(posterior)
  library(priorsense)
}

# Random helper functions ------------------------------------------------------

#' Function to define a Normal mixture
#' 
#' This was inspired from rBesT.
#' But wanted to avoid depending on rBesT as it calls stanHeaders (rstan).
#' stanHeaders occasionaly caused unresolvable conflicts with cmdstanR.
#' 
#' @param w numeric vector containing the weight of each mixture component
#' @param m numeric vector containing the means of each mixture component
#' @param sd numeric vector containing the standard deviation of
#'        each mixture component
#' @return An object of class `norm_mix` which is a dataframe of 
#' where each column is a normal mixture characterized by a weight (row1),
#' a mean (row2) and a standard deviation (row3).
norm_mix <- function(w, m, sd) {
  
  ### Argument checks
  assert_numeric(w)
  assert_numeric(m)
  assert_numeric(sd, low = 1e-6)
  if (n_distinct(sapply(list(w, m, sd), length)) != 1) {
    stop("w, m, sd must be of equal length")
  }
  if (sum(w) != 1) {
    stop("weights (w) of mixture should sum to one")
  }
  
  # Return
  out <- t(data.frame("w" = w, "m" = m, "s" = sd))
  structure(out, class = "norm_mix")
}

# Helper functions for cmdstanR ------------------------------------------------

#' Function to automatically assess convergence from cmdstanr output
#'
#' Non convergence defined if one rhat > 1.05 or one divergent transition
#' @param obj a CmdStanMCMC object
#' @return logical, equal to true if convergence achieved wrt above criteria
assess_convg <- function(obj) {
  assert_class(obj, "CmdStanMCMC")
  
  # Extracting convergence summaries of interest
  r_hats <- obj$summary(variables = NULL, "rhat")$rhat
  divgs <- subset_draws(obj$sampler_diagnostics(), variable = "divergent__")
  n_divg <- sum(apply(divgs, 2, sum))
  
  # Convergence as specified above
  all(r_hats < 1.05) & (n_divg == 0)
}

#' Function to run a cmdstanR model until convergence or timeout
#' 
#' @param model a CmdStanModel object.
#' @param timeout numeric. Seconds after which model stops running
#' (i.e. stan is timed-out)
#' @param max_try integerish. Maximum number of attempts to run the model.
#' Note: failure for non-convergence or timeout are not distinguished atm.
#' @return list. list made of
#'  \itemize{
#' \item `fit`: a cmdstanR fit, for posterior analysis.
#' Note this returns NA if failure reason is timeout.
#' \item `conv`: a logical, denoting whether convergence is achieved.
#' Note that convergence is asssessed using `assess_convg` fun
#' \item `timeout`: a logical, denoting whether the stan run was timed-out.
#' }
run_cmdstan_convg <- function(model,
                              timeout, 
                              max_try = 4,
                             ...) {
  assert_class(model, "CmdStanModel")
  assert_numeric(timeout, len = 1)
  assert_integerish(max_try, len = 1, lower = 2)
  
  ### First run, with timeout
  dots <- list(...)
  try_i <- 1
  
  # https://github.com/stan-dev/cmdstanr/issues/474
  # for inspiration on cmdstanr timeout
  fit <- try(
    R.utils::withTimeout(
      do.call(model[["sample"]], dots),
      timeout = timeout,
      cpu = Inf,
      silent = T
    )
  )
  
  # If the fit has timeout or not converged re-run
  if (is(fit, "try-error")) {
    keep_run <- T
  } else {
    conv <- assess_convg(fit)
    keep_run <- !conv
  }
  
  while(keep_run) {
    
    # Update number of attempts
    try_i <- try_i + 1
    
    # Try a new seed
    dots$seed <- dots$seed + 1001
    
    # Re-run
    fit <- try(
      R.utils::withTimeout(
        do.call(model[["sample"]], dots),
        timeout = timeout,
        cpu = Inf,
        silent = T
      )
    )
    
    # If the fit has timeout or not converged re-run
    if (is(fit, "try-error")) {
      keep_run <- T
    } else {
      conv <- assess_convg(fit)
      keep_run <- !conv
    }
    
    # Stop running when we reached maximum number of tries
    if (try_i == max_try) {
      keep_run <- F
    }
  }
  
  # Returned object depends on timeout
  time_flag <- is(fit, "try-error")
  
  if (time_flag) {
    out <- NA
    conv <- NA
  } else {
    out <- fit
  }
  
  ### Return
  list(fit = out, conv = conv, timeout = time_flag)
}


# Functions to analyse Bayesian trial ------------------------------------------

#' Function to calculate prior / data commensurability (delta)
#' based on the metric by Gravestock and Held (2017)
#' assuming that the prior is normally distributed.
#' 
#' @param data data-frame. Should contains the endpoint and covariate data
#' from the trial, but only for the control arm.
#' @param mu_h numeric. Contains mean of the prior.
#' @param sd_h numeric. Contains standard deviation of the prior.
#' @param zs logical. If TRUE (default) the delta calculation is rounded,
#' this can avoid numerical issues in algorithm for calculation.
#' @return A list of two objects:
#' - `delta` the commensurability parameter
#' - `sd` the "adapted" standard deviation of the power prior 
#' based on commensurability
calc_commens <- function(data, mu_h, sd_h, zs = T) {
  assert_data_frame(data)
  must_nm <- c("Y", "X_VA0")
  assert_names(colnames(data), must.include = must_nm)
  assert_numeric(mu_h, len = 1)
  assert_numeric(sd_h, len = 1, lower = 0)
  assert_logical(zs, len = 1)
  
  # Fit regression on current data
  reg_cur <- summary(lm(Y ~ 1 + X_VA0, data = data))
  
  # Extract maximum likelihood based estimates of treatment effect and sd
  mu_c <- reg_cur$coeff[1, 1]
  sd_c <- reg_cur$coeff[1, 2]
  
  # Calculate delta, commensurability parameter
  # https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fpst.1814&file=pst1814-sup-0001-Supplementary.pdf
  # See Section A.3
  # delta = 1 max commensurability
  # delta -> 0 low commensurability
  delta <- sd_h^2 / (max((mu_h - mu_c)^2 , (sd_h^2) + (sd_c^2)) - sd_c^2)
  
  # Approximate to avoid numerical error issues
  if (zs) {
    delta <- zapsmall(delta)
  }
  
  # Standard deviation of prior, as function of commensurability value
  sd_pri <- switch(
    as.character(delta == 1),
    "TRUE"  = sd_h,
    "FALSE" = sqrt((mu_h - mu_c)^2  - sd_c^2)
  )
  
  # Return: 1) delta val 2) sd for pow pri 
  list("delta" = delta, "sd_pri" = sd_pri)
}

#' Function to run Bayesian analysis of a trial
#' 
#' @param data_tr (data frame) containing trial data.
#'        Must have cols `ID`, `X_VA0`, `Y`, `ARM`
#' @param true_trt_eff true treatment effect of `data_tr`
#'        Needed to get BIAS, SD and RMSD of treatment effect
#' @param ni_margin (numeric) chosen non-inferiority margin
#' @param alpha (numeric in 0-1) the chosen significance level
#' @param mod_stan (CmdStanModel object) containing a stan model
#' @param soc_pri (norm-mix object) containing the definition of the prior
#'        for the standard-of-care treatment effect.
#' @param trt_pri (norm-mix object) containing the definition of the prior
#'        for the experimental treatment effect.
#' @param cov_va0_pri (norm-mix object) containing the definition of the prior
#'        for the baseline covariate effect
#' @param sd_pri (norm-mix object) containing the definition of the prior
#'        for the standard deviation of the residual error
#' @param power (logical) if set to true the `soc_pri` prior is robustified,
#'        based on Gravestock and Held (2019)
#' @param args_cmdstan (list) parameters that specify the Stan's HMC run. 
#'        For full details see \link[cmdstanr]{sample} 
#' @param seed (integerish) Seed for the current analysis.
#' @return A list of three elements.
#'  \itemize{
#'   \item CONV - Contains information on whether convergence 
#'   (no rhat > 1.05, no divergent transitions) has been achieved and whether
#'   timeout has been reached
#'   \item RES - Contains all the important quantities of
#'    the posteriors of interest (eg mean, 90% intervals...) as well as
#'    decision made (ie. accept / reject non-inferiority)
#'   \item SENS - Contains output from power-sensitivity analysis.
#' }
analyse_trial <- function(data_tr, 
                          true_trt_eff,
                          ni_margin = 3.9,
                          alpha = 0.025,
                          mod_stan,
                          soc_pri,
                          trt_pri,
                          cov_va0_pri,
                          sd_pri,
                          power = F,
                          args_cmdstan,
                          seed = 1825) {
  ### Arguments checks
  assert_data_frame(data_tr)
  assert_names(colnames(data_tr), must.include = c("ID", "X_VA0", "Y", "ARM"))
  assert_numeric(ni_margin, len = 1)
  assert_numeric(alpha, lower = 0, upper = 1, len = 1)
  assert_class(mod_stan, "CmdStanModel")
  assert_class(soc_pri, "norm_mix")
  assert_class(trt_pri, "norm_mix")
  assert_class(cov_va0_pri, "norm_mix")
  assert_class(sd_pri, "norm_mix")
  if (ncol(soc_pri) > 2) {
    stop("soc_pri can only be a norm_mix object with at most two components")
  }
  if (ncol(trt_pri) != 1 | ncol(cov_va0_pri) != 1 | ncol(sd_pri) != 1) {
    txt <- paste(
      "trt_pri, cov_va0_pri, sd_pri can only be norm_mix object",
      "with at most one component"
    )
    stop(txt)
  }
  assert_list(args_cmdstan)
  must_cmdstan <- c(
    "chains", "parallel_chains", "refresh", "iter_warmup", "iter_sampling", 
    "adapt_delta", "max_treedepth", "show_messages", "timeout", "max_try"
  )
  assert_names(names(args_cmdstan), must.include = must_cmdstan)
  assert_logical(power, len = 1)
  assert_integerish(seed, len = 1)
  
  ### Dealing with MAP priors
  n_mix_soc_pri <- ncol(soc_pri)
  map_pri_flag <- n_mix_soc_pri == 2
  if (map_pri_flag) {
    lam <- soc_pri[1, 1]
    mu_soc <- soc_pri["m", ]
    sd_soc <- soc_pri["s", ]
  } else {
    lam <- 1L
    # IMPORTANT note: 999 is a dummy number which is later ignored by stan
    # (needed so that priors are always specified as two component mixture,
    # also when there is only one component)
    mu_soc <- c(soc_pri["m", 1], 999)
    sd_soc <- c(soc_pri["s", 1], 999)
  }
  
  ### Dealing with power priors
  if (power) {
    if (map_pri_flag) {
      stop("power prior does not work for mixture (MAP) priors")
    }
    # Calculate delta parameter from Gravenstock and Held for normal dist
    comm <- calc_commens(
      data = data_tr %>% filter(ARM == 0),
      mu_h = soc_pri["m", 1], 
      sd_h = soc_pri["s", 1], 
      zs = T    
    )
    sd_soc <- c(comm$sd_pri, 999)
  }
  
  ### Stan data
  data_stan <- list(
    # actual trial data
    N = n_distinct(data_tr$ID),
    Y = data_tr$Y,
    X = as.matrix(data_tr[, c("ARM", "X_VA0")]),
    # Prior parameters
    N_mix_pri = n_mix_soc_pri,
    mix_w = lam,
    mu_soc = mu_soc,
    mu_trt = trt_pri["m", 1],
    mu_va0_cov = cov_va0_pri["m", 1],
    mu_sig = sd_pri["m", 1],
    sd_soc = sd_soc,
    sd_trt = trt_pri["s", 1],
    sd_va0_cov = cov_va0_pri["s", 1],
    sd_sig = sd_pri["s", 1]
  )
  
  ### Run Stan HMC
  # (until convergence or timeout)
  cmdstan_out <- run_cmdstan_convg(
    model = mod_stan,
    max_try = args_cmdstan$max_try,
    timeout = args_cmdstan$timeout,
    data = data_stan,
    seed = seed,
    chains = args_cmdstan$chains,
    parallel_chains = args_cmdstan$parallel_chains,
    refresh = args_cmdstan$refresh,
    iter_warmup = args_cmdstan$iter_warmup,
    iter_sampling = args_cmdstan$iter_sampling,
    adapt_delta = args_cmdstan$adapt_delta,
    max_treedepth = args_cmdstan$max_treedepth,
    show_messages = args_cmdstan$show_messages
  )
  
  # Define data frames for time-out cases
  conv <- with(cmdstan_out, tibble("CONV" = conv, "TIMEOUT" = timeout))
  
  ### Output 
  # If timeout cmdstan_out has no fit object -> return NA df
  if (cmdstan_out$timeout) {
    res <- tibble(
      "DEC" = NA, "NI_PROB" = NA,
      "TRT_mean" = NA, "SOC_mean" = NA, "SDRES_mean" = NA, "VA0_mean" = NA, 
      "TRT_sd" = NA, "SOC_sd" = NA, "SDRES_sd" = NA, "VA0_sd" = NA, 
      "TRT_Q025" = NA, "SOC_Q025" = NA, "SDRES_Q025" = NA, "VA0_Q025" = NA, 
      "TRT_Q975"= NA, "SOC_Q975" = NA, "SDRES_Q975", "VA0_Q975" = NA,
      "POW_DELTA" = NA, "POW_SD_SOC" = NA
    )
    
    pss_summ <- tibble(
      "variable" = NA, "mean" = NA, "median" = NA, "sd" = NA, "alpha" = NA,
      "n_eff" = NA, "pareto_k" = NA, "cjs_dist" = NA, "component" = NA
    )
  
  } else {
    
    fit <- cmdstan_out$fit
    
    ### Extracting posterior draws
    fit_draws <- as_draws_df(
      subset_draws(
        fit$draws(), 
        variable = c("soc_eff", "trt_eff", "va0_eff", "sigma")
      )
    )
    
    ### Running power sensitivity analysis
    pss <- powerscale_sequence(fit)
    pss_summ <- summarise_draws(pss)[[1]] %>% 
      filter(variable == "trt_eff") %>% 
      select(variable, mean, median, sd, alpha, n_eff, pareto_k, cjs_dist, component) %>% 
      as_tibble()
    
    ### Bayesian decision making
    # Calculate posterior probability
    # i.e. how much of the mass is below the non inferiority margin
    p_prob <- mean(fit_draws$trt_eff < -ni_margin)
    # if more mass than specified by alpha below NI margin
    #  -> can NOT claim that new treatment is NI to SOC (dec = F)
    # if less mass than specified by alpha below NI margin
    #  -> can claim that new treatment is NI to SOC (dec = T)
    dec <- (p_prob < alpha)
    
    ### Treatment effect bias, and RMSD
    bias_trt <- mean(fit_draws$trt_eff - true_trt_eff)
    rmsd_trt <- sqrt(mean(((fit_draws$trt_eff - true_trt_eff) ^ 2)))
    
    ### Output processing
    res <- left_join(
      summarize_draws(fit_draws, "mean", "sd"),
      summarize_draws(fit_draws,  ~ quantile(., probs = c(0.025, 0.975))),
      by = "variable"
    ) %>% 
      # Nice variable naming
      mutate(variable = case_when(
        variable == "trt_eff" ~ "TRT",
        variable == "va0_eff" ~ "VA0",
        variable == "soc_eff" ~ "SOC",
        variable == "sigma" ~ "SDRES"
      )) %>% 
      rename(Q025 = "2.5%", Q975 = "97.5%") %>% 
      # Output in one line (useful for many simulations)
      pivot_wider(
        names_from = variable, 
        values_from = c("mean", "sd", "Q025", "Q975"),
        names_glue = "{variable}_{.value}"
      ) %>% 
      # adding RMSD and BIAS
      mutate(TRT_BIAS = bias_trt) %>% 
      mutate(TRT_RMSD = rmsd_trt) %>% 
      # Bayesian decision
      mutate(NI_PROB = p_prob) %>% 
      mutate(DEC = dec) %>% 
      # Power prior output
      mutate(POW_DELTA = ifelse(power, comm$delta, NA_real_)) %>% 
      mutate(POW_SD_SOC = ifelse(power, comm$sd_pri, NA_real_)) %>% 
      # Nice order
      select(
        DEC, NI_PROB,
        TRT_mean, SOC_mean, SDRES_mean, VA0_mean,
        TRT_sd, SOC_sd, SDRES_sd, VA0_sd,
        TRT_Q025, SOC_Q025, SDRES_Q025, VA0_Q025,
        TRT_Q975, SOC_Q975, SDRES_Q975, VA0_Q975,
        TRT_BIAS, TRT_RMSD,
        POW_DELTA, POW_SD_SOC
      )
  }
  
  ### Return
  # Three dataframes
  # CONV: information on convergence and timeout
  # RES: posterior summaries from trial
  # SENS: prior sensitivity analysis
  list(CONV = conv, RES = res, SENS = pss_summ)
}

#' Wrapper function to run a trial with 4 different priors
#' (non-informative, informative, informative + power, mixture)
run_trial_i <- function(data_i, 
                        true_trt_eff,
                        ni_margin = 3.9,
                        alpha = 0.025,
                        mod_stan,
                        soc_pri_non_info,
                        soc_pri_info,
                        soc_pri_mix,
                        trt_pri,
                        cov_va0_pri,
                        sd_pri,
                        args_cmdstan,
                        seed = 1825) {
  
  ### Initialization
  res_l <- vector("list", 4)
  names(res_l) <- c("NON_INFO", "INFO", "POWER", "MIX")

  ### Non-info prior
  res_l[[1]] <- analyse_trial(
    data_tr = data_i, 
    true_trt_eff = true_trt_eff,
    ni_margin = n_i_m,
    alpha = alpha,
    mod_stan = mod,
    soc_pri = soc_pri_non_info,
    trt_pri = trt_pri,
    cov_va0_pri = cov_va0_pri,
    sd_pri = sd_pri,
    power = F,
    args_cmdstan = cmdstan_args,
    seed = seed
  )
  
  ### Info prior
  res_l[[2]] <- analyse_trial(
    data_tr = data_i, 
    true_trt_eff = true_trt_eff,
    ni_margin = n_i_m,
    alpha = alpha,
    mod_stan = mod,
    soc_pri = soc_pri_info,
    trt_pri = trt_pri,
    cov_va0_pri = cov_va0_pri,
    sd_pri = sd_pri,
    power = F,
    args_cmdstan = cmdstan_args,
    seed = seed
  )
  
  
  ### Info power prior
  res_l[[3]] <- analyse_trial(
    data_tr = data_i, 
    true_trt_eff = true_trt_eff,
    ni_margin = n_i_m,
    alpha = alpha,
    mod_stan = mod,
    soc_pri = soc_pri_info,
    trt_pri = trt_pri,
    cov_va0_pri = cov_va0_pri,
    sd_pri = sd_pri,
    power = T,
    args_cmdstan = cmdstan_args,
    seed = seed
  )
  
  ### MAP prior
  res_l[[4]] <- analyse_trial(
    data_tr = data_i, 
    true_trt_eff = true_trt_eff,
    ni_margin = n_i_m,
    alpha = alpha,
    mod_stan = mod,
    soc_pri = soc_pri_mix,
    trt_pri = trt_pri,
    cov_va0_pri = cov_va0_pri,
    sd_pri = sd_pri,
    power = F,
    args_cmdstan = cmdstan_args,
    seed = seed
  )
  

  ### Output
  # Was convergence achieved
  # Append this to res and sens_df for easy results analysis
  conv_df <- lapply(res_l, function(x) x$CONV) %>% 
    bind_rows(.id = "PRIOR") 
  
  # MCMC results
  res_df <- lapply(res_l, function(x) x$RES) %>% 
    bind_rows(.id = "PRIOR") %>% 
    left_join(conv_df, .)
  
  # Sensitivity analysis dataframe
  sens_df <- lapply(res_l, function(x) x$SENS) %>% 
    bind_rows(.id = "PRIOR") %>% 
    left_join(conv_df, .)
 
  ### Return
  list(RES = res_df, SENS = sens_df)
}

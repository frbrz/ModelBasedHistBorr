## AUTHOR: Francesco Brizzi
## AIM: Set of functions to simulate individuals based on the PK/BCVA and PK/BCVA-TRIAL models.
##      simulate_new_inds simulates all the longitudinal profile.
##      simulate_trial reduces the longitudinal profile to a cross-sectional endpoint.

# Setup ------------------------------------------------------------------------

library(tidyverse)
library(lixoftConnectors)
library(RsSimulx)
library(checkmate)

# Launch simulx
initializeLixoftConnectors(
  software = "simulx",
  path = "C:/ProgramData/Lixoft/MonolixSuite2020R1/",
  force = T
)

# Function to simulate individuals from model ----------------------------------
# i.e. simulate BCVA and CVA given baseline characteristics 
# using either the one EMAX or the trial EMAX model
# (assuming q4w treatment)


#' Function to simulate the clinical enpoint over time
#' 
#' Simulate the clinical endpoints (BCVA and CVA),
#' given a set of baseline charactersitics
#' and parameters estimated from  PK/BCVA or trial-specific PK/BCVA model.
#'
#' Assumptions:
#' - Q4w 0.5mg dosing for all individuals
#' 
#' @param model character. Must be one of:
#'  - `one_emax`: PK/BCVA model
#'  - `trial_emax`": trial-specific PK/BCVA model
#' @param unc_pop_pars logical.
#'  If set equal to true the uncertainty in population parameters
#'  is taken into account in the simulations. Otherwise not.
#' @param base_df data-frame. Contains baseline characteristics and 
#' for pre-treated patients the number of pre-trial baseline injections
#' @param n_sim integerish. Number of simulations 
#' @param mon_tri integerish. 
#' Trial length in months (or visits which are monthly) 
#' @param mon_out vector of integerish.
#' If NULL (default) endpoint values at all months (visits) are outputed. 
#' Otherwise, for computational efficiency, only a subset of visits / months
#' can be specified
#' @param naive logical.
#' If TRUE no pre-treatment is assumed.
#' If FALSE individuals are assumed to be pre-treated with number of 
#' injections as specified in base_df$PRE_VIS
#' @param seed seed
#' @return A long data-frame with columns:
#' - `rep`: simulation number
#' - `id`: simulation ID
#' - `time`: time of observation
#' - `sim_VA_TRANS`: simulated BCVA on transformed scale
#' - `X_VA0`: baseline BCVA
#' - `COV_AGE`: age covariate
#' - `X_AGE`: age regressor (same as age-covariate)
#' - `CAT_STD`: study covariate for trial specific model
simulate_new_inds <- function(model,
                              unc_pop_pars = T,
                              base_df,
                              n_sim,
                              mon_tri = 12,
                              mon_out = NULL,
                              naive = F,
                              seed) {
  ### Arguments checks
  assert_character(model, len = 1)
  assert_choice(model, c("one_emax", "trial_emax"))
  assert_logical(unc_pop_pars, len = 1)
  assert_data_frame(base_df)
  must_nm <- c("id", "X_VA0", "COV_AGE", "PRE_VIS", "X_AGE")
  assert_names(colnames(base_df), must.include = must_nm)
  assert_integerish(n_sim, len = 1, lower = 1)
  assert_integerish(seed, len = 1, lower = 1)
  if (! length(unique(base_df$id)) == nrow(base_df)) {
    stop("Non unique IDs in base_df")
  }
  assert_integerish(mon_tri, len = 1)
  assert_integerish(mon_out, lower = 0, upper = mon_tri, null.ok = T)
  assert_logical(naive, len = 1)
  
  ### Monolix specific
  
  # Temporary directory for monolix stuff
  mlx_temp_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/Temp"
  # Create temporary folder if needed
  if (!dir.exists(mlx_temp_wd)) {
    dir.create(mlx_temp_wd)
  }
  
  # Read mlx project from simulx
  mlx_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes"
  mod_wd <- switch (
    model,
    "one_emax" = file.path(mlx_wd, "pkbcva_trans2"),
    "trial_emax" = file.path(mlx_wd, "pkbcva_trans2_trial"),
    stop("Unknown model")
  )
  
  ### Simulate uncertainty in population parameters
  # (must be done before launching simulx to avoid conflicts)
  if (unc_pop_pars) {
    file_pop <- file.path(mlx_temp_wd, "pop.csv")
    set.seed(seed)
    sim_pop_pars <- simpopmlx(
      n = n_sim, 
      project = paste0(mod_wd, ".mlxtran"),
      outputFilename = file_pop
    )
  }
  
  ### Simulx specific
  
  # Launch simulix
  importMonolixProject(paste0(mod_wd, ".mlxtran"))
  
  # Temporary directory where smlx stores stuff
  temp_wd <- dirname(getTreatmentElements()$mlx_Adm1$file)
  
  # Further argument checks on base_df and trt_df
  must_covs <- read.csv(
    file.path(temp_wd, "mlx_Cov.txt"), nrows = 1, header = T
  ) %>% 
    select(-id) %>% 
    colnames()
  must_reg <- read.csv(
    file.path(temp_wd, "mlx_Reg.txt"), nrows = 1, header = T
  ) %>% 
    select(-id, -time) %>% 
    colnames()
  must_cols <- c("id", must_covs, must_reg)
  if (any(must_cols == "CAT_STD")) {
    # Dirty trick, as CAT_STD is not supplied
    # when we use the trial model, so remove from must_covs
    must_cols <- must_cols[- which(must_cols == "CAT_STD")]
  }
  assert_names(colnames(base_df), must.include = must_cols)
  
  ### Population parameters 
  if (unc_pop_pars) {
    # Use file of simulated uncertainty for pop-parameters
    definePopulationElement("pop_par", element = file_pop)
    unlink(file_pop)
  } else{
    # use estimated population parameters
    pop_pars <- getPopulationElements()$mlx_Pop$data[-1]
    definePopulationElement("pop_par", element = pop_pars)
  }
  
  # If using the trial model, replicate base_df
  # so that n_id individuals are simulated for each trial
  if (model == "trial_emax") {
    n_id <- nrow(base_df)
    stds <- unique(getCovariateElements()$mlx_Cov$data$CAT_STD)
    base_df_tmp <- as.data.frame(sapply(base_df, rep.int, length(stds)))
    base_df_tmp$id <- 1:(n_id * length(stds))
    base_df_tmp$CAT_STD <- rep(stds, each = n_id)
    base_df <- base_df_tmp
  }
  
  ### Covariates 
  cov_df <- base_df %>% 
    select(id, all_of(must_covs))
  file_cov <- file.path(mlx_temp_wd, "cov.csv")
  write_csv(cov_df, file_cov)
  defineCovariateElement("cov", element = file_cov)
  unlink(file_cov)
  
  ### Treatment
  
  # Simulation time depends on wheter naive
  mon_max <- ifelse(naive, mon_tri, mon_tri + max(base_df$PRE_VIS) + 1)
  
  # Defining time of treatment (q4w IVT RBZ)
  t_trt <- seq(0, mon_max * 28, by = 28) + 0.042
  t_ids <- unique(cov_df$id)
  trt_df <- data.frame(
    id = rep(t_ids, each = length(t_trt)),
    time = rep(t_trt, length(t_ids)),
    amount = 500
  )
  file_trt <- file.path(mlx_temp_wd, "trt.csv")
  write_csv(trt_df, file_trt)
  defineTreatmentElement("trt", element = list(data = file_trt))
  unlink(file_trt)
  
  ### Regressor
  reg_df <- left_join(
    trt_df %>% select(-amount), 
    base_df %>% 
      select(id, all_of(must_reg)),
    by = "id"
  ) 
  file_reg <- file.path(mlx_temp_wd, "reg.csv")
  write_csv(reg_df, file_reg)
  defineRegressorElement("reg", element = file_reg)
  unlink(file_reg)
  
  ### Output
  # Maximum month depends on wheter naive population
  mon_max <- ifelse(naive, mon_tri, mon_tri + max(base_df$PRE_VIS) + 1)

  # Full output visit
  t_ids <- unique(cov_df$id)
  
  # Reducing output based on mon_vis
  if (naive) {
    # When population naive output does not depend on pre baselien visits
    # i.e. not id-dependent
    if (is_null(mon_out)) {
      t_obs <- seq(0, mon_max * 28, by = 28) 
    } else {
      t_obs <- mon_out * 28
    }
    obs_df <- data.frame(
      id = rep(t_ids, each = length(t_obs)),
      time = rep(t_obs, length(t_ids))
    )
  } else {
    # When population pre-loaded, output is dependent on id
    if (is.null(mon_out)) {
      t_obs <- seq2(
        from = (base_df$PRE_VIS + 1) * 28, 
        to   = (base_df$PRE_VIS + 1 + mon_tri) * 28, 
        by   = 28
      )
    } else {
      # NOTE: 0 is to always display baseline visit
      t_obs <- lapply(base_df$PRE_VIS + 1, function(x) (c(0, mon_out) + x) * 28) 
    }
    
    obs_df <- data.frame(
      id = rep(t_ids, each = length(t_obs[[1]])),
      time = unlist(t_obs)
    )
  }
  
  # Final observation data-frame
  out_df <- obs_df %>% 
    # do not want to simulate baseline as we know it
    filter(time != 0)
  file_out <- file.path(mlx_temp_wd, "out.csv")
  write_csv(out_df, file_out)
  defineOutputElement(
    name = "out",
    element = list(data = file_out, output = "Y")
  )
  unlink(file_out)
  
  ### Simulx simulation
  
  # Remove old simulation groups (otherwise screws up)
  grp_mlx <- sapply(getGroups(), function(x) x$name)
  lapply(grp_mlx, function(x) removeGroup(x))
  
  el <- c("pop_par", "trt", "reg", "cov", "out")
  
  # Seed (otherwise same residual error for each sim_par)
  setProjectSettings("seed" = seed)
  
  # Define groups needed for simulation
  addGroup("Sim")
  setGroupElement(group = "Sim", elements = el)
  setGroupSize(group = "Sim", size = n_distinct(out_df$id))
  setNbReplicates(n_sim)
  
  # Run simulation
  runSimulation()
  res_read <- getSimulationResults()$res$Y
  
  # Remove temporary directories used for simulation
  unlink(mlx_temp_wd, recursive = T)
  
  ### Output processing
  if (naive) {
    # Baseline output dataframe
    out_base_df <- tibble(
      id = rep(base_df$id, n_sim),
      rep = rep(1:n_sim, each = nrow(base_df)),
      time = 0,
      Y = rep(base_df$X_VA0, n_sim)
    )
    res <- bind_rows(out_base_df, res_read) %>% 
      # Merge baseline covariates
      left_join(base_df, by = "id") 
  } else {
    res <- res_read %>% 
      # Merge baseline covariates
      left_join(base_df, by = "id") %>% 
      # Re-define baseline, as first on trial-observation
      group_by(rep, id) %>%
      mutate(
        X_VA0 = Y[1],
        time = time - time[1]
      ) %>% 
      ungroup()
  }
  
  res <- res %>% 
    rename(sim_VA_TRANS = Y) %>% 
    # Endpoint on untransformed scale
    mutate(
      sim_VA = inv_trans_fun(sim_VA_TRANS),
      X_VA0 = inv_trans_fun(X_VA0),
      sim_CVA = sim_VA - X_VA0
    ) %>% 
    # Order in reasonable way
    arrange(rep, id, time) 
  
  ### Return
  res
}

# Function to simulate trial ---------------------------------------------------

#' Function to simulate a model-based two trial arms
#' This is generates a population of n_trt + n_soc patients based on the
#' simulate_new_inds function
#' 
#' @param n_trt number of patients in the new treatment arm.
#' @param n_soc number of patients in the SOC arm.
#' @param trt_eff treatment effect of the new drug.
#' @param hist_cont conflict between current and historical control arms.
#' @param ... arguments passed to function simulate_new_inds function.
#' @return A long data-frame ready for statistical analysis, with columns:
#' - `rep`: simulation number
#' - `id`: simulation ID
#' - `time`: time of observation
#' - `sim_CVA`: simulated change from baseline BCVA 
#' - `X_VA0`: centered (at 50) baseline BCVA 
simulate_trial <- function(n_trt, n_soc,trt_eff, hist_conf, base_df, ...) {
  
  ### Argument checks
  assert_integerish(n_trt, len = 1, lower = 1, upper = nrow(base_df))
  assert_integerish(n_soc, len = 1, lower = 1, upper = nrow(base_df))
  assert_numeric(trt_eff, len = 1)
  assert_numeric(hist_conf, len = 1)
  
  ### Jointly simulate SOC and TRT arms
  # (this is to speed up simulations)
  n_trial <- n_trt + n_soc
  
  if (nrow(base_df) != n_trial) {
    stop("base_df has different number of rows than n_trt + n_soc.")
  }

  ### Run simulation
  sim_ids <- simulate_new_inds(base_df = base_df, ...)
  
  ### Simulations cleaning
  res_tmp <- sim_ids %>% 
    # remove baseline
    filter(time > 0) %>% 
    # Only keep important cols for stat analysis
    # to reduce memory footprint
    select(rep, id, sim_CVA, X_VA0) %>% 
    # declare arm: 1 new trt, 0 SOC
    mutate(ARM = ifelse(id <= n_trt, 1L, 0L))
  
  # Monolix is not perfect, across simulation mean diff btwn arms shall be 0
  # but there is little empirical bias, which we must correct for
  emp_bias <- res_tmp %>%
    group_by(rep, ARM) %>%
    summarise(MU = mean(sim_CVA)) %>% 
    pivot_wider(names_from = ARM, values_from = MU) %>% 
    mutate(DIFF = `1` - `0`) %>%
    ungroup() %>%
    summarise(a = mean(DIFF)) %>% 
    pull(a)
  
  res <- res_tmp %>%
    # Adding historical conflict
    # (which shifts both SOC and new trt arms)
    mutate(sim_CVA = sim_CVA + hist_conf) %>% 
    # Adding treatment effect + empirical bias
    mutate(sim_CVA = ifelse(
      ARM == 1L, sim_CVA + trt_eff, sim_CVA + emp_bias
    )) %>% 
    # Centering VA0
    mutate(X_VA0 = X_VA0 - 50)
  
  ### Return
  res
}


# Helper function for quick frequentist trial analysis -------------------------
# Not planning to do those, just act as a safety check 
# as well as are needed to derive a prior of interest


#' Function to do a frequentist approximate power analysis
#' @param data simulated model-based data
#' @return average power for simulated data
run_pow_analysis <- function(data) {
  dec <- rep(NA, max(data$rep))
  
  for(sim in 1:max(data$rep)) {
    
    # Progress
    if (sim %% 100 == 0) {print (sim)}
    
    # Data
    dat <- data[data$rep == sim, ]
    
    # Run regression  
    reg <- lm("sim_CVA ~ 1 + ARM + X_VA0", data = dat)
    
    # Nice output via broom
    # Note: alpha is for one sided NI decision,
    # so need CI for 2*alpha
    reg_summ <- broom::tidy(reg, conf.int = T, conf.lev = 1 - 2 * 0.025)
    
    # Non inferiority decision
    # if lower bound of CI > - ni_margin
    dec[sim] <- reg_summ[reg_summ$term == "ARM", "conf.low", drop = T] > - 3.9
  }
  
  # Power (correct decision) averaged over simulations
  sum(dec) / length(dec)
}



#' Function to generate prior 
#' for n supplementary patients, based on model-based simulated data
#' @param data simulated data from NLME model
#' @param n size of patients in prior
#' @param seed for reproducibility
#' @return Distribution of SOC effect parameters, 
#' one for each simulated dataset
generate_prior <- function(data, n, seed = 123) {
  
  # Argument checks
  assert_data_frame(data)
  must_nm <- c("rep", "id", "sim_CVA", "X_VA0", "ARM")
  assert_names(colnames(data), must.include = must_nm)
  assert_integerish(n, len = 1, lower = 1, upper = n_distinct(data$id))
  
  if (! all(data$ARM == 0)) {
    stop("Prior shall only consider SOC patients")
  }
  
  n_sim <- max(data$rep)
  res <- rep(NA, n_sim)
  for(sim in 1:n_sim) {
    
    # Progress
    if (sim %% 100 == 0) {print (sim)}
    
    # Data
    # reducing full dataset to specific simulation
    dat <- data[data$rep == sim, ]
    
    # Subsampling IDs
    set.seed(seed)
    ids <- sample(dat$id, n, replace = F)
    dat <- dat[match(ids, dat$id), ]
    
    # Run regression  
    reg <- lm("sim_CVA ~ 1 + X_VA0", data = dat)
    
    # Nice output via broom
    # Note: alpha is for one sided NI decision,
    # so need CI for 2*alpha
    reg_summ <- broom::tidy(reg, conf.int = T, conf.lev = 1 - 2 * 0.025)
    
    # Non inferiority decision
    # if lower bound of CI > - ni_margin
    res[sim] <- reg_summ[reg_summ$term == "(Intercept)", "estimate", drop = T]
  }
  
  # Return
  res
} 

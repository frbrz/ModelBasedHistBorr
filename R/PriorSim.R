## AUTHOR: Francesco Brizzi
## AIM: Based on the baseline characteristics + dosing info (in MasterTrialSim.R)
##      simulate longitudinal BCVA + CBCVA profiles for individuals
##      in treatment-naive and pre-treated trials.

# Setup ------------------------------------------------------------------------
rm(list = ls())

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

# Useful functions -------------------------------------------------------------

source("W:/Users/brizzif/VCPaper/R/FunHelp.R")

# Base characteristics of individuals in new trial -----------------------------

# Defined in master_sim.R
base_df <- data.table::fread("S:/VCPaper/TrialSimData/base_df.csv")

# Global options for simulations -----------------------------------------------

n_id <- nrow(base_df)
n_sim <- 1000
sd <- 1807
mon_sim <- 16

# Working directories of mlx models --------------------------------------------

trt_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/pkbcva_trans2/"
std_emax_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/pkbcva_trans2_trial/"
save_wd <- "S:/VCPaper/TrialSimData"

# Function to make simulations for prior ---------------------------------------
# i.e. simulate BCVA and CVA given baseline characteristics 
# using either the one EMAX or the trial EMAX model
# (we assume q4w treatment by default)


#' Function to simulate the enpoint over time
#' 
#' Simulate the clinical endpoints (BCVA and CVA),
#' given a set of baseline charactersitics
#' and parameters estimated from  PK/BCVA or trial-specific PK/BCVA model.
#'
#' Assumptions:
#' - Q4w 0.5mg dosing for all individuals
#' 
#' @param model character: must be one of:
#'  - `one_emax`: PK/BCVA model
#'  - `trial_emax`": trial-specific PK/BCVA model
#' @param unc_pop_pars logical.
#'  If set equal to true the uncertainty in population parameters
#'  is taken into account in the simulations. Otherwise not.
#' @param base_df supplied data-frame of baseline characteristics
#' @param n_sim number of simulations to be done
#' @param mon_max number of visits (months) to be longitudinally simulated
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
                              mon_max = 12,
                              seed) {
  # Arguments checks
  assert_character(model, len = 1)
  assert_choice(model, c("one_emax", "trial_emax"))
  assert_logical(unc_pop_pars, len = 1)
  assert_data_frame(base_df)
  assert_names(colnames(base_df), must.include = c("id"))
  assert_integerish(n_sim, len = 1, lower = 1)
  assert_integerish(seed, len = 1, lower = 1)
  if (! length(unique(base_df$id)) == nrow(base_df)) {
    stop("Non unique IDs in base_df")
  }
  
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
  
  # Simulated uncertainty in population parameters
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
  
  # Launch simulx
  importMonolixProject(paste0(mod_wd, ".mlxtran"))
  
  # Temporary directory where smlx stores stuff
  temp_wd <- dirname(getTreatmentElements()$mlx_Adm1$file)
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
  
  # Further checks on base_df and trt_df
  assert_names(colnames(base_df), must.include = must_cols)
  
  # Population parameters 
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
  
  # Covariates 
  cov_df <- base_df %>% 
    select(id, all_of(must_covs))
  file_cov <- file.path(mlx_temp_wd, "cov.csv")
  write_csv(cov_df, file_cov)
  defineCovariateElement("cov", element = file_cov)
  unlink(file_cov)
  
  # Treatment
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
  
  # Regressor
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
  
  # Output
  out_df <- trt_df %>% 
    select(-amount) %>% 
    mutate(time = time - 0.042) %>% 
    # do not want to simulate baseline as we know it
    filter(time != 0)
  file_out <- file.path(mlx_temp_wd, "out.csv")
  write_csv(out_df, file_out)
  defineOutputElement(
    name = "out",
    element = list(data = file_out, output = "Y")
  )
  unlink(file_out)
  
  # Creating a baseline output dataframe
  out_base_df <- tibble(
    id = rep(base_df$id, n_sim),
    rep = rep(1:n_sim, each = nrow(base_df)),
    time = 0,
    Y = rep(base_df$X_VA0, n_sim)
  )
  
  # Remove old simulation groups (otherwise screws up)
  grp_mlx <- sapply(getGroups(), function(x) x$name)
  lapply(grp_mlx, function(x) removeGroup(x))
  
  el <- c("pop_par", "trt", "reg", "cov", "out")
  
  # Seed (otherwise same residual error for each sim_par)
  setProjectSettings("seed" = seed)
  addGroup("Sim")
  setGroupElement(group = "Sim", elements = el)
  setGroupSize(group = "Sim", size = n_distinct(out_df$id))
  setNbReplicates(n_sim)
  
  # Simulation
  tictoc::tic()
  runSimulation()
  tictoc::toc()  
  res_read <- getSimulationResults()$res$Y
  
  # Remove temporary directories and clean mem
  unlink(mlx_temp_wd, recursive = T)
#  gc()
  
  # Minimal data-cleaning
  res <- bind_rows(out_base_df, res_read) %>% 
    rename(sim_VA_TRANS = Y) %>% 
    left_join(base_df, by = "id") %>% 
    arrange(rep, id, time)
  
  # Return
  res
}

# Simulating 1000 ids, one-emax model ------------------------------------------


# Monolix is stupid and cannot deal with simulations that are too big...
# Therefore instead of doing all n_sim at once
# do n_sim / chunks simulations at once...

n_chunks <- 5
n_sim_ch <- n_sim / n_chunks
if (n_sim %% n_chunks != 0) {
  warning ("Not splitting simulations in equal size!!!")
  n_sim_ch <- floor(n_sim_ch)
}

one_emax_unc_l <- vector("list", n_chunks)
for (i in 1:n_chunks) {
  print(i)
  
  one_emax_unc_l[[i]] <- simulate_new_inds(
    model = "one_emax",
    unc_pop_pars = T,
    base_df = base_df,
    n_sim = n_sim_ch,
    mon_max = mon_sim,
    seed = sd + i - 1
  ) %>% 
    # correct simulation index
    mutate(rep = rep + n_sim_ch * (i - 1))
}
one_emax_unc <- bind_rows(one_emax_unc_l)

data.table::fwrite(
  one_emax_unc, file = file.path(save_wd, "pri_one_emax_unc.csv")
)

# Simulating 1000 ids, trial-emax model ----------------------------------------

# Splitting by ID is more efficient
n_chunks <- 20
n_ids_ch <- nrow(base_df) / n_chunks
n_std <- 6
if (nrow(base_df) %% n_chunks != 0) {
  warning ("Not splitting simulations in equal size!!!")
  n_ids_ch <- floor(n_ids_ch)
}

for (i in 1:n_chunks) {
  
  # breaking n_sim by 2 to avoid out of memory of monolix
  print(i)
  ids <- (1 + (i - 1) * n_ids_ch):(i * n_ids_ch)
  ids1 <- (1 + (i - 1) * n_ids_ch * n_std):(i * n_ids_ch  * n_std)
  ids_df <- data.frame(id = 1:(n_ids_ch * n_std), ID = ids1)
  
  save_txt <- paste0("pri_trial_emax_unc_ch", i, ".csv")
  save_str <- file.path(save_wd, save_txt)
  
  trial_emax_unc_tmp_1 <- simulate_new_inds(
    model = "trial_emax",
    unc_pop_pars = T,
    base_df = base_df[ids, ],
    n_sim = n_sim / 2,
    mon_max = mon_sim,
    seed = sd
  ) 
  
  a <- trial_emax_unc_tmp_1 %>% 
    left_join(ids_df, by = "id") %>% 
    select(-id) %>% 
    rename(id = ID) %>% 
    select(id, everything())
  
  rep_df <- data.frame(rep = 1:(n_sim/2), REP = ((n_sim/2 + 1):n_sim))
  
  trial_emax_unc_tmp_2 <- simulate_new_inds(
    model = "trial_emax",
    unc_pop_pars = T,
    base_df = base_df[ids, ],
    n_sim = n_sim / 2,
    mon_max = mon_sim,
    seed = sd
  )
  
  b <- trial_emax_unc_tmp_2 %>% 
    left_join(ids_df, by = "id") %>%
    select(-id) %>%
    rename(id = ID) %>%
    left_join(rep_df, by = "rep") %>% 
    select(-rep) %>% 
    rename(rep = REP) %>% 
    select(id, rep, everything())
  
  trial_emax_unc <- bind_rows(a, b)
  
  # Storing every single simulation chunk
  # bc Monolix often goes out of memory
  data.table::fwrite(trial_emax_unc, file = save_str)
}



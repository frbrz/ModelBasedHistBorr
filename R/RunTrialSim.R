## AUTHOR: Francesco Brizzi
## AIM: Run trial simulations
##      Analyses many datasets (one true treatment effect, one sample size, one historical conflict)
##      one at the times using different priors
##      (informative, weakly informative, power, mixture) priors
##      Notes:
##      * This file externally called on the cluster by the SubmitLSF.sh file.
##      * All datasets with CHUNK = ch (specified by array of SubmitLSF)
##        are analysed sequentially in this script 

# Setup ------------------------------------------------------------------------

rm(list = ls())

# Using fread and fwrite, without exporting whole data.table
fread <- data.table::fread
fwrite <- data.table::fwrite

# Other packages in sourced file with functions

# Cluster setup ----------------------------------------------------------------

on_cluster <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

# Grepping simulation index ----------------------------------------------------

if (on_cluster) {
  arg <- commandArgs(trailingOnly = TRUE)
  ch <- as.numeric(arg[1])
} else {
  ch <- 1
}

# Working directories and loading data -----------------------------------------

if (on_cluster) {
  base_wd <- "/home/brizzif/global/VCPaper"
  scratch_wd <- "/scratch/site/u/brizzif/VCPaper"
  cmdstan_wd <- "/home/brizzif/.cmdstan/cmdstan-2.27.0"
} else {
  base_wd <- "W:/Users/brizzif/VCPaper"
  scratch_wd <- "S:/VCPaper"
  cmdstan_wd <- "C:/Users/brizzif/Documents/.cmdstanr/cmdstan-2.27.0/"
}

stan_wd <- file.path(base_wd, "Stan")
dat_wd <- file.path(scratch_wd, "TrialSimData")
save_wd <- file.path(scratch_wd, "TrialSimRes1")

# Useful functions -------------------------------------------------------------

source(file.path(base_wd, "R", "FunTrialAnalysis.R"))

# Loading simulation master ----------------------------------------------------

master_sim_df <- fread(file.path(dat_wd, "sim_master_df.csv"))

# Focusing on specific simulation components -----------------------------------

# Subset of simulations to be run
sim_df <- master_sim_df %>% 
  filter(CHUNK == ch)

if (n_distinct(sim_df$DATA_SIM_IND) > 1) {
  stop("This code only work for a unique DATA_SIM_IND")
}
if (n_distinct(sim_df$DATA_IND) > 1) {
  stop("This code only work for a unique DATA_IND")
}
if (n_distinct(sim_df$N_PBO) > 1) {
  stop("This code only work for a unique N_PBO")
}
n_pbo <- unique(sim_df$N_PBO)
trt_eff_true <- unique(sim_df$TRT)

# Load data, based on DATA_SIM_IND ---------------------------------------------

data_txt <- paste0("Data_", unique(sim_df$DATA_IND), ".csv")
data_sc_read <- fread(file.path(dat_wd, data_txt)) %>% 
  # renaming for consistency with analyse_trial
  rename(Y = sim_CVA, ID = id) %>%
  # tibble easier format to handle
  as_tibble()

# In case number of placebo smaller than 90
# 1) use a subset of the data
# 2) recorrect for empirical bias
if (n_pbo < 90) {
  data_sc <- data_sc_read %>% 
    filter(ID %in% c(1:90, 91:(90 + n_pbo)))
  
  emp_bias <- data_sc %>% 
    group_by(rep, ARM) %>% 
    summarise(MU = mean(Y)) %>% 
    pivot_wider(names_from = ARM, values_from = MU) %>% 
    mutate(DIFF = `1` - `0`) %>% 
    pull(DIFF) %>% 
    mean() %>% 
    {. - trt_eff_true}
  
  data_sc <- data_sc %>% 
    mutate(Y = ifelse(ARM == 1, Y, Y + emp_bias))
} else {
  data_sc <- data_sc_read
}
rm(data_sc_read)

# Compile stan file for analysis -----------------------------------------------

# Safety measure, at times cmdstan path not found in R
set_cmdstan_path(cmdstan_wd)
mod <- cmdstan_model(stan_file = file.path(stan_wd, "trial_anal.stan"))

# Default simulation arguments -------------------------------------------------
# These are common to all simulations

cmdstan_args <- list(
  chains = 3,
  parallel_chains = 3,
  refresh = 500,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.99,
  max_treedepth = 12,
  max_try = 5,
  timeout = 80, # in seconds, avoids stuff getting stuck on cluster
  show_messages	= T
)

if (on_cluster) {
  # Reducing memory footprint of output files
  cmdstan_args$refresh <- 0
  cmdstan_args$show_messages <- F
}

# Trial design non inferiority margin and alpha, and endpoint timing
n_i_m <- 3.9
alpha <- 0.025

# Seed for simulations
seed <- 1825

# Define priors -----------------------------------------------------------------

# SOC priors
pri_mean <- 1.52
pri_sd_info <- 2.32 # Corresponds to 30 patients
pri_sd_non_info <- 15

pri_soc_non_info <- norm_mix(w = 1, m = pri_mean, sd = pri_sd_non_info)
pri_soc_info <- norm_mix(w = 1, m = pri_mean, sd = pri_sd_info)
pri_soc_mix <- norm_mix(
  w = c(0.7, 0.3),
  m = c(pri_mean, pri_mean),
  sd = c(pri_sd_info, pri_sd_non_info)
)

# Non informative priors for other parameters
trt_pri <- norm_mix(w = 1, m = 0, sd = 15)
cov_va0_pri <- norm_mix(w = 1, m = 0, sd = 15)
sd_pri <- norm_mix(w = 1, m = 0, sd = 15)

# Run MCMC via stan ------------------------------------------------------------

# indices of data to be run
data_ind <- unique(sim_df$DATA)

# Looping over specific dataset
res_l <- setNames(vector("list", length = length(data_ind)), data_ind)
for (i in seq_along(data_ind)) {
  
  tictoc::tic()
  
  print(paste("Simulation i", i, "out of", length(data_ind)))
  
  ### Data specific for this simulation
  data_tr <- data_sc[data_sc$rep == data_ind[i], ] 
  
  ### Running trial with different priors
  res_l[[i]] <- run_trial_i(
    data_i = data_tr, 
    true_trt_eff = trt_eff_true,
    ni_margin = n_i_m,
    alpha = alpha,
    mod_stan = mod,
    soc_pri_non_info = pri_soc_non_info,
    soc_pri_info = pri_soc_info,
    soc_pri_mix = pri_soc_mix,
    trt_pri = trt_pri,
    cov_va0_pri = cov_va0_pri,
    sd_pri = sd_pri,
    args_cmdstan = cmdstan_args,
    seed = seed
  )
  tictoc::toc()
}

# Clean up results -------------------------------------------------------------

# Results as dataframe / csv
res_mcmc <- bind_rows(map(res_l, ~ .$RES), .id = "DATA") %>% 
  mutate(DATA = as.integer(DATA)) %>% 
  left_join(sim_df, ., by = "DATA")
res_sens <- bind_rows(map(res_l, ~ .$SENS), .id = "DATA") %>% 
  mutate(DATA = as.integer(DATA)) %>% 
  left_join(sim_df, ., by = "DATA")

# Store ------------------------------------------------------------------------

save_txt_res <- paste0("ResChunk_", unique(sim_df$CHUNK), ".csv")
save_txt_sens <- paste0("SensChunk_", unique(sim_df$CHUNK), ".csv")
fwrite(res_mcmc, file.path(save_wd, save_txt_res))
fwrite(res_sens, file.path(save_wd, save_txt_sens))

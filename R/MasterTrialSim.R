## AUTHOR: Francesco Brizzi
## AIM: Define all simulation scenario to be investigated
##      1) define baseline characteristics and dosing information 
##      2) Define sample size, true treatment effect and number of replicate datasets
##      These info are use to generate simulated data in SimData.R
##      (in chunks to avoid memory issues)

# Setup ------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)

# Helpful function -------------------------------------------------------------

trans_fun <- function(x) {2 - log10(100 - x)}

# Simulating baseline characteristics ------------------------------------------

# Maximum number of individuals considered
n_id <- 1000

# Baseline characteristic assumptions
cov_sim_dist <- list(
  mean_age = 77,
  median_age = 78,
  sd_age = 8,
  mean_va0 = 55,
  sd_va0 = 12
)

# Distribution of pre-baseline characteristics
# when patients are not-naive
n_vis_pre <- c(4, 5, 6)
w_vis_pre <- c(1/2, 1/4, 1/4)

set.seed(1807)
va0 <- with(cov_sim_dist, rnorm(n_id, mean_va0, sd_va0))
age <- with(cov_sim_dist, rnorm(n_id, mean_age, sd_age) - median_age)
pre_vis <- sample(n_vis_pre, n_id, prob = w_vis_pre, replace = T)

# In data-frame format
base_df <- data.frame(
  id = rep(1:n_id), 
  X_VA0 = va0, 
  COV_AGE = age,
  PRE_VIS = pre_vis
) %>% 
  mutate(X_AGE = COV_AGE) %>% 
  mutate(X_VA0 = trans_fun(X_VA0))

# Definition for simulation ----------------------------------------------------

# n_pbo assumed 90 and 60
# based on power calculations with frequentist design, see SimData.R
# number of placebo individuals
n_pbo <- c(90, 60)
# number of simulations for each scenario
n_rep <- 2000

# chunk_size is the size of how many datasets
# should be run within a parallel simulation on the cluster
# FYI: running a dataset with all (5) priors and frequentist desings
#      takes ~ 7.5 sec on the cluster
#      so a chunk size of 125 results in a ~15 min run time
chunk_size <- 250

# Defining conflict
c_seq <- seq(-6, 6, by = 0.25)

# Defining treatment effect
trt_eff <- c(-3.9, 0)

# # Grid with all combo ---------------------------------------------------------

# All simulations that have to be done
sim_master_df <- expand_grid(
  N_PBO = n_pbo,
  TRT = trt_eff,
  CONFLICT = c_seq,
  DATA = 1:n_rep
) %>%
  # Needed to preserve that DATA_SIM_IND = 1, for n = 90
  mutate(N_PBO = factor(N_PBO, levels = c(90, 60))) %>% 
  # DATA_IND refers to unique combinations of TRT, CONFLICT
  # To speed up datasimulation process N_PBO < 90 is a subset of simulated data
  mutate(DATA_IND = as.numeric(interaction(TRT, CONFLICT))) %>% 
  # DATA_SIM_IND refers to unique combinations of TRT, CONFLICT and N_PBO
  # this is the minimum simulation scenario
  mutate(DATA_SIM_IND = as.numeric(interaction(TRT, CONFLICT, N_PBO))) %>%
  arrange(DATA_SIM_IND) %>%
  # Within each DATA_SIM_IND (i.e. simulation scenario)
  # create a chunk index breaking the 1:n_rep data
  # each (except the last) of size chunk_size
  # CHUNK is a slurm simulation index
  # so if slurm gets CHUNK = 12
  # all simulations with CHUNK = 12 will be run serially (not in parallel)
  # whereas different chunks will be run in parallel
  # (ideally chunk_size = 1 but this exceeds cluster quota)
  group_by(DATA_SIM_IND) %>%
  mutate(CHUNK_MAX = ceiling(n() / chunk_size)) %>%
  mutate(CHUNK_INIT = CHUNK_MAX * (DATA_SIM_IND - 1)) %>%
  mutate(CHUNK_GRP = ((DATA - 1) %/% chunk_size) + 1) %>%
  mutate(CHUNK = CHUNK_GRP + CHUNK_INIT) %>%
  ungroup() %>%
  select(-CHUNK_MAX, -CHUNK_GRP, -CHUNK_INIT)

# Safety check
# CHUNK is a continuous sequence from 1:max(CHUNK)
if (! all(1:max(sim_master_df$CHUNK) %in% sim_master_df$CHUNK)) {
  stop ("CHUNK is NOT a continuous sequence")
}

# Storing ----------------------------------------------------------------------

store_wd <- "S:/VCPaper/TrialSimData/"

data.table::fwrite(base_df, file.path(store_wd, "base_df.csv"))
data.table::fwrite(sim_master_df, file.path(store_wd, "sim_master_df.csv"))

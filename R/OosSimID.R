## AUTHOR: Francesco Brizzi
## AIM: Run out of sample model simulations (i.e. prior-predictive distributions)
##      predicting the trial that was left out from the estimation procedure
##      based on the PK/BCVA and PK/BCVA-trial models.

# Setup ------------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(lixoftConnectors)
library(checkmate)

# Simulx
initializeLixoftConnectors(
  software = "simulx",
  path = "C:/ProgramData/Lixoft/MonolixSuite2020R1/",
  force = T
)

# Monolix results for full (not-cv) model
trt_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/pkbcva_trans2"
# Create a temporary directory to store intermediate file / res
mlx_temp_wd <- "W:/Users/brizzif/MIDD_VA2DF/Monolix/ModelRes/Temp"
# Create temporary folder if needed
if (! dir.exists(mlx_temp_wd)) {
  dir.create(mlx_temp_wd)
}


# Helpful plotting functions
source("W:/Users/brizzif/MIDD_VA2DF/R/fun_help.R")

# Full dataset
full_data <- read_data_project(trt_wd, ret = "data")

# Saving options ---------------------------------------------------------------
save_wd <- "S:/VCPaper/OOSSim" 
save_flag <- T

# Function to simulate one individual ------------------------------------------

sim_one_ind <- function(id,
                        full_data,
                        pop_pars,
                        rep_by_id) {
  
  # Argument checks to be implemented
  tictoc::tic()
  
  # Data for only the specific ID
  data_id <- full_data %>% 
    filter(ID == id) %>% 
    # Do not care here about what happens after v9
    # saving computing time!
    filter(VIS < 10)
  
  # In which repetitions (sim) is the ID involved
  reps <- rep_by_id[[which(names(rep_by_id) == id)]] %>% 
    as.numeric() %>% 
    sort()
  n_reps <- length(reps)

  # Define population parameters
  # pop pars for the specific reps in which ID involved
  pop_p <- pop_pars_df[reps, ] %>% 
    # STD and pop clash with formatting of Mlx save
    select(-STD, -pop) 
  file_pop <- file.path(mlx_temp_wd, "pop_pars.csv")
  write_csv(pop_p, file_pop)
  definePopulationElement("pop_par", element = file_pop)
  # unlink(file_pop)
  
  # Define regressor (via temporary file)
  base_reg <- data_id %>% 
    select(ID, TIME, starts_with("X_")) %>% 
    rename(id = ID, time = TIME) %>% 
    as.data.frame()
  file_reg <- file.path(mlx_temp_wd, "reg.csv")
  write_csv(base_reg, file_reg)
  defineRegressorElement("reg", element = file_reg)
  
  # Define treatment
  tr_df <- data_id %>% 
    filter(EVID == 1) %>% 
    select(ID, TIME, AMT) %>% 
    rename(id = ID, time = TIME, amount = AMT) %>% 
    as.data.frame()
  file_trt <- file.path(mlx_temp_wd, "trt.csv")
  write_csv(tr_df, file_trt)
  defineTreatmentElement(
    "trt", 
    element = list(admID = 1, probaMissDose = 0, data = file_trt)
  )
  
  # Output - only at visit 9
  out_df <- data_id %>% 
    filter(EVID == 0, VIS == 9) %>% 
    # Do not want to simulate 0 as baseline
    # filter(TIME > 0) %>% 
    select(ID, TIME) %>% 
    rename(id = ID, time = TIME) %>% 
    as.data.frame()
  file_out <- file.path(mlx_wd, "out.csv")
  write_csv(out_df, file_out)
  defineOutputElement(
    name = "out",
    element = list(data = file_out, output = "Y")
  )
  
  # Define covariates
  all_covs <- getCovariateElements()$mlx_Cov$data %>% 
    select(-id) %>% 
    colnames()
  cov_df <- data_id %>% 
    filter(TIME == 0) %>% 
    select(ID, all_of(all_covs)) %>% 
    as.data.frame()
  file_cov <- file.path(mlx_wd, "cov.csv")
  write_csv(cov_df, file_cov)
  defineCovariateElement("cov", element = file_cov)
  
  # Simulation groups
  # Remove old simulation groups (safety measure)
  grp_mlx <- sapply(getGroups(), function(x) x$name)
  lapply(grp_mlx, function(x) removeGroup(x))
  
  addGroup("Sim")
  setGroupElement(
    group = "Sim", elements =  c("pop_par", "reg", "cov", "out", "trt")
  )
  setGroupSize(group = "Sim", size = 1)
  setNbReplicates(n_reps)
  
  # Run simulation
  runSimulation()
  gc()
  tictoc::toc()
  
  # Output
  # (consistent with one of OOS sim)
  res_smx <- getSimulationResults()$res$Y
  res_smx$rep <- reps
  res_smx$ID <- id
  
  res_smx <- res_smx %>% 
    rename(sim_VA_TRANS = Y) %>% 
    rename(TIME = time) %>% 
    select(rep, ID, TIME, sim_VA_TRANS) %>% 
    mutate(STD = "POOL") %>% 
    as_tibble()
}


# Global options for simulations ---------------------------------------

# Number of simulations
n_sim <- 500
# (Max) number of IDS (N_ess) for simulations
n_ess <- 683

# Simulate params with uncertainty + boot covs for each study --------

all_std <- data.frame(
  SHORT = c("ANC", "AVE", "HAR", "MAR", "STA"),
  LONG = c("ANCHOR", "AVENUE", "HARBOR", "MARINA", "STAIRWAY")
) 

ids_boot_l <- list()
sim_pop_pars_l <- list()
for(i in seq_along(all_std$SHORT)) {
  print(i)
  std <- all_std$SHORT[i]
  str_mod <- paste0("pkbcva_", "trans2", "_cv")
  
  # Monolix and temporary monolix folder
  res_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/" 
  mlx_wd <- file.path(res_wd, paste0(str_mod, "_", std))
  
  # Simulate population params with uncertainty for each study
  pop_txt <- paste0("pop_", std)
  file_pop <- file.path(mlx_temp_wd, paste0(pop_txt, ".csv"))
  
  sim_pop_pars_l[[i]] <- RsSimulx::simpopmlx(
    n = n_sim, 
    project = paste0(mlx_wd, ".mlxtran"),
    outputFilename = file_pop
  ) %>% 
    mutate(STD = std)
  
  # IDs for specific study
  ids <- full_data %>% 
    # Ensure that we have readout at time 9
    filter(VIS == 9) %>% 
    filter(str_detect(CAT_STD, std)) %>% 
    filter(CAT_TRT == "0.5MG_q4w") %>% 
    pull(ID) %>% 
    unique()
  
  # Bootstrap ids of interest for n_ess
  ids_boot_l[[i]] <- replicate(
    n_sim,
    sample(ids, n_ess, replace = T),
    simplify = F
  ) %>% 
    set_names(1:n_sim) %>% 
    bind_rows() %>% 
    pivot_longer(
      cols = everything(),
      values_to = "ID",
      names_to = "rep"
    ) %>% 
    arrange(rep) %>% 
    mutate(STD = std) 
}

# IDs included in each bootstrap
ids_boot_df <- bind_rows(ids_boot_l)
# All population pars with unc by study
pop_pars_df <- bind_rows(sim_pop_pars_l)

# List containing repetition index for each id
rep_by_id_l <- list()
for(i in 1:n_distinct(ids_boot_df$ID)) {
  id <- unique(ids_boot_df$ID)[i]
  rep_by_id_l[[i]] <- ids_boot_df[ids_boot_df$ID == id, "rep", drop = T]
} 
names(rep_by_id_l) <- unique(ids_boot_df$ID)

# Running simulations -------------------------------------------------------

### LOADING SAME (non CV proj for all)
# fine as same struct mod...
# and passing all pars by hand
importMonolixProject(paste0(trt_wd, ".mlxtran"))


for (s in 1:length(all_std$SHORT)) {
  
  ids_std <- ids_boot_df %>% 
    filter(STD == all_std$SHORT[s]) %>% 
    pull(ID) %>% 
    unique()

  # If list not reinitialised, copy from prev-trial  
  res_l <- list()
  for (i in 1:length(ids_std)) {
    print(paste("std", s, "/", 5, "id:", i, "/", length(ids_std)))
  
    res_l[[i]] <- sim_one_ind(
      id = ids_std[i],
      full_data = full_data,
      pop_pars = pop_pars_df,
      rep_by_id = rep_by_id_l
    ) 
  }
  
  res <- bind_rows(res_l) 
  if (save_flag) {
    save_txt <- paste0("OOS_PRI_ID_", all_std$SHORT[s], ".csv")
    data.table::fwrite(res, file.path(save_wd, save_txt))
  }
}

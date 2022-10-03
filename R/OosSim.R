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
trt_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/pkbcva_trans2/"

# Helpful plotting functions
source("W:/Users/brizzif/MIDD_VA2DF/R/fun_help.R")

# Full dataset
full_data <- read_data_project(trt_wd, ret = "data")

# Saving options ---------------------------------------------------------------
save_wd <- "S:/VCPaper/OOSSim" 
save_flag <- T

# Function for out of sample predictions ---------------------------------------

#' Function to simulate out of sample studies
#' 
#' Based on the dosing regimen and baseline characteristics, 
#' longitudinal endpoint outcomes are simulated.
#' Based on the name of a out of sample study (i.e. `std`) 
#' and the name of a monolix model in a certain location (i.e. `model`)
#' this function carries on simulation based on simulx
#' 
#' @param std string. The study to be predicted out of sample.
#' Can be one of c("ANC", "AVE", "HAR", "MAR", "STA")
#' @param model string. The monolix name of the model. Note the model must be
#' in: "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/" 
#' @param type string. One of:
#' - `sim`: VPC type of simulation (with res error + between subj var)
#' - `pred`: expectation based on covariates (no res error + no BSV)
#' @param by_trial logical. 
#' - If `TRUE` then the trial-specific model is used for simulation, so that 
#'   the same individual is predicted for n_trial based on trial-specific
#'   vision asymptotyc steady states estimated via the trial-specific model.
#' - If `FALSE` the pooled model is used
#' @param full_data dataframe. The full data including all studies
#' (i.e. also the oos study)
#' @param n_sim integerish. The number of simulations
#' @param unc_pop_pars logical. 
#' If `TRUE` uncertainty in the population parameters is included in simulation
#' @param seed integerish. Seed
#' @return A data-frame made of
#' - `IO`: the ID of the out of sample study
#' - `rep`: the simulation replicate (only for type == `sim`)
#' - `sim_VA_trans`: the simulated transformed VA score 
#' - `STD`: either trial-specific membership for the study-specific model
#'    or simply equal to `pool` for the pooled model.
simulate_oos_std <- function(std,
                             model = "trans2", 
                             type = "sim", 
                             by_trial = F,
                             n_sim = 500,
                             full_data = full_data,
                             unc_pop_pars = T,
                             seed = 1825) {
  
  # Argument checks
  assert_character(std, len = 1)
  assert_character(model, len = 1)
  assert_character(type, len = 1)
  assert_logical(by_trial, len = 1)
  assert_integerish(n_sim, len = 1, lower = 1)
  assert_choice(type, c("sim", "pred"))
  if (type == "pred" & (n_sim != 1)) {
    txt <- paste(
      "With type = `pred`, it does not make sense to run more than one sim. \n",
      "n_sim is internally set to 1"
    )
    warning(txt)
    n_sim <- 1
  }
  
  if (type == "pred" & unc_pop_pars == T) {
    txt <- paste(
      "type = pred only makes sense without uncertainty in pop pars",
      "unc_pop_pars internally set to F"
    )
    unc_pop_pars <- F
  }
  
  # Working directory from where to access results
  res_wd <- "W:/Users/brizzif/VCPaper/Monolix/ModelsRes/" 
  trial_txt <- ifelse(by_trial, "_trial", "")
  str_mod <- paste0("pkbcva_", model, trial_txt, "_cv")
  all_files <- list.files(res_wd)
  cv_files <- all_files[grepl(".mlxtran", all_files)]
  cv_files <- cv_files[grepl(str_mod, cv_files)]
  # Checking that supplied std has a corresponding model
  ok_std_tmp <- strsplit(str_replace_all(cv_files, ".mlxtran", ""), "_")
  ok_std <- sapply(ok_std_tmp, function(x) x[length(x)])
  assert_choice(std, ok_std)
  
  # Monolix and temporary monolix folder
  mlx_wd <- file.path(res_wd, paste0(str_mod, "_", std))
  mlx_temp_wd <- "W:/Users/brizzif/MIDD_VA2DF/Monolix/ModelRes/Temp"
  
  # Create temporary folder if needed
  if (! dir.exists(mlx_temp_wd)) {
    dir.create(mlx_temp_wd)
  }
  
  # Getting data for the study out of sample
  std_oos <- tail(strsplit(basename(mlx_wd), "_")[[1]], 1)
  data_oos <- full_data %>% 
    filter(str_detect(CAT_STD, std_oos))
  
  # Monolix is stupid and it is easier to work with 1:n IDs
  orig_ids <- data_oos$ID
  id_dic_base <- data.frame(
    ID = unique(data_oos$ID),
    id = 1:n_distinct(unique(data_oos$ID))
  )
  data_oos <- data_oos %>% 
    left_join(id_dic_base, by = "ID") %>% 
    select(-ID) %>% 
    rename(ID = id) %>% 
    select(ID, everything())
  
  # If by_trial = T 
  # Augment oos data to create n_oos patinets
  # for each trial - 1 (the left out one)
  if (by_trial) {
    std_all <- unique(full_data$CAT_STD)
    std_aug <- std_all[!str_detect(std_all, std_oos)]
    l <- length(std_aug)
    std_aug_num <- as.numeric(factor(std_aug, levels = std_all))
    n_r_orig <- nrow(data_oos)
    
    id_dic <- data.frame(
      id = rep(data_oos$ID, l),
      id_orig = rep(orig_ids, l),
      CAT_STD = rep(std_aug, each = n_r_orig),
      CAT_STD_NUM = rep(std_aug_num, each = n_r_orig)
    ) %>% 
      mutate(ID = as.numeric(interaction(id, CAT_STD_NUM))) %>% 
      select(-CAT_STD_NUM)
    
    data_oos <- bind_rows(replicate(l, data_oos, simplify = FALSE))
    data_oos$CAT_STD <- id_dic$CAT_STD
    data_oos$ID <- id_dic$ID
  }
  
  # Simulating uncertainty in population parameters
  # (must be done before launching simulx to avoid conflicts)
  if (unc_pop_pars) {
    file_pop <- file.path(mlx_temp_wd, "pop.csv")
    set.seed(seed)
    sim_pop_pars <- RsSimulx::simpopmlx(
      n = n_sim, 
      project = paste0(mlx_wd, ".mlxtran"),
      outputFilename = file_pop
    )
  }
  
  ### Simulation stuff 
  # Launch simulx
  importMonolixProject(paste0(mlx_wd, ".mlxtran"))
  
  # Simulation output
  out_str <- switch(type, "sim" = "Y", "pred" = "VA")
  out_nm <-  switch(type, "sim" = "sim_VA_TRANS", "pred" = "VA_TRANS_PRED")
  
  # Population parameters 
  if (unc_pop_pars) {
    # Use file of simulated uncertainty for pop-parameters
    definePopulationElement("pop_par", element = file_pop)
    unlink(file_pop)
  } else{
    # use estimated population parameters
    pop_pars <- getPopulationElements()$mlx_Pop$data[-1]
    if (type == "pred") {
      # No random effects when type == "pred"
      pop_pars[1, grep("omega", colnames(pop_pars))] <- 0
    }
    definePopulationElement("pop_par", element = pop_pars)
  }
  
  # Define regressor (via temporary file)
  base_reg <- data_oos %>% 
    select(ID, TIME, starts_with("X_")) %>% 
    rename(id = ID, time = TIME) %>% 
    as.data.frame()
  file_reg <- file.path(mlx_temp_wd, "reg.csv")
  write_csv(base_reg, file_reg)
  defineRegressorElement("reg", element = file_reg)
  
  # Define covariates
  cov_data <- getCovariateElements()$mlx_Cov$data
  
  if (!is.null(cov_data)) {
    cov_inc <- TRUE
    all_covs <- getCovariateElements()$mlx_Cov$data %>% 
      select(-id) %>% 
      colnames()
    cov_df <- data_oos %>% 
      group_by(ID) %>% 
      slice(1) %>% 
      ungroup() %>% 
      select(ID, all_of(all_covs)) %>% 
      rename(id = ID) %>% 
      as.data.frame()
    file_cov <- file.path(mlx_temp_wd, "cov.csv")
    write_csv(cov_df, file_cov)
    defineCovariateElement("cov", element = file_cov)
  } else {
    cov_inc <- FALSE
  }
  
  # Define treatment
  tr_df <- data_oos %>% 
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
  
  # Output
  out_df <- data_oos %>% 
    filter(EVID == 0) %>% 
    # Do not want to simulate 0 as baseline
    # filter(TIME > 0) %>% 
    select(ID, TIME) %>% 
    rename(id = ID, time = TIME) %>% 
    as.data.frame()
  file_out <- file.path(mlx_wd, "out.csv")
  write_csv(out_df, file_out)
  defineOutputElement(
    name = "out",
    element = list(data = file_out, output = out_str)
  )
  
  # Baseline data-frame
  # to be collated after simulations, as baseline no simulated
  base_df_tmp <- data_oos %>% 
    filter(EVID == 0) %>% 
    # Do not want to simulate 0 as baseline
    filter(TIME == 0)
  ids <- base_df_tmp$ID
  
  base_df <- data.frame(
    id = rep(ids, n_sim),
    rep = rep(1:n_sim, each = length(ids)),
    TIME = 0,
    A = rep(base_df_tmp$X_VA0, n_sim)
  ) %>% 
    rename(!!out_nm := A) 
  if (type == "pred") {
    # if type = "pred" there is no "rep" in output... remove
    base_df <- base_df %>% 
      select(-rep)
  }
  
  # Cleanup
  # Remove temporary dir and files
  unlink(file_reg)
  unlink(file_trt)
  unlink(file_out)
  unlink(mlx_temp_wd, recursive = T)
  # Remove old groups
  grp_mlx <- sapply(getGroups(), function(x) x$name)
  lapply(grp_mlx, function(x) removeGroup(x))
  
  # Simulation elements
  els <- c("pop_par", "trt", "reg", "out")
  if (cov_inc) {
    els <- c(els, "cov")
    unlink(file_cov)
  }
  
  # Actual simulation
  addGroup("Sim")
  setGroupElement(group = "Sim", elements = els)
  setGroupSize(group = "Sim", size = n_distinct(data_oos$ID))
  setNbReplicates(n_sim)
  
  # Run simulation
  setProjectSettings(seed = seed)
  tictoc::tic()
  runSimulation()
  tictoc::toc()
  
  # Collect results
  res_smx_tmp <- getSimulationResults()$res
  res_smx_tmp <- res_smx_tmp[[out_str]]
  
  # Minimal cleanup of results
  # for efficiency most data cleaning after reading results
  # so that no need to read in large files
  res_smx <- res_smx_tmp %>% 
    # remove TIME == 0 simulations
    # because baseline is known
    # (could have done this defining out, but screws up IDs
    # in case individuals with only baseline)
    filter(time > 0) %>% 
    rename(!!out_nm := all_of(out_str)) %>% 
    rename(TIME = time) %>% 
    # re-add baseline = 0 simulations
    bind_rows(base_df)
  
  if (type == "pred") {
    res_smx <- res_smx %>% 
      arrange(id, TIME)
  } else {
    res_smx <- res_smx %>% 
      arrange(rep, id, TIME)
  }
  
  # Convert stupid default simulx IDs back to original ones
  if (by_trial) {
    id_dic1 <- id_dic %>% 
      select(-id) %>% 
      rename(id = ID) %>% 
      unique()
    
    res_smx <- res_smx %>% 
      left_join(id_dic1, by = "id") %>% 
      select(-id) %>%
      rename(ID = id_orig) %>% 
      rename(STD = CAT_STD)
    
  } else {
    res_smx <- res_smx %>% 
      left_join(id_dic_base, by = "id") %>% 
      select(-id) %>% 
      mutate(STD = "POOL")
  }
  
  if (type == "sim") {
    res_smx <- res_smx %>% 
      select(rep, ID, TIME, everything())
  } else {
    res_smx <- res_smx %>% 
      select(ID, TIME, everything())
  }
  
  # Return
  res_smx
}

# Running out of sample predictions --------------------------------------------

# Bc of bug in lixoftConnectors: with n > 200 doesnot work for HAR and MAR
# goes out of memory
# n_sim <- if_else(std %in% c("HAR", "MAR"), 200, 600)
# therefore sample at most n = 200, by using a "chunk" strategy

all_std <- c("ANC", "AVE", "HAR", "MAR", "STA")
n_sim <- 800
n_ch <- 4
n_sim_ch <- floor(n_sim / n_ch)

# Running for pooled model
for (s in 1:length(all_std)) {
  # Selecting one out of sample study 
  print(all_std[s])
  std <- all_std[s]
  
  # Storing texts
  vpc_trans_txt <- paste0("OOS_TRANS_VPC_", std)
  pred_trans_txt <- paste0("OOS_TRANS_PRED_", std, ".csv")
  
  # Out of sample predictions (VPC + PRED) for trans model
  print("VPC_trans")
  
  vpc_tmp <- list()
  
  for (i in 1:n_ch) {
    ch_txt <- paste0("_CH", i, ".csv")
    print(paste("Chunk", i))
    
    vpc_tmp <- simulate_oos_std(
      std = std,
      model = "trans2", 
      type = "sim", 
      n_sim = n_sim_ch,
      full_data = full_data,
      unc_pop_pars = T,
      by_trial = F,
      seed = 1825 + i - 1 
    )  %>% 
      # correct simulation index
      mutate(rep = rep + n_sim_ch * (i - 1))
    
    if (save_flag) {
      save_txt <- paste0(vpc_trans_txt, ch_txt)
      data.table::fwrite(vpc_tmp, file.path(save_wd, save_txt))
    }
  }
  
  print("PRED_trans")
  pred_trans <- simulate_oos_std(
    std = std,
    model = "trans2", 
    type = "pred", 
    n_sim = 1,
    full_data = full_data,
    unc_pop_pars = F,
    by_trial = F,
    seed = 1825
  )
  
  if (save_flag) {
    data.table::fwrite(pred_trans, file.path(save_wd, pred_trans_txt))
  }
}


# Running for pooled model
for (s in 1:length(all_std)) {
  # Selecting one out of sample study 
  print(all_std[s])
  std <- all_std[s]

  # Storing texts
  vpc_trans_txt <- paste0("OOS_TRANS_TRIAL_VPC_", std)

  # Selecting number of chunks (bigger trials have more chunks)
  n_ch <- switch(
    std,
    "MAR" = 20,
    "HAR" = 20,
    "AVE" = 10,
    "ANC" = 10,
    "STA" = 4
  )
  n_sim_ch <- floor(n_sim / n_ch)
  
  # Out of sample predictions (VPC + PRED) for trans model
  print("VPC_trans")
  
  vpc_l <- list()
  
  for (i in 1:n_ch) {
    print(paste("Chunk", i))
    
    # Storing text
    ch_txt <- paste0("_CH", i, ".csv")
    
    # Running simulations
    vpc_tmp <- simulate_oos_std(
      std = std,
      model = "trans2", 
      type = "sim", 
      n_sim = n_sim_ch,
      full_data = full_data,
      unc_pop_pars = T,
      by_trial = T,
      seed = 1825 + i - 1 
    )  %>% 
      # correct simulation index
      mutate(rep = rep + n_sim_ch * (i - 1))
    
    if (save_flag) {
      save_txt <- paste0(vpc_trans_txt, ch_txt)
      data.table::fwrite(vpc_tmp, file.path(save_wd, save_txt))
    }
  }
}
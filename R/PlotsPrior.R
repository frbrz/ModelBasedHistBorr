## AUTHOR: Francesco Brizzi
## AIM: Creating plots for model-informed and other priors 
##      in the context of treatment-naive and pre-treated trial populations.
##      Yields Figure 4 (Naive + PreTreated) prior comparison in manuscript.

# Setup ------------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(checkmate)
library(ggdist)
library(data.table)
library(viridis)

# Saving outputs
save_flag <- F

# Cluster setup ----------------------------------------------------------------

on_cluster <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

# Working directories and loading data -----------------------------------------

if (on_cluster) {
  base_wd <- "~/global/VCPaper"
  scratch_wd <- "~/scratch/VCPaper"
  save_wd <- "/pstore/scratch/u/brizzif/VCPaper/Figures"
  full_data <- read_csv("~/MIDD_VA2DF/Data/reg_trans_data.csv") 
} else {
  base_wd <- "W:/Users/brizzif/VCPaper"
  scratch_wd <- "S:/VCPaper"
  save_wd <- "S:/VCPaper/Figures"
  proj_wd <- getwd()
  trt_wd <- file.path(proj_wd, "Monolix/ModelsRes/pkbcva_trans2") 
  full_data <- read_data_project(trt_wd, ret = "data")
}
load_wd <- file.path(scratch_wd, "TrialSimData")
oos_wd <- file.path(scratch_wd, "OOSSim")
fig_wd <- file.path(scratch_wd, "Figures")

# Useful functions -------------------------------------------------------------

source(file.path(base_wd, "R/FunHelp.R"))

# Function to check simulation output ------------------------------------------

#' Function to check simulation data obtained via simulx
#' 
#' Idea: checking data is computationally expensive due to their large size
#' so it is better to only do it once, when data are loaded
#' 
#' @param data is a simulation data-frame, obtained via the `tidy_res`function
#' @return data as input, but with extra (dataOK) class
#' if all safety checks have been satisfied
check_sim_data <- function(data) {
  # Safety checks on data
  must_col <- c("ID", "rep", "TIME", "VIS", "X_VA0", "PRE_VIS", "SIM_VA", "SIM_CVA")
  assert_data_frame(data)
  assert_names(colnames(data), must.include = must_col)
  
  # Use data.table for efficiency reasons
  if (!is.data.table(data)) {
    data <- data.table(data)
  }
  
  by_trial <- ifelse("CAT_STD" %in% colnames(data), T, F)
  
  if (by_trial) {
    n_std <- length(unique(data$CAT_STD))
    n_ids <- length(unique(data$ID))
    
    n_ids_tr <- data[, length(unique(ID)), by = CAT_STD]
    n_id <- n_ids_tr$V1[1]
    
    if (length(unique(n_ids_tr$V1)) != 1) {
      warning("The supplied data have unequal IDs per study")
    }
    if (n_ids != (n_id * n_std)) {
      stop("Expect n_std * n_ids IDs in this function... STOP")
    }
  } else {
    n_id <- length(unique(data$ID))
  }
  
  # Checking if IDs are consistent across reps
  id_un <- data[, length(unique(ID)), by = rep]
  if (length(unique(id_un$V1)) != 1) {
    stop("Non consisten IDs across rep")
  }
  
  # Append dataOK class, meaning checks were succesfull
  class(data) <- append("dataOK", class(data))
  
  # Return 
  data
}

# Function to modify simulation output -----------------------------------------

#' Function to reshape simulations obtained via simulx
#' 
#' The following ordered is implicit in the function:
#' - 1) First a subset of ID is selected if needed
#' (via id_sub or n_sub)
#' - 2) The population is made treatment naive if needed
#' - 3) Visits are averaged if requested
#' - 4) Only visits of interest are returned
#' 
#' Because of the order vis_avg and vis_sub refer to
#' the visits post-trial baseline (rather than from diagnosis)
#' for pre-treated population
#' 
#' @param data is a data-frame, obtained via the `tidy_res`function
#' @param n_sub specifying a subset of n_ids, rather than using all ids,
#' by using first `n_sub` ids (for each study).
#' @param id_sub specifying a subset of ids, rather than using all ids,
#' by passing the specific `ID`s to be subset
#' @param make_pre_trt (logical).
#' If FALSE considers that patients are treatment naive
#' If TRUE considers that patients have received monthly pre-dosing
#' for a number of month specified in the `PRE_VIS` column of `data`
#' @param vis_avg averages BCVA of specified `vis_avg` visits.
#' Useful in case BCVA averaged in statistical endpoint
#' @param vis_sub subset of visits of interest,
#' rather than outputting entire data-frame
#' @param center_va0 (numeric, between 1 and 100) 
#' If not NULL create a column X_VA0_CENTER where X_VA0 is centered
#' by the supplied center_va0 value
#' @return An updated `data` object
sim_data_reshape <- function(data,
                             n_sub = NULL,
                             id_sub = NULL,
                             make_pre_trt = F,
                             vis_avg = NULL,
                             vis_sub = NULL,
                             center_va0 = NULL) {
  
  # Safety checks
  assert_class(data, "dataOK")
  assert_logical(make_pre_trt, len = 1)
  assert_integerish(vis_avg, null.ok = T)
  assert_integerish(vis_sub, null.ok = T)
  assert_integerish(center_va0, null.ok = T, len = 1, lower = 1, upper = 100)
  if ((!is.null(n_sub)) && (!is.null(id_sub))) {
    stop("At most one of n_sub and id_sub can be specified")
  }
  
  # Getting data-spec
  # base_data safe to use (for efficiency) after safety check
  base_data <- data[data$rep == 1, ]
  by_trial <- ifelse("CAT_STD" %in% colnames(data), T, F)
  n_std <- ifelse(by_trial, length(unique(data$CAT_STD)), 1L)
  n_id <- ifelse(
    by_trial, 
    length(unique(data$ID)) / n_std, 
    length(unique(data$ID))
  )
  
  # Subsampling first n_sub IDs
  if (!is.null(n_sub)) {
    if (n_sub > n_id) {
      stop("n_sub greater than number of ID in data")
    }
    if (by_trial) {
      id_s <- tapply(
        base_data$ID,
        base_data$CAT_STD,
        function(x) head(unique(x), n_sub)
      ) %>% 
        unlist()
    } else {
      id_s <- head(unique(base_data$ID), n_sub)
    }
    
    data <- data[which(data$ID %in% id_s), ]
  }
  
  # Subsampling by supplied id_sub
  if (!is.null(id_sub)) {
    data <- data[which(data$ID %in% id_sub), ]
    
    n_id_std <- data[, length(unique(ID)), by = CAT_STD]
    if (! n_distinct(n_id_std$V1) == 1) {
      txt <- paste(
        "After subsampling with supplied id_sub",
        "have unequal number of ID per study"
      )
      warning(txt)
    }
  }
  
  # Make pre-treated data from naive
  if (make_pre_trt) {
    # Old code, replaced with more efficient data.table version
    # data <- data %>% 
    #   mutate(VIS_BASE = PRE_VIS + 1) %>% 
    #   filter(VIS >= VIS_BASE) %>% 
    #   mutate(VIS = VIS - VIS_BASE) %>% 
    #   group_by(ID, rep) %>% 
    #   mutate(TIME = TIME - TIME[1]) %>% 
    #   mutate(X_VA0 = SIM_VA[1]) %>% 
    #   mutate(SIM_CVA = SIM_VA - X_VA0) %>% 
    #   ungroup()
    
    data <- data[VIS >= PRE_VIS + 1][
      , VIS := VIS - (PRE_VIS + 1)][
        , c("TIME", "X_VA0") := .(TIME - TIME[1], SIM_VA[1]), keyby = .(rep, ID)
      ][
        , SIM_CVA := SIM_VA - X_VA0
      ][]
  }
  
  # If averaging BCVA value across visits
  # helpful for stat analysis
  if (! is.null(vis_avg)) {
    vis_ok <- unique(data$VIS)
    if (! all(vis_avg %in% vis_ok)) {
      stop("Trying to average inexistent visits")
    }
    
    # Data points which do not need averaging
    data0 <- filter(data, ! VIS %in% vis_avg)
    
    # Data points to be averaged
    data1 <- data %>% 
      filter(VIS %in% vis_avg) %>% 
      group_by(ID, rep) %>% 
      mutate(
        SIM_VA = mean(SIM_VA),
        SIM_CVA = SIM_VA - X_VA0,
      ) %>% 
      ungroup()
    
    # Merging the two sub dataframes
    data <- bind_rows(data0, data1) %>% 
      arrange(rep, ID, VIS)
  }
  
  # Reducing visits of interest
  if (!is.null(vis_sub)) {
    vis_ok <- unique(data$VIS)
    if (! all(vis_sub %in% vis_ok)) {
      stop("Trying to select inexistent visits")
    }
    
    data <- data[which(data$VIS %in% vis_sub), ]
  }
  
  if (!is.null(center_va0)) {
    data <- data %>% 
      mutate(X_VA0_CENTER = X_VA0 - center_va0)
  }
  
  # Return
  data
}

# Function to calculate the statistical parameter ------------------------------

#' Function to use model simulations to construct statistical prior
#' @param data is a data-frame, obtained via the `sim_data_reshape`function
#' @param fm (formula) formula for statistical analysis,
#' which is fit to the simulated data.
#' @param sample_std (logical) should a subsample of n_id patients be taken when
#' we use patients from all trials (have n_id * n_trial patients in total).
#' @param std_weights (vector) probability to be sampled from a given study.
#' Only applies when sample_std = T
#' @return A tibble with two columns:
#' - `rep` the simulation index
#' - `PAR` the estimated statistical parameter
calc_stat_par <- function(data,
                          fm,
                          sample_std = F,
                          std_weights = NULL,
                          seed = 123) {
  
  assert_data_frame(data)
  assert_class(fm, "formula")
  must_col <- c("rep", "ID", all.vars(fm))
  assert_names(colnames(data), must.include = must_col)
  assert_integerish(seed)
  
  by_trial <- "CAT_STD" %in% colnames(data)
  if (! by_trial) {
    if (sample_std) {
      stop("sample_std does not make sense without by_trial (CAT_STD) data")
    }
    if (! is.null(std_weights)) {
      stop("std_weights do not make sense without by_trial (CAT_STD) data")
    }
  }
  
  # Subsampling data if needed
  if (sample_std) {
    if (! by_trial) {
      stop("sample_std only works with by_trial (CAT_STD) data")
    }
    if (is.null(std_weights)) {
      stop("Supply std_weights, with sample_std")
    }
    assert_data_frame(std_weights, nrow = n_distinct(data$CAT_STD))
    assert_names(names(std_weights), permutation.of = c("CAT_STD", "WEIGHT"))
    if (!all(std_weights$CAT_STD %in% unique(data$CAT_STD))) {
      stop("names in std_weights$CAT_STD not matching data$CAT_STD")
    }
    if (all.equal(sum(std_weights$WEIGHT), 1) != TRUE) {
      stop("Supplied weights do not sum to one")
    }
    
    n_std <- tapply(data$ID, data$CAT_STD, function(x) length(unique(x))) 
    if (n_distinct(n_std) != 1) {
      txt <- paste(
        "Uneven number of IDs per study...",
        "this function was not made for this."
      )
      stop(txt)
    }
    # Assuming number of IDs in each trial same
    n_sam <- n_std[1]
    
    unids <- unique(data$UNID)
    n_unids <- length(unids)
    reps <- unique(data$rep)
    n_reps <- length(reps)
    set.seed(seed)
    
    # Have one participant simulated for a number of studies
    # Sample study, so that one participant sampled from each study
    sam_std <- sample(std_weights$CAT_STD, n_reps * n_sam, replace = T, prob = std_weights$WEIGHT)
    
    sam_df <- data.frame(
      UNID = rep(unids, n_reps), 
      rep = rep(reps, each = n_sam), 
      CAT_STD = sam_std
    ) 
    
    data <- left_join(sam_df, data, by = c("UNID", "CAT_STD", "rep")) %>% 
      arrange(rep, ID)
  }
  
  # Index over which to apply the regression
  # if pooled model -> REP
  # if trial model + subsampling -> REP
  # if trail model and no subsampling  -> interaction(REP, STD)
  if (by_trial && (!sample_std)) {
    # levels preserves order of CAT_STD
    stds <- unique(data$CAT_STD)
    data <- data %>% 
      mutate(CAT_STD_NUM = as.numeric(factor(CAT_STD, levels = stds))) %>% 
      mutate(IND = as.numeric(interaction(rep, CAT_STD_NUM)))
  } else {
    data <- data %>% 
      mutate(IND = rep)
  }
  
  # Run regression for each index
  n_ind <- max(data$IND)
  coef_reg_l <- lapply(
    1:n_ind,
    function(i) {
      # display progress
      if (i %% 100 == 0) print(paste(i, "/", n_ind))
      
      # Run regression and return coefficients
      dat <- data[data$IND == i, ]
      reg <- lm(fm, data = dat)
      c(IND = i, coef(reg), SIGMA = summary(reg)$sigma)
    }
  )
  
  out <- bind_rows(coef_reg_l) 
  if (by_trial && (!sample_std)) {
    out <- out %>% 
      left_join(data %>% select(IND, CAT_STD) %>% unique())
  }
  
  # Return
  out
}

# Function to run all priors in one go -----------------------------------------

#' Function to construct priors in many different ways in one go
#' 
#' Basically this is a wrapper around calc_stat_par
#' just allowing flexibility for endpoint and baseline analysis months
#' as well as running the prior for the first n_sub patients only.
#' 
#' @inheritParams calc_stat_par
#' @param n_sub: size of the subsample of patients to be used.
#' Note that for each study IDs 1:n_sub are used.
#' @return A tibble fit for plotting via the `plot_all_priors` function
run_all_priors <- function(df_pool,
                           df_trial,
                           month_anal = c(8,9),
                           naive = T,
                           n_sub = NULL,
                           fm,
                           std_w,
                           seed = 1807) {
  
  if (length(month_anal) == 1) {
    v_avg <- NULL
    v_sub <- month_anal
  } else {
    v_avg <- month_anal
    v_sub <- month_anal[length(month_anal)]
  }
  
  ### Preparing data
  
  data_one  <- sim_data_reshape(
    data = df_pool,
    make_pre_trt = !naive,
    vis_avg = v_avg,
    vis_sub = v_sub,
    n_sub = n_sub,
    center_va0 = 50
  )
  data_trial <- sim_data_reshape(
    data = df_trial,
    make_pre_trt = !naive,
    vis_avg = v_avg,
    vis_sub = v_sub,
    n_sub = n_sub,
    center_va0 = 50
  )
  
  # Pooled model - with and without pop par uncertainty
  print("Priors for approach 1")
  app1 <- calc_stat_par(
    data = data_one,
    sample_std = F,
    std_weights = NULL,
    fm = fm,
    seed = seed
  )
  
  # Trial model, no pooling - with and without pop par uncertainty
  print("Priors for approach 2")
  app2 <- calc_stat_par(
    data = data_trial,
    sample_std = F,
    std_weights = NULL,
    fm = fm,
    seed = seed
  )
  
  # # Trial model, pooling with equal weights - with and without pop par unc
  print("Priors for approach 3")
  app3 <- calc_stat_par(
    data = data_trial,
    sample_std = T,
    std_weights = std_w,
    fm = fm,
    seed = seed
  )
  
  ### Results data-frame, format suitable for plotting
  res_df <- bind_rows(
    app1 %>%
      mutate(APPROACH = "NO BTV"),
    app2 %>%
      mutate(APPROACH = "TRIAL-SS"),
    app3 %>%
      mutate(APPROACH = "BTV")
  ) %>% 
    mutate(CAT_STD = if_else(is.na(CAT_STD), "POOL", CAT_STD))
  
  ### Return
  res_df
}

# Quick efficient function to calculate mean CBCVA over time
# avoids overheard of calculate_quantile function
calc_mean_time <- function(df) {
  by_trial <- ifelse("CAT_STD" %in% colnames(df), T, F)
  out <- NA
  if (by_trial) {
    out <- df %>% 
      group_by(TIME, CAT_STD, rep) %>% 
      summarise(MEAN = mean(SIM_CVA, na.rm = T)) %>% 
      ungroup() %>% 
      group_by(TIME, CAT_STD) %>% 
      summarise(
        MEAN_LWR = quantile(MEAN, 0.05),
        MEAN_UPR = quantile(MEAN, 0.95)
      )
  } else {
    out <- df %>% 
      group_by(TIME, rep) %>% 
      summarise(MEAN = mean(SIM_CVA, na.rm = T)) %>% 
      ungroup() %>% 
      group_by(TIME) %>% 
      summarise(
        MEAN_LWR = quantile(MEAN, 0.05),
        MEAN_UPR = quantile(MEAN, 0.95)
      )
  }
  out
}

# Reading and cleaning up simulations results ----------------------------------

# NOTE: In paper considering N_ESS = 683 at max
# so only loading 683 patients (683 * 6 trials, in trial model) for comp efficacy
one_emax_unc <- fread(file.path(load_wd, "pri_one_emax_unc.csv")) %>% 
  tidy_res() %>% 
  filter(ID %in% 1:683) %>% 
  # Minimal data-cleaning to help memory mgmt
  select(-SIM_VA_TRANS, -SIM_CVA_TRANS, -COV_AGE, -X_AGE)

# Loading all chunks of trial_emax
files <- list.files(load_wd)
files_ch <- files[str_detect(files, "pri_trial_emax_unc_ch")]
files_path <- file.path(load_wd, files_ch) %>% 
  str_sort(numeric = T) %>% 
  # Only select 1000 sims for first 700 ids
  # (otherwise memory explodes)
  .[1:14]

trial_emax_unc <- rbindlist(
  lapply(files_path, function(x) {
    print(x)
    tidy_res(fread(x)) %>% 
      # Minimal data-cleaning to help memory mgmt
      select(-SIM_VA_TRANS, -SIM_CVA_TRANS, -COV_AGE, -X_AGE) %>% 
      filter(rep %in% 1:750)
  }),
  use.names = TRUE
) 

# Cross validation of priors
files <- list.files(oos_wd)
files_cv <- files[str_detect(files, "OOS_PRI_ID")]

cv_one_emax_unc <- rbindlist(
  lapply(files_cv, function(x) {
    print(x)
    fread(file.path(oos_wd, x)) %>% 
      left_join(full_data %>% select(ID, VIS, TIME, CAT_STD, X_VA0)) %>% 
      mutate(
        SIM_VA = inv_trans_fun(sim_VA_TRANS),
        X_VA0 = inv_trans_fun(X_VA0),
        SIM_CVA = SIM_VA - X_VA0
      )
  }),
  use.names = TRUE
) 

# Checking that obtained prior is in good format -------------------------------

one_emax_unc <- check_sim_data(one_emax_unc)
trial_emax_unc <- check_sim_data(trial_emax_unc)

# Adding a unique identifier to trial_emax_unc
# ID = {1, 51, 101, ...} have same baseline characteristics, 
# but different trial steady state
# UNID characterizes patients with same IDs but different trials
un_id <- trial_emax_unc %>% pull(X_VA0) %>% {match(., unique(.))}
trial_emax_unc$UNID <- un_id

# General "options" for obtaining priors ---------------------------------------

# Equal weights
std_w_df <-  tribble(
  ~CAT_STD, ~WEIGHT,
  "ANCHOR", 1/6,
  "AVENUE", 1/6,
  "HARBOR", 1/6,
  "MARINA", 1/6,
  "STAIRWAY", 1/6,
  "PIER", 1/6
)

# Setup colors
shades_grey <- colorRampPalette(c("grey30", "grey80"))
cols_df <- tibble(
  BRK = c("POOL", "ANCHOR", "AVENUE", "HARBOR", "MARINA", "PIER", "STAIRWAY"),
  VAL = c("black", viridis(6)),
  LTY = c("2262", "22", "42", "44", "13", "F282", "12223242"),
  PCH = c(NA, 0:5)
)

# Calculating mean profiles: naive and pre-trt ---------------------------------
# For efficiency reasons, do not store sim_data_reshape res

# Treatment naive
vpc_naive_const_one_emax_unc_683 <- calc_mean_time(
  sim_data_reshape(one_emax_unc, n_sub = 683, make_pre_trt = F)
)
vpc_naive_const_trial_emax_unc_683 <- calc_mean_time(
  sim_data_reshape(trial_emax_unc, n_sub = 683, make_pre_trt = F)
)
vpc_naive_const_one_emax_unc_30 <- calc_mean_time(
  sim_data_reshape(one_emax_unc, n_sub = 30, make_pre_trt = F)
)
vpc_naive_const_trial_emax_unc_30 <- calc_mean_time(
  sim_data_reshape(trial_emax_unc, n_sub = 30, make_pre_trt = F)
)

gc()

# Pre-treated
vpc_pre_treat_const_one_emax_unc_683 <- calc_mean_time(
  sim_data_reshape(one_emax_unc, n_sub = 683, make_pre_trt = T)
)
vpc_pre_treat_const_trial_emax_unc_683 <- calc_mean_time(
  sim_data_reshape(trial_emax_unc, n_sub = 683, make_pre_trt = T)
)
vpc_pre_treat_const_one_emax_unc_30 <- calc_mean_time(
  sim_data_reshape(one_emax_unc, n_sub = 30, make_pre_trt = T)
)
vpc_pre_treat_const_trial_emax_unc_30 <- calc_mean_time(
  sim_data_reshape(trial_emax_unc, n_sub = 30, make_pre_trt = T)
)

# Plot mean longitudinal profiles ----------------------------------------------

plot_df <- bind_rows(
  vpc_naive_const_one_emax_unc_683 %>% 
    mutate(
      MODEL = "PK/BCVA", 
      CAT_STD = "POOL",
      N_ESS = "683",
      POP = "NAIVE"
    ),
  vpc_naive_const_trial_emax_unc_683 %>% 
    mutate(
      MODEL = "PK/BCVA-TRIAL",
      N_ESS = "683",
      POP = "NAIVE"
    ),
  vpc_naive_const_one_emax_unc_30 %>% 
    mutate(
      MODEL = "PK/BCVA", 
      CAT_STD = "POOL",
      N_ESS = "30",
      POP = "NAIVE"
    ),
  vpc_naive_const_trial_emax_unc_30 %>% 
    mutate(
      MODEL = "PK/BCVA-TRIAL",
      N_ESS = "30",
      POP = "NAIVE"
    ),
  vpc_pre_treat_const_one_emax_unc_683 %>% 
    mutate(
      MODEL = "PK/BCVA", 
      CAT_STD = "POOL",
      N_ESS = "683",
      POP = "PRE-TRT"
    ),
  vpc_pre_treat_const_trial_emax_unc_683 %>% 
    mutate(
      MODEL = "PK/BCVA-TRIAL",
      N_ESS = "683",
      POP = "PRE-TRT"
    ),
  vpc_pre_treat_const_one_emax_unc_30 %>% 
    mutate(
      MODEL = "PK/BCVA", 
      CAT_STD = "POOL",
      N_ESS = "30",
      POP = "PRE-TRT"
    ),
  vpc_pre_treat_const_trial_emax_unc_30 %>% 
    mutate(
      MODEL = "PK/BCVA-TRIAL",
      N_ESS = "30",
      POP = "PRE-TRT"
    )
) %>% 
  mutate(N_ESS = as.factor(N_ESS)) %>% 
  mutate(MODEL = as.factor(MODEL)) %>% 
  mutate(POP = as.factor(POP)) %>% 
  mutate(MONTH = TIME / 28) %>% 
  rename(STUDY = CAT_STD) %>% 
  filter(MONTH < 10)

levels(plot_df$N_ESS) <- c(expression("N[ESS] == 30"), expression("N[ESS] == 683"))
# levels(pri_data$APPROACH) <- c(expression("NO~BTV"), expression("BTV"))

# Define one common graph for both
plot_long_sim_data <- function(df) {
  df %>% 
    ggplot() +
    geom_line(aes(x = MONTH, y = MEAN_LWR, col = STUDY, linetype = STUDY), lwd = 0.9) +
    geom_line(aes(x = MONTH, y = MEAN_UPR, col = STUDY, linetype = STUDY), lwd = 0.9) +
    facet_grid(N_ESS ~ MODEL, labeller = label_parsed) +
    theme_bw(base_size = 22) +
    scale_x_continuous(breaks = scales::pretty_breaks()) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) + 
    scale_color_manual(breaks = cols_df$BRK, values = cols_df$VAL) +
    # BW version (for main doc), cols are OK in supp mat
    # scale_color_manual(breaks = cols_df$BRK, values = cols_df$VAL) +
    scale_linetype_manual(breaks = cols_df$BRK, values = cols_df$LTY, name = "STUDY") +
    theme(legend.position = "bottom") +
    ylab("CBCVA [ETDRS letters]") +
    xlab("Month") 
}
  
# Plot naive and pre-treated
p1 <- plot_df %>% 
  filter(POP == "NAIVE") %>% 
  plot_long_sim_data()

p2 <- plot_df %>% 
  filter(POP == "PRE-TRT") %>% 
  plot_long_sim_data()

if (save_flag) {
  ggsave(file.path(fig_wd, "CVALongNaive.eps"), p1, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "CVALongPreTrt.eps"), p2, width = 12, height = 9, dpi = "retina")
}

# Getting stat priors -----------------------------------------------------------

naive_pri_683 <- run_all_priors(
  df_pool = one_emax_unc,
  df_trial = trial_emax_unc,
  month_anal = 9,
  naive = T,
  fm = formula("SIM_CVA ~ 1 + X_VA0_CENTER"),
  std_w = std_w_df,
  n_sub = 683
)

gc() 

naive_pri_30 <- run_all_priors(
  df_pool = one_emax_unc,
  df_trial = trial_emax_unc,
  month_anal = 9,
  naive = T,
  fm = formula("SIM_CVA ~ 1 + X_VA0_CENTER"),
  std_w = std_w_df,
  n_sub = 30
)

gc() 

pre_trt_pri_683 <- run_all_priors(
  df_pool = one_emax_unc,
  df_trial = trial_emax_unc,
  month_anal = 9,
  naive = F,
  fm = formula("SIM_CVA ~ 1 + X_VA0_CENTER"),
  std_w = std_w_df,
  n_sub = 683
)

gc() 

pre_trt_pri_30 <- run_all_priors(
  df_pool = one_emax_unc,
  df_trial = trial_emax_unc,
  month_anal = 9,
  naive = F,
  fm = formula("SIM_CVA ~ 1 + X_VA0_CENTER"),
  std_w = std_w_df,
  n_sub = 30
)

gc() 

pri_data <- bind_rows(
  naive_pri_683 %>% 
    filter(APPROACH != "TRIAL-SS") %>% 
    rename(PAR = "(Intercept)") %>% 
    mutate(N_ESS = "683") %>% 
    mutate(POP = "NAIVE"),
  naive_pri_30 %>%
    filter(APPROACH != "TRIAL-SS") %>% 
    rename(PAR = "(Intercept)") %>% 
    mutate(N_ESS = "30") %>% 
    mutate(POP = "NAIVE"),
  pre_trt_pri_683 %>% 
    filter(APPROACH != "TRIAL-SS") %>% 
    rename(PAR = "(Intercept)") %>% 
    mutate(N_ESS = "683") %>% 
    mutate(POP = "PRE-TRT"),
  pre_trt_pri_30 %>%
    filter(APPROACH != "TRIAL-SS") %>% 
    rename(PAR = "(Intercept)") %>% 
    mutate(N_ESS = "30") %>% 
    mutate(POP = "PRE-TRT"),
) %>% 
  mutate(N_ESS = as.factor(N_ESS)) %>% 
  mutate(APPROACH = factor(APPROACH, levels = c("NO BTV", "BTV")))

levels(pri_data$N_ESS) <- c(expression("N[ESS] == 30"), expression("N[ESS] == 683"))
levels(pri_data$APPROACH) <- c(expression("NO~BTV"), expression("BTV"))

# Summary statistics (mean and sd) data-frame
stat_df <- pri_data %>%
  group_by(APPROACH, N_ESS, POP) %>%
  mutate(
    MEAN = round(mean(PAR), 2),
    SD = round(sd(PAR), 2)
  ) %>%
  ungroup() %>%
  filter(IND == 1) %>%
  select(POP, APPROACH, N_ESS, MEAN, SD) %>%
  unique()

# Common part of plot for pre-trt and trt-naive
plot_dens <- function(pri_df) {
  pri_df %>% 
    ggplot(aes(x = PAR)) +
    geom_density(aes()) +
    stat_pointinterval(.width = c(.66, .95)) +
    facet_grid(N_ESS ~ APPROACH, labeller = label_parsed) +
    xlab(expression(tau[C])) +
    ylab("Density") +
    theme_bw(base_size = 18) +
    scale_x_continuous(breaks = scales::pretty_breaks()) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) + 
    theme(legend.position = "bottom") 
}

p1 <- plot_dens(
  pri_df = pri_data %>% filter(POP == "NAIVE")
) + 
  coord_cartesian(ylim = c(0, 0.5)) +
  geom_text(
    aes(label = paste("MEAN =", MEAN, "\n", "SD =", SD)),
    x = 5,
    y = 0.35,
    size = 4.5,
    data = stat_df %>% filter(POP == "NAIVE")
  ) 
p2 <- plot_dens(pri_data %>% filter(POP == "PRE-TRT")) + 
  coord_cartesian(ylim = c(0, 0.8)) +
  geom_text(
    aes(label = paste("MEAN =", MEAN, "\n", "SD =", SD)),
    x = -1,
    y = 0.6,
    size = 4.5,
    data = stat_df %>% filter(POP == "PRE-TRT")
  ) 

if (save_flag) {
  ggsave(file.path(fig_wd, "StatParDensNaive.eps"), p1, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "StatParDensPreTrt.eps"), p2, width = 12, height = 9, dpi = "retina")
}

# Comparing all priors ---------------------------------------------------------

library(RBesT)
library(cmdstanr)

# STAN power prior with 
pp_mod <- cmdstan_model(stan_file = file.path(getwd(), "Stan/pow_pri.stan"))

# Linear priors
data_read <- read_csv("~/MIDD_VA2DF/Data/full_data.csv") 

# Note: Some studies only know SE because have access to data
meta_df <- tribble(
  ~STUDY, ~TRT, ~TIME, ~CBCVA, ~SE,
  "MARINA","0.5mgq4w", 12, 6.6, 0.812, 
  "MARINA","0.3mgq4w", 12, 5.4, 0.976,
  "ANCHOR","0.5mgq4w", 12, 11.3, 1.24, 
  "ANCHOR","0.3mgq4w", 12, 8.3, 1.28,
  "PIER", "0.3mgq12w", 12, -0.2, 2.08,
  "PIER", "0.5mgq12w", 12, -1.6, 1.74,
  "SAILOR", "0.3mgPRN", 12, 0.5, 1, # SE HERE IS BALLPARK FROM PUBLI
  "SAILOR", "0.5mgPRN", 12, 2.3, 1, # SE HERE IS BALLPARK FROM PUBLI
  "SUSTAIN", "0.3/0.5mgPRN", 12, 3.6, 0.61,
  "EXCITE", "0.3mgq4w", 12, 8.3, 1.56,
  "EXCITE", "0.3mgq12w", 12, 4.9, 1.74,
  "EXCITE", "0.5mgq12w", 12, 3.8, 1.62,
  "CATT", "0.5mgq4w", 12, 8.5, 0.8,
  "CATT", "0.5mgPRN", 12, 6.8, 0.8,
  "IVAN", "0.5mgq4w + PRN", 24, 4.9, 0.92,
  "VIEW1", "0.5mgq4w", 24, 8.1, 0.87,
  "VIEW2", "0.5mgq4w", 24, 9.4, 0.80,
  "HARBOR", "0.5mgq4w", 12, 10.1, 0.80,
  "HARBOR", "0.5mgPRN", 12, 8.2, 0.80,
  "HARBOR", "0.3mgq4w", 12, 9.2, 0.88,
  "HARBOR", "0.3mgPRN", 12, 8.6, 0.83,
  "MANTA", "0.5mgPRN", 12, 4.3, 1.02, # 6.3 - 4.3 / 1.96 = SE, BALLPARK FROM PUB
  "GEFAL", "0.5mgPRN", 12, 2.93, 1.11,
  "LUCAS", "0.5mgT&E", 12, 7.6, 1.10,
  "CEDAR+SEQUIOIA", "0.5mgq4w", 12, 9.5, 0.6,
  "AVENUE", "0.5mgq4w", 9, 8.5, 1.30,
  "STAIRWAY", "0.5mgq4w", 12, 9.6, 1.28
) %>% 
  mutate(ARM = paste(STUDY, TRT, sep = "_"))


# Deriving linear priors
data <- data_read %>% 
  select(ID, VIS, TIME, Y, CVA, X_VA0, starts_with("CAT_"), starts_with("COV_")) %>% 
  filter(!is.na(Y)) %>% 
  rename(
    STUDY = CAT_STD,
    TRT = CAT_TRT,
    SEX = CAT_SEX,
    CNV = CAT_CNV,
    AGE = COV_AGE,
    LEAKAGE = COV_LEA,
    BBCVA  = X_VA0
  ) 

data_soc <- data %>% 
  # Select only M9
  filter(VIS == 9) %>% 
  # Select only SOC
  filter(TRT == "0.5MG_q4w") %>% 
  # Obvious as only dealing with soc
  select(-TRT)

sig <- 12
map_mcmc_all <- gMAP(
  cbind(CBCVA, SE) ~ 1 | ARM, 
  data = meta_df,
  family = "gaussian",
  beta.prior = cbind(0, sig),
  tau.dist = "HalfNormal",
  tau.prior = cbind(0, sig / 2)
)
print(map_mcmc_all)

map_mcmc_soc <- gMAP(
  cbind(CBCVA, SE) ~ 1 | ARM, 
  data = meta_df %>% filter(TRT == "0.5mgq4w"),
  family = "gaussian",
  beta.prior = cbind(0, sig),
  tau.dist = "HalfNormal",
  tau.prior = cbind(0, sig / 2)
)
print(map_mcmc_soc)

# Power prior, alpha = 0.2
dat_pp <- list(
  N = nrow(data_soc),
  M = 2, #-2 covs: intercept + BBCVA
  alpha = 0.2,
  y = data_soc$CVA,
  X = cbind(1, data_soc$BBCVA - 50) # centering, important!
)
pp_02_stan <- pp_mod$sample(
  data = dat_pp,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.9
) 

# Power prior, alpha = 1
dat_pp$alpha <- 1
pp_1_stan <- pp_mod$sample(
  data = dat_pp,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.9
) 

# Plotting priors for new trial

# Naive model based priors
plot_naive_df_nlme <- bind_rows(
  naive_pri_30 %>% 
    filter(APPROACH != "TRIAL-SS") %>% 
    mutate(APPROACH = ifelse(APPROACH == "NO BTV", "NOBTV", APPROACH)) %>% 
    mutate(PRIOR = paste0("MI_", APPROACH, "_30")),
  naive_pri_683 %>%
    filter(APPROACH != "TRIAL-SS") %>% 
    mutate(APPROACH = ifelse(APPROACH == "NO BTV", "NOBTV", APPROACH)) %>% 
    mutate(PRIOR = paste0("MI_", APPROACH, "_683"))
) %>% 
  rename(PAR = "(Intercept)") %>% 
  select(APPROACH, PAR, PRIOR)

# Pre-treated model based priors
plot_pre_trt_df_nlme <- bind_rows(
  pre_trt_pri_30 %>% 
    filter(APPROACH != "TRIAL-SS") %>% 
    mutate(APPROACH = ifelse(APPROACH == "NO BTV", "NOBTV", APPROACH)) %>% 
    mutate(PRIOR = paste0("MI_", APPROACH, "_30")),
  pre_trt_pri_683 %>%
    filter(APPROACH != "TRIAL-SS") %>% 
    mutate(APPROACH = ifelse(APPROACH == "NO BTV", "NOBTV", APPROACH)) %>% 
    mutate(PRIOR = paste0("MI_", APPROACH, "_683"))
) %>% 
  rename(PAR = "(Intercept)") %>% 
  select(APPROACH, PAR, PRIOR)

plot_new_df_stan <- data.frame(
  MAP_ALL = as.data.frame(map_mcmc_all$fit)$theta_resp,
  MAP_SOC = as.data.frame(map_mcmc_soc$fit)$theta_resp,
  PP_02 = as.vector(pp_02_stan$draws(variables = "beta[1]")),
  PP_1 = as.vector(pp_1_stan$draws(variables = "beta[1]"))
) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "PRIOR",
    values_to = "PAR"
  ) 

# Plot without BTV, as main text -----------------------------------------------

plot_new_df <- bind_rows(
  plot_naive_df_nlme %>% mutate(POP = "NAIVE"), 
  plot_pre_trt_df_nlme %>% mutate(POP = "PRE-TREATED") , 
  plot_new_df_stan %>% mutate(POP = "NAIVE")
) %>% 
  filter(PRIOR %in% c("MI_NOBTV_30", "MI_NOBTV_683", "MAP_ALL", "MAP_SOC", "PP_02", "PP_1")) %>% 
  # renaming priors
  mutate(PRIOR = ifelse(PRIOR == "MI_NOBTV_30", "MI_NC/3", PRIOR)) %>% 
  mutate(PRIOR = ifelse(PRIOR == "MI_NOBTV_683", "MI_NH", PRIOR)) %>% 
  mutate(PRIOR = factor(
    PRIOR,
    levels = c("MAP_ALL", "MAP_SOC", "PP_1", "PP_02", "MI_NH", "MI_NC/3"),
    labels = c("MAP_ALL", "MAP_SOC", "PP_1", "PP_02", "MI_NH", "MI_NC/3")
  ))

ann_text_df <- data.frame(
  x = c(4, NA),
  y = c(0, NA),
  TXT = rep(paste("NH =", nrow(data_soc), " NC = 90"), 2),
  POP = c("NAIVE", "PRE-TREATED")
)

p <- plot_new_df %>% 
  ggplot() +
  aes(x = PRIOR, y = PAR, group = PRIOR) +
  theme_classic(base_size = 20) +
  stat_interval(position = "dodge", .width = c(.5, .66, .95)) +
  geom_text(aes(x = x, y = y, label = TXT), data = ann_text_df, inherit.aes = F) +
  ylab(expression(beta[0])) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
  ) + 
  scale_color_viridis(discrete = T) +
  facet_wrap(~ POP, scales = "free") +
  labs(col = "Interval")

if (save_flag) {
  ggsave(file.path(fig_wd, "PriNewStd.eps"), p, width = 12, height = 9, dpi = "retina")
}

# Plot with BTV, as in supp mat ------------------------------------------------
plot_new_df <- bind_rows(
  plot_naive_df_nlme %>% mutate(POP = "NAIVE"), 
  plot_pre_trt_df_nlme %>% mutate(POP = "PRE-TREATED"), 
  plot_new_df_stan %>% mutate(POP = "NAIVE")
) %>% 
  mutate(PRIOR = factor(
    PRIOR,
    levels = c("MAP_ALL", "MAP_SOC", "PP_02", "PP_1", "MI_NOBTV_30", "MI_BTV_30", "MI_NOBTV_683", "MI_BTV_683"),
    labels = c("MAP_ALL", "MAP_SOC", "PP_02", "PP_1", "MI_NOBTV_NC/3", "MI_BTV_NC/3", "MI_NOBTV_NH", "MI_BTV_NH")
  ))

p <- plot_new_df %>% 
  ggplot() +
  aes(x = PRIOR, y = PAR, group = PRIOR) +
  theme_classic(base_size = 18) +
  stat_interval(position = "dodge", .width = c(.5, .66, .95)) +
  ylab(expression(beta[0])) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) + 
  scale_color_viridis(discrete = T) +
  facet_wrap(~ POP, scales = "free") +
  labs(col = "Quantile") 

if (save_flag) {
  ggsave(file.path(fig_wd, "PriNewStd_supp.eps"), p, width = 12, height = 9, dpi = "retina")
}

# Priors cross validation ------------------------------------------------------


stds <- c("ANCHOR", "AVENUE", "HARBOR", "MARINA", "STAIRWAY")
res_l <- list()
for (i in seq_along(stds)) {
  
  print(i)
  s <- stds[i]
  
  # All data, common part of interest
  dat_tmp <- full_data %>% 
    # Select one study
    filter(VIS == 9, EVID == 0, CAT_TRT == "0.5MG_q4w") %>% 
    mutate(
      X_VA0 = inv_trans_fun(X_VA0),
      Y = inv_trans_fun(Y),
      CVA = Y - X_VA0
    ) %>% 
    select(ID, VIS, CVA, X_VA0, CAT_STD) %>% 
    mutate(X_VA0 = X_VA0 - 50)
  
  dat_obs <- dat_tmp %>% filter(CAT_STD == s)
  dat_fit <- dat_tmp %>% filter(CAT_STD != s)
  
  # True effect
  true_beta0 <- coef(lm("CVA ~ 1 + X_VA0", data = dat_obs))[1]
  
  # MI priors
  n_pp <- n_distinct(dat_fit$ID)
  n_std_div3 <- ceiling(n_distinct(dat_obs$ID) / 3)
  if (s == "STAIRWAY") {
    n_std_div3 <- ceiling(n_distinct(dat_obs$ID) / 2) 
  }
  
  dat_mi <- cv_one_emax_unc %>% 
    mutate(X_VA0 = X_VA0 - 50) %>% 
    filter(CAT_STD == s) %>%  
    group_by(rep) 
    
  mi_pri1 <- dat_mi %>% 
    slice_sample(n = n_pp, replace = T) %>% 
    summarise(BETA0 = coef(lm(SIM_CVA ~ 1 + X_VA0))[1])
  n_iter <- length(mi_pri1$BETA0)
  mi_pri2 <- dat_mi %>% 
    slice_sample(n = n_std_div3, replace = T) %>% 
    summarise(BETA0 = coef(lm(SIM_CVA ~ 1 + X_VA0))[1])

  # Power prior alpha = 1 
  dat_pp <- list(
    N = nrow(dat_fit),
    M = 2, #-2 covs: intercept + BBCVA
    alpha = 1,
    y = dat_fit$CVA,
    X = cbind(1, dat_fit$X_VA0) 
  )
  
  pp_1_stan <- pp_mod$sample(
    data = dat_pp,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.9
  ) 
  
  # Power prior alpha = 0.2 
  dat_pp <- list(
    N = nrow(dat_fit),
    M = 2, #-2 covs: intercept + BBCVA
    alpha = 0.2,
    y = dat_fit$CVA,
    X = cbind(1, dat_fit$X_VA0) 
  )
  
  pp_02_stan <- pp_mod$sample(
    data = dat_pp,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.9
  ) 
  
  # Meta-analytic prior
  sig <- 12
  map_mcmc_all <- gMAP(
    cbind(CBCVA, SE) ~ 1 | ARM, 
    data = meta_df %>% filter(STUDY != s),
    family = "gaussian",
    beta.prior = cbind(0, sig),
    tau.dist = "HalfNormal",
    tau.prior = cbind(0, sig / 2)
  )
  print(map_mcmc_all)
  
  map_mcmc_soc <- gMAP(
    cbind(CBCVA, SE) ~ 1 | ARM, 
    data = meta_df %>% filter(TRT == "0.5mgq4w") %>% filter(STUDY != s),
    family = "gaussian",
    beta.prior = cbind(0, sig),
    tau.dist = "HalfNormal",
    tau.prior = cbind(0, sig / 2)
  )
  print(map_mcmc_soc)
  
  
  # Results
  res_l[[i]] <- data.frame(
    STD = s,
    N_FIT = nrow(dat_fit),
    N_STD = nrow(dat_obs),
    BETA0_OBS = true_beta0,
    MI_NH = mi_pri1$BETA0,
    MI_2 = mi_pri2$BETA0,
    PP_1 = as.vector(pp_1_stan$draws(variables = "beta[1]"))[1:n_iter],
    PP_02 = as.vector(pp_02_stan$draws(variables = "beta[1]"))[1:n_iter],
    MAP_ALL = as.data.frame(map_mcmc_all$fit)$theta_resp[1:n_iter],
    MAP_SOC = as.data.frame(map_mcmc_soc$fit)$theta_resp[1:n_iter]
  )
}

res_df <- bind_rows(res_l) %>% 
  rename("MI_NC/3" = MI_2)

obs_df <- res_df %>% 
  group_by(STD) %>% 
  slice(1) %>% 
  select(BETA0_OBS, N_FIT, N_STD) %>% 
  mutate(
    N_FIT = paste("NH =", N_FIT),
    N_STD = paste("NC =", N_STD)
  )

p <- res_df %>% 
  select(-BETA0_OBS, -N_FIT, -N_STD) %>% 
  pivot_longer(
    cols = -STD,
    values_to = "Y",
    names_to = "PRI"
  ) %>% 
  mutate(PRI = factor(
    PRI,
    levels = c("MAP_ALL", "MAP_SOC", "PP_02", "PP_1", "MI_NH", "MI_NC/3")
  )) %>% 
  ggplot() +
  aes(x = PRI, y = Y) + 
  stat_interval(position = "dodge", .width = c(.5, .66, .95)) + 
  geom_hline(aes(yintercept = BETA0_OBS), data = obs_df) +
  geom_text(aes(label = N_FIT), x = 5, y = 0, data = obs_df) +
  geom_text(aes(label = N_STD), x = 5, y = -2.5, data = obs_df) +
  facet_wrap(~ factor(
    STD,
    levels = c("HARBOR", "MARINA", "ANCHOR", "AVENUE", "STAIRWAY")
  )) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)
  ) + 
  scale_color_viridis(discrete = T) +
  labs(col = "Interval") +
  ylab(expression(beta[0])) +
  xlab("PRIOR")


if (save_flag) {
  ggsave(file.path(fig_wd, "PriCV.eps"), p, width = 12, height = 9, dpi = "retina")
}


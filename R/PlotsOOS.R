## AUTHOR: Francesco Brizzi
## AIM: Longitudinal BCVA and CBCVA profiles were simulated (i.e. prior predictive distributions)
##      based upon the PK/BCVA and PK/BCVA trial models (in OosSim.R) for out-of-sample trials.
##      Here, we read the results and create plots from the results, namely:
##      - Out of sample visual predictive checks for the the PK/BCVA and PK/BCVA trial models.
##      - Out of sample comparison at month 9 for standard of care ranibizumab regimens
##        for the the power prior, PK/BCVA, and PK/BCVA trial models.

# Setup ------------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(checkmate)
library(data.table)
library(cmdstanr)

save_flag <- T

# Cluster setup ----------------------------------------------------------------

on_cluster <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

# Working directories and loading data -----------------------------------------

if (on_cluster) {
  base_wd <- "/pstore/home/brizzif/VCPaper"
  scratch_wd <- "/pstore/scratch/u/brizzif/VCPaper"
} else {
  base_wd <- "W:/Users/brizzif/VCPaper"
  scratch_wd <- "S:/VCPaper"
}

vpc_sim_wd <- file.path(scratch_wd, "OOSSim")
fig_wd <- file.path(scratch_wd, "Figures")

# Useful functions -------------------------------------------------------------

source(file.path(base_wd, "R/FunHelp.R"))

# Function to subsample VPC simulations ----------------------------------------

subsample_vpc_sim_by_std <- function(vpc_sim_df, prob_df, seed = 1824) {
  assert_data_frame(prob_df)
  assert_names(colnames(prob_df), must.include = c("STD", "PROB"))
  assert_data_frame(vpc_sim_df)
  assert_names(colnames(vpc_sim_df), must.include = c("ID", "rep", "STD"))
  
  n_rep <- length(unique(vpc_sim_df$rep))
  n_ids <- length(unique(vpc_sim_df$ID))
  stds <- unique(vpc_sim_df$STD)
  
  if (! all(stds %in% prob_df$STD)) {
    stop("Supplied study names do not correspond to those of vpc_sim_df$STD")
  }
  
  p <- prob_df %>% 
    filter(STD %in% stds) %>% 
    pull(PROB)
  if (sum(p) != 1) {
    stop("probs for CV scheme do not sum to 1")
  }
  
  set.seed(seed)
  sam_std <- sample(stds, n_ids * n_rep, replace = T, prob = p)
  df <- data.frame(
    ID = rep(unique(vpc_sim_df$ID), n_rep), 
    rep = rep(unique(vpc_sim_df$rep), each = n_ids), 
    STD = sam_std
  )
  
  left_join(df, vpc_sim_df, by = c("ID", "rep", "STD"))
}

# Function to clean VPC and PRED res -------------------------------------------

#' Function to tidy VPC simulations
tidy_vpc <- function(res, data) {
  
  # Argument checks
  assert_data_frame(res)
  must_cols <- c("rep", "ID", "TIME", "sim_VA_TRANS", "STD")
  assert_names(colnames(res), must.include = must_cols)
  assert_data_frame(data)
  
  # Get nice results
  out <- res %>% 
    left_join(
      data %>% 
        filter(EVID == 0) %>% 
        select(
          ID, TIME, VIS, Y, CVA, 
          starts_with("X_"), starts_with("CAT"), starts_with("COV")
        ),
      by = c("ID", "TIME")
    ) %>% 
    # Create consistent cols and names so that can use plot vpc function
    rename(
      VA_TRANS = Y,
      CVA_TRANS = CVA,
      SIM_VA_TRANS = sim_VA_TRANS
    ) %>%
    mutate(
      SIM_CVA_TRANS = SIM_VA_TRANS - X_VA0,
      SIM_VA = inv_trans_fun(SIM_VA_TRANS),
      SIM_CVA = SIM_VA - inv_trans_fun(X_VA0),
      VA = inv_trans_fun(VA_TRANS),
      CVA = VA - inv_trans_fun(X_VA0)
    ) %>%
    mutate(CAT_ARM = paste(CAT_STD, CAT_TRT)) %>% 
    # Return tibble
    as_tibble()
  
   out <- out %>% 
     # Nice order
     select(
       rep, ID, TIME, VIS, starts_with("X_"),
       VA, CVA, VA_TRANS, CVA_TRANS,
       SIM_VA, SIM_CVA, SIM_VA_TRANS, SIM_CVA_TRANS,
       starts_with("CAT_"), starts_with("COV_")
     ) 
  
  # Return
  out
}

# Function to create PRED from linear model ------------------------------------

get_pop_pred_pp <- function(data, fit, seed = 123) {
  
  assert_class(fit, "CmdStanFit")
  assert_data_frame(data$pred)
  assert_names(colnames(data$pred), must.include = c("ID", "BBCVA"))
  
  X <- cbind(1, data$pred$BBCVA - 50)
  betas <- cbind(
    as.vector(fit$draws(variables = "beta[1]")),
    as.vector(fit$draws(variables = "beta[2]"))
  )
  sig <- as.vector(fit$draws(variables = "sigma"))
  
  n_id <- nrow(X)
  n_sim <- length(sig)
  
  set.seed(seed)
  
  res_l <- lapply(
    1:n_sim,
    function(i) {
      pred_df <- data.frame(
        ID = data$pred$ID,
        BBCVA = data$pred$BBCVA,
        REP = i,
        SIM_CVA = X %*% betas[i, ] + rnorm(n_id, 0, sig[i])
      ) %>% 
        mutate(SIM_VA = BBCVA + SIM_CVA)
    }
  )
  
  bind_rows(res_l)
}


# Function for CV scheme -------------------------------------------------------

prepare_cv_data <- function(data, oos_trial) {
  # Checks
  must_cols <- c("STUDY", "CVA", "BBCVA", "ID")
  assert_data_frame(data)
  assert_names(colnames(data), must.include = must_cols)
  ok_std <- unique(data$STUDY)
  assert_character(oos_trial, len = 1)
  assert_choice(oos_trial, ok_std)
  
  # Data, excluding out of sample study
  data_std <- data %>% 
    filter(STUDY != oos_trial) %>% 
    select(ID, CVA, BBCVA, STUDY) %>% 
    mutate(STUDY_CAT = as.numeric(as.factor(STUDY)))
  
  # Data to predict
  data_pred <- data %>% 
    filter(STUDY == oos_trial) %>% 
    select(ID, CVA, BBCVA) 
  
  # Aggregate data for RBest
  data_agg <- data_std %>% 
    group_by(STUDY) %>% 
    summarise(
      SIGMA = sd(CVA),
      N = n(),
      CVA_MEAN = mean(CVA),
      CVA_SE = SIGMA / sqrt(N)
    ) 
  
  # Return individual and aggregate data
  out <- list(
    "ind" = data_std, 
    "agg" = data_agg,
    "pred" = data_pred
  )
}

# Compiling stan models --------------------------------------------------------

pp_mod <- cmdstan_model(stan_file = file.path(getwd(), "Stan/pow_pri.stan"))

# Loading full data ------------------------------------------------------------

if (on_cluster) {
  full_data <- read_csv("~/MIDD_VA2DF/Data/reg_trans_data.csv") 
} else {
  proj_wd <- getwd()
  trt_wd <- file.path(proj_wd, "Monolix/ModelsRes/pkbcva_trans2") 
  full_data <- read_data_project(trt_wd, ret = "data")
}

# Global options for reading VPC simulations -----------------------------------

all_std <- c("ANC", "AVE", "HAR", "MAR", "STA")
files <- list.files(vpc_sim_wd)

# VPC results from trial model, with subsampling
# NOTE: PROB = 1/5 rather than 1/6 bc one trial used for model fitting
p_df <- data.frame(
  STD = c("ANCHOR", "AVENUE", "HARBOR", "MARINA", "PIER", "STAIRWAY"),
  PROB = rep(1/5, 6)
)

# Loading VPC results ----------------------------------------------------------

# VPC results from pooled model
vpc_pool_res <- list()
for(i in 1:length(all_std)) {
  std <- all_std[i]
  print(std)
  txt <- paste0("OOS_TRANS_VPC_", std)
  
  files_std <- files[str_detect(files, txt)]
  
  vpc_pool_res[[i]] <- rbindlist(
    lapply(files_std, function(x) {
      print(x)
      fread(file.path(vpc_sim_wd, x)) %>%
        tidy_vpc(data = full_data) %>% 
        # Memory management
        select(
          -X_AGE, -VA_TRANS, -CVA_TRANS, -SIM_VA_TRANS, - SIM_CVA_TRANS,
          -CAT_SEX, -CAT_CNV, -COV_AGE, -COV_LEA, -COV_VA0
        )
    }),
    use.names = TRUE
  )
}

vpc_trial_res <- list()
for(i in 1:length(all_std)) {
  std <- all_std[i]
  print(std)
  txt <- paste0("OOS_TRANS_TRIAL_VPC_", std)
  
  files_std <- files[str_detect(files, txt)]
  
  vpc_trial_res[[i]] <- rbindlist(
    lapply(files_std, function(x) {
      print(x)
      fread(file.path(vpc_sim_wd, x)) %>%
        subsample_vpc_sim_by_std(prob_df = p_df) %>% 
        tidy_vpc(data = full_data) %>% 
        # Memory management
        select(
          -X_AGE, -VA_TRANS, -CVA_TRANS, -SIM_VA_TRANS, - SIM_CVA_TRANS,
          -CAT_SEX, -CAT_CNV, -COV_AGE, -COV_LEA, -COV_VA0
        )
    }),
    use.names = TRUE
  )
}

# Minimal data cleaning to VPC output -----------------------------------------

clean_vpc <- . %>% 
  bind_rows() %>% 
  mutate(
    CAT_ARM = factor(
      CAT_ARM, 
      levels = c(
        "ANCHOR 0.3MG_q4w", "ANCHOR 0.5MG_q4w", "HARBOR 0.5MG_q4w", 
        "HARBOR 2MG_q4w", "MARINA 0.3MG_q4w", "MARINA 0.5MG_q4w", 
        "AVENUE 0.5MG_q4w", "STAIRWAY 0.5MG_q4w"
      )
    )
  ) %>% 
  mutate(CAT_ARM = str_replace_all(CAT_ARM, "_", " ")) %>% 
  select(-CAT_TRT, -CAT_STD) # %>% 
  # Step below helpful, but uneeded and slow
  # arrange(rep, ID, TIME)

vpc_pool_res <- clean_vpc(vpc_pool_res)
vpc_trial_res <- clean_vpc(vpc_trial_res)

# Cross-validation loop to obtain linear model predictions ---------------------

data_soc <- full_data %>%
  # Select only observations
  filter(EVID == 0) %>% 
  # Select only M9
  filter(VIS == 9) %>% 
  # Select only SOC
  filter(CAT_TRT == "0.5MG_q4w") %>% 
  # Obvious as only dealing with soc
  select(-CAT_TRT) %>% 
  # To keep consistent with existing functions
  mutate(BBCVA = COV_VA0 * 100) %>% 
  # Bringing BCVA back to untransformed scale
  mutate(
    Y = inv_trans_fun(Y),
    X_VA0 = inv_trans_fun(X_VA0),
    CVA = round(Y - X_VA0)
  ) %>% 
  rename(STUDY = CAT_STD, BCVA = Y) %>% 
  select(ID, BCVA, CVA, BBCVA, STUDY)

stds <- unique(data_soc$STUDY)
res_pp <- list()
for (i in 1:length(stds)) {
  
  print(paste("i =", i, "/", length(stds)))
  
  # OOS study
  oos_std <- stds[i]
  
  # Eliminate OOS study from data-fitting process
  dat <- prepare_cv_data(data = data_soc, oos_trial = oos_std)

  dat_pp <- list(
    N = nrow(dat$ind),
    M = 2, #-2 covs: intercept + BBCVA
    alpha = 1,
    y = dat$ind$CVA,
    X = cbind(1, dat$ind$BBCVA - 50) # centering, important!
  )
  
  # Power prior, alpha = 1
  pp_1_stan <- pp_mod$sample(
    data = dat_pp,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.9
  ) 
  
  as.vector(pp_1_stan$draws(variables = "beta[1]"))
  
  # Create tibble with results
  res_pp[[i]] <- data.frame(
    get_pop_pred_pp(data = dat, fit = pp_1_stan),
    STUDY = oos_std
  )
}

# Getting VPC constructors -------------------------------------------------

vpc_const_pool_va <- calculate_vpc_quant(
  sim_df = vpc_pool_res %>% select(-VA) %>% rename(VA = SIM_VA),
  data_df = vpc_pool_res %>% filter(rep == 1),
  out_col = "VA",
  rep_col = "rep",
  time_col = "VIS",
  strat_cols = c("CAT_ARM")
) 

vpc_const_pool_cva <- calculate_vpc_quant(
  sim_df = vpc_pool_res %>% select(-CVA) %>% rename(CVA = SIM_CVA),
  data_df = vpc_pool_res %>% filter(rep == 1),
  out_col = "CVA",
  rep_col = "rep",
  time_col = "VIS",
  strat_cols = c("CAT_ARM")
) 

vpc_const_trial_va <- calculate_vpc_quant(
  sim_df = vpc_trial_res %>% select(-VA) %>% rename(VA = SIM_VA),
  data_df = vpc_trial_res %>% filter(rep == 1),
  out_col = "VA",
  rep_col = "rep",
  time_col = "VIS",
  strat_cols = c("CAT_ARM")
) 

vpc_const_trial_cva <- calculate_vpc_quant(
  sim_df = vpc_trial_res %>% select(-CVA) %>% rename(CVA = SIM_CVA),
  data_df = vpc_trial_res %>% filter(rep == 1),
  out_col = "CVA",
  rep_col = "rep",
  time_col = "VIS",
  strat_cols = c("CAT_ARM")
) 

# Data manipulation for plotting ---------------------------------------

reshape_const <- function(vpc_const) {
  assert_class(vpc_const, "vpc_const")
  
  # Put simulations in long format for ggplot
  vpc_const$sim %>% 
    pivot_longer(
      cols = c(
        "MEAN_LWR", "MEAN_UPR", "MEDIAN_LWR", "MEDIAN_UPR",
        "Q10_LWR", "Q10_UPR", "Q90_LWR", "Q90_UPR"
      ),
      names_to = "TMP",
      values_to = "Y"
    ) %>% 
    separate(TMP, into = c("QUANTILE", "BOUND")) %>% 
    pivot_wider(
      names_from = BOUND,
      values_from = Y
    ) %>% 
    mutate(
      QUANTILE = factor(QUANTILE, levels = c("Q90", "MEDIAN", "MEAN", "Q10"))
    ) 
} 

# VPC plot ---------------------------------------------------------------------

sim_df <- bind_rows(
  vpc_const_pool_cva %>%
    reshape_const() %>%
    mutate(TYPE = "PK/BCVA (NO BTV)") %>% 
    mutate(OUT = "CVA"),
  vpc_const_trial_cva %>% 
    reshape_const() %>%
    mutate(TYPE = "PK/BCVA-TRIAL (BTV)") %>% 
    mutate(OUT = "CVA"),
  vpc_const_pool_va %>%
    reshape_const() %>%
    mutate(TYPE = "PK/BCVA (NO BTV)") %>% 
    mutate(OUT = "VA"),
  vpc_const_trial_va %>% 
    reshape_const() %>%
    mutate(TYPE = "PK/BCVA-TRIAL (BTV)") %>% 
    mutate(OUT = "VA")
) 
dat_df <- bind_rows(
  vpc_const_pool_cva$data %>% mutate(OUT = "CVA"),
  vpc_const_pool_va$data %>% mutate(OUT = "VA")
)
  
shades_grey <- colorRampPalette(c("grey30", "grey80"))
cols <- shades_grey(4)

vpc_plot <- function(simul_df, data_df) {
  
  simul_df %>% 
    ggplot() +
    aes(x = VIS) +
    geom_ribbon(aes(ymin = LWR, ymax = UPR), lwd = 1.2, fill = "grey90") +
    geom_point(aes(y = Y), data = data_df, lwd = 1.3) +
    geom_line(aes(y = Y), data = data_df, lwd = 0.7) +
    facet_grid(
      fct_relevel(QUANTILE,'Q90','MEDIAN','MEAN','Q10') ~ CAT_ARM, 
      scales = "free_y",
      labeller = label_wrap_gen(width = 5)
    ) +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.title = element_blank()
    ) +
    xlab("Months")
}

p_cva <- vpc_plot(
  simul_df = sim_df %>% filter(OUT == "CVA") %>% filter(TYPE == "PK/BCVA (NO BTV)"),
  data_df = dat_df %>% filter(OUT == "CVA")
) +
  ylab("CBCVA [ETDRS Letters]") +
  theme(legend.position = "none")

p_va <- vpc_plot(
  simul_df = sim_df %>% filter(OUT == "VA") %>% filter(TYPE == "PK/BCVA (NO BTV)"),
  data_df = dat_df %>% filter(OUT == "VA")
) +
  ylab("BCVA [ETDRS Letters]") +
  theme(legend.position = "none")

p_cva_supp <- vpc_plot(
  simul_df = sim_df %>% filter(OUT == "CVA"),
  data_df = dat_df %>% filter(OUT == "CVA")
) +
  ylab("CBCVA [ETDRS Letters]") 

p_va_supp <- vpc_plot(
  simul_df = sim_df %>% filter(OUT == "VA"),
  data_df = dat_df %>% filter(OUT == "VA")
) +
  ylab("BCVA [ETDRS Letters]") 

if (save_flag) {
  ggsave(file.path(fig_wd, "OOSVPC_CVA.eps"), p_cva, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "OOSVPC_BCVA.eps"), p_va, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "OOSVPC_CVA_supp.eps"), p_cva_supp, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "OOSVPC_BCVA_supp.eps"), p_va_supp, width = 12, height = 9, dpi = "retina")
}


# Comparing OOS prediction at month 9 ------------------------------------------

lin_pred <- bind_rows(res_pp) %>% 
  pivot_longer(
    cols = c("SIM_VA", "SIM_CVA"),
    names_to = "OUT",
    values_to = "Y"
  ) %>% 
  group_by(REP, STUDY, OUT) %>% 
  summarise(
    Q10 = quantile(Y, 0.1),
    Q90 = quantile(Y, 0.9),
    MEAN = mean(Y),
    MEDIAN = median(Y)
  ) %>% 
  pivot_longer(
    cols = c("Q10", "Q90", "MEDIAN", "MEAN"),
    names_to = "STAT",
    values_to = "Y"
  ) %>% 
  group_by(STAT, STUDY, OUT) %>% 
  summarise(
    Q05 = quantile(Y, 0.05),
    Q95 = quantile(Y, 0.95)
  ) %>% 
  mutate(CAT_ARM = paste(STUDY, "0.5MG q4w")) %>% 
  rename(LWR = Q05, UPR = Q95) %>% 
  mutate(MODEL = "TRIAL") %>% 
  mutate(OUT = str_remove_all(OUT, "SIM_"))

plot_df <- bind_rows(
  lin_pred %>%
    ungroup() %>%
    select(-STUDY),
  sim_df %>% 
    ungroup() %>%
    rename(STAT = QUANTILE, MODEL = TYPE) %>%
    filter(VIS == 9) %>% 
    filter(str_detect(CAT_ARM, "0.5MG q4w")) %>% 
    select(-VIS)
) %>% 
  mutate(
    XMIN = case_when(
      MODEL == "TRIAL" ~ 0.5,
      MODEL == "PK/BCVA (NO BTV)" ~ 1.5,
      MODEL == "PK/BCVA-TRIAL (BTV)" ~ 2.5,
      TRUE ~ NA_real_
    )
  ) %>% 
  mutate(XMAX = XMIN + 1)  %>% 
  mutate(MODEL = case_when(
    MODEL ==  "TRIAL" ~ "TRIAL",
    MODEL == "PK/BCVA (NO BTV)" ~ "PK/BCVA",
    MODEL == "PK/BCVA-TRIAL (BTV)" ~ "PK/BCVA-TRIAL",
    TRUE ~ NA_character_
  )) %>%
  mutate(OUT = case_when(
    OUT ==  "CVA" ~ "CBCVA",
    OUT ==  "VA" ~ "BCVA",
    TRUE ~ NA_character_
  )) 

data_soc_summ <- data_soc %>% 
  rename(VA = BCVA) %>% 
  pivot_longer(
    cols = c("VA", "CVA"),
    names_to = "OUT",
    values_to = "Y"
  ) %>% 
  group_by(STUDY, OUT) %>% 
  summarise(
    Q10 = quantile(Y, 0.1),
    Q90 = quantile(Y, 0.9),
    MEAN = mean(Y),
    MEDIAN = median(Y)
  ) %>% 
  pivot_longer(
    cols = c("Q10", "Q90", "MEDIAN", "MEAN"),
    names_to = "STAT",
    values_to = "Y"
  ) %>% 
  mutate(CAT_ARM = paste(STUDY, "0.5MG q4w")) %>%
  mutate(OUT = case_when(
    OUT ==  "CVA" ~ "CBCVA",
    OUT ==  "VA" ~ "BCVA",
    TRUE ~ NA_character_
  )) 

plot_m9_preds <- function(pred_df, dat_df) {
  
  mods <- unique(pred_df$MODEL)
  
  # Deriving center of intervals
  xs <- pred_df %>% 
    mutate(X = XMIN + ((XMAX - XMIN)/2)) %>%
    pull(X)
  
  # Creating a copy of dat_df for each X
  dat_df <- expand_grid(dat_df, X = xs)
  
  # Plot
  pred_df %>% 
    ggplot() +
    aes(group = STAT, col = STAT) +
    geom_rect(
      aes(xmin = XMIN, xmax = XMAX, ymin = LWR, ymax = UPR), 
      fill = NA, 
      lty = 2
    ) + 
    geom_point(aes(y = Y, x = X, col = STAT, group = STAT), data = dat_df) + 
    geom_point(aes(y = Y, x = X, col = STAT, group = STAT), data = dat_df) + 
    facet_grid(OUT ~ CAT_ARM, labeller = label_wrap_gen(width = 5), scales = "free") +
    scale_x_continuous(breaks = 1:n_distinct(mods), labels = unique(mods)) +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.title = element_blank()
    ) +
    scale_color_viridis(discrete = T) + 
    ylab("ETDRS Letters") +
    xlab("MODEL")
}

p1 <- plot_m9_preds(
  pred_df = plot_df %>% filter(!MODEL %in% "PK/BCVA-TRIAL"),
  dat_df = data_soc_summ
)
  
p2 <- plot_m9_preds(
  pred_df = plot_df,
  dat_df = data_soc_summ
)

if (save_flag) {
  ggsave(file.path(fig_wd, "OOSSOC9M.eps"), p1, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "OOSSOC9M_supp.eps"), p2, width = 12, height = 9, dpi = "retina")
}

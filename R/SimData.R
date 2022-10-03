## AUTHOR: Francesco Brizzi
## AIM: Generate simulated data for clinical trial simulations.
##      Baseline characteristics and dosing regimen defined in MasterTrialSim.R
##      Longitudinal profiles for each simulation are simulated based on PK/BCVA model
##      only cross-sectional (9 months) CBCVA outcomes are retained
##      and linear conflict can be addded to introduce time bias
##      (i.e. kappa parameter in paper)

# Setup ------------------------------------------------------------------------

rm(list = ls())
save_flag <- T

# Working directories ----------------------------------------------------------

base_wd <- "W:/Users/brizzif/VCPaper"
scratch_wd <- "S:/VCPaper"
fig_wd <- file.path(base_wd, "Figures")
r_wd <- file.path(base_wd, "R")
store_wd <- file.path(scratch_wd, "TrialSimData")
fig_wd <- file.path(base_wd, "Figures")

# Loading useful functions -----------------------------------------------------

source(file.path(r_wd, "FunHelp.R"))
source(file.path(r_wd, "FunDataSim.R"))

# Loading simulation parameters ------------------------------------------------

# Simulation master
master_sim_df <- data.table::fread(file.path(store_wd, "sim_master_df.csv"))
# Baseline characteristics
base_read_df <- data.table::fread(file.path(store_wd, "base_df.csv"))

# Parameters for data-generation -----------------------------------------------

# Scenarios (TRT_EFF, CONFLICT) to be simulated
# Scenarios with fewer individuals will just be considered 
# as a subsample of full dataset
scenario_data_sim_df <- master_sim_df %>% 
  select(TRT, CONFLICT, DATA_IND) %>% 
  unique() %>% 
  as_tibble()

# Number of datasets per scenario
n_rep <- max(master_sim_df$DATA)

# Assumed trial lenght (months)
n_mon_tri <- 9

# Sample size choice -----------------------------------------------------------
# for 1:1 design

# Data simulation
dat_pow_an <- simulate_trial(
  n_trt = 100, 
  n_soc = 100,
  trt_eff = 0,
  hist_conf = 0,
  base_df = cbind(
    id = 1:200, rbind(base_read_df[1:100, -1], base_read_df[1:100, -1])
  ),
  model = "one_emax",
  unc_pop_pars = T,
  n_sim = n_rep,
  mon_tri = n_mon_tri,
  mon_out = n_mon_tri,
  seed = 1825 + 42, # +42 to keep consistent with wider data simulation
  naive = F
)

# Running power analyses with incrementally fewer patients
run_pow_analysis(data = dat_pow_an)
# 100:100 -> ~ 95% power
run_pow_analysis(data = dat_pow_an %>% filter(id %in% c(1:90, 101:190)))
# 90:90  -> ~ 91% power
run_pow_analysis(data = dat_pow_an %>% filter(id %in% c(1:80, 101:180)))
# 80:80 ->  88% power

run_pow_analysis(data = dat_pow_an %>% filter(id %in% c(1:90, 101:170)))
# 90:70 -> 88% power
run_pow_analysis(data = dat_pow_an %>% filter(id %in% c(1:90, 101:160)))
# 90:60  -> ~ 85% power

# Posterior sigma distribution plot --------------------------------------------

# Varying N of trial
n_seq <- seq(10, 100, 5)
res_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(res_df) <- c("Mean", "Q05", "Q95")

for (i in 1:length(n_seq)) {
  n <- n_seq[i]
  print(n)
  
  # Subsampling patients to obtain smaller n
  dat <- dat_pow_an %>% 
    filter(ARM == 0) %>% 
    group_by(rep) %>% 
    sample_n(n) %>%  
    ungroup()
  
  # Calculating standard deviation
  sigs <- sapply(
    1:max(dat$rep), 
    function(x) {
      summary(lm(sim_CVA ~ 1 + X_VA0, data = dat[dat$rep == x, ]))$sigma
    }
  )
  
  # Return
  res_df[i, ] <- c(mean(sigs), quantile(sigs, probs = c(0.05, 0.95)))
}

p <- res_df %>% 
  mutate(N = n_seq) %>% 
  pivot_longer(-N, names_to = "Quantities", values_to = "Y") %>% 
  ggplot() +
  geom_smooth(aes(x = N, y = Y, group = Quantities), se = FALSE, col = "black") +
  theme_bw(base_size = 18) +
  ylab(expression(sigma)) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) 

ggsave(file.path(fig_wd, "SigmaPredDist.eps"), p, width = 12, height = 9, dpi = 300)

# Running data-generation process ----------------------------------------------

# Based on above results, choose 90:90 as base scenario
# + 90:60 for reduced sample size scenario
n_soc <- 90

# Redefining baseline characteristics dataframe
# Assuming SOC and TRT have same baseline characteristics
base_df <- rbind(base_read_df[1:n_soc, ], base_read_df[1:n_soc, ])
base_df$id <- 1:(2 * n_soc)

for (i in 1:nrow(scenario_data_sim_df)) {
  
  print(i)
  
  trt <- scenario_data_sim_df[i, "TRT", drop = T]
  conf <- scenario_data_sim_df[i, "CONFLICT", drop = T]
  data_ind <- scenario_data_sim_df[i, "DATA_IND", drop = T]
  
  # Data simulation
  # IMPORTANT: if i == 17 -> add nrow(scenario_data_sim_df) + 1 to random seed
  # as seed 1825 gives very weird results
  # (likely caused by simulx bug)
  tictoc::tic()
  sim_dat <- simulate_trial(
    n_trt = n_soc, 
    n_soc = n_soc,
    trt_eff = trt,
    hist_conf = conf,
    base_df = base_df,
    model = "one_emax",
    unc_pop_pars = T,
    n_sim = n_rep,
    mon_tri = n_mon_tri,
    mon_out = n_mon_tri,
    seed = 1825 + i + ifelse(i == 17, nrow(scenario_data_sim_df) + 1, 0),
    naive = F
  )
  tictoc::toc()
  
  # Saving if needed
  if (save_flag) {
    save_txt <- paste0("Data_", data_ind, ".csv")
    data.table::fwrite(sim_dat, file.path(store_wd, save_txt))
  }  
}

# Checks for i = 42 (no conflict, no trt eff) ----------------------------------

i <- scenario_data_sim_df %>%
  filter(TRT == 0 & CONFLICT == 0) %>% 
  pull(DATA_IND)
trt <- scenario_data_sim_df[i, "TRT", drop = T]
conf <- scenario_data_sim_df[i, "CONFLICT", drop = T]
data_ind <- scenario_data_sim_df[i, "DATA_IND", drop = T]

# Data simulation
tictoc::tic()
sim_dat <- simulate_trial(
  n_trt = n_soc, 
  n_soc = n_soc,
  trt_eff = trt,
  hist_conf = conf,
  base_df = base_df,
  model = "one_emax",
  unc_pop_pars = T,
  n_sim = n_rep,
  mon_tri = n_mon_tri,
  mon_out = n_mon_tri,
  seed = 1825 + i,
  naive = F
)
tictoc::toc()

run_pow_analysis(dat = sim_dat)
# Power is ~90%

sim_dat %>% 
  group_by(rep, ARM) %>% 
  summarise(a = mean(sim_CVA)) %>% 
  ungroup() %>% 
  group_by(ARM) %>% 
  summarise(a = mean(a))

# Prior for full study (n = 90)
# ---> need for prior mean
pri_all <- generate_prior(data = sim_dat %>% filter(ARM == 0), n = 90)
summary(pri_all)
# ---> prior mean shall be centered around 1.52
sd(pri_all)

# Prior for 30 additional patients
pri_30 <- generate_prior(data = sim_dat %>% filter(ARM == 0), n = 30)
mean(pri_30)
sd(pri_30)
# ---> prior has sd of 2.32

# Prior for 60 additional patients
pri_60 <- generate_prior(data = sim_dat %>% filter(ARM == 0), n = 60)
mean(pri_60)
sd(pri_60)
# ---> prior has sd of 1.52
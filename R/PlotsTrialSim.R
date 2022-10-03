## AUTHOR: Francesco Brizzi
## AIM: Plots results from trial simulations 
##      i.e. Fig 6 in manuscript

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(RColorBrewer)
library(distributional)
library(ggdist)
library(checkmate)
library(viridis)

# Using fread and fwrite, without exporting whole data.table
fread <- data.table::fread
fwrite <- data.table::fwrite

save_flag <- FALSE

# Cluster setup ----------------------------------------------------------------

on_cluster <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

# Working directories and loading data -----------------------------------------

if (on_cluster) {
  base_wd <- "/pstore/home/brizzif/VCPaper"
  scratch_wd <- "/pstore/scratch/u/brizzif/VCPaper"
  cmdstan_wd <- "/pstore/home/brizzif/.cmdstanr/cmdstan-2.27.0/"
  save_wd <- "/pstore/scratch/u/brizzif/VCPaper/Figures"
} else {
  base_wd <- "W:/Users/brizzif/VCPaper"
  scratch_wd <- "S:/VCPaper"
  cmdstan_wd <- "C:/Users/brizzif/Documents/.cmdstanr/cmdstan-2.27.0/"
  save_wd <- "S:/VCPaper/Figures"
  
}

res_wd <- file.path(scratch_wd, "TrialSimRes")
data_wd <- file.path(scratch_wd, "TrialSimData")
fig_wd <- file.path(scratch_wd, "Figures")

# Read in simulation master ----------------------------------------------------

master_sim_df <- fread(file.path(data_wd, "sim_master_df.csv"))

# Reading in results ----------------------------------------------------------

res_files <- list.files(res_wd, pattern = "Chunk")

# Check do I have all results?
all(unique(master_sim_df$CHUNK) %in% parse_number(res_files))

res_l <- list()
for (i in 1:length(res_files)) {
  if (i %% 50 == 0) print(i)
  res_l[[i]] <- fread(file.path(res_wd, res_files[i])) 
}

# Make results a data.frame
res <- bind_rows(res_l) %>% 
  # Better names / types for plotting
  mutate(
    PRIOR = case_when(
      PRIOR == "NON_INFO" ~ "WI",
      PRIOR == "INFO" ~ "MI",
      PRIOR == "POWER" ~ "MIP",
      PRIOR == "MIX" ~ "MIM",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(PRIOR = as.factor(as.character(PRIOR)))

# Cleaning memory
rm(res_l)

# Checking results for convergence and timeout -----------------------------------

# Check, and remove if needed, not-converged simulations
if (any(res$CONV == FALSE)) {
  # Count non converged simulations
  count(res, CONV)
  # Remove simulations which have not converged
  res <- filter(res, CONV == FALSE)
}

# Check, and remove if needed, timed-out simulations
if (any(res$TIMEOUT == TRUE)) {
  # Count non converged simulations
  count(res, TIMEOUT)
  # Remove simulations which have not converged
  res <- filter(res, TIMEOUT == TRUE)
}

# Obtain number of converged iteration for each scenario
res <- res %>% 
  group_by(DATA_SIM_IND, PRIOR) %>% 
  mutate(N_CONV = sum(CONV)) %>% 
  ungroup() 

# Simulated data plots ---------------------------------------------------------

data_files <- list.files(data_wd, pattern = "Data")

data_l <- list()
for (i in 1:length(data_files)) {
  if (i %% 50 == 0) print(i)
  data_ind <- parse_number(data_files[i])
  data_l[[i]] <- fread(file.path(data_wd, data_files[i])) %>% 
    mutate(DATA_SIM_IND = data_ind) 
}

# Adding useful info
data_df <- bind_rows(data_l) %>% 
  rename(DATA = rep) %>% 
  left_join(
    res %>% 
      select(DATA_SIM_IND, DATA, TRT, CONFLICT) %>% 
      unique()
  )

# Cleaning memory
# rm(data_l)

# Plotting standard deviations
# RECALL: DATA_SIM_IND = interaction(TRT, CONFLICT)
sd_df <- data_df %>% 
  group_by(DATA, DATA_SIM_IND) %>% 
  summarise(SD = sd(sim_CVA)) 

sd_df %>% 
  ggplot() +
  geom_histogram(aes(x = SD)) 


# Plotting effects
eff_data_df <- data_df %>% 
  group_by(DATA, DATA_SIM_IND, ARM) %>% 
  summarise(MU = mean(sim_CVA)) %>% 
  pivot_wider(names_from = ARM, values_from = MU) %>% 
  mutate(DIFF = `1` - `0`) 

# Examine that simulated effects are reasonable
# Important because some seeds fail with Monolix/Simulx
n_row_fig <- 5
n_col_fig <- 5
n_pg_fig <- n_row_fig * n_col_fig
n_pg <- ceiling(max(eff_data_df$DATA_SIM_IND) / n_pg_fig)

for (i in 1:n_pg) {
  ind <- (n_pg_fig * (i - 1) + 1):(n_pg_fig * i)   
    
  p <- eff_data_df %>% 
    filter(DATA_SIM_IND %in% ind) %>% 
    ggplot() +
    geom_histogram(aes(x = DIFF)) +
    facet_wrap(~ DATA_SIM_IND, nrow = n_row_fig, ncol = n_col_fig)
  
  print(p)
}

# Global plotting definitions --------------------------------------------

# Defining colors
cols <- viridis(4)
graph_df <- tribble(
  ~ PRIOR, ~ COL, ~ LTY, ~SHAPE,
  "WI", cols[4], 1L, 21,
  "MI", cols[3], 22L, 22,
  "MIM", cols[2], 42L, 23,
  "MIP", cols[1], 44L, 24,
)

# Datasets for plotting --------------------------------------------------

# Type I error case, when treatment effect assumed to be -3.9
t1_df <- res %>% 
  # Type I error data only when TRT = -3.9
  filter(TRT == -3.9) %>% 
  # group by sim scenario and prior
  group_by(DATA_SIM_IND, PRIOR, N_CONV, CONFLICT, N_PBO) %>% 
  # Scale results by number of converged
  summarise(TOT_ERR = sum(DEC)) %>% 
  mutate(T1_ERROR = TOT_ERR / N_CONV) %>% 
  ungroup()

# Power case, when treatment effect assumed to be 0
pow_df <- res %>% 
  filter(TRT > -3.9) %>% 
  # group by sim scenario and prior
  group_by(DATA_SIM_IND, PRIOR, N_CONV, TRT, CONFLICT, N_PBO) %>% 
  # Scale results by number of converged
  summarise(TOT_GOOD = sum(DEC)) %>% 
  mutate(POWER = TOT_GOOD / N_CONV) 

op_char_df <- res %>% 
  filter(TRT > -3.9) %>% 
  # group by sim scenario and prior
  group_by(DATA_SIM_IND, PRIOR, N_CONV, TRT, CONFLICT, N_PBO) %>% 
  summarise(
    BIAS = mean(TRT_BIAS),
    SD = mean(TRT_sd),
    RMSD = mean(TRT_RMSD),
    DEC = mean(TRT_Q025)
  )

commens_df <- res %>%
  filter(PRIOR %in% c("MIP")) %>% 
  group_by(PRIOR, CONFLICT, TRT, N_PBO) %>%  
  summarise(
    DELTA = mean(POW_DELTA),
    SD_PRI = mean(POW_SD_SOC)
  )


# Type I error + bias plots ------------------------------------------------------

t1_df %>%
  ggplot() +
  aes(x = CONFLICT, y = T1_ERROR, group = PRIOR, col = PRIOR) + 
  geom_point() +
  geom_smooth(se = FALSE, span = 0.35, lwd = 1.2) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  theme_bw() +
  facet_wrap(~ N_PBO, labeller = label_bquote(N[PBO]==.(N_PBO))) +
  theme(legend.position = "bottom") +
  ylab("Type I error") +
  xlab(expression(kappa))


pow_df %>% 
  ggplot() +
  aes(x = CONFLICT, y = POWER, group = PRIOR, col = PRIOR, lty = PRIOR) +
  geom_smooth(se = FALSE, span = 0.4, lwd = 1.2) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(
    ~ N_PBO, 
    labeller = label_bquote(cols = N[PBO]==.(N_PBO))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  xlab(expression(kappa))


p1 <- bind_rows(
  t1_df %>% mutate(TYPE = "TYPE I ERROR") %>% rename(Y = T1_ERROR) %>% mutate(TRT = -3.9),
  pow_df %>% mutate(TYPE = "POWER") %>% select(-TOT_GOOD) %>% rename(Y = POWER)
) %>%  
  ggplot() +
  aes(x = CONFLICT, y = Y, group = PRIOR, col = PRIOR, lty = PRIOR) +
  geom_smooth(se = FALSE, span = 0.5, lwd = 1.5) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  facet_grid(
    TYPE ~ N_PBO, 
    labeller = label_bquote(cols = N[SOC]==.(N_PBO)),
    scales = "free_y"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  xlab(expression(kappa)) +
  ylab("")

# Other operating characteristics plots --------------------------------------------------

op_char_df %>% 
  ggplot() +
  aes(x = CONFLICT, y = BIAS, group = PRIOR, col = PRIOR, lty = PRIOR) +
  geom_point() + 
  geom_smooth(se = FALSE, span = 0.4, lwd = 1.2) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(
    ~ N_PBO, 
    labeller = label_bquote(cols = N[PBO]==.(N_PBO))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  xlab(expression(kappa))

op_char_df %>% 
  ggplot() +
  aes(x = CONFLICT, y = SD, group = PRIOR, col = PRIOR, lty = PRIOR) +
  geom_smooth(se = FALSE, span = 0.4, lwd = 1.2) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(
    ~ N_PBO, 
    labeller = label_bquote(cols = N[PBO]==.(N_PBO))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  xlab(expression(kappa))

op_char_df %>% 
  ggplot() +
  aes(x = CONFLICT, y = RMSD, group = PRIOR, col = PRIOR, lty = PRIOR) +
  geom_smooth(se = FALSE, span = 0.4, lwd = 1.2) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(
    ~ N_PBO, 
    labeller = label_bquote(cols = N[PBO]==.(N_PBO))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  xlab(expression(kappa))

op_char_df %>% 
  ggplot() +
  aes(x = CONFLICT, y = DEC, group = PRIOR, col = PRIOR, lty = PRIOR) +
  geom_point() + 
  geom_smooth(se = FALSE, span = 0.4, lwd = 1.2) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(
    ~ N_PBO, 
    labeller = label_bquote(cols = N[PBO]==.(N_PBO))
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  xlab(expression(kappa))

p2 <- op_char_df %>% 
  pivot_longer(
    cols = c("BIAS", "SD", "RMSD"),
    names_to = "TYPE",
    values_to = "Y"
  ) %>%  
  mutate(TYPE = factor(TYPE, levels = c("BIAS", "SD", "RMSD"))) %>% 
  ggplot() +
  aes(x = CONFLICT, y = Y, group = PRIOR, col = PRIOR, lty = PRIOR) +
  geom_smooth(se = FALSE, span = 0.5, lwd = 1.5) +
  scale_color_manual(breaks = graph_df$PRIOR, values = graph_df$COL, name = "PRIOR") +
  scale_linetype_manual(breaks = graph_df$PRIOR, values = graph_df$LTY, name = "PRIOR") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  facet_grid(
    TYPE ~ N_PBO, 
    labeller = label_bquote(cols = N[SOC]==.(N_PBO)),
    scales = "free_y"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  xlab(expression(kappa)) +
  ylab("")

# Delta (commensurability plot) --------------------------------------------

plt_df <- commens_df %>% 
  mutate(N_PBO = as.factor(N_PBO)) %>% 
  pivot_longer(cols = c("DELTA", "SD_PRI"), names_to = "TYPE", values_to = "Y") %>% 
  mutate(TYPE = as.factor(TYPE)) %>% 
  mutate(TRT = as.factor(TRT))

levels(plt_df$TYPE) <- c(expression(hat(delta)), expression(sigma[H]/hat(delta)))
levels(plt_df$TRT) <- c(expression(beta^"*" == -3.9), expression(beta^"*" == 0))

p3 <- plt_df %>% 
  ggplot() +
  aes(x = CONFLICT, y = Y, group = N_PBO, col = N_PBO) +
  geom_point(lwd = 1.2) +
  # geom_smooth(se = FALSE, span = 0.4, lwd = 1.2) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  facet_grid(
    TYPE ~ TRT, 
    labeller = label_parsed,
    scales = "free"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  ylab("") +
  xlab(expression(kappa)) +
  scale_color_viridis(name = expression(N[SOC]), discrete = T) 

# Prior expectations -------------------------------------------------------

calc_mean_oc <- function(oc, pri, n_pbo, n_sam = 1e5) {
  
  assert_character(oc, len = 1)
  assert_choice(oc, c("pow", "t1e"))
  assert_character(pri, len = 1)
  assert_choice(pri,  c("WI", "MI", "MIP", "MIM"))
  assert_integerish(n_pbo, len = 1)
  assert_choice(n_pbo,  c(60, 90))
  assert_integerish(n_sam, len = 1, lower = 0)
  
  # Sample from expected distribution with 90 individuals
  set.seed(1546)
  sd_sam <- switch(
    as.character(n_pbo), 
    "60" = 1.49,
    "90" = 1.26
  )
  pri_sam <- rnorm(n_sam, 0, sd_sam)
  
  df <- switch(
    oc,
    "pow" = pow_df %>% rename(OC = POWER),
    "t1e" = t1_df %>% rename(OC = T1_ERROR),
    stop("Unknown OC")
  )
  
  sub_df <- df %>% 
    filter(PRIOR == pri, N_PBO == n_pbo) 

  confs <- sub_df$CONFLICT
  oc <- sub_df$OC
  
  mean(sapply(pri_sam, function(x) oc[which.min(abs(x-confs))]))
}

des_df <- expand.grid(PRI = unique(t1_df$PRIOR), N_PBO = unique(t1_df$N_PBO)) %>%
  mutate(DES = paste(PRI, N_PBO, sep = "_")) 
des <- des_df$DES
tab_df <- data.frame(matrix(nrow = length(des), ncol = 2))
colnames(tab_df) <- c("TYPE I ERROR", "POWER")

for (i in 1:nrow(des_df)) {
  pri <- as.character(des_df[i, "PRI"])
  n_pbo <- as.numeric(des_df[i, "N_PBO"])
  tab_df[i, "POWER"] <- calc_mean_oc(oc = "pow", pri = pri, n_pbo = n_pbo)
  tab_df[i, "TYPE I ERROR"] <- calc_mean_oc(oc = "t1e", pri = pri, n_pbo = n_pbo)
}

tab <- cbind(des_df[, 1:2], tab_df[, 1:2])
colnames(tab)[1] <- "PRIOR"
colnames(tab)[2] <- "N_SOC"

xtable::xtable(tab, digits = 4)

# Expected powers
calc_mean_oc(oc = "pow", pri = "WI", n_pbo = 90)
calc_mean_oc(oc = "pow", pri = "MIP", n_pbo = 90)
calc_mean_oc(oc = "pow", pri = "MI", n_pbo = 60)
calc_mean_oc(oc = "pow", pri = "MIM", n_pbo = 60)

# Expected t1 errors
calc_mean_oc(oc = "t1e", pri = "WI", n_pbo = 90)
calc_mean_oc(oc = "t1e", pri = "MIP", n_pbo = 90)
calc_mean_oc(oc = "t1e", pri = "MI", n_pbo = 60)
calc_mean_oc(oc = "t1e", pri = "MIM", n_pbo = 60)

# Storing ------------------------------------------------------------------

ggsave(file = file.path(fig_wd, "BiasPowSim.eps"), p1, width = 8, height = 8, dpi = "retina")
ggsave(file = file.path(fig_wd, "OtherOpCharSim.eps"), p2, width = 8, height = 8, dpi = "retina")
ggsave(file = file.path(fig_wd, "PowPriComm.eps"), p3, width = 12, height = 9, dpi = "retina")


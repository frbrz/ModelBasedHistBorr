## AUTHOR: Francesco Brizzi
## AIM: Goodness of fit plots for sham, PK/BCVA, and PK/BCVA-TRIAL models.

# Setup ------------------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(checkmate)
library(cowplot)
library(ggdist)
library(lixoftConnectors)
library(viridis)

# Simulx
initializeLixoftConnectors(
  software = "simulx",
  path = "C:/ProgramData/Lixoft/MonolixSuite2020R1/",
  force = T
)

# Saving options
save_flag <- FALSE

proj_wd <- getwd()
fig_wd <- file.path(proj_wd, "Figures")

# Helpful plotting functions
source(file.path(proj_wd, "R", "FunHelp.R"))

# Monlix result directories ----------------------------------------------------

## FOR FULL REPRODUCIBILITY, PBO MODEL (AND THEREFORE ALL MODELS) SHALL BE RERUN
pbo_wd <- "W:/Users/brizzif/MIDD_VA2DF/Monolix/ModelRes/sham_trans"
trt_wd <- file.path(proj_wd, "Monolix/ModelsRes/pkbcva_trans2") 
# std_emax_wd <- "W:/Users/brizzif/MIDD_VA2DF/Monolix/ModelRes/pkbcva_trans_trial/"
# base_wd <- "W:/Users/brizzif/MIDD_VA2DF/Monolix/ModelRes/pkbcva_raw"

# Running data reader / constructors -------------------------------------------

pbo_mlx <- read_mlx(pbo_wd)
trt_mlx <- read_mlx(trt_wd)
# std_emax_mlx <- read_mlx(trt_wd)

# Placebo plots ----------------------------------------------------------------

p1 <- plot_dv_pred(pbo_mlx, out = "VA", trans = T, type = "IPRED", subsamp = 0.8) + 
  ylim(0, 1) +
  xlim(0, 1)
p2 <- plot_dv_pred(pbo_mlx, out = "VA", trans = F, type = "IPRED", subsamp = 0.8) +
  xlim(0, 100) +
  ylim(0, 100)
p3 <- plot_resid_scatter(pbo_mlx, type = "IWRES", x_ax = "PRED", subsamp = 0.8)
p4 <- plot_resid_scatter(pbo_mlx, type = "IWRES", x_ax = "TIME", subsamp = 0.8) + 
  scale_x_continuous(
    breaks = seq(0, 800, 60),
    labels = seq(0, 800 / 30, 2)
  ) +
  xlab("TIME [Months]")
p5 <- plot_resid_qq(pbo_mlx, type = "IWRES", subsamp = 0.8)
p6 <- plot_resid_hist(pbo_mlx, type = "IWRES", subsamp = 0.8)
p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2) 
if (save_flag) {
  ggsave(file.path(fig_wd, "IndDiagSham.eps"), p, width = 14, height = 14, dpi = "print")
}

p1 <- plot_dv_pred(pbo_mlx, out = "VA", trans = T, type = "PRED", subsamp = 0.8) + 
  ylim(0, 1) +
  xlim(0, 1)
p2 <- plot_dv_pred(pbo_mlx, out = "VA", trans = F, type = "PRED", subsamp = 0.8) +
  xlim(0, 100) +
  ylim(0, 100)
p3 <- plot_resid_scatter(pbo_mlx, type = "PWRES", x_ax = "PRED", subsamp = 0.8)
p4 <- plot_resid_scatter(pbo_mlx, type = "PWRES", x_ax = "TIME", subsamp = 0.8) +
  scale_x_continuous(
    breaks = seq(0, 800, 60),
    labels = seq(0, 800 / 30, 2)
  ) +
  xlab("TIME [Months]")
p5 <- plot_resid_qq(pbo_mlx, type = "PWRES", subsamp = 0.8)
p6 <- plot_resid_hist(pbo_mlx, type = "PWRES", subsamp = 0.8)
p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
if (save_flag) {
  ggsave(file.path(fig_wd, "PopDiagSham.eps"), p, width = 14, height = 14, dpi = "print")
}

p1 <- plot_mlx_vpc(pbo_mlx, out = "VA", strat = "CAT_STD", trans = F, bw = T) +
  xlab("Month") +
  ylab("BCVA [ETDRS Letters]") 
p2 <- plot_mlx_vpc(pbo_mlx, out = "VA", strat = "CAT_VA0", trans = F, bw = T) +
  xlab("Month") +
  ylab("BCVA [ETDRS Letters]") 
if (save_flag) {
  ggsave(file.path(fig_wd, "VPCArmSham.eps"), p1, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "VPCVA0Sham.eps"), p2, width = 12, height = 9, dpi = "retina")
}

# PK/BCVA plots ----------------------------------------------------------------

p1 <- plot_dv_pred(trt_mlx, out = "VA", trans = T, type = "IPRED", subsamp = 0.75) + 
  ylim(0, 2) +
  xlim(0, 2)
p2 <- plot_dv_pred(trt_mlx, out = "VA", trans = F, type = "IPRED", subsamp = 0.75) +
  xlim(0, 100) +
  ylim(0, 100)
p3 <- plot_resid_scatter(trt_mlx, type = "IWRES", x_ax = "PRED", subsamp = 0.75)
p4 <- plot_resid_scatter(trt_mlx, type = "IWRES", x_ax = "TIME", subsamp = 0.75) +
  scale_x_continuous(
    breaks = seq(0, 800, 60),
    labels = seq(0, 800 / 30, 2)
  ) +
  xlab("TIME [Months]")
p5 <- plot_resid_qq(trt_mlx, type = "IWRES", subsamp = 0.75)
p6 <- plot_resid_hist(trt_mlx, type = "IWRES", subsamp = 0.75) +
  xlim(-5, 5)
p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2) 
if (save_flag) {
  ggsave(file.path(fig_wd, "IndDiagPKBCVA.eps"), p, width = 14, height = 14, dpi = "print")
}

p1 <- plot_dv_pred(trt_mlx, out = "VA", trans = T, type = "PRED", subsamp = 0.75) + 
  ylim(0, 2) +
  xlim(0, 1)
p2 <- plot_dv_pred(trt_mlx, out = "VA", trans = F, type = "PRED", subsamp = 0.75) +
  xlim(0, 100) +
  ylim(0, 100)
p3 <- plot_resid_scatter(trt_mlx, type = "PWRES", x_ax = "PRED", subsamp = 0.75)
p4 <- plot_resid_scatter(trt_mlx, type = "PWRES", x_ax = "TIME", subsamp = 0.75) +
  scale_x_continuous(
    breaks = seq(0, 800, 60),
    labels = seq(0, 800 / 30, 2)
  ) +
  xlab("TIME [Months]")
p5 <- plot_resid_qq(trt_mlx, type = "PWRES", subsamp = 0.75)
p6 <- plot_resid_hist(trt_mlx, type = "PWRES", subsamp = 0.75) +
  xlim(-5, 5)
p <- plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)
if (save_flag) {
  ggsave(file.path(fig_wd, "PopDiagPKBCVA.eps"), p, width = 14, height = 14, dpi = "print")
}

p1 <- plot_mlx_vpc(trt_mlx, "VA", strat = c("CAT_STD", "CAT_TRT"), trans = F, bw = T) +
  xlab("Month") +
  ylab("BCVA [ETDRS Letters]") 
p2 <- plot_mlx_vpc(trt_mlx, "VA", strat = "CAT_VA0", trans = F, bw = T)  +
  xlab("Month") +
  ylab("BCVA [ETDRS Letters]") +
  facet_grid(
    fct_relevel(QUANTILE,'Q90','MEDIAN','MEAN','Q10') ~ STRAT, 
    scales = "free_y",
    labeller = label_wrap_gen(width = 5)
  ) 
if (save_flag) {
  ggsave(file.path(fig_wd, "VPCArmPKBCVA.eps"), p1, width = 12, height = 9, dpi = "retina")
  ggsave(file.path(fig_wd, "VPCVA0PKBCVA.eps"), p2, width = 12, height = 9, dpi = "retina")
}

# Dose response ----------------------------------------------------------------

# Importing Monolix project
importMonolixProject(paste0(trt_wd, ".mlxtran"))

all_trt_df <- data.frame(
  amt = c(2000, 500, 300, 500, 300),
  q = c(30, 30, 30, 90, 90),
  nm = c("2000_q4w", "500_q4w", "300_q4w", "500_q12w", "300_q12w")
)

# PRED-like dose response
res <- list()
for (i in 1:nrow(all_trt_df)) {
  trt_i <- all_trt_df[i, ]
  res[[i]] <- sim_one_ind(
    time_sim = seq(0, 1000, by = 2),
    cov_reg = data.frame(X_VA0 = trans_fun(50), X_AGE = 0, COV_AGE = 0, CAT_STD = "A"),
    trt_df = data.frame(
      time = c(0, 30, 60, seq(60 + trt_i$q, 1000, by = trt_i$q)),
      amount = trt_i$amt
    ),
    reff = F,
    res_err = F,
    n_sim = 1,
    group_nm = trt_i$nm
  )
}
res <- bind_rows(res)

p <- res %>%
  rename(TREATMENT = group) %>% 
  group_by(TREATMENT, time) %>% 
  mutate(TREATMENT = str_replace_all(TREATMENT, "_", " ")) %>% 
  summarise(VA = median(VA)) %>%
  mutate(REGIMEN = factor(
      TREATMENT, 
      levels = c("300 q4w", "500 q4w", "2000 q4w", "300 q12w", "500 q12w"),
      labels = c("300mg q4w", "500mg q4w", "2000mg q4w", "300mg q12w", "500mg q12w"),
    )
  ) %>% 
  ggplot() + 
  geom_line(aes(x = time, y = VA, group = REGIMEN, col = REGIMEN), lwd = 1) +
  theme_bw(base_size = 22) +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = scales::pretty_breaks()) + 
  ylab("BCVA [ETDRS letters]") + 
  scale_x_continuous(
    breaks = seq(0, 1000, 60),
    labels = seq(0, 1000 / 30, 2)
  ) +
  xlab("Month")+
  scale_color_viridis(discrete = T)

if (save_flag) {
  ggsave(file.path(fig_wd, "DoseRespFinMod.eps"), p, width = 12, height = 9, dpi = "retina")
}

# Simulate one individual ------------------------------------------------------

ind_va30 <- sim_one_ind(
  time_sim = seq(0, 1000, by = 2),
  cov_reg = data.frame(X_VA0 = trans_fun(30), X_AGE = 0, COV_AGE = 0, CAT_STD = "A"),
  trt_df = data.frame(
    time = c(0, 30, 60, seq(60 + 30, 1000, by = 30)),
    amount = 500
  ),
  reff = T,
  res_err = T,
  n_sim = 500, 
  group_nm = "30"
)

ind_va50 <- sim_one_ind(
  time_sim = seq(0, 1000, by = 2),
  cov_reg = data.frame(X_VA0 = trans_fun(50), X_AGE = 0, COV_AGE = 0, CAT_STD = "A"),
  trt_df = data.frame(
    time = c(0, 30, 60, seq(60 + 30, 1000, by = 30)),
    amount = 500
  ),
  reff = T,
  res_err = T,
  n_sim = 500, 
  group_nm = "50"
)

ind_va70 <- sim_one_ind(
  time_sim = seq(0, 1000, by = 2),
  cov_reg = data.frame(X_VA0 = trans_fun(70), X_AGE = 0, COV_AGE = 0, CAT_STD = "A"),
  trt_df = data.frame(
    time = c(0, 30, 60, seq(60 + 30, 1000, by = 30)),
    amount = 500
  ),
  reff = T,
  res_err = T,
  n_sim = 500, 
  group_nm = "70"
)

tmp_df <- bind_rows(ind_va30, ind_va50, ind_va70) %>% 
  mutate(group = paste("BBCVA =", group)) %>% 
  filter(time %in% c(0:10, seq(10, 1000, 5))) %>% 
  group_by(time, group)

mean_med_df <- tmp_df %>% 
  summarise(MEAN = mean(VA), MEDIAN = median(VA)) %>% 
  group_by(group) %>% 
  mutate(
    MEAN = predict(loess(MEAN ~ time)),
    MEDIAN = predict(loess(MEDIAN ~ time))
  ) %>% 
  pivot_longer(
    cols = c("MEAN", "MEDIAN"),
    names_to = "STAT",
    values_to = "VA"
  )

p <- tmp_df %>% 
  median_qi(VA, .width = c(.5, .66, .9)) %>%
  # Smoothing the edges of the prediction interval via loess
  group_by(.width, group) %>% 
  mutate(.lower = predict(loess(.lower ~ time))) %>% 
  mutate(.upper = predict(loess(.upper ~ time))) %>% 
  mutate(VA = predict(loess(VA ~ time))) %>% 
  mutate(.width = paste(.width * 100, "%")) %>% 
  ggplot() +
  # automatically uses aes(fill = fct_rev(ordered(.width)))
  geom_interval(aes(x = time, y = VA, ymin = .lower, ymax = .upper)) +
  geom_line(
    aes(x = time, y = VA, group = STAT, lty = STAT),
    lwd = 1.2, data = mean_med_df
  ) +
  scale_color_viridis(name = "PREDICTION INTERVAL", discrete = T) +
  scale_linetype_discrete(name = "") +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  scale_x_continuous(
    breaks = seq(0, 1000, 90),
    labels = seq(0, 1000 / 30, 3)
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
  ) +
  facet_wrap(~ group) +
  ylab("BCVA [ETDRS letters]") + 
  xlab("Month") 


if (save_flag) {
  ggsave(file.path(fig_wd, "IndsFinMod.eps"), p, width = 12, height = 9, dpi = "retina")
}
## AUTHOR: Francesco Brizzi
## AIM: Various preliminary plots illustrating ranibizumab data

# Setup ------------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(meta)
library(virdis)

# do we want to save plots / figs
save_flag <- F

# Cluster setup ----------------------------------------------------------------

on_cluster <- ifelse(.Platform$OS.type == "windows", FALSE, TRUE)

# Working directories setup ----------------------------------------------------

if (on_cluster) {
  data_wd <- "/pstore/home/brizzif/MIDD_VA2DF/Data/"
  fig_wd <- "/pstore/scratch/u/brizzif/VCPaper/Figures/"
} else {
  data_wd <- "W:/Users/brizzif/MIDD_VA2DF/Data/"
  fig_wd <- "W:/Users/brizzif/VCPaper/Figures/"
}

# Data loading -----------------------------------------------------------------

data_read <- read_csv(file.path(data_wd, "full_data.csv")) 

# Data manipulation ------------------------------------------------------------

data <- data_read %>% 
  select(ID, X_VA0, starts_with("CAT_"), starts_with("COV_")) %>% 
  rename(
    STD = CAT_STD,
    TRT = CAT_TRT,
    SEX = CAT_SEX,
    CNV = CAT_CNV,
    AGE = COV_AGE,
    LEAKAGE = COV_LEA,
    VA  = X_VA0
  ) %>% 
  mutate(
    CNV = case_when(
      CNV == 1 ~ "Occult",
      CNV == 2 ~ "Minimally Classic",
      CNV == 3 ~ "Predominantly Classic"
    )
  ) %>%
  mutate(SEX = if_else(SEX == "F", "Female", "Male")) %>% 
  filter(!(TRT %in% c("0.5MG_PRN", "2MG_PRN")))


data_plot <- data_read %>% 
  select(ID, VIS , Y, CVA, X_VA0, starts_with("CAT_"), starts_with("COV_")) %>% 
  filter(!is.na(Y)) %>% 
  rename(
    STUDY = CAT_STD,
    TRT = CAT_TRT,
    SEX = CAT_SEX,
    CNV = CAT_CNV,
    AGE = COV_AGE,
    LEAKAGE = COV_LEA,
    BBCVA  = X_VA0
  ) %>% 
  mutate(
    CNV = case_when(
      CNV == 1 ~ "Occult",
      CNV == 2 ~ "Minimally Classic",
      CNV == 3 ~ "Predominantly Classic"
    )
  ) %>%
  mutate(SEX = if_else(SEX == "F", "Female", "Male")) %>% 
  rename(VA = Y)

# Table of baseline characteristics --------------------------------------------

base_tab <- table1(
  ~ SEX + CNV + VA + AGE + LEAKAGE  | STD * TRT, 
  data = data, 
  overall = F
)

# Mapping html format to R data frame
tab <- as.data.frame(read_html(base_tab) %>% html_table(fill = TRUE))

# Plotting clinical outcomes ---------------------------------------------------

p <- data_plot %>% 
  filter(!TRT  %in% c("SHAM", "2MG_PRN", "0.5MG_PRN")) %>% 
  mutate(TRT = str_replace_all(TRT, "_", " ")) %>% 
  mutate(
    TRT = factor(
      TRT, 
      levels = c("0.3MG q4w", "0.5MG q4w", "2MG q4w", "0.3MG q12w", "0.5MG q12w"),
      labels = c("0.3MG q4w", "0.5MG q4w (SoC)", "2MG q4w", "0.3MG q12w", "0.5MG q12w")
    )
  ) %>% 
  mutate(
    STUDY = factor(
      STUDY,
      levels = c("ANCHOR", "MARINA", "PIER", "HARBOR", "AVENUE", "STAIRWAY")
    )
  ) %>% 
  rename(BCVA = VA, CBCVA = CVA) %>% 
  group_by(VIS, STUDY, TRT) %>% 
  summarise(
    N = n(),
    SE = sd(CBCVA) / sqrt(N),
    CBCVA = mean(CBCVA),
  ) %>% 
  ungroup() %>% 
  mutate(VIS = as.numeric(VIS)) %>% 
  ggplot() +
  aes(x = VIS, y = CBCVA, group = STUDY, col = STUDY) +
  geom_line(lwd = 1.3) +
  geom_errorbar(
    aes(ymin = CBCVA - SE, ymax = CBCVA + SE), 
    width = 0,
    position = position_dodge(0.75),
    lwd = 0.7
  ) +
  scale_x_continuous(breaks = seq(0, 26, 4)) + 
  scale_color_viridis(discrete = TRUE) +
  facet_wrap(~ TRT) + 
  theme_minimal(base_size = 22) +
  theme(legend.position = "bottom") + 
  xlab("Month") +
  ylab("CBCVA [ETDRS Letters]") 

if (save_flag) {
  ggsave(file.path(fig_wd, "Fig1.eps"), p, width = 12, height = 9, dpi = "retina")  
}

# CBCVA results ----------------------------------------------------------------

data_plot %>%
  filter(VIS == 12 & STUDY %in% c("ANCHOR", "MARINA", "PIER")) %>%
  group_by(TRT, STUDY) %>%
  summarise(N = n(), SE = sd(CVA) / sqrt(N))

# Note: Some studies only know SE because have access to data
cbcva_rbz_df <- tribble(
  ~STUDY, ~TRT, ~PHASE,  ~TIME, ~N, ~CBCVA, ~SE,
  "MARINA","0.5mgq4w", "III", 24, 239, 6.6, 0.812, 
  "MARINA","0.3mgq4w", "III", 24, 238, 5.4, 0.976,
  "ANCHOR","0.5mgq4w", "III", 24, 140, 11.3, 1.24, 
  "ANCHOR","0.3mgq4w", "III", 24, 140, 8.3, 1.28,
  "PIER", "0.3mgq12w", "IIIb", 12, 59, -0.2, 2.08,
  "PIER", "0.5mgq12w", "IIIb", 12, 61, -1.6, 1.74,
  "SAILOR", "0.3mgPRN", "IIIb", 12, 462, 0.5, 1, # SE HERE IS BALLPARK FROM PUB
  "SAILOR", "0.5mgPRN", "IIIb", 12, 490, 2.3, 1, # SE HERE IS BALLPARK FROM PUB
  "SUSTAIN", "0.3/0.5mgPRN", "IIIb", 12, 513, 3.6, 0.61,
  "EXCITE", "0.3mgq4w", "IIIb", 12, 115, 8.3, 1.56,
  "EXCITE", "0.3mgq12w", "IIIb", 12, 120, 4.9, 1.74,
  "EXCITE", "0.5mgq12w", "IIIb", 12, 118, 3.8, 1.62,
  "CATT", "0.5mgq4w", "IV", 24, 301, 8.5, 0.8,
  "CATT", "0.5mgPRN", "IV", 24, 298, 6.8, 0.8,
  "IVAN", "0.5mgq4w + PRN", "IV", 24, 271, 4.9, 0.92,
  "VIEW1", "0.5mgq4w", "III", 24, 306, 8.1, 0.87,
  "VIEW2", "0.5mgq4w", "III", 24, 303, 9.4, 0.80,
  "HARBOR", "0.5mgq4w", "IIIb", 12, 274, 10.1, 0.80,
  "HARBOR", "0.5mgPRN", "IIIb", 12, 275, 8.2, 0.80,
  "HARBOR", "0.3mgq4w", "IIIb", 12, 274, 9.2, 0.88,
  "HARBOR", "0.3mgPRN", "IIIb", 12, 274, 8.6, 0.83,
  "MANTA", "0.5mgPRN", "IV", 12, 163, 4.3, 1.02, # 6.3 - 4.3 / 1.96 = SE, BALLPARK FROM PUB
  "GEFAL", "0.5mgPRN", "IV", 12, 246, 2.93, 1.11,
  "LUCAS", "0.5mgT&E", "IV", 12, 221, 7.6, 1.10,
  "CEDAR+SEQUIOIA", "0.5mgq4w", "III", 12, 630, 9.5, 0.6,
  "AVENUE", "0.5mgq4w", "II", 9, 67, 8.5, 1.30,
  "STAIRWAY", "0.5mgq4w", "II", 12, 16, 9.6, 1.28
) %>% 
  as.data.frame()

cbcva_soc_df <- cbcva_rbz_df %>% 
  filter(TRT == "0.5mgq4w")

# Meta-analysis res
data_plt_tmp <- cbcva_rbz_df %>% 
  mutate(labels = paste(STUDY, TRT)) %>% 
  # Specify if subject level data available
  mutate(
    SUBLEVDAT = ifelse(
      STUDY %in% c("ANCHOR", "AVENUE", "HARBOR", "MARINA", "PIER", "STAIRWAY"),
      "Yes",
      "No")
  )

data_plt_soc_tmp <- cbcva_rbz_df %>%
  filter(TRT == "0.5mgq4w") %>%
  mutate(labels = paste(STUDY, TRT))

meta_all <- meta::metagen(
  CBCVA,
  SE,
  data = data_plt_tmp,
  studlab = data_plt_tmp$labels,
  comb.fixed = FALSE,
  comb.random = TRUE,
  method.tau = "SJ",
  hakn = TRUE,
  prediction = TRUE
) 

meta_soc <- meta::metagen(
  CBCVA,
  SE,
  data = data_plt_soc_tmp,
  studlab = data_plt_soc_tmp$labels,
  comb.fixed = FALSE,
  comb.random = TRUE,
  method.tau = "SJ",
  hakn = TRUE,
  prediction = TRUE
) 

dat_res <- as.data.frame(meta_all) %>% 
  select(studlab, TE, lower, upper) %>% 
  rename(labels = studlab)

# data frame for results + prediction
dat_pred <- data.frame(
  labels = c(
    "RANDOM EFFECTS MODEL (ALL)", "RANDOM EFFECTS MODEL (q4w 0.5mg)",
    "PREDICTION INTERVAL (ALL)", "PREDICTION INTERVAL (q4w 0.5mg)"
  ),
  TE = c(
    meta_all$TE.random, meta_soc$TE.random,
    meta_all$TE.random, meta_soc$TE.random
  ),
  lower = c(
    meta_all$lower.random, meta_soc$lower.random,
    meta_all$lower.predict, meta_soc$lower.predict 
  ),
  upper = c(
    meta_all$upper.random, meta_soc$upper.random,
    meta_all$upper.predict, meta_soc$upper.predict
  )
)

vir_col <-  viridis(2)

dat_plot <- bind_rows(dat_res %>% add_row(), dat_pred) %>% 
  left_join(data_plt_tmp, by = "labels") %>% 
  mutate(STUDY = coalesce(STUDY, labels)) %>% 
  mutate(ID_TMP = seq_along(labels)) %>% 
  # making q4w 0.5mg arms grey
  mutate(COL = ifelse(TRT == "0.5mgq4w", "gray95", "white")) %>% 
  mutate(COL = ifelse(str_detect(STUDY, "q4w 0.5mg"), "gray95", COL)) %>% 
  # defining color of result bars
  mutate(COL1 = case_when(
    STUDY == "RANDOM EFFECTS MODEL (q4w 0.5mg)" ~ vir_col[1],
    STUDY == "PREDICTION INTERVAL (q4w 0.5mg)" ~ vir_col[1],
    STUDY == "RANDOM EFFECTS MODEL (ALL)" ~ vir_col[2],
    STUDY == "PREDICTION INTERVAL (ALL)" ~ vir_col[2],
    TRUE ~ "black"
  )) %>% 
  # To reverse the order of ID (i.e. put pred int at bottom)
  mutate(ID = max(ID_TMP) - ID_TMP + 1) %>% 
  mutate(ID = as.factor(ID))

tit_df <- data.frame(
  ID = as.numeric(dat_plot$ID[1]) + 1,
  COL = dat_plot$COL[2],
  PHASE = "PHASE",
  N = "N",
  TIME = "DUR",
  TRT = "REGIMEN",
  STUDY = "STUDY",
  SUBLEVDAT = "SUBJ DAT"
)

geom_text_size <- 5.2

p1 <- ggplot(dat_plot, aes(x = TE, y = ID, xmin = lower, xmax = upper)) +
  geom_hline(aes(yintercept = ID, colour = COL), size = 7) + 
  geom_pointrange(shape = 22, aes(colour = COL1), size = 1.1) +
  geom_vline(xintercept = dat_pred$TE[1], color = vir_col[2], lty = 3, lwd = 1.3) +
  geom_vline(xintercept = dat_pred$TE[2], color = vir_col[1], lty = 3, lwd = 1.3) +
  geom_hline(yintercept = 5, lty = 2) +
  theme_classic(base_size = 18) +
  scale_colour_identity() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
  xlab("Mean CBCVA") +
  xlim(-6, 20) +
  coord_cartesian(ylim = c(1, tit_df$ID))

p2 <- ggplot(data = dat_plot, aes(y = ID)) +
  geom_hline(aes(yintercept = ID, colour = COL), size = 7) +
  geom_text(aes(x = 0, label = STUDY), size = geom_text_size, hjust = 0) +
  geom_text(aes(x = 3.5, label = TRT),size = geom_text_size - 0.5,  hjust = 0) +
  geom_text(aes(x = 7, label = N), size = geom_text_size) +
  geom_text(aes(x = 8.5, label = PHASE), size = geom_text_size) +
  geom_text(aes(x = 10, label = TIME), size = geom_text_size, hjust = 1) +
  geom_text(aes(x = 11.5, label = SUBLEVDAT), size = geom_text_size, hjust = 1) +
  geom_text(aes(x = 0, label = STUDY), size = geom_text_size, hjust = 0, fontface = "bold", data = tit_df) +
  geom_text(aes(x = 3.5, label = TRT), size = geom_text_size, hjust = 0, fontface = "bold", data = tit_df) +
  geom_text(aes(x = 7, label = N), size = geom_text_size, fontface = "bold", data = tit_df) +
  geom_text(aes(x = 8.5, label = PHASE), size = geom_text_size, fontface = "bold", data = tit_df) +
  geom_text(aes(x = 10, label = TIME), size = geom_text_size, hjust = 1, fontface = "bold", data = tit_df) +
  geom_text(aes(x = 12, label = SUBLEVDAT), size = geom_text_size, hjust = 1, fontface = "bold", data = tit_df) +
  geom_hline(aes(yintercept = tit_df$ID - 0.5)) +
  scale_colour_identity() +
  theme_void(base_size = 18) + 
  theme(plot.margin = margin(6, 0, 53, 0))  +
  geom_hline(yintercept = 5, lty = 2) +
  coord_cartesian(ylim = c(1, tit_df$ID))

p <- gridExtra::grid.arrange(p2, p1, ncol = 2, widths = c(1.3, 1))

if (save_flag) {
  ggsave(file.path(fig_wd, "Meta.eps"), p, width = 13, height = 9, dpi = "retina")
}

# Baseline BCVA for nAMD studies -----------------------------------------------

# Values taken from OCEAN publication when possible 
# va0_df <- tribble(
#   ~STUDY, ~TRT, ~VA0, ~SD, ~N,
#   "MARINA","0.5mgq4w", 53.7, 12.8, 240,
#   "MARINA","0.3mgq4w", 53.1, 12.9, 238,
#   "ANCHOR","0.5mgq4w", 47.1, 13.0, 140,
#   "ANCHOR","0.3mgq4w", 47, 13.1, 140,
#   "PIER", "0.3mgq12w", 55.8, 12.2, 60, 
#   "PIER", "0.5mgq12w", 53.7, 15.5, 61,
#   "SAILOR", "0.3mgPRN", 55, 12.5, 462, ### Only trt naive
#   "SAILOR", "0.5mgPRN", 48.9, 13.8, 490, ### Only trt naive
#   "SUSTAIN", "0.3/0.5mgPRN", 56.1, 12.19, 512,
#   "EXCITE", "0.3mgq4w", 56.5, 12.9, 115,
#   "EXCITE", "0.3mgq12w", 55.8, 11.8, 120,
#   "EXCITE", "0.5mgq12w", 57.7, 13.06, 118,
#   "CATT", "0.5mgq4w", 59.9, 14.2, 146,
#   "CATT", "0.5mgPRN", 61.6, 13.1, 287,
#   "IVAN", "0.5mgq4w + PRN", 61.8, 15.0, 314,
#   "VIEW1", "0.5mgq4w", 54, 13.4, 304,
#   "VIEW2", "0.5mgq4w", 53.8, 13.5, 291,
#   "HARBOR", "0.5mgq4w", 54.2, 13.3, 275,
#   "HARBOR", "0.5mgPRN", 54.5, 11.7, 275,
#   "HARBOR", "0.3mgq4w", 53.5, 13.1, 274,
#   "HARBOR", "0.3mgPRN", 53.5, 13.2, 273,
#   "MANTA", "0.5mgPRN", 56.4, 13.5, 163,
#   "GEFAL", "0.5mgPRN", 
#   "LUCAS", "0.5mgT&E", 
#   "CEDAR+SEQUIOIA", "0.5mgq4w", 
#   "AVENUE", "0.5mgq4w", 
#   "STAIRWAY", "0.5mgq4w", 
# )

base_df <- tribble(
  ~STUDY,  ~N, ~VA0_MEAN, ~VA0_SD, ~AGE_MEAN, ~AGE_SD,
  "MARINA", 716, 53.5, 13.5, 77, 8,
  "ANCHOR", 423, 46.5, 13.2, 77, 8,
  "PIER", 184, 55, 14, 78.5, 7.5,
  "SAILOR", 2378, 52, 14, 78.7, 8, ## TRT NAIVE ONLY
  "SUSTAIN", 513, 56.1, 12.9, 75.1, 8.06,
  "EXCITE", 353, 56.7, 12.4, 75.3, 7.56,
  "CATT", 1185, 60.5, 13.5, 79, 7.5,
  "IVAN", 561, 61.4, 15.3, 77, 7.4,
  "VIEW1", 1217, 55, 13.25, 78, 8,
  "VIEW2", 1240, 53, 14, 74, 8.5,
  "HARBOR", 1097, 54, 13, 78.5, 8.3,
  "MANTA", 317, 56.75, 13, 77, 8,
  "GEFAL", 501, 55, 14, 79, 7,
  "LUCAS", 441, 61, 13.5, 78.5, 8,
  "HAWK", 1087, 60.6, 13.71, 76.5, 8.68, 
  "HARRIER", 739, 61.2, 12.76, 75.1, 8.24,
  "CEDAR+SEQUIOIA", 1888, 57, 12.5, 75.5, 8.25,
  "AVENUE", 263, 56, 12, 78, 8.5,
  "STAIRWAY", 71, 58.4, 11, 78.5, 8.47
)

# Plot
p1 <- base_df %>% 
  mutate(LOWER = AGE_MEAN - AGE_SD) %>% 
  mutate(UPPER = AGE_MEAN + AGE_SD) %>% 
  rename(AGE = AGE_MEAN) %>% 
  mutate(STUDY = factor(STUDY, levels = unique(STUDY))) %>% 
  ggplot() +
  aes(x = fct_rev(STUDY), y = AGE, ymin = LOWER, ymax = UPPER) +
  geom_pointrange() + 
  #  geom_hline(yintercept = mean(cbcva_rbz_df$CBCVA), lty = 2) +  # add a dotted line at mean
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Arm") +
  ylab("Age +/- SD") +
  theme_bw()  

p2 <- base_df %>% 
  mutate(LOWER = VA0_MEAN - VA0_SD) %>% 
  mutate(UPPER = VA0_MEAN + VA0_SD) %>% 
  rename(VA0 = VA0_MEAN) %>% 
  mutate(STUDY = factor(STUDY, levels = unique(STUDY))) %>% 
  ggplot() +
  aes(x = fct_rev(STUDY), y = VA0, ymin = LOWER, ymax = UPPER) +
  geom_pointrange() + 
  #  geom_hline(yintercept = mean(cVA0_rbz_df$CVA0), lty = 2) +  # add a dotted line at mean
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Arm") +
  ylab("BCVA +/- SD") +
  theme_bw()  

cowplot::plot_grid(p1, p2)


low_df <- base_df %>% 
  mutate(AGE_LOWER = AGE_MEAN - AGE_SD) %>% 
  mutate(VA0_LOWER = VA0_MEAN - VA0_SD) %>% 
  select(-N, -VA0_SD, -AGE_SD, -VA0_MEAN, -AGE_MEAN) %>% 
  pivot_longer(cols = c(AGE_LOWER, VA0_LOWER), values_to = "LOWER", names_to = "COV") %>% 
  mutate(COV = str_remove_all(COV, "_LOWER"))

upp_df <- base_df %>% 
  mutate(AGE_UPPER = AGE_MEAN + AGE_SD) %>% 
  mutate(VA0_UPPER = VA0_MEAN + VA0_SD) %>% 
  select(-N, -VA0_SD, -AGE_SD, -VA0_MEAN, -AGE_MEAN) %>% 
  pivot_longer(cols = c(AGE_UPPER, VA0_UPPER), values_to = "UPPER", names_to = "COV") %>% 
  mutate(COV = str_remove_all(COV, "_UPPER"))

mean_df <-  base_df %>% 
  select(STUDY, VA0_MEAN, AGE_MEAN) %>% 
  pivot_longer(cols = c(VA0_MEAN, AGE_MEAN), values_to = "MEAN", names_to = "COV") %>% 
  mutate(COV = str_remove_all(COV, "_MEAN"))

plot_df <- reduce(list(low_df, upp_df, mean_df), left_join)

plot_df %>% 
  mutate(STUDY = factor(STUDY, levels = unique(STUDY))) %>% 
  ggplot() +
  aes(x = fct_rev(STUDY), y = MEAN, ymin = LOWER, ymax = UPPER) +
  geom_pointrange() + 
  #  geom_hline(yintercept = mean(y), lty = 2) +  # add a dotted line at mean
  coord_flip() +  # flip coordinates (puts labels on y axis) 
  facet_wrap(~ COV, scales = "free_x") +
  xlab("Arm") +
  ylab("MEAN +/- SD") +
  theme_bw(base_size = 18)  

# Model schematic --------------------------------------------------------------

base <- 0.5
gain_trt <- 0.1
ss_trt <- base + gain_trt
ss_sh <- 0.2
t_h_sh <- 125
t_h_trt <- 30


cols <- viridis(2)

mod_df <- tibble(TIME = 0:(365 * 2)) %>% 
  mutate(SHAM = base + (ss_sh - base) * (1 - exp(- log(2) / t_h_sh * TIME))) %>% 
  mutate(TRT  = base + gain_trt * (1 - exp(- log(2) / t_h_trt * TIME)))

seg_df <- tibble(
  ARM = c("SHAM", "TRT"),
  XBEG = c(t_h_sh, t_h_trt),
  XEND = c(t_h_sh, t_h_trt),
  YBEG = c(0, 0),
  YEND = c(
    mod_df[mod_df$TIME == t_h_sh, "SHAM", drop = T],
    mod_df[mod_df$TIME == t_h_trt, "TRT", drop = T]
  )
)

p <- mod_df %>% 
  pivot_longer(cols = c(SHAM, TRT), names_to = "ARM", values_to = "BCVA") %>% 
  ggplot() +
  geom_line(aes(x = TIME, y = BCVA, group = ARM, col = ARM), lwd = 1) +
  geom_segment(
    aes(x = XBEG, xend = XEND, y = YBEG, yend = YEND, group = ARM, col = ARM),
    lty = 2,
    lwd = 1,
    data = seg_df
  ) + 
  annotate(
    "segment", 
    x = t_h_sh,
    xend = t_h_trt, 
    y = 0.1,
    yend = 0.1, 
    arrow = arrow(), 
  ) +
  annotate(
    "segment", 
    x = 700,
    xend = 700, 
    y = mod_df[mod_df$TIME == 700, "SHAM", drop = T],
    yend = mod_df[mod_df$TIME == 700, "TRT", drop = T], 
    arrow = arrow()
  ) +
  annotate(
    "text",
    x = 700,
    y = ss_trt + 0.03,
    label = "SS[i*',T']",
    fontface = 2,
    col = cols[2],
    size = 8.5,
    parse = T
  ) +
  annotate(
    "text",
    x = 700,
    y = ss_sh - 0.03,
    label = "SS[i*',S']",
    fontface = 2,
    col = cols[1],
    size = 8.5,
    parse = T
  ) +
  annotate(
    "text",
    x = 700,
    y = (ss_sh + ss_trt) / 2,
    label = expression(atop("Increasing", paste(~C[i], "(t)"))),
    fontface = 2,
    size = 7.5
  ) +
  annotate(
    "text",
    x = (t_h_sh + t_h_trt) / 2,
    y = 0.1,
    label = expression(atop("Increasing", paste(~C[i], "(t)"))),
    fontface = 2,
    size = 6.5
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.title = element_blank()
  ) +
  ylab(expression(VTS[i])) +
  xlab("Month") +
  scale_x_continuous(
    breaks = c(0, t_h_trt, t_h_sh, seq(200, 800, 100)),
    labels = c(
      NA, 
      expression(over(log(2), paste(k[i*',S'], "(1 -", k[i*',T'], ")"))), 
      expression(over(log(2), k[i*',S'])),
      seq(3, 16, 2)
    )
  ) +
  scale_color_viridis(discrete = T)

if (save_flag) {
  ggsave(file.path(fig_wd, "ModSch.eps"), p, width = 12, height = 9, dpi = "retina")
}

# Plot mapping BCVA to transformed VTS ----------------------------------------

bcva_df <- tibble(BCVA = 0:99) %>% 
  mutate(BCVA_TRANS = 2 - log10(100 - BCVA))

p <- bcva_df %>% 
  ggplot() +
  geom_line(aes(y = BCVA, x = BCVA_TRANS)) +
  xlab(expression(paste("VTS"))) +
  theme_bw(base_size = 24) +
  scale_x_continuous(breaks = seq(0, 2, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) 

if (save_flag) {
  ggsave(file.path(fig_wd, "TransMap.eps"), p, width = 12, height = 9, dpi = "retina")
}


# Pre and post transformation longitudinal plot --------------------------------

p <- data_plot %>%
  filter(!(TRT %in% c("0.5MG_PRN", "2MG_PRN", "SHAM"))) %>% 
  mutate(ARM = paste(STUDY, TRT, sep = " ")) %>% 
  mutate(ARM = str_replace_all(ARM, "_", " ")) %>% 
  group_by(STUDY, TRT, VIS, ARM) %>% 
  mutate(BCVA_TRANS = 2 - log10(100 - VA)) %>%
  rename(BCVA = VA) %>% 
  summarise(
    across(
      c("BCVA", "BCVA_TRANS"),
      ~ c(mean(.), quantile(., probs = c(0.5, 0.1, 0.9)))
    ),
    QUANTILE = c("MEAN", "MEDIAN", "Q10", "Q90")
  ) %>% 
  rename("VTS" = BCVA_TRANS) %>% 
  pivot_longer(
    cols = c("BCVA", "VTS"),
    names_to = "TRANS",
    values_to = "Y"
  ) %>% 
  ungroup() %>% 
  filter(QUANTILE %in% c("MEAN", "MEDIAN", "Q10", "Q90")) %>% 
  mutate(VIS = as.integer(VIS)) %>% 
  mutate(QUANTILE = factor(QUANTILE, levels = c("Q10", "MEDIAN", "MEAN", "Q90"))) %>%
  ggplot() +
  geom_line(aes(x = VIS, y = Y, col = QUANTILE, group = QUANTILE)) +
  geom_point(aes(x = VIS, y = Y, col = QUANTILE, group = QUANTILE)) +
  facet_grid(
    TRANS ~ ARM, 
    scales = "free", 
    labeller = label_wrap_gen(width = 5),
    switch = "y"
  ) + 
  theme_minimal(base_size = 19) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 12, angle = 90),
    strip.placement = "outside"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks()) + 
  labs(col = "Summary statistic") + 
  # ggtitle("Endpoint distribution, after transformations") +
  xlab("Month") +
  ylab("") +
  scale_color_viridis(discrete = T)

if (save_flag) {
  ggsave(file.path(fig_wd, "TransLong.eps"), p, width = 12, height = 9, dpi = "retina")
}

library(checkmate)

# Function to recreate default ggplot colors -----------------------------------

gg_col_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Vectorized sequence function -------------------------------------------------

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"), SIMPLIFY = F)

# BCVA to VTS and inverse transformation functions -----------------------------

trans_fun <- function(x) {2 - log10(100 - x)}
inv_trans_fun <- function(x) {100 - 10 ^ (2 - x)}

# Function to return Monolix data from proj directory --------------------------

# Returns monolix data (or directory)
# from the project location
read_data_project <- function(wd_proj, ret) {
  # Arguments check
  assert_character(wd_proj, len = 1, null.ok = F)
  assert_character(ret, len = 1, null.ok = F)
  assert_choice(ret, c("data", "dir"))
  
  # Grepping mlxtran file
  files <- list.files(file.path(wd_proj, ".Internals"))
  mod_file <- files[grepl("mlxtran", files)]
  mlx_wd <- file.path(wd_proj, ".Internals")
  # reading in mlxtran file
  mlx_file <- readLines(file.path(mlx_wd, mod_file))
  # reducing mlx_file to [FILEINFO] section where datainfo stored
  mlx_begin <- which(str_detect(mlx_file, "\\[FILEINFO\\]"))
  mlx_sec <- which(str_detect(mlx_file, "\\["))
  mlx_end <- mlx_sec[which(mlx_sec == mlx_begin) + 1] - 1
  mlx_file_info <- mlx_file[mlx_begin:mlx_end]
  # Extracting line where data directory stored
  file_str <- mlx_file_info[str_which(mlx_file_info, "^file")]
  # Only keeping part of line which are related to working directory
  data_wd <- str_remove_all(file_str, "file")
  data_wd <- str_remove_all(data_wd, "=")
  data_wd <- str_remove_all(data_wd, "\\'")
  data_wd <- str_remove_all(data_wd, "[:blank:]")
  dat_wd <- file.path(mlx_wd, data_wd)
  
  out <- NA
  if (ret == "data") {
    out <- data.table::fread(dat_wd)
  } else {
    out <- dat_wd
  }
  out
}

# Function to read pop-pars ----------------------------------------------------

mlx_pop_pars <- function(wd_proj) {
  p <- file.path(wd_proj, "PopulationParameters.txt")
  read.csv(p)
}

# Function to return important results file from Monolix -----------------------
# (these can be used in R for plotting)

# Return all helpful datasets for doing diagnostic monolix plots
read_mlx <- function(wd_proj, read_vpc = T, trans = T) {
  assert_character(wd_proj, len = 1, null.ok = F)
  assert_logical(read_vpc, len = 1, null.ok = F)
  
  res <- list()
  
  # Working directories
  preds_wd <- file.path(wd_proj, "predictions.txt")
  resid_wd <- file.path(wd_proj, "ChartsData/ScatterPlotOfTheResiduals/Y_residuals.txt")
  vpc_wd <- file.path(wd_proj, "ChartsData/VisualPredictiveCheck/Y_simulations.txt")
  
  # data
  print("Reading data")
  data_read <- read_data_project(wd_proj, ret = "data")
  data <- data_read %>% 
    filter(!is.na(Y)) %>% 
    select(ID, VIS, TIME, Y, CVA, X_VA0, starts_with("CAT"), starts_with("COV")) %>% 
    mutate(COV_VA0 = COV_VA0 * 100) %>% 
    mutate(CAT_VA0 = cut(COV_VA0, breaks = c(-1, 40, 50, 60, 70, 101))) 
  
  # IPRED and PRED, PWRES and IWRES
  print("Reading residuals")
  if (file.exists(resid_wd)) {
    imp_cols <- c("ID", "time", "prediction_pwRes", "pwRes", "prediction_iwRes_mode", "iwRes_mode")
    resid_read <- data.table::fread(
      resid_wd, 
      select = imp_cols, 
      colClasses = list(numeric = 1:length(imp_cols)),
      sep =','
    )
    
    resid_tmp <- resid_read %>% 
      rename(
        TIME = time,
        PWRES = pwRes,
        IWRES = iwRes_mode,
      ) %>% 
      left_join(data, by = c("ID", "TIME"))
    
    if (trans) {
      resid <- resid_tmp %>% 
        rename(
          VA_TRANS_PRED = prediction_pwRes,
          VA_TRANS_IPRED = prediction_iwRes_mode,
          VA_TRANS = Y,
          CVA_TRANS = CVA
        ) %>% 
        mutate(
          CVA_TRANS_PRED = VA_TRANS_PRED - X_VA0,
          CVA_TRANS_IPRED = VA_TRANS_IPRED - X_VA0,
          VA_PRED = inv_trans_fun(VA_TRANS_PRED),
          VA_IPRED = inv_trans_fun(VA_TRANS_IPRED),
          CVA_PRED = VA_PRED - inv_trans_fun(X_VA0),
          CVA_IPRED = VA_IPRED - inv_trans_fun(X_VA0),
          VA = inv_trans_fun(VA_TRANS),
          CVA = VA - inv_trans_fun(X_VA0)
        ) %>% 
        select(
          ID, TIME, VIS, X_VA0, VA, CVA, VA_TRANS, CVA_TRANS,
          starts_with("CAT_"), starts_with("COV_"),
          VA_PRED, VA_TRANS_PRED, CVA_PRED, CVA_TRANS_PRED,
          VA_IPRED, VA_TRANS_IPRED, CVA_IPRED, CVA_TRANS_IPRED,
          PWRES, IWRES
        )
    } else {
      resid <- resid_tmp %>% rename(
        VA_PRED = prediction_pwRes,
        VA_IPRED = prediction_iwRes_mode,
        VA = Y
      ) %>% 
        mutate(
          CVA_PRED = VA_PRED - X_VA0,
          CVA_IPRED = VA_IPRED - X_VA0
        ) %>% 
        select(
          ID, TIME, VIS, X_VA0, VA, CVA, 
          starts_with("CAT_"), starts_with("COV_"),
          VA_PRED, CVA_PRED, VA_IPRED, CVA_IPRED, PWRES, IWRES
        )
    }
    res$resid <- resid
  } else {
    warning("residual file not found... did you save charts?")
    res$resid <- NULL
  }
  
  # VPC
  print("Reading VPC")
  if (file.exists(resid_wd) & read_vpc) {
    imp_cols <- c("rep", "ID", "time", "sim_Y")
    vpc_read <- data.table::fread(
      vpc_wd, 
      select = imp_cols, 
      colClasses = list(numeric = 1:length(imp_cols)),
      sep =','
    )  
    
    vpc_tmp <- vpc_read %>% 
      rename(TIME = time) %>% 
      left_join(data, by = c("ID", "TIME"))
    
    if (trans) {
      vpc <- vpc_tmp %>% 
        rename(SIM_VA_TRANS = sim_Y) %>% 
        mutate(SIM_CVA_TRANS = SIM_VA_TRANS - X_VA0) %>% 
        rename(VA_TRANS = Y, CVA_TRANS = CVA) %>% 
        mutate(
          SIM_VA = inv_trans_fun(SIM_VA_TRANS),
          SIM_CVA = SIM_VA - inv_trans_fun(X_VA0),
          VA = inv_trans_fun(VA_TRANS),
          CVA = VA - inv_trans_fun(X_VA0)
        )
    } else {
      vpc <- vpc_tmp %>% 
        rename(SIM_VA = sim_Y) %>% 
        mutate(SIM_CVA = SIM_VA - X_VA0) %>% 
        rename(VA = Y) 
    }
    res$vpc <- vpc
  } else {
    warning("VPC file not found... did you export VPC data?")
    res$vpc <- NULL
  }
  
  class(res) <- "read_mlx"
  res
}


# Goodness-of-fit plots --------------------------------------------------------

plot_dv_pred <- function(obj, 
                         out = "VA", 
                         type = "PRED", 
                         trans = F,
                         vis_filter = NULL,
                         subsamp = NULL) {
  
  # Argument checks
  assert_class(obj, "read_mlx")
  assert_character(type, len = 1, null.ok = F)
  assert_logical(trans, len = 1, null.ok = F)
  assert_choice(out, c("VA", "CVA"))
  assert_choice(type, c("PRED", "IPRED"))
  assert_double(subsamp, null.ok = T, lower = 0.01, upper = 0.99)
  
  # x and y variables based on arguments
  trans_txt <- ifelse(trans == T, "_TRANS", "")
  obs_var <- sym(paste0(out, trans_txt))
  pred_var <- sym(paste0(out, trans_txt, "_", type))
  
  # ylabel
  ylab <- ifelse(out == "VA", "BCVA", out)
  ylab <- ifelse(trans, "VTS", ylab)
  
  plot_df <- obj$resid
  if (!is.null(vis_filter)) {
    all_vis <- unique(plot_df$VIS)
    if (! all(vis_filter %in% all_vis)) {
      txt <- paste(
        "vis_filter can only be a subset of: ,\n" ,
        paste(all_vis, collapse = ", ")
      )
      stop(txt)
    }
  }
  
  # Subsampling, can be useful to reduce size when too many points
  if (!is.null(subsamp)) {
    plot_df <- sample_frac(plot_df, subsamp)
  }
  
  # Plot
  plot_df %>% 
    ggplot() + 
    geom_point(aes(x = {{pred_var}}, y = {{obs_var}}), col = "gray50") +
    geom_abline(intercept = 0, slope = 1, lty = 2, lwd = 1.3) +
    geom_smooth(aes(x = {{pred_var}}, y = {{obs_var}}), col = "red", se = FALSE) +
    ggpmisc::stat_poly_eq(
      formula = "y ~ x",
      aes(x = {{pred_var}}, y = {{obs_var}}, label = paste(..rr.label.., sep = "~~~")),
      parse = TRUE
    ) +
    theme_bw(base_size = 20) +
    xlab(paste0(type, " [", ylab,"]")) +
    ylab(ylab)
}

plot_resid_scatter <- function(obj, type = "PWRES", x_ax = "PRED", subsamp = NULL) {
  
  # Argument checks
  assert_class(obj, "read_mlx")
  assert_character(type, len = 1, null.ok = F)
  assert_character(x_ax, len = 1, null.ok = F)
  assert_choice(type, c("PWRES", "IWRES"))
  x_ax_ok <- c(
    "PRED", 
    "TIME", 
    obj$resid %>% select(starts_with("CAT"), starts_with("COV")) %>% colnames()
  )
  assert_choice(x_ax, x_ax_ok)
  assert_double(subsamp, null.ok = T, lower = 0.01, upper = 0.99)
  
  # x and y variables based on arguments
  x_ax_arg <- x_ax
  if (x_ax == "PRED") {
    x_ax <- ifelse(type == "PWRES", "VA_TRANS_PRED", "VA_TRANS_IPRED")
    x_lab_tmp <- ifelse(type == "PWRES", "PRED [VTS]", "IPRED [VTS]")
  }
  x_ax_sym <- sym(x_ax)
  y_sym <- sym(type)
  
  # x label
  x_lab <- ifelse(x_ax_arg == "PRED", x_lab_tmp, paste(x_ax, "[DAYS]"))
  
  plot_df <- obj$resid
  # Subsampling, can be useful to reduce size when too many points
  if (!is.null(subsamp)) {
    plot_df <- sample_frac(plot_df, subsamp)
  }
  
  # Plot
  plot_df %>% 
    ggplot() +
    geom_point(aes(x = {{x_ax_sym}}, y = {{y_sym}}), col = "gray50") +
    geom_smooth(aes(x = {{x_ax_sym}}, y = {{y_sym}}), col = "red", se = FALSE) +
    geom_abline(intercept = 0, slope = 0, lty = 2, lwd = 1.4) +
    geom_hline(yintercept = 2, lty = 2, lwd = 1.05) +
    geom_hline(yintercept = -2, lty = 2, lwd = 1.05) +
    theme_bw(base_size = 20) +
    xlab(x_lab)
}

plot_resid_qq <- function(obj, type = "PWRES", subsamp = NULL) {
  
  # Argument checks
  assert_class(obj, "read_mlx")
  assert_character(type, len = 1, null.ok = F)
  assert_choice(type, c("PWRES", "IWRES"))
  assert_double(subsamp, null.ok = T, lower = 0.01, upper = 0.99)
  
  type_sym <- sym(type)
  
  plot_df <- obj$resid
  # Subsampling, can be useful to reduce size when too many points
  if (!is.null(subsamp)) {
    plot_df <- sample_frac(plot_df, subsamp)
  }
  
  # Plot
  plot_df %>% 
    ggplot(aes(sample = !!type_sym)) + 
    stat_qq() + 
    stat_qq_line() +
    theme_bw(base_size = 20) +
    xlab("Theoretical Normal Scores") +
    ylab(type)
}

plot_resid_hist <- function(obj, type = "PWRES", subsamp = NULL) {
  
  # Argument checks
  assert_class(obj, "read_mlx")
  assert_character(type, len = 1, null.ok = F)
  assert_choice(type, c("PWRES", "IWRES")) 
  assert_double(subsamp, null.ok = T, lower = 0.01, upper = 0.99)
  
  type_sym <- sym(type)
  
  plot_df <- obj$resid
  # Subsampling, can be useful to reduce size when too many points
  if (!is.null(subsamp)) {
    plot_df <- sample_frac(plot_df, subsamp)
  }
  
  # Plot
  plot_df %>% 
    ggplot(aes(x = !!type_sym)) + 
    geom_histogram() +
    theme_bw(base_size = 20) +
    ylab("Frequency")
}

# VPC (based on general VPC plotting) 
plot_mlx_vpc <- function(obj, out, trans = F, strat = NULL, bw = F) {
  assert_class(obj, "read_mlx")
  # strat is safety checked later in calculate vpc_quant
  assert_logical(trans, len = 1)
  assert_character(out, len = 1)
  assert_choice(out, c("VA", "CVA"))
  
  trans_txt <- ifelse(trans == T, "_TRANS", "")
  y_var <- paste0(out, trans_txt)
  sim_y_var <- paste0("SIM_", y_var)
  
  sim_cols <- c("rep", "ID", "VIS", sim_y_var, strat)
  dat_cols <- c("ID", "VIS", y_var, strat)
  
  vpc_const <- calculate_vpc_quant(
    sim_df = obj$vpc %>% 
      select(all_of(sim_cols)) %>%
      rename({{y_var}} := {{sim_y_var}}),
    data_df = obj$vpc %>% 
      filter(rep == 1) %>% 
      select(all_of(dat_cols)),
    out_col = y_var,
    rep_col = "rep", 
    time_col = "VIS", 
    strat_cols = strat
  ) 
  
  # Plotting
  p <- NULL
  if (bw) {
    p <- plot_vpc_const_bw1(vpc_const)  
  } else {
    p <- plot_vpc_const(vpc_const)  
  }
  p
}

# General VPC plotting ---------------------------------------------------------

# Would be nice to use data.table here for efficiency reasons
calculate_vpc_quant <- function(sim_df,
                                data_df,
                                out_col,
                                rep_col = "rep", 
                                time_col, 
                                pi_size = 0.9,
                                strat_cols = NULL) {
  # Arguments checks
  assert_character(rep_col, len = 1)
  assert_character(time_col, len = 1)
  assert_character(strat_cols, max.len = 2, null.ok = T)
  cols_sim_df <- c(out_col, rep_col, time_col, strat_cols)
  cols_data_df <- c(out_col, time_col, strat_cols)
  assert_data_frame(sim_df, min.cols = length(cols_sim_df))
  assert_data_frame(data_df, min.cols = length(cols_data_df), null.ok = T)
  assert_names(colnames(sim_df), must.include = cols_sim_df)
  sim_df <- as.data.frame(sim_df)
  class_strat_sim_df <- sapply(sim_df, class)[strat_cols]
  ok_class <- c("character", "integer", "factor")
  if (! all(class_strat_sim_df %in% ok_class)) {
    txt <- paste(
      "Class of columns:", 
      paste(class_strat_sim_df, collapse = ", "),
      "in sim_df shall be one of:",
      paste(ok_class, collapse = ", ")
    )
    stop(txt)
  }
  # Checks on data_df only if data_df supplied
  if (!is.null(data_df)) {
    assert_names(colnames(data_df), must.include = cols_data_df)
    data_df <- as.data.frame(data_df)
    class_strat_data_df <- sapply(data_df[, strat_cols], class)
    if (! all(class_strat_data_df %in% ok_class)) {
      txt <- paste(
        "Class of columns:", 
        paste(class_strat_data_df, collapse = ", "),
        "in data_df shall be one of:",
        paste(ok_class, collapse = ", ")
      )
      stop(txt)
    }
  }
  
  # Function to calculate 90 quantiles
  q_low <- (1 - pi_size) / 2
  q_high <- 1 - q_low
  q_l <- list(LWR = ~ quantile(.x, q_low), UPR = ~ quantile(.x, q_high))
  
  # Creating stratification factors for VPC
  all_cols <- c(out_col, rep_col, time_col, strat_cols)
  grp_one <- syms(all_cols[-1])
  grp_two <- syms(all_cols[-c(1, 2)])
  out_sym <- sym(out_col)
  
  # Simulations
  vpc_simu <- sim_df %>% 
    group_by(!!!grp_one) %>% 
    summarise(
      MEAN = mean(!!out_sym, na.rm = T),
      MEDIAN = median(!!out_sym, na.rm = T),
      Q10 = quantile(!!out_sym, 0.1, na.rm = T),
      Q90 = quantile(!!out_sym, 0.9, na.rm = T)
    ) %>% 
    ungroup() %>% 
    group_by(!!!grp_two) %>% 
    summarise(
      across(
        .cols = c(MEAN, MEDIAN, Q10, Q90),
        .fns = q_l,
        .names = "{.col}_{.fn}"
      )
    )
  
  # Actual observations
  if (is.null(data_df)) {
    vpc_dat <- NULL
  } else {
    vpc_dat <- data_df %>% 
      group_by(!!!grp_two) %>% 
      summarise(
        MEAN = mean(!!out_sym, na.rm = T),
        MEDIAN = median(!!out_sym, na.rm = T),
        Q10 = quantile(!!out_sym, 0.1, na.rm = T),
        Q90 = quantile(!!out_sym, 0.9, na.rm = T)
      ) %>% 
      pivot_longer(
        cols = c("MEAN", "MEDIAN", "Q10", "Q90"),
        names_to = "QUANTILE",
        values_to = "Y"
      )
  }
  
  # Output and return
  out <- list(
    sim = vpc_simu, 
    data = vpc_dat,
    out_col = out_col,
    time_col = time_col,
    strat_cols = strat_cols
  )
  class(out) <- c("vpc_const")
  out
}

plot_vpc_const <- function(obj) {
  assert_class(obj, "vpc_const")
  
  # Declaring symbols
  time_sym <- sym(obj$time_col)
  out_sym <- sym(obj$out_col)
  
  cols <- gg_col_hue(4)
  
  # Plot
  p <- obj$sim %>% 
    mutate(!!time_sym := as.integer(!!time_sym)) %>% 
    ggplot() +
    geom_ribbon(
      aes(x = !!time_sym, ymin = MEAN_LWR, ymax = MEAN_UPR), 
      fill = cols[1], 
      alpha = 0.4
    ) +
    geom_ribbon(
      aes(x = !!time_sym, ymin = MEDIAN_LWR, ymax = MEDIAN_UPR), 
      fill = cols[2], 
      alpha = 0.4
    ) +
    geom_ribbon(
      aes(x = !!time_sym, ymin = Q10_LWR, ymax = Q10_UPR),
      fill = cols[3],
      alpha = 0.4
    ) +
    geom_ribbon(
      aes(x = !!time_sym, ymin = Q90_LWR, ymax = Q90_UPR), 
      fill = cols[4],
      alpha = 0.4
    ) +
    ylab(obj$out_col) +
    scale_x_continuous(breaks = scales::pretty_breaks()) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) + 
    theme_bw(base_size = 20) +
    theme(legend.position = "bottom")
  
  if (!is.null(obj$data)) {
    # Data only plotted if supplied via constructor
    p <- p + 
      geom_point(
        aes(x = !!time_sym, y = Y, group = QUANTILE, col = QUANTILE), 
        lwd = 2, 
        data = obj$data
      ) +
      geom_line(
        aes(x = !!time_sym, y = Y, group = QUANTILE, col = QUANTILE),
        lwd = 1.1,
        data = obj$data
      ) 
  }
  
  if (!is.null(obj$strat_cols)) {
    # Stratification factors only plotted if supplied via constructor
    strat <- obj$strat_cols
    f_facet <- switch(
      as.character(length(strat)),
      "1" = as.formula(paste("~", strat)),
      "2" = as.formula(paste(strat[1], "~", strat[2])),
      stop("Can not handle more than two stratification variables")
    )
    p <- p + facet_wrap(f_facet, scales = "free") 
  }
  p
}


plot_vpc_const_bw <- function(obj) {
  assert_class(obj, "vpc_const")
  
  # Declaring symbols
  time_sym <- sym(obj$time_col)
  out_sym <- sym(obj$out_col)
  
  cols <- gray.colors(4, start = 0, end = 0.5, gamma = 2.2, alpha = 0.8, rev = FALSE)
  llwd <- 1.2
  
  # Plot
  p <- obj$sim %>% 
    mutate(!!time_sym := as.integer(!!time_sym)) %>% 
    ggplot() +
    # geom_line(aes(x = !!time_sym, y = MEAN_LWR), lty = 2, col = cols[1], lwd = llwd) +
    # geom_line(aes(x = !!time_sym, y = MEAN_UPR), lty = 2, col = cols[1], lwd = llwd) +
    # geom_line(aes(x = !!time_sym, y = MEDIAN_LWR), lty = 2, col = cols[4], lwd = llwd) +
    # geom_line(aes(x = !!time_sym, y = MEDIAN_UPR), lty = 2, col = cols[4], lwd = llwd) +
    # geom_line(aes(x = !!time_sym, y = Q10_LWR), lty = 2, col = cols[1], lwd = llwd) +
    # geom_line(aes(x = !!time_sym, y = Q10_UPR), lty = 2, col = cols[1], lwd = llwd) +
    # geom_line(aes(x = !!time_sym, y = Q90_LWR), lty = 2, col = cols[1], lwd = llwd) +
    # geom_line(aes(x = !!time_sym, y = Q90_UPR), lty = 2, col = cols[1], lwd = llwd) +
    geom_ribbon(
      aes(x = !!time_sym, ymin = MEAN_LWR, ymax = MEAN_UPR),
      fill = cols[1],
      alpha = 0.3
    ) +
    geom_ribbon(
      aes(x = !!time_sym, ymin = MEDIAN_LWR, ymax = MEDIAN_UPR),
      fill = cols[4],
      alpha = 0.3
    ) +
    geom_ribbon(
      aes(x = !!time_sym, ymin = Q10_LWR, ymax = Q10_UPR),
      fill = cols[2],
      alpha = 0.3
    ) +
    geom_ribbon(
      aes(x = !!time_sym, ymin = Q90_LWR, ymax = Q90_UPR),
      fill = cols[2],
      alpha = 0.3
    ) +
    ylab(obj$out_col) +
    scale_x_continuous(breaks = scales::pretty_breaks()) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) + 
    theme_bw(base_size = 18) +
    theme(legend.position = "bottom")
  
  if (!is.null(obj$data)) {
    # Data only plotted if supplied via constructor
    p <- p + 
      geom_point(
        aes(x = !!time_sym, y = Y), lwd = 2, alpha = 0.6, pch = 1,
        data = obj$data %>% filter(QUANTILE == "MEAN"),
      )  +
      geom_point(
        aes(x = !!time_sym, y = Y), lwd = 2, alpha = 0.6, pch = 2, 
        data = obj$data %>% filter(QUANTILE == "MEDIAN"),
      )  +
      geom_point(
        aes(x = !!time_sym, y = Y), lwd = 2, alpha = 0.6, pch = 1,
        data = obj$data %>% filter(QUANTILE == "Q10"),
      )  +
      geom_point(
        aes(x = !!time_sym, y = Y), lwd = 2, alpha = 0.6, pch = 1,
        data = obj$data %>% filter(QUANTILE == "Q90"),
      )  +
      scale_color_manual(values = cols)
  }
  
  if (!is.null(obj$strat_cols)) {
    # Stratification factors only plotted if supplied via constructor
    strat <- obj$strat_cols
    f_facet <- switch(
      as.character(length(strat)),
      "1" = as.formula(paste("~", strat)),
      "2" = as.formula(paste(strat[1], "~", strat[2])),
      stop("Can not handle more than two stratification variables")
    )
    p <- p + facet_wrap(f_facet, scales = "free") 
  }
  p
}

plot_vpc_const_bw1 <- function(obj) {
  assert_class(obj, "vpc_const")
  
  # Declaring symbols
  time_sym <- sym(obj$time_col)
  out_sym <- sym(obj$out_col)
  need_strat <- !is.null(obj$strat_cols)
  
  # Data handling
  sim_df <- obj$sim %>% 
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
  
  dat_df <- obj$data 
  
  # Preparing columns to stratify plot with
  if (need_strat) {
   strat_l <- length(obj$strat_cols)
   strat_fm <- switch(
     as.character(strat_l),
     "1" = as.character(obj$strat_cols[1]),
     "2" = paste0("paste(", obj$strat_cols[1], ", ", obj$strat_cols[2], ")"),
     stop("Can not handle more than two stratification variables")
   )
   
   sim_df <- sim_df %>% 
     mutate(STRAT = !!rlang::parse_expr(strat_fm)) %>% 
     mutate(STRAT = str_replace_all(STRAT, "_", " "))
   dat_df <- dat_df %>% 
     mutate(STRAT = !!rlang::parse_expr(strat_fm)) %>% 
     mutate(STRAT = str_replace_all(STRAT, "_", " "))
  }

  # Plot
  p <- sim_df %>% 
    mutate(!!time_sym := as.integer(!!time_sym)) %>% 
    ggplot() +
    aes(x = !!time_sym) +
    geom_ribbon(aes(ymin = LWR, ymax = UPR), fill = "gray80") +
    facet_grid(
      fct_relevel(QUANTILE,'Q90','MEDIAN','MEAN','Q10') ~ STRAT, 
      scales = "free_y",
      labeller = label_wrap_gen(width = 5)
    ) +
    geom_point(aes(y = Y), data = dat_df) +
    geom_line(aes(y = Y), data = dat_df) +
    scale_x_continuous(breaks = scales::pretty_breaks()) + 
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ylab(obj$out_col) 
    
  p
}


# Function to simulate one individual ------------------------------------------

# Function to simulate one individual, given:
# - simulation times (time_sim)
# - baseline covariates and regressor (reg cannot be time dept)
# - treatment data-frame containing two cols:
#  `time` at which dose administered
#  `amount`  at which dose administered
# - reff if random effect variability to be included
# - res_err if to include residual error
# - n_sim number of simulations
# - group_nm: name of group (useful for plotting)
sim_one_ind <- function(time_sim,
                        cov_reg,
                        trt_df,
                        reff,
                        res_err,
                        n_sim = 1, 
                        group_nm = "Sim") {
  
  # Temporary directory where smlx stores stuff
  temp_wd <- "C:/Users/brizzif/lixoft/simulx/simulx2020R1/tmp/Untitled\\ExternalFiles\\"
  must_covs <- read.csv(file.path(temp_wd, "mlx_Cov.txt"), nrows = 1, header = T) %>% 
    select(-id) %>% 
    colnames()
  must_reg <- read.csv(file.path(temp_wd, "mlx_Reg.txt"), nrows = 1, header = T) %>% 
    select(-id, -time) %>% 
    colnames()
  must_trt_covs <- c(must_covs, must_reg)
  
  # Argument checks
  assert_numeric(time_sim)
  assert_data_frame(cov_reg, ncols = length(must_trt_covs))
  assert_data_frame(trt_df, ncols = 2)
  assert_names(colnames(cov_reg), permutation.of = must_trt_covs)
  assert_names(colnames(trt_df), permutation.of = c("time", "amount"))
  assert_logical(reff, len = 1, null.ok = F)
  assert_integerish(n_sim, len = 1, lower = 1)
  assert_character(group_nm, len = 1)
  
  # Population parameters
  par_mlx <- getPopulationElements()
  # Set random effects to zero if needed
  if (reff == FALSE) {
    par_mlx$mlx_Pop$data[1, grep("omega", colnames(par_mlx$mlx_Pop$data))] <- 0
  }
  # Pass pop pars to smlx
  definePopulationElement(
    name = "par_pop",
    element = par_mlx$mlx_Pop$data[, -1, drop = F]
  )
  
  # Simulation output
  # Include res error if needed
  out_str <- ifelse(res_err == T, "Y", "VA")
  t_max <- tail(time_sim, 1)
  out_df <- data.frame(time = time_sim)
  # Pass to smlx
  defineOutputElement(
    name = "out",
    element = list(data = out_df, output = out_str)
  )
  
  # Treatment (read from trt_df)
  defineTreatmentElement(
    "trt", 
    element = list(admID = 1, probaMissDose = 0, data = trt_df)
  )
  
  # Covariates
  defineCovariateElement(name = "cov", element = cov_reg[, must_covs, drop = F])
  
  # Regressor
  defineRegressorElement(
    name = "reg", 
    element = data.frame(time = time_sim, cov_reg[, must_reg, drop = F])
  )
  
  # Simulation groups
  # Remove old simulation groups (safety measure)
  grp_mlx <- sapply(getGroups(), function(x) x$name)
  lapply(grp_mlx, function(x) removeGroup(x))
  
  addGroup(group_nm)
  setGroupElement(
    group = group_nm, elements =  c("par_pop", "reg", "cov", "out", "trt")
  )
  setGroupSize(group = group_nm, size = 1)
  setNbReplicates(n_sim)
  
  # Run simulation
  tictoc::tic()
  runSimulation()
  tictoc::toc()
  gc()
  
  # Output
  out_txt <- paste0("getSimulationResults()$res$", out_str)
  res_smx <- eval(parse(text = out_txt)) %>% 
    rename(VA_TRANS = all_of(out_str)) %>% 
    mutate(CVA_TRANS = VA_TRANS - cov_reg[1, "X_VA0"]) %>% 
    mutate(VA = inv_trans_fun(VA_TRANS)) %>% 
    mutate(CVA = VA - inv_trans_fun(cov_reg[1, "X_VA0"])) %>% 
    mutate(group = group_nm) %>% 
    as_tibble()
}

# Function to clean-up simulation results --------------------------------------

#' Function to tidy-up results 
#' 
#' @param results object from simulations obtained via the `simulate_new_inds`
#' function in `prior_sim.R`
tidy_res <- function(res, make_naive = F) {
  
  # Get nice results
  out <- res %>% 
    rename(
      ID = id,
      TIME = time,
      SIM_VA_TRANS = sim_VA_TRANS
    ) %>%
    mutate(
      VIS = TIME / 28,
      SIM_CVA_TRANS = SIM_VA_TRANS - X_VA0,
      SIM_VA = inv_trans_fun(SIM_VA_TRANS),
      X_VA0 = inv_trans_fun(X_VA0),
      SIM_CVA = SIM_VA - X_VA0,
    ) %>%
    #    mutate(CAT_ARM = paste(CAT_STD, CAT_TRT)) %>% 
    # Return tibble
    as_tibble() 
  
  if (make_naive) {
    out <- out %>% 
      mutate(VIS_BASE = PRE_VIS + 1) %>% 
      filter(VIS >= VIS_BASE) %>% 
      mutate(VIS = VIS - VIS_BASE) %>% 
      group_by(ID, rep) %>% 
      mutate(TIME = TIME - TIME[1]) %>% 
      mutate(X_VA0 = SIM_VA[1]) %>% 
      mutate(SIM_CVA = SIM_VA - X_VA0) %>% 
      ungroup()
  }
  
  # Return
  out
}
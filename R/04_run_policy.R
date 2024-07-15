##### Load Libraries #####
library(doParallel)
library(doRNG)
library(doSNOW)
library(foreach)
library(dtplyr)
library(dplyr)
library(readr)
library(rootSolve)
library(tidyverse)
library(APCtools)
library(data.table)

source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

##### Functions #####
foi_transition_matrix_assort <- function(v_betas, w = w2) {
  v_betas_1 <- v_betas[1:8]
  v_betas_2 <- v_betas[9:16]
  v_betas_3 <- v_betas[17:24]
  v_betas_4 <- v_betas[25:32]

  beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                   breaks = waifw_breaks,
                                   w = w,
                                   upper = TRUE)
  beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                   breaks = waifw_breaks,
                                   w = w,
                                   upper = TRUE)
  beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                   breaks = waifw_breaks,
                                   w = w,
                                   upper = TRUE)
  beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                   breaks = waifw_breaks,
                                   w = w,
                                   upper = TRUE)
  beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                   breaks = waifw_breaks,
                                   w = w,
                                   upper = TRUE)
  Beta <- rbind(cbind(beta1, beta0, beta0, beta0),
                cbind(beta0, beta2, beta0, beta0),
                cbind(beta0, beta0, beta3, beta0),
                cbind(beta0, beta0, beta0, beta4))
  return(Beta)
}

foi_transition_matrix_rp <- function(v_betas, w = w2) {
  v_betas_1  <- v_betas[1:8]
  v_betas_2  <- v_betas[9:16]
  v_betas_3  <- v_betas[17:24]
  v_betas_4  <- v_betas[25:32]
  v_betas_5  <- v_betas[33:40]
  v_betas_6  <- v_betas[41:48]
  v_betas_7  <- v_betas[49:56]
  v_betas_8  <- v_betas[57:64]
  v_betas_9  <- v_betas[65:72]
  v_betas_10 <- v_betas[73:80]

  beta1  <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta2  <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta3  <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta4  <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta5  <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta6  <- get_transmission_matrix(betas = v_betas_6,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta7  <- get_transmission_matrix(betas = v_betas_7,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta8  <- get_transmission_matrix(betas = v_betas_8,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta9  <- get_transmission_matrix(betas = v_betas_9,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta10 <- get_transmission_matrix(betas = v_betas_10,  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)

  Beta <- rbind(cbind(beta1, beta5, beta6,  beta7),
                cbind(beta5, beta2, beta8,  beta9),
                cbind(beta6, beta8, beta3,  beta10),
                cbind(beta7, beta9, beta10, beta4))
  return(Beta)
}

alpha_matrix_within <- function(v_alpha, n_ages = 80) {
  # Assign alpha values
  alpha_11  <- v_alpha[1]
  alpha_22  <- v_alpha[2]
  alpha_33  <- v_alpha[3]
  alpha_44  <- v_alpha[4]

  # Create matrix with alpha values
  m_alpha_11 <- matrix(data = rep(alpha_11, n_ages^2), nrow = n_ages)
  m_alpha_22 <- matrix(data = rep(alpha_22, n_ages^2), nrow = n_ages)
  m_alpha_33 <- matrix(data = rep(alpha_33, n_ages^2), nrow = n_ages)
  m_alpha_44 <- matrix(data = rep(alpha_44, n_ages^2), nrow = n_ages)
  m_zero     <- matrix(data = rep(0, n_ages^2), nrow = n_ages)

  # Create matrix with alpha values
  m_alpha_within <- rbind(cbind(m_alpha_11, m_zero, m_zero, m_zero),
                          cbind(m_zero, m_alpha_22, m_zero, m_zero),
                          cbind(m_zero, m_zero, m_alpha_33, m_zero),
                          cbind(m_zero, m_zero, m_zero, m_alpha_44))
  return(m_alpha_within)
}

alpha_matrix_between <- function(v_alpha, n_ages = 80) {
  # Assign alpha values
  alpha_12  <- v_alpha[1]
  alpha_13  <- v_alpha[2]
  alpha_14  <- v_alpha[3]
  alpha_21  <- v_alpha[4]
  alpha_23  <- v_alpha[5]
  alpha_24  <- v_alpha[6]
  alpha_31  <- v_alpha[7]
  alpha_32  <- v_alpha[8]
  alpha_34  <- v_alpha[9]
  alpha_41  <- v_alpha[10]
  alpha_42  <- v_alpha[11]
  alpha_43  <- v_alpha[12]

  # Create matrix with alpha values
  m_alpha_12 <- matrix(data = rep(alpha_12, n_ages^2), nrow = n_ages)
  m_alpha_13 <- matrix(data = rep(alpha_13, n_ages^2), nrow = n_ages)
  m_alpha_14 <- matrix(data = rep(alpha_14, n_ages^2), nrow = n_ages)
  m_alpha_21 <- matrix(data = rep(alpha_21, n_ages^2), nrow = n_ages)
  m_alpha_23 <- matrix(data = rep(alpha_23, n_ages^2), nrow = n_ages)
  m_alpha_24 <- matrix(data = rep(alpha_24, n_ages^2), nrow = n_ages)
  m_alpha_31 <- matrix(data = rep(alpha_31, n_ages^2), nrow = n_ages)
  m_alpha_32 <- matrix(data = rep(alpha_32, n_ages^2), nrow = n_ages)
  m_alpha_34 <- matrix(data = rep(alpha_34, n_ages^2), nrow = n_ages)
  m_alpha_41 <- matrix(data = rep(alpha_41, n_ages^2), nrow = n_ages)
  m_alpha_42 <- matrix(data = rep(alpha_42, n_ages^2), nrow = n_ages)
  m_alpha_43 <- matrix(data = rep(alpha_43, n_ages^2), nrow = n_ages)
  m_zero     <- matrix(data = rep(0, n_ages^2), nrow = n_ages)

  # Create matrix with alpha values
  m_alpha_between <- rbind(cbind(m_zero, m_alpha_12, m_alpha_13, m_alpha_14),
                           cbind(m_alpha_21, m_zero, m_alpha_23, m_alpha_24),
                           cbind(m_alpha_31, m_alpha_32, m_zero, m_alpha_34),
                           cbind(m_alpha_41, m_alpha_42, m_alpha_43, m_zero))
  return(m_alpha_between)
}

# Function to get prevalence by age
get_prev <- function(df, race) {
  v_totals <- paste0(c("i0w_", "i0r_", "i1r_"), race)
  prev_strain <- dtplyr::lazy_dt(df) %>%
    dplyr::select(., lubridate::intersect(starts_with("i"),
                                          dplyr::ends_with(race))) %>%
    dplyr::select(-one_of(v_totals)) %>%
    dplyr::as_tibble()

  l_prev <- purrr::map(purrr::set_names(c("I0w", "I0r", "I1r")),
                       ~dplyr::select(prev_strain, dplyr::starts_with(.x)))
  l_prev <- lapply(l_prev, as.matrix)

  m_prev <- l_prev[[1]] + l_prev[[2]] + l_prev[[3]]
  return(m_prev)
}

get_foi <- function(m_prev, m_waifw, scen_n, p_trans) {
  m_foi <- matrix(nrow = nrow(m_prev), ncol = ncol(m_prev))
  for (j in 1:nrow(m_prev)) {
    m_foi[j,] <- (p_trans * m_waifw) %*% m_prev[j,]
  }
  df_foi <- dtplyr::lazy_dt(m_foi) %>%
    dplyr::mutate(time = row_number() - 1) %>%
    dplyr::mutate(scen = scen_n) %>%
    dplyr::mutate(id = i) %>%
    dplyr::as_tibble()

  return(df_foi)
}

run_scenarios <- function(v_params) {
  # Make copy of v_parameter and update for burn-in
  v_parameter_burn         <- v_params
  v_parameter_burn$sigma   <- 0
  v_parameter_burn$v_alpha <- rep(0, length(v_parameter_burn$v_alpha))
  v_parameter_burn$v_psi   <- 0
  v_parameter_burn$v_eta   <- rep(0, length(v_parameter_burn$v_eta))
  v_parameter_burn$burn    <- TRUE

  # Run burn-in period and store starting states
  burn_results <- runsteady(y = v_parameter_burn$v_state,
                            times = c(0, 1E5), func = sis_abr_model,
                            parms = v_parameter_burn,
                            verbose = FALSE)
  v_burn_state <- burn_results$y

  # Simulate forward to 2024 from 1991
  v_parameter_2024            <- v_params
  v_parameter_2024$v_state    <- v_burn_state
  v_parameter_2024$v_eta      <- rep(0, length(v_parameter_2024$v_eta))
  v_parameter_2024$to_present <- TRUE
  v_parameter_2024$v_time     <- seq(0, 33, by = 1)

  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_2024,
    v_race_names = v_race_names
  )

  v_current_state <- as.numeric(tmp_results[nrow(tmp_results),
                                            2:1601])

  v_results <- list()
  # Scenario 1: No test and treat
  v_parameter_1         <- v_params
  v_parameter_1$v_time  <- seq(0, 50, by = 1)
  v_parameter_1$v_state <- v_current_state
  v_parameter_1$v_eta   <- rep(0, length(v_parameter_1$v_eta))

  scen_1_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_1,
    v_race_names = v_race_names
  )
  v_results[[1]] <- scen_1_results

  # Scenario 2: 18+, test & treat adults
  v_parameter_2          <- v_params
  v_parameter_2$v_state <- v_current_state
  v_parameter_2$v_time   <- seq(0, 50, by = 1)

  scen_2_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_2,
    v_race_names = v_race_names
  )
  v_results[[2]] <- scen_2_results

  # Scenario 3: Targeted 18+ test & treat Hispanic and NH Black
  v_parameter_3          <- v_params
  v_parameter_3$v_state  <- v_current_state
  v_parameter_3$v_eta    <- c(c(rep(0, 17), rep(1, 80 - 17)),
                              rep(0, 80),
                              c(rep(0, 17), rep(1, 80 - 17)),
                              rep(0, 80))
  v_parameter_3$v_time   <- seq(0, 50, by = 1)

  scen_3_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_3,
    v_race_names = v_race_names
  )
  v_results[[3]] <- scen_3_results

  # Scenario 4: 40+, test & treat adults
  v_parameter_4          <- v_params
  v_parameter_4$v_state  <- v_current_state
  v_parameter_4$v_eta    <- c(c(rep(0, 39), rep(1, 80 - 39)),
                              c(rep(0, 39), rep(1, 80 - 39)),
                              c(rep(0, 39), rep(1, 80 - 39)),
                              c(rep(0, 39), rep(1, 80 - 39)))
  v_parameter_4$v_time   <- seq(0, 50, by = 1)

  scen_4_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_4,
    v_race_names = v_race_names
  )
  v_results[[4]] <- scen_4_results

  # Scenario 5: Targeted 40+, test & treat adults
  v_parameter_5          <- v_params
  v_parameter_5$v_state  <- v_current_state
  v_parameter_5$v_eta    <- c(c(rep(0, 39), rep(1, 80 - 39)),
                              rep(0, 80),
                              c(rep(0, 39), rep(1, 80 - 39)),
                              rep(0, 80))
  v_parameter_5$v_time   <- seq(0, 50, by = 1)

  scen_5_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_5,
    v_race_names = v_race_names
  )
  v_results[[5]] <- scen_5_results

  v_results <- mapply(function(x, y) "[<-"(x, "scenario", value = y),
                      v_results, seq(1, 5), SIMPLIFY = FALSE)
  df_all <- dplyr::bind_rows(v_results)
  m_hisp <- df_all %>%
    dplyr::select(c(time, s0_hisp, i0w_hisp, i0r_hisp,
                    s1_hisp, i1r_hisp, scenario)) %>%
    dplyr::mutate(
      race = "Hispanic"
    ) %>%
    dplyr::rename(
      s0  = s0_hisp,
      i0w = i0w_hisp,
      i0r = i0r_hisp,
      s1  = s1_hisp,
      i1r = i1r_hisp
    )
  m_black <- df_all %>%
    dplyr::select(c(time, s0_black, i0w_black, i0r_black,
                    s1_black, i1r_black, scenario)) %>%
    dplyr::mutate(
      race = "NH Black"
    ) %>%
    dplyr::rename(
      s0  = s0_black,
      i0w = i0w_black,
      i0r = i0r_black,
      s1  = s1_black,
      i1r = i1r_black
    )
  m_white <- df_all %>%
    dplyr::select(c(time, s0_white, i0w_white, i0r_white,
                    s1_white, i1r_white, scenario)) %>%
    dplyr::mutate(
      race = "NH White"
    ) %>%
    dplyr::rename(
      s0  = s0_white,
      i0w = i0w_white,
      i0r = i0r_white,
      s1  = s1_white,
      i1r = i1r_white
    )
  m_other <- df_all %>%
    dplyr::select(c(time, s0_other, i0w_other, i0r_other,
                    s1_other, i1r_other, scenario)) %>%
    dplyr::mutate(
      race = "Other"
    ) %>%
    dplyr::rename(
      s0  = s0_other,
      i0w = i0w_other,
      i0r = i0r_other,
      s1  = s1_other,
      i1r = i1r_other
    )
  df_all_long <- rbind(m_hisp, m_black, m_white, m_other) %>%
    dplyr::mutate(I = i0w + i0r + i1r) %>%
    dplyr::mutate(Ir = i0r + i1r) %>%
    dplyr::mutate(Iw = i0w) %>%
    dplyr::mutate(scen = ifelse(scenario > 1, 1, 0)) %>%
    dplyr::mutate(scen = factor(scen, levels = c(0, 1))) %>%
    dplyr::mutate(scenario = factor(scenario, levels = c(1,2,3,4,5),
                                    labels = c("No Treatment",
                                               "18+: Test and Treat",
                                               "Targeted 18+: Test and Treat",
                                               "40+: Test and Treat",
                                               "Targeted 40+: Test and Treat")))
  df_all_long$id <- i

  # FOI by age, race, time, scenario
  # Scenario 1
  m_prev_hisp_1  <- get_prev(scen_1_results, "hisp")
  m_prev_white_1 <- get_prev(scen_1_results, "white")
  m_prev_black_1 <- get_prev(scen_1_results, "black")
  m_prev_other_1 <- get_prev(scen_1_results, "other")
  m_prev_all_1   <- as.matrix(cbind(m_prev_hisp_1, m_prev_white_1,
                                    m_prev_black_1, m_prev_other_1))
  df_foi_1 <- get_foi(m_prev = m_prev_all_1, m_waifw = v_params$m_waifw,
                      scen_n = 1, p_trans)
  # Scenario 2
  m_prev_hisp_2  <- get_prev(scen_2_results, "hisp")
  m_prev_white_2 <- get_prev(scen_2_results, "white")
  m_prev_black_2 <- get_prev(scen_2_results, "black")
  m_prev_other_2 <- get_prev(scen_2_results, "other")
  m_prev_all_2 <- cbind(m_prev_hisp_2, m_prev_white_2,
                        m_prev_black_2, m_prev_other_2)
  df_foi_2 <- get_foi(m_prev = m_prev_all_2, m_waifw = v_params$m_waifw,
                      scen_n = 2, p_trans)
  # Scenario 3
  m_prev_hisp_3  <- get_prev(scen_3_results, "hisp")
  m_prev_white_3 <- get_prev(scen_3_results, "white")
  m_prev_black_3 <- get_prev(scen_3_results, "black")
  m_prev_other_3 <- get_prev(scen_3_results, "other")
  m_prev_all_3 <- cbind(m_prev_hisp_3, m_prev_white_3,
                        m_prev_black_3, m_prev_other_3)
  df_foi_3 <- get_foi(m_prev = m_prev_all_3, m_waifw = v_params$m_waifw,
                      scen_n = 3, p_trans)

  # Scenario 4
  m_prev_hisp_4  <- get_prev(scen_4_results, "hisp")
  m_prev_white_4 <- get_prev(scen_4_results, "white")
  m_prev_black_4 <- get_prev(scen_4_results, "black")
  m_prev_other_4 <- get_prev(scen_4_results, "other")
  m_prev_all_4 <- cbind(m_prev_hisp_4, m_prev_white_4,
                        m_prev_black_4, m_prev_other_4)
  df_foi_4 <- get_foi(m_prev = m_prev_all_4, m_waifw = v_params$m_waifw,
                      scen_n = 4, p_trans)

  # Scenario 5
  m_prev_hisp_5  <- get_prev(scen_5_results, "hisp")
  m_prev_white_5 <- get_prev(scen_5_results, "white")
  m_prev_black_5 <- get_prev(scen_5_results, "black")
  m_prev_other_5 <- get_prev(scen_5_results, "other")
  m_prev_all_5 <- cbind(m_prev_hisp_5, m_prev_white_5,
                        m_prev_black_5, m_prev_other_5)
  df_foi_5 <- get_foi(m_prev = m_prev_all_5, m_waifw = v_params$m_waifw,
                      scen_n = 5, p_trans)
  # Combine
  df_foi <- rbind(df_foi_1, df_foi_2, df_foi_3, df_foi_4, df_foi_5)

  return(list(v_burn_state = v_burn_state,
              v_current_state = v_current_state,
              df_foi = df_foi,
              df_all_long = df_all_long,
              v_parameter = v_parameter))
}

##### Get Parameter Sets #####
# Assortative WAIFW parameters
v_m_resamp_assort <- c("imis_resample_1_assort_w2", "imis_resample_2_assort_w2",
                       "imis_resample_3_assort_w2", "imis_resample_4_assort_w2",
                       "imis_resample_5_assort_w2", "imis_resample_6_assort_w2",
                       "imis_resample_7_assort_w2", "imis_resample_8_assort_w2",
                       "imis_resample_9_assort_w2", "imis_resample_10_assort_w2")

for (i in seq_along(v_m_resamp_assort)) {
  assign(v_m_resamp_assort[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_assort[i], ".csv"))))
}

l_m_resamp_assort <- list(imis_resample_1_assort_w2, imis_resample_2_assort_w2,
                          imis_resample_3_assort_w2, imis_resample_4_assort_w2,
                          imis_resample_5_assort_w2, imis_resample_6_assort_w2,
                          imis_resample_7_assort_w2, imis_resample_8_assort_w2,
                          imis_resample_9_assort_w2, imis_resample_10_assort_w2)

# RP WAIFW parameters
v_m_resamp_rp <- c("imis_resample_1_rp_w2", "imis_resample_2_rp_w2",
                   "imis_resample_3_rp_w2", "imis_resample_4_rp_w2",
                   "imis_resample_5_rp_w2", "imis_resample_6_rp_w2",
                   "imis_resample_7_rp_w2", "imis_resample_8_rp_w2",
                   "imis_resample_9_rp_w2", "imis_resample_10_rp_w2")

for (i in seq_along(v_m_resamp_rp)) {
  assign(v_m_resamp_rp[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_rp[i], ".csv"))))
}

l_m_resamp_rp <- list(imis_resample_1_rp_w2, imis_resample_2_rp_w2,
                      imis_resample_3_rp_w2, imis_resample_4_rp_w2,
                      imis_resample_5_rp_w2, imis_resample_6_rp_w2,
                      imis_resample_7_rp_w2, imis_resample_8_rp_w2,
                      imis_resample_9_rp_w2, imis_resample_10_rp_w2)

# Alpha parameters
v_m_resamp_alphas <- c("imis_resample_1_alphas_w2", "imis_resample_2_alphas_w2",
                       "imis_resample_3_alphas_w2", "imis_resample_4_alphas_w2",
                       "imis_resample_5_alphas_w2", "imis_resample_6_alphas_w2",
                       "imis_resample_7_alphas_w2", "imis_resample_8_alphas_w2",
                       "imis_resample_9_alphas_w2", "imis_resample_10_alphas_w2")

for (i in seq_along(v_m_resamp_alphas)) {
  assign(v_m_resamp_alphas[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_alphas[i], ".csv"))))
}

l_m_resamp_alphas <- list(imis_resample_1_alphas_w2, imis_resample_2_alphas_w2,
                          imis_resample_3_alphas_w2, imis_resample_4_alphas_w2,
                          imis_resample_5_alphas_w2, imis_resample_6_alphas_w2,
                          imis_resample_7_alphas_w2, imis_resample_8_alphas_w2,
                          imis_resample_9_alphas_w2, imis_resample_10_alphas_w2)

assort_resamp <- imis_resample_1_assort_w2
rp_resamp <- imis_resample_1_rp_w2
alphas_resamp <- imis_resample_1_alphas_w2

w1 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)
w2 <- matrix(c(1, 1, 1, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")
v_break <- waifw_breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
groups  <- c(1:80)
n_ages  <- length(groups)
v_p <- seq(from = 0.005, to = 0.05, by = 0.005)

v_race       <- c("Hispanic", "NH White", "NH Black", "Other")
v_race_names <- c("hisp", "white", "black", "other")

v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race,
  waifw = NULL,
  trt_year = 50000000,
  end_t = 1000,
  ages = groups,
  prob = FALSE,
  p_trans = v_p[1]
)

v_race_prop <- v_parameter$v_pop / sum(v_parameter$v_pop)
v_race_prop <- rep(v_race_prop, each = length(groups))
v_race_prop <- rep(v_race_prop, 5)

v_parameter$v_state <- v_parameter$v_state * v_race_prop
v_parameter$b <- v_parameter$b * (v_parameter$v_pop / sum(v_parameter$v_pop))

##### Parallel Code #####
comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}

n <- 1000
n_cores <- 80

my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

doSNOW::registerDoSNOW(my_cluster)
doRNG::registerDoRNG(seed = 12345)

pb <- txtProgressBar(max = n, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
                   .options.snow = opts,
                   .packages = c("deSolve", "dplyr",
                                 "rootSolve", "dtplyr",
                                 "purrr", "lubridate")) %dopar% {
                                   # Outputs:
                                   # 1) FOI by age and time
                                   # 2)
                                   v_race       <- c("Hispanic", "NH White", "NH Black", "Other")
                                   v_race_names <- c("hisp", "white", "black", "other")
                                   p_trans <- v_parameter$p_trans


                                   # W1
                                   # Generate the 3 mixing pattern matrices (A, RP, PM)
                                   waifw_assort <- foi_transition_matrix_assort(assort_resamp[i, ], w = w1)
                                   waifw_rp     <- foi_transition_matrix_rp(rp_resamp[i, ], w = w1)
                                   v_alphas_within     <- as.numeric(alphas_resamp[i, c(1, 6, 11, 16)])
                                   v_alphas_within_off <- as.numeric(1 - alphas_resamp[i, c(1, 6, 11, 16)])
                                   v_alphas_between    <- as.numeric(alphas_resamp[i, c(2:5, 7:10, 12:15)])
                                   m_within      <- alpha_matrix_within(v_alphas_within)
                                   m_within_off  <- alpha_matrix_within(v_alphas_within_off)
                                   m_between     <- alpha_matrix_between(v_alphas_between)
                                   waifw_w1  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
                                     (m_between * waifw_rp)

                                   v_parameter_w1 <- v_parameter
                                   v_parameter_w1$m_waifw <- waifw_w1

                                   # W2
                                   # Generate the 3 mixing pattern matrices (A, RP, PM)
                                   waifw_assort <- foi_transition_matrix_assort(assort_resamp[i, ], w = w2)
                                   waifw_rp     <- foi_transition_matrix_rp(rp_resamp[i, ], w = w2)
                                   v_alphas_within     <- as.numeric(alphas_resamp[i, c(1, 6, 11, 16)])
                                   v_alphas_within_off <- as.numeric(1 - alphas_resamp[i, c(1, 6, 11, 16)])
                                   v_alphas_between    <- as.numeric(alphas_resamp[i, c(2:5, 7:10, 12:15)])
                                   m_within      <- alpha_matrix_within(v_alphas_within)
                                   m_within_off  <- alpha_matrix_within(v_alphas_within_off)
                                   m_between     <- alpha_matrix_between(v_alphas_between)
                                   waifw_w2  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
                                     (m_between * waifw_rp)

                                   v_parameter_w2 <- v_parameter
                                   v_parameter_w2$m_waifw <- waifw_w2

                                   results_w1     <- run_scenarios(v_parameter_w1)
                                   results_w2     <- run_scenarios(v_parameter_w2)


                                   list(df_foi_w1 = results_w1$df_foi,
                                        df_all_long_w1 = results_w1$df_all_long,
                                        df_foi_w2 = results_w2$df_foi,
                                        df_all_long_w2 = results_w2$df_all_long)
                                 }
close(pb)
parallel::stopCluster(cl = my_cluster)

## Save results
library(data.table)

fwrite(results$df_foi_w1,
       file = "results/paper/df_policy_foi_w1.csv", row.names = FALSE)
fwrite(results$df_all_long_w1,
       file = "results/paper/df_policy_long_w1.csv", row.names = FALSE)
fwrite(results$df_foi_w2,
       file = "results/paper/df_policy_foi_w2.csv", row.names = FALSE)
fwrite(results$df_all_long_w2,
       file = "results/paper/df_policy_long_w2.csv", row.names = FALSE)

## Load Results
library(data.table)
df_foi_w1          <- fread(file = "results/paper/df_policy_foi_w1.csv")
df_all_long_w1     <- fread(file = "results/paper/df_policy_long_w1.csv")
df_foi_w2          <- fread(file = "results/paper/df_policy_foi_w2.csv")
df_all_long_w2     <- fread(file = "results/paper/df_policy_long_w2.csv")

# FOI Figures
plot_hexamap_foi <- function(df_foi_long, w) {
  # Output: Hexamap for each race/ethnicity
  colnames(df_foi_long) <- c(paste0(c(rep("hisp", 80),
                                      rep("white", 80),
                                      rep("black", 80),
                                      rep("other", 80)), "_", rep(1:80, 4)),
                             "time", "scen", "id")
  df_foi_long <- df_foi_long %>%
    gather(age_race, foi, hisp_1:other_80)

  df_foi_long[c("race", "age")] <- str_split_fixed(df_foi_long$age_race, "_", 2)
  df_foi_long$period <- df_foi_long$time + 2024
  df_foi_long$age <- as.numeric(df_foi_long$age)
  df_foi_long$foi <- df_foi_long$foi * 0.05  # delete after policy is rerun

  l_df_foi <- split(df_foi_long, df_foi_long$race)
  v_scen <- unique(df_foi_long$scen)
  # y_max <- max(df_foi_long$foi)


  # for (i in 1:length(v_scen)) {
  #   data_tmp <- l_df_foi$hisp %>% filter(scen == v_scen[i])
  #   par(mar = c(2, 2, 1, 1))
  #   # pdf(paste0("results/paper/plot_", race, "_foi_", i, "_", w, ".pdf"))
  #   plot_APChexamap(dat = data_tmp, y_var = "foi")
  # }
  for (race in c("hisp", "white", "black", "other")) {
    if (race == "hisp") {
      data <- l_df_foi$hisp
      l_df_scen <- split(data, data$scen)
    } else if (race == "white") {
      data <- l_df_foi$white
      l_df_scen <- split(data, data$scen)
    } else if (race == "black") {
      data <- l_df_foi$black
      l_df_scen <- split(data, data$scen)
    } else if (race == "other") {
      data <- l_df_foi$other
      l_df_scen <- split(data, data$scen)
    }
    pdf(paste0("results/paper/hexamap_", race, "_foi_", w, ".pdf"))
    # m <- matrix(c(1, 2, 3,
    #               4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)
    # layout(mat = m, heights = c(0.5, 0.5))
    par(mfrow = c(1,5))
    for (i in 1:length(v_scen)) {
      plot_APChexamap(dat = l_df_scen[[i]], y_var = "foi")
      # title(main = paste0("Policy ", i))
    }
    dev.off()
  }
}

colnames(df_foi_w1) <- c(paste0(c(rep("hisp", 80),
                                    rep("white", 80),
                                    rep("black", 80),
                                    rep("other", 80)), "_", rep(1:80, 4)),
                           "time", "scen", "id")
df_foi_long <- df_foi_w1 %>%
  gather(age_race, foi, hisp_1:other_80)

df_foi_long[c("race", "age")] <- str_split_fixed(df_foi_long$age_race, "_", 2)
df_foi_long$period <- df_foi_long$time + 2024
df_foi_long$age <- as.numeric(df_foi_long$age)
df_foi_long$foi <- df_foi_long$foi * 0.05  # delete after policy is rerun

l_df_foi_w1 <- split(df_foi_long, df_foi_long$race)
l_df_scen <- split(l_df_foi_w1$hisp, l_df_foi_w1$hisp$scen)
# png(paste0("results/paper/hexamap_hisp_foi_w1.png"))
m <- matrix(c(1, 2, 3,
              4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(0.5, 0.5))
plot_APChexamap(dat = l_df_scen[[1]], y_var = "foi")
frame()
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
plot_APChexamap(dat = l_df_scen[[2]], y_var = "foi")
popViewport(3)
plot_APChexamap(dat = l_df_scen[[3]], y_var = "foi")
plot_APChexamap(dat = l_df_scen[[4]], y_var = "foi")
plot_APChexamap(dat = l_df_scen[[5]], y_var = "foi")
dev.off()

plot_hexamap_foi(df_foi_w1, "w1")
plot_hexamap_foi(df_foi_w2, "w2")

library(ggplot2)
library(ggsci)
library(ggpubr)


# Prevalence mean 95% Credible Interval
graph_results_grid <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  df_notrt <- df_all %>%
    mutate(I = ifelse(race == "Hispanic", I / (608.095 / 3297.255),
                      ifelse(race == "NH White", I / (1960.185 / 3297.255),
                             ifelse(race == "NH Black", I / (402.225 / 3297.255),
                                    I / (326.750 / 3297.255))))) %>%
    mutate(Ir = ifelse(race == "Hispanic", Ir / (608.095 / 3297.255),
                       ifelse(race == "NH White", Ir / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Ir / (402.225 / 3297.255),
                                     Ir / (326.750 / 3297.255))))) %>%
    mutate(Iw = ifelse(race == "Hispanic", Iw / (608.095 / 3297.255),
                       ifelse(race == "NH White", Iw / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Iw / (402.225 / 3297.255),
                                     Iw / (326.750 / 3297.255))))) %>%
    gather(strain, prev, I:Iw) %>%
    group_by(time, scen, race, strain, scenario) %>%
    summarize(
      mean_prev = mean(prev),
      lower_prev = quantile(prev, c(0.025)),
      upper_prev = quantile(prev, c(0.975))
      )
  df_notrt <- df_notrt %>%
    # filter(strain == "I") %>%
    mutate(strain = ifelse(strain == "I", "All Strains",
                           ifelse(strain == "Ir", "Treatment Resistant",
                                  "Treatment Sensitive"))) %>%
    mutate(time = time + 2024)

  plot_1 <- ggplot(data = df_notrt,
                   aes(x = time, color = as.factor(scenario))) +
    geom_ribbon(mapping = aes(ymin = lower_prev, ymax = upper_prev,
                              alpha = 0.05, fill = scenario)) +
    geom_line(mapping = aes(y = mean_prev)) +
    geom_vline(xintercept = 2025, linetype = "dashed") +
    facet_grid(race ~ strain,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_fill_nejm() +
    scale_color_nejm() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.8) ,
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(2024, 2074)) +
    xlab("Year") + ylab("Prevalence") +
    guides(fill = guide_legend(title = "Scenario:", nrow = 2, byrow = TRUE),
           color = "none",
           alpha = "none") +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          strip.text = element_text(size = 15),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))

  return(plot_1)
}

# Model trace by scenario (5)
df_all_long_w1 <- df_all_long_w1  %>%
  dplyr::mutate(scenario = ifelse(scenario == "No Treatment", 1,
                                  ifelse(scenario == "18+: Test and Treat", 2,
                                         ifelse(scenario == "Targeted 18+: Test and Treat", 3,
                                                ifelse(scenario == "40+: Test and Treat", 4, 5))))) %>%
  dplyr::mutate(scenario = factor(scenario, levels = c(1,2,3,4,5),
                                  labels = c("No Treatment",
                                             "18+: Test and Treat",
                                             "Targeted 18+: Test and Treat",
                                             "40+: Test and Treat",
                                             "Targeted 40+: Test and Treat")))

df_all_long_w2 <- df_all_long_w2  %>%
  dplyr::mutate(scenario = ifelse(scenario == "No Treatment", 1,
                                  ifelse(scenario == "18+: Test and Treat", 2,
                                         ifelse(scenario == "Targeted 18+: Test and Treat", 3,
                                                ifelse(scenario == "40+: Test and Treat", 4, 5))))) %>%
  dplyr::mutate(scenario = factor(scenario, levels = c(1,2,3,4,5),
                                  labels = c("No Treatment",
                                             "18+: Test and Treat",
                                             "Targeted 18+: Test and Treat",
                                             "40+: Test and Treat",
                                             "Targeted 40+: Test and Treat")))


plot_prev_w1 <- graph_results_grid(df_all_long_w1) +
  ggtitle("Prevalence by Strain and Race/Ethnicity") +
  theme(plot.title = element_text(size = 20))
plot_prev_w2 <- graph_results_grid(df_all_long_w2) +
  ggtitle("Prevalence by Strain and Race/Ethnicity") +
  theme(plot.title = element_text(size = 20))

ggsave("results/paper/plot_prev_w1.png", plot_prev_w1,
       width = 10, height = 10)
ggsave("results/paper/plot_prev_w2.png", plot_prev_w2,
       width = 10, height = 10)


## Box plots
pct_impact <- function(df_all) {
  require(ggplot2)
  require(dplyr)
  require(ggsci)

  df_all <- df_all %>%
    filter(time == 0 | time == 50)

  df_wide_I <- df_all %>%
    select(c(time, scenario, waifw, race, id, I)) %>%
    mutate(I = ifelse(race == "Hispanic", I / (608.095 / 3297.255),
                       ifelse(race == "NH White", I / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", I / (402.225 / 3297.255),
                                     I / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = I)

  df_wide_Ir <- df_all %>%
    select(c(time, scenario, waifw, race, id, Ir)) %>%
    mutate(Ir = ifelse(race == "Hispanic", Ir / (608.095 / 3297.255),
                       ifelse(race == "NH White", Ir / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Ir / (402.225 / 3297.255),
                                     Ir / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = Ir)

  df_wide_Iw <- df_all %>%
    select(c(time, scenario, waifw, race, id, Iw)) %>%
    mutate(Iw = ifelse(race == "Hispanic", Iw / (608.095 / 3297.255),
                       ifelse(race == "NH White", Iw / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Iw / (402.225 / 3297.255),
                                     Iw / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = Iw)

  df_wide_I$strain  <- "All Strains"
  df_wide_Ir$strain <- "Treatment Resistant"
  df_wide_Iw$strain <- "Treatment Sensitive"

  df_wide <- rbind(df_wide_I, df_wide_Ir, df_wide_Iw)

  colnames(df_wide) <- c("waifw", "race", "id", "t_0_sq",
                         "t_50_sq", "t_0_p1", "t_50_p1",
                         "t_0_p2", "t_50_p2",
                         "t_0_p3", "t_50_p3",
                         "t_0_p4", "t_50_p4", "strain")

  df_wide <- df_wide %>%
    dplyr::mutate(net_diff_p1 = (t_50_p1 - t_50_sq)) %>%
    dplyr::mutate(net_diff_p2 = (t_50_p2 - t_50_sq)) %>%
    dplyr::mutate(net_diff_p3 = (t_50_p3 - t_50_sq)) %>%
    dplyr::mutate(net_diff_p4 = (t_50_p4 - t_50_sq))

  df_long <- df_wide %>%
    pivot_longer(cols = starts_with("net_diff"),
                 names_to = "policy",
                 names_prefix = "policy_",
                 values_to = "net_diff")
  df_long <- df_long %>%
    dplyr::mutate(policy = ifelse(policy == "net_diff_p1", 1,
                                    ifelse(policy == "net_diff_p2", 2,
                                           ifelse(policy == "net_diff_p3", 3, 4)))) %>%
    dplyr::mutate(policy = factor(policy, levels = c(1,2,3,4),
                                    labels = c("18+: Test and Treat",
                                               "Targeted 18+: Test and Treat",
                                               "40+: Test and Treat",
                                               "Targeted 40+: Test and Treat")))

  plot1 <- ggplot(data = df_long, aes(x = strain, y = net_diff, fill = waifw)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_grid(race ~ policy,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Strain") + ylab("Net Change in Prevalence (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    guides(fill = guide_legend(title = "Mixing Pattern")) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 14),
          legend.position = "right")

  return(plot1)
}

df_all_long_assort$waifw <- "Assortative Mixing"
df_all_long_rp$waifw     <- "Random Mixing"
df_all_long_pm$waifw     <- "Proportional Mixing"

scenario_results_all <- as.data.frame(rbind(df_all_long_assort,
                                            df_all_long_rp,
                                            df_all_long_pm))

pct_plots <- pct_impact(scenario_results_all)
pct_plots

pct_plots <- pct_plots +
  ggtitle("Impact of Test-and-Treat Policies: Year 50")

ggsave(plot = pct_plots, filename = "results/plot_w1_net.png",
       height = 11, width = 16)

### Regression analysis
scale_pct <- function(x) (x*100)

split_data_clean <- function(df_all) {
  df_all <- df_all %>%
    filter(time == 0 | time == 50)

  df_wide_I <- df_all %>%
    select(c(time, scenario, waifw, race, id, I)) %>%
    mutate(I = ifelse(race == "Hispanic", I / (608.095 / 3297.255),
                      ifelse(race == "NH White", I / (1960.185 / 3297.255),
                             ifelse(race == "NH Black", I / (402.225 / 3297.255),
                                    I / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = I) %>%
    mutate(assort = ifelse(waifw == "Assortative Mixing", 1, 0)) %>%
    mutate(rp = ifelse(waifw == "Random Mixing", 1, 0)) %>%
    mutate(pm = ifelse(waifw == "Proportional Mixing", 1, 0))

  df_wide_Ir <- df_all %>%
    select(c(time, scenario, waifw, race, id, Ir)) %>%
    mutate(Ir = ifelse(race == "Hispanic", Ir / (608.095 / 3297.255),
                       ifelse(race == "NH White", Ir / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Ir / (402.225 / 3297.255),
                                     Ir / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = Ir) %>%
    mutate(assort = ifelse(waifw == "Assortative Mixing", 1, 0)) %>%
    mutate(rp = ifelse(waifw == "Random Mixing", 1, 0)) %>%
    mutate(pm = ifelse(waifw == "Proportional Mixing", 1, 0))

  df_wide_Iw <- df_all %>%
    select(c(time, scenario, waifw, race, id, Iw)) %>%
    mutate(Iw = ifelse(race == "Hispanic", Iw / (608.095 / 3297.255),
                       ifelse(race == "NH White", Iw / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Iw / (402.225 / 3297.255),
                                     Iw / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = Iw) %>%
    mutate(assort = ifelse(waifw == "Assortative Mixing", 1, 0)) %>%
    mutate(rp = ifelse(waifw == "Random Mixing", 1, 0)) %>%
    mutate(pm = ifelse(waifw == "Proportional Mixing", 1, 0))

  df_wide_I$strain  <- "All Strains"
  df_wide_Ir$strain <- "Treatment Resistant"
  df_wide_Iw$strain <- "Treatment Sensitive"

  colnames(df_wide_I) <- c("waifw", "race", "id", "t_0_sq",
                           "t_50_sq", "t_0_p1", "t_50_p1",
                           "t_0_p2", "t_50_p2", "t_0_p3", "t_50_p3",
                           "t_0_p4", "t_50_p4", "assort", "rp", "pm", "strain")
  colnames(df_wide_Ir) <- c("waifw", "race", "id", "t_0_sq",
                            "t_50_sq", "t_0_p1", "t_50_p1",
                            "t_0_p2", "t_50_p2", "t_0_p3", "t_50_p3",
                            "t_0_p4", "t_50_p4", "assort", "rp", "pm", "strain")
  colnames(df_wide_Iw) <- c("waifw", "race", "id", "t_0_sq",
                            "t_50_sq", "t_0_p1", "t_50_p1",
                            "t_0_p2", "t_50_p2", "t_0_p3", "t_50_p3",
                            "t_0_p4", "t_50_p4", "assort", "rp", "pm", "strain")

  # df_wide_I <- df_wide_I %>%
  #   mutate(across(c("t_0_sq", "t_50_sq", "t_0_p1", "t_50_p1", "t_0_p2",
  #                   "t_50_p2", "t_0_p3", "t_50_p3", "t_0_p4", "t_50_p4"),
  #                 scale_pct))
  # df_wide_Ir <- df_wide_Ir %>%
  #   mutate(across(c("t_0_sq", "t_50_sq", "t_0_p1", "t_50_p1", "t_0_p2",
  #                   "t_50_p2", "t_0_p3", "t_50_p3", "t_0_p4", "t_50_p4"),
  #                 scale_pct))
  # df_wide_Iw <- df_wide_Iw %>%
  #   mutate(across(c("t_0_sq", "t_50_sq", "t_0_p1", "t_50_p1", "t_0_p2",
  #                   "t_50_p2", "t_0_p3", "t_50_p3", "t_0_p4", "t_50_p4"),
  #                 scale_pct))


  df_wide_I <- df_wide_I %>%
    mutate(net_diff_p1 = (t_50_p1 - t_50_sq)) %>%
    mutate(net_diff_p2 = (t_50_p2 - t_50_sq)) %>%
    mutate(net_diff_p3 = (t_50_p3 - t_50_sq)) %>%
    mutate(net_diff_p4 = (t_50_p4 - t_50_sq))

  df_wide_Ir <- df_wide_Ir %>%
    mutate(net_diff_p1 = (t_50_p1 - t_50_sq)) %>%
    mutate(net_diff_p2 = (t_50_p2 - t_50_sq)) %>%
    mutate(net_diff_p3 = (t_50_p3 - t_50_sq)) %>%
    mutate(net_diff_p4 = (t_50_p4 - t_50_sq))

  df_wide_Iw <- df_wide_Iw %>%
    mutate(net_diff_p1 = (t_50_p1 - t_50_sq)) %>%
    mutate(net_diff_p2 = (t_50_p2 - t_50_sq)) %>%
    mutate(net_diff_p3 = (t_50_p3 - t_50_sq)) %>%
    mutate(net_diff_p4 = (t_50_p4 - t_50_sq))

  l_wide_I  <- split(df_wide_I, as.factor(df_wide_I$race))
  l_wide_Ir <- split(df_wide_Ir, as.factor(df_wide_Ir$race))
  l_wide_Iw <- split(df_wide_Iw, as.factor(df_wide_Iw$race))

  l_df_all <- list(l_wide_I[[1]], l_wide_I[[2]], l_wide_I[[3]],
                   l_wide_I[[4]], l_wide_Ir[[1]], l_wide_Ir[[2]],
                   l_wide_Ir[[3]], l_wide_Ir[[4]], l_wide_Iw[[1]],
                   l_wide_Iw[[2]], l_wide_Iw[[3]], l_wide_Iw[[4]])

  return(l_df_all)
}

# reg_and_tex <- function(l_df) {
#   library(texreg)
#   l_df_all <- list
#   l_tex <- list()
#   for (i in 1:12) {
#     tmp_reg <- lm(lm(net_diff_p1 ~ t_0_sq + rp + pm, data = l_df[[i]]))
#     tmp_tex <- createTexreg(coef.names = c("Intercept", "Prev. SQ - T0",
#                                            "Random Partnership",
#                                            "Proportional Mixing"),
#                             coef = tmp_reg$coefficients,
#                             se = summary(tmp_reg)$coefficients[, 2],
#                             pvalues = summary(tmp_reg)$coefficients[, 4])
#     l_tex[[i]] <- tmp_tex
#   }
#   return(l_tex)
# }
#
# l_df_all      <- split_data_clean(scenario_results_all)
#
# l_tex_all  <- reg_and_tex(l_df_all)
#
# all_latex <- texreg(l_tex_all, digits = 4, stars = c(0.001, 0.01, 0.05),
#                     booktabs = TRUE,
#                     custom.model.names = c("I_hisp", "I_black", "I_white",
#                                            "I_other", "Ir_hisp", "Ir_black",
#                                            "Ir_white", "Ir_other", "Iw_hisp",
#                                            "Iw_black", "Iw_white", "Iw_other"))
# write.table(all_latex, file = "results/reg_latex_all.txt")

library(modelsummary)
library(ggplot2)
library(ggsci)

l_df_all      <- split_data_clean(scenario_results_all)

results_18 <- lapply(1:12, function(x) lm(net_diff_p1 ~ t_0_sq + rp + assort, data = l_df_all[[x]])) |>
  modelplot(draw = FALSE) |>
  transform("Policy" = "18+ Test and Treat")
results_18_target <- lapply(1:12, function(x) lm(net_diff_p2 ~ t_0_sq + rp + assort, data = l_df_all[[x]])) |>
  modelplot(draw = FALSE) |>
  transform("Policy" = "Targeted 18+ Test and Treat")
results_40 <- lapply(1:12, function(x) lm(net_diff_p3 ~ t_0_sq + rp + assort, data = l_df_all[[x]])) |>
  modelplot(draw = FALSE) |>
  transform("Policy" = "40+ Test and Treat")
results_40_target <- lapply(1:12, function(x) lm(net_diff_p4 ~ t_0_sq + rp + assort, data = l_df_all[[x]])) |>
  modelplot(draw = FALSE) |>
  transform("Policy" = "Targeted 40+ Test and Treat")


results <- rbind(results_18, results_18_target, results_40, results_40_target) %>%
  filter(term == "rp" | term == "assort") %>%
  mutate(term = ifelse(term == "rp", "Random Partnership", "Assortative Mixing" )) %>%
  filter(term == "Assortative Mixing") %>%
  dplyr::mutate(Policy = ifelse(Policy == "18+ Test and Treat", 1,
                                ifelse(Policy == "Targeted 18+ Test and Treat", 2,
                                       ifelse(Policy == "40+ Test and Treat", 3, 4)))) %>%
  dplyr::mutate(Policy = factor(Policy, levels = c(1,2,3,4),
                                labels = c("18+: Test and Treat",
                                           "Targeted 18+ Test and Treat",
                                           "40+ Test and Treat",
                                           "Targeted 40+ Test and Treat")))

results$model <- factor(rep(c("Hispanic - All Strains", "NH Black - All Strains",
                              "NH White - All Strains", "NH Other - All Strains",
                              "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                              "NH White - Resistant Strains", "NH Other - Resistant Strains",
                              "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                              "NH White - Sensitive Strains", "NH Other - Sensitive Strains"), 2),
                        levels = c("Hispanic - All Strains", "NH Black - All Strains",
                                   "NH White - All Strains", "NH Other - All Strains",
                                   "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                                   "NH White - Resistant Strains", "NH Other - Resistant Strains",
                                   "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                                   "NH White - Sensitive Strains", "NH Other - Sensitive Strains"))

coef_plot <- ggplot(results, aes(y = estimate, x = Policy, xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(aes(color = Policy), position = position_dodge(width = .5)) +
  facet_wrap(~model,
             labeller = labeller(model = label_wrap_gen(20))) +
  theme_bw() +
  scale_color_nejm() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.02, 0.02)) +
  # scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
  ylab("Percentage Point Difference in Prevalence") +
  xlab(element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "bottom")
coef_plot
ggsave("results/w1_coef_plot_assort.png", coef_plot, height = 7, width = 9)

results <- rbind(results_18, results_18_target, results_40, results_40_target) %>%
  filter(term == "rp" | term == "assort") %>%
  mutate(term = ifelse(term == "rp", "Random Partnership", "Assortative Mixing" )) %>%
  filter(term == "Random Partnership") %>%
  dplyr::mutate(Policy = ifelse(Policy == "18+ Test and Treat", 1,
                                ifelse(Policy == "Targeted 18+ Test and Treat", 2,
                                       ifelse(Policy == "40+ Test and Treat", 3, 4)))) %>%
  dplyr::mutate(Policy = factor(Policy, levels = c(1,2,3,4),
                                labels = c("18+: Test and Treat",
                                           "Targeted 18+ Test and Treat",
                                           "40+ Test and Treat",
                                           "Targeted 40+ Test and Treat")))

results$model <- factor(rep(c("Hispanic - All Strains", "NH Black - All Strains",
                              "NH White - All Strains", "NH Other - All Strains",
                              "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                              "NH White - Resistant Strains", "NH Other - Resistant Strains",
                              "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                              "NH White - Sensitive Strains", "NH Other - Sensitive Strains"), 2),
                        levels = c("Hispanic - All Strains", "NH Black - All Strains",
                                   "NH White - All Strains", "NH Other - All Strains",
                                   "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                                   "NH White - Resistant Strains", "NH Other - Resistant Strains",
                                   "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                                   "NH White - Sensitive Strains", "NH Other - Sensitive Strains"))

coef_plot <- ggplot(results, aes(y = estimate, x = Policy, xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(aes(color = Policy), position = position_dodge(width = .5)) +
  facet_wrap(~model,
             labeller = labeller(model = label_wrap_gen(20))) +
  theme_bw() +
  scale_color_nejm() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     # breaks = c(-0.08, -0.04, 0, 0.04, 0.08),
                     limits = c(-0.15, 0.15)) +
  # scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
  ylab("Percentage Point Difference in Prevalence") +
  xlab(element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "bottom")
coef_plot
ggsave("results/w1_coef_plot_rp.png", coef_plot, height = 7, width = 9)

### Only look at proportional mixing results
l_df_pm <- lapply(l_df_all, function(x) filter(x, waifw == "Proportional Mixing"))
l_df_pm <- lapply(l_df_pm, function(x) pivot_longer(x,
                                                    cols = starts_with("net_diff"),
                                                    names_to = "policy",
                                                    names_prefix = "policy_",
                                                    values_to = "net_diff"))
l_df_pm_hisp  <- rbind(l_df_pm[[1]], l_df_pm[[5]], l_df_pm[[9]]) %>%
  pivot_longer(cols = starts_with("net_diff"),
               names_to = "policy",
               names_prefix = "policy_",
               values_to = "net_diff")
l_df_pm_white <- rbind(l_df_pm[[2]], l_df_pm[[6]], l_df_pm[[10]]) %>%
  pivot_longer(cols = starts_with("net_diff"),
               names_to = "policy",
               names_prefix = "policy_",
               values_to = "net_diff")
l_df_pm_black <- rbind(l_df_pm[[3]], l_df_pm[[7]], l_df_pm[[11]]) %>%
  pivot_longer(cols = starts_with("net_diff"),
               names_to = "policy",
               names_prefix = "policy_",
               values_to = "net_diff")
l_df_pm_other <- rbind(l_df_pm[[4]], l_df_pm[[8]], l_df_pm[[12]]) %>%
  pivot_longer(cols = starts_with("net_diff"),
               names_to = "policy",
               names_prefix = "policy_",
               values_to = "net_diff")
l_df_pm_all_race <- list(l_df_pm_hisp, l_df_pm_white,
                         l_df_pm_black, l_df_pm_other)

results_pm_all <- lapply(1:12, function(x) lm(net_diff ~ t_0_sq + as.factor(policy), data = l_df_pm[[x]])) |>
  modelplot(draw = FALSE) %>%
  filter(term != "(Intercept)") %>%
  filter(term != "t_0_sq")
results_pm_all$model <- factor(rep(c("Hispanic - All Strains", "NH Black - All Strains",
                              "NH White - All Strains", "NH Other - All Strains",
                              "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                              "NH White - Resistant Strains", "NH Other - Resistant Strains",
                              "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                              "NH White - Sensitive Strains", "NH Other - Sensitive Strains"), 3),
                        levels = c("Hispanic - All Strains", "NH Black - All Strains",
                                   "NH White - All Strains", "NH Other - All Strains",
                                   "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                                   "NH White - Resistant Strains", "NH Other - Resistant Strains",
                                   "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                                   "NH White - Sensitive Strains", "NH Other - Sensitive Strains"))
results_pm_all$term <- factor(c(rep("Targeted 18+ Test and Treat", 12),
                                rep("40+ Test and Treat", 12),
                                rep("Targeted 40+ Test and Treat", 12)),
                              levels = c("Targeted 18+ Test and Treat",
                                         "40+ Test and Treat",
                                         "Targeted 40+ Test and Treat"))

coef_plot_test <- ggplot(results_pm_all, aes(y = estimate, x = term, xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(aes(color = term), position = position_dodge(width = .5)) +
  facet_wrap(~model,
             labeller = labeller(model = label_wrap_gen(20))) +
  theme_bw() +
  scale_color_nejm() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     # breaks = c(-0.08, -0.04, 0, 0.04, 0.08),
                     limits = c(-0.15, 0.3)) +
  # scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
  ylab("Percentage Point Difference in Prevalence") +
  xlab(element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  guides(color = guide_legend(title = "Scenario:")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "bottom")
coef_plot_test
ggsave("results/w1_coef_plot_pm.png", coef_plot_test, height = 7, width = 9)


### Ratio of resistance to susceptible over time
graph_strain_ratio <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  df_trt_pct <- df_all %>%
    mutate(I = ifelse(race == "Hispanic", I / (608.095 / 3297.255),
                      ifelse(race == "NH White", I / (1960.185 / 3297.255),
                             ifelse(race == "NH Black", I / (402.225 / 3297.255),
                                    I / (326.750 / 3297.255))))) %>%
    mutate(Ir = ifelse(race == "Hispanic", Ir / (608.095 / 3297.255),
                       ifelse(race == "NH White", Ir / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Ir / (402.225 / 3297.255),
                                     Ir / (326.750 / 3297.255))))) %>%
    mutate(Iw = ifelse(race == "Hispanic", Iw / (608.095 / 3297.255),
                       ifelse(race == "NH White", Iw / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Iw / (402.225 / 3297.255),
                                     Iw / (326.750 / 3297.255))))) %>%
    mutate(trt_sus_pct = Iw / I) %>%
    mutate(trt_res_pct = Ir / I) %>%
    mutate(race = ifelse(race == "Other", "NH Other", race)) %>%
    group_by(time, race, scenario) %>%
    summarize(
      mean_trt_res_pct = mean(trt_res_pct),
      lower_trt_res_pct = quantile(trt_res_pct, c(0.025)),
      upper_trt_res_pct = quantile(trt_res_pct, c(0.975))
    ) %>%
    mutate(time = time + 2024)

  plot_1 <- ggplot(data = df_trt_pct,
                   aes(x = time, color = as.factor(scenario))) +
    geom_ribbon(mapping = aes(ymin = lower_trt_res_pct,
                              ymax = upper_trt_res_pct,
                              alpha = 0.05, fill = scenario)) +
    geom_line(mapping = aes(y = mean_trt_res_pct)) +
    geom_vline(xintercept = 2025, linetype = "dashed") +
    facet_wrap(race ~ .,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_fill_nejm() +
    scale_color_nejm() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 1) ,
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    # scale_y_continuous() +
    xlab("Year") + ylab("Percent of Cases Treatment Resistant") +
    guides(fill = guide_legend(title = "Scenario:", nrow = 2, byrow = TRUE),
           color = "none",
           alpha = "none") +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          strip.text = element_text(size = 15),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))

  return(plot_1)
}
pct_treatment_res_plot <- graph_strain_ratio(df_all_long_w1)
ggsave("results/paper/pct_res_plot_w1.png", pct_treatment_res_plot, height = 8, width = 12)

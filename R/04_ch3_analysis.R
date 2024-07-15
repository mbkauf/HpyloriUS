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

get_foi <- function(m_prev, m_waifw, scen_n) {
  m_foi <- matrix(nrow = nrow(m_prev), ncol = ncol(m_prev))
  for (j in 1:nrow(m_prev)) {
    m_foi[j,] <- m_waifw %*% m_prev[j,]
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
  # # Scenario 2: One-time treat everyone
  # v_parameter_2                 <- v_params
  # v_parameter_2$v_state         <- v_current_state
  # v_parameter_2$v_time          <- seq(0, 50, by = 1)
  # v_parameter_2$trt_year        <- 1
  # v_parameter_2$one_time_policy <- TRUE
  # v_parameter_2$trt_sus_policy  <- TRUE
  # v_parameter_2$v_eta           <- rep(0, length(v_parameter_2$v_eta))

  # scen_2_results <- get_sis_abr_model_results_all(
  #   v_params = v_parameter_2,
  #   v_race_names = v_race_names
  # )
  # v_results[[2]] <- scen_2_results

  # Scenario 2: 18+, test & treat adults
  v_parameter_2          <- v_params
  v_parameter_2$v_state <- v_current_state
  v_parameter_2$v_time   <- seq(0, 50, by = 1)

  scen_2_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_2,
    v_race_names = v_race_names
  )
  v_results[[2]] <- scen_2_results

  v_results <- mapply(function(x, y) "[<-"(x, "scenario", value = y),
                      v_results, seq(1, 2), SIMPLIFY = FALSE)
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
    dplyr::mutate(scenario = factor(scenario, levels = c(1, 2),
                                    labels = c("No Treatment",
                                               # "18+: One-time Mass Treatment",
                                               "18+: Test and Treat")))
  df_all_long$id <- i

  # FOI by age, race, time, scenario
  # Scenario 1
  m_prev_hisp_1  <- get_prev(scen_1_results, "hisp")
  m_prev_white_1 <- get_prev(scen_1_results, "white")
  m_prev_black_1 <- get_prev(scen_1_results, "black")
  m_prev_other_1 <- get_prev(scen_1_results, "other")
  m_prev_all_1 <- cbind(m_prev_hisp_1, m_prev_white_1,
                        m_prev_black_1, m_prev_other_1)
  df_foi_1 <- get_foi(m_prev = m_prev_all_1, m_waifw = v_params$m_waifw,
                      scen_n = 1)
  # Scenario 2
  m_prev_hisp_2  <- get_prev(scen_2_results, "hisp")
  m_prev_white_2 <- get_prev(scen_2_results, "white")
  m_prev_black_2 <- get_prev(scen_2_results, "black")
  m_prev_other_2 <- get_prev(scen_2_results, "other")
  m_prev_all_2 <- cbind(m_prev_hisp_2, m_prev_white_2,
                        m_prev_black_2, m_prev_other_2)
  df_foi_2 <- get_foi(m_prev = m_prev_all_2, m_waifw = v_params$m_waifw,
                      scen_n = 2)
  # # Scenario 3
  # m_prev_hisp_3  <- get_prev(scen_3_results, "hisp")
  # m_prev_white_3 <- get_prev(scen_3_results, "white")
  # m_prev_black_3 <- get_prev(scen_3_results, "black")
  # m_prev_other_3 <- get_prev(scen_3_results, "other")
  # m_prev_all_3 <- cbind(m_prev_hisp_3, m_prev_white_3,
  #                       m_prev_black_3, m_prev_other_3)
  # df_foi_3 <- get_foi(m_prev = m_prev_all_3, m_waifw = v_params$waifw,
  #                     scen_n = 3)

  # Combine
  df_foi <- rbind(df_foi_1, df_foi_2)  #, df_foi_3
  return(list(df_foi = df_foi,
              df_all_long = df_all_long,
              v_parameter = v_parameter))
}
run_scenarios_race <- function(v_params) {
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
  # # Scenario 2: One-time treat everyone
  # v_parameter_2                 <- v_params
  # v_parameter_2$v_state         <- v_current_state
  # v_parameter_2$v_time          <- seq(0, 50, by = 1)
  # v_parameter_2$trt_year        <- 1
  # v_parameter_2$one_time_policy <- TRUE
  # v_parameter_2$trt_sus_policy  <- TRUE
  # v_parameter_2$v_eta           <- rep(0, length(v_parameter_2$v_eta))

  # scen_2_results <- get_sis_abr_model_results_all(
  #   v_params = v_parameter_2,
  #   v_race_names = v_race_names
  # )
  # v_results[[2]] <- scen_2_results

  # Scenario 2: 18+, test & treat adults
  v_parameter_2          <- v_params
  v_parameter_2$v_state  <- v_current_state
  v_parameter_2$v_eta    <- c(c(rep(0, 17), rep(1, 80 - 17)),
                              rep(0, 80),
                              c(rep(0, 17), rep(1, 80 - 17)),
                              rep(0, 80))
  v_parameter_2$v_time   <- seq(0, 50, by = 1)

  scen_2_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_2,
    v_race_names = v_race_names
  )
  v_results[[2]] <- scen_2_results

  v_results <- mapply(function(x, y) "[<-"(x, "scenario", value = y),
                      v_results, seq(1, 2), SIMPLIFY = FALSE)
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
    dplyr::mutate(scenario = factor(scenario, levels = c(1, 2),
                                    labels = c("No Treatment",
                                               # "18+: One-time Mass Treatment",
                                               "18+: Test and Treat")))
  df_all_long$id <- i

  # FOI by age, race, time, scenario
  # Scenario 1
  m_prev_hisp_1  <- get_prev(scen_1_results, "hisp")
  m_prev_white_1 <- get_prev(scen_1_results, "white")
  m_prev_black_1 <- get_prev(scen_1_results, "black")
  m_prev_other_1 <- get_prev(scen_1_results, "other")
  m_prev_all_1 <- cbind(m_prev_hisp_1, m_prev_white_1,
                        m_prev_black_1, m_prev_other_1)
  df_foi_1 <- get_foi(m_prev = m_prev_all_1, m_waifw = v_params$m_waifw,
                      scen_n = 1)
  # Scenario 2
  m_prev_hisp_2  <- get_prev(scen_2_results, "hisp")
  m_prev_white_2 <- get_prev(scen_2_results, "white")
  m_prev_black_2 <- get_prev(scen_2_results, "black")
  m_prev_other_2 <- get_prev(scen_2_results, "other")
  m_prev_all_2 <- cbind(m_prev_hisp_2, m_prev_white_2,
                        m_prev_black_2, m_prev_other_2)
  df_foi_2 <- get_foi(m_prev = m_prev_all_2, m_waifw = v_params$m_waifw,
                      scen_n = 2)
  # # Scenario 3
  # m_prev_hisp_3  <- get_prev(scen_3_results, "hisp")
  # m_prev_white_3 <- get_prev(scen_3_results, "white")
  # m_prev_black_3 <- get_prev(scen_3_results, "black")
  # m_prev_other_3 <- get_prev(scen_3_results, "other")
  # m_prev_all_3 <- cbind(m_prev_hisp_3, m_prev_white_3,
  #                       m_prev_black_3, m_prev_other_3)
  # df_foi_3 <- get_foi(m_prev = m_prev_all_3, m_waifw = v_params$waifw,
  #                     scen_n = 3)

  # Combine
  df_foi <- rbind(df_foi_1, df_foi_2)  #, df_foi_3
  return(list(df_foi = df_foi,
              df_all_long = df_all_long,
              v_parameter = v_parameter))
}

##### Get Parameter Sets #####
# Assortative WAIFW parameters
v_m_resamp_assort <- c("imis_resample_1_assort", "imis_resample_2_assort",
                       "imis_resample_3_assort", "imis_resample_4_assort",
                       "imis_resample_5_assort", "imis_resample_6_assort",
                       "imis_resample_7_assort", "imis_resample_8_assort",
                       "imis_resample_9_assort", "imis_resample_10_assort")

for (i in seq_along(v_m_resamp_assort)) {
  assign(v_m_resamp_assort[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_assort[i], ".csv"))))
}

l_m_resamp_assort <- list(imis_resample_1_assort, imis_resample_2_assort,
                          imis_resample_3_assort, imis_resample_4_assort,
                          imis_resample_5_assort, imis_resample_6_assort,
                          imis_resample_7_assort, imis_resample_8_assort,
                          imis_resample_9_assort, imis_resample_10_assort)

# RP WAIFW parameters
v_m_resamp_rp <- c("imis_resample_1_rp", "imis_resample_2_rp",
                   "imis_resample_3_rp", "imis_resample_4_rp",
                   "imis_resample_5_rp", "imis_resample_6_rp",
                   "imis_resample_7_rp", "imis_resample_8_rp",
                   "imis_resample_9_rp", "imis_resample_10_rp")

for (i in seq_along(v_m_resamp_rp)) {
  assign(v_m_resamp_rp[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_rp[i], ".csv"))))
}

l_m_resamp_rp <- list(imis_resample_1_rp, imis_resample_2_rp,
                      imis_resample_3_rp, imis_resample_4_rp,
                      imis_resample_5_rp, imis_resample_6_rp,
                      imis_resample_7_rp, imis_resample_8_rp,
                      imis_resample_9_rp, imis_resample_10_rp)

# Alpha parameters
v_m_resamp_alphas <- c("imis_resample_1_alphas", "imis_resample_2_alphas",
                       "imis_resample_3_alphas", "imis_resample_4_alphas",
                       "imis_resample_5_alphas", "imis_resample_6_alphas",
                       "imis_resample_7_alphas", "imis_resample_8_alphas",
                       "imis_resample_9_alphas", "imis_resample_10_alphas")

for (i in seq_along(v_m_resamp_alphas)) {
  assign(v_m_resamp_alphas[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_alphas[i], ".csv"))))
}

l_m_resamp_alphas <- list(imis_resample_1_alphas, imis_resample_2_alphas,
                          imis_resample_3_alphas, imis_resample_4_alphas,
                          imis_resample_5_alphas, imis_resample_6_alphas,
                          imis_resample_7_alphas, imis_resample_8_alphas,
                          imis_resample_9_alphas, imis_resample_10_alphas)

assort_resamp <- imis_resample_1_assort
rp_resamp <- imis_resample_1_rp
alphas_resamp <- imis_resample_1_alphas

w2 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
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

                                   # Generate the 3 mixing pattern matrices (A, RP, PM)
                                   waifw_assort <- foi_transition_matrix_assort(assort_resamp[i, ])
                                   waifw_rp     <- foi_transition_matrix_rp(rp_resamp[i, ])
                                   v_alphas_within     <- as.numeric(alphas_resamp[i, c(1, 6, 11, 16)])
                                   v_alphas_within_off <- as.numeric(1 - alphas_resamp[i, c(1, 6, 11, 16)])
                                   v_alphas_between    <- as.numeric(alphas_resamp[i, c(2:5, 7:10, 12:15)])
                                   m_within      <- alpha_matrix_within(v_alphas_within)
                                   m_within_off  <- alpha_matrix_within(v_alphas_within_off)
                                   m_between     <- alpha_matrix_between(v_alphas_between)
                                   waifw_pm  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
                                     (m_between * waifw_rp)

                                   v_parameter_assort <- v_parameter
                                   v_parameter_assort$m_waifw <- waifw_assort

                                   v_parameter_rp <- v_parameter
                                   v_parameter_rp$m_waifw <- waifw_rp

                                   v_parameter_pm <- v_parameter
                                   v_parameter_pm$m_waifw <- waifw_pm

                                   results_assort <- run_scenarios(v_parameter_assort)
                                   results_rp     <- run_scenarios(v_parameter_rp)
                                   results_pm     <- run_scenarios(v_parameter_pm)

                                   results_assort_race <- run_scenarios_race(v_parameter_assort)
                                   results_rp_race     <- run_scenarios_race(v_parameter_rp)
                                   results_pm_race     <- run_scenarios_race(v_parameter_pm)

                                   list(df_foi_assort = results_assort$df_foi,
                                        df_all_long_assort = results_assort$df_all_long,
                                        df_foi_rp = results_rp$df_foi,
                                        df_all_long_rp = results_rp$df_all_long,
                                        df_foi_pm = results_pm$df_foi,
                                        df_all_long_pm = results_pm$df_all_long,
                                        df_foi_assort_race = results_assort_race$df_foi,
                                        df_all_long_assort_race = results_assort_race$df_all_long,
                                        df_foi_rp_race = results_rp_race$df_foi,
                                        df_all_long_rp_race = results_rp_race$df_all_long,
                                        df_foi_pm_race = results_pm_race$df_foi,
                                        df_all_long_pm_race = results_pm_race$df_all_long)
                                 }
close(pb)
parallel::stopCluster(cl = my_cluster)



## Save results
library(data.table)
fwrite(results$df_foi_assort,
       file = "results/ch3/abr/df_method_policy_foi_a.csv", row.names = FALSE)
fwrite(results$df_all_long_assort,
       file = "results/ch3/abr/df_method_policy_long_a.csv", row.names = FALSE)
fwrite(results$df_foi_assort_race,
       file = "results/ch3/abr/df_method_policy_foi_race_a.csv", row.names = FALSE)
fwrite(results$df_all_long_assort_race,
       file = "results/ch3/abr/df_method_policy_long_race_a.csv", row.names = FALSE)

fwrite(results$df_foi_rp,
       file = "results/ch3/abr/df_method_policy_foi_r.csv", row.names = FALSE)
fwrite(results$df_all_long_rp,
       file = "results/ch3/abr/df_method_policy_long_r.csv", row.names = FALSE)
fwrite(results$df_foi_rp_race,
       file = "results/ch3/abr/df_method_policy_foi_race_r.csv", row.names = FALSE)
fwrite(results$df_all_long_rp_race,
       file = "results/ch3/abr/df_method_policy_long_race_r.csv", row.names = FALSE)

fwrite(results$df_foi_pm,
       file = "results/ch3/abr/df_method_policy_foi_p.csv", row.names = FALSE)
fwrite(results$df_all_long_pm,
       file = "results/ch3/abr/df_method_policy_long_p.csv", row.names = FALSE)
fwrite(results$df_foi_pm_race,
       file = "results/ch3/abr/df_method_policy_foi_race_p.csv", row.names = FALSE)
fwrite(results$df_all_long_pm_race,
       file = "results/ch3/abr/df_method_policy_long_race_p.csv", row.names = FALSE)

## Load Results
library(data.table)
df_foi_assort      <- fread(file = "results/ch3/abr/df_method_policy_foi_a.csv")
df_all_long_assort <- fread(file = "results/ch3/abr/df_method_policy_long_a.csv")
df_foi_rp          <- fread(file = "results/ch3/abr/df_method_policy_foi_r.csv")
df_all_long_rp     <- fread(file = "results/ch3/abr/df_method_policy_long_r.csv")
df_foi_pm          <- fread(file = "results/ch3/abr/df_method_policy_foi_p.csv")
df_all_long_pm     <- fread(file = "results/ch3/abr/df_method_policy_long_p.csv")

df_foi_assort_race      <- fread(file = "results/ch3/abr/df_method_policy_foi_race_a.csv")
df_all_long_assort_race <- fread(file = "results/ch3/abr/df_method_policy_long_race_a.csv")
df_foi_rp_race          <- fread(file = "results/ch3/abr/df_method_policy_foi_race_r.csv")
df_all_long_rp_race     <- fread(file = "results/ch3/abr/df_method_policy_long_race_r.csv")
df_foi_pm_race          <- fread(file = "results/ch3/abr/df_method_policy_foi_race_p.csv")
df_all_long_pm_race     <- fread(file = "results/ch3/abr/df_method_policy_long_race_p.csv")

# FOI Figures (Assort)
plot_hexamap_foi <- function(df_foi_long, waifw_type = c("a", "r", "p")) {
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

  df_foi_long_1 <- df_foi_long %>%
    filter(scen == 1)
  l_df_foi <- split(df_foi_long_1, df_foi_long_1$race)

  if (waifw_type == "a") {
    waifw_text <- "assort"
  } else if (waifw_type == "r") {
    waifw_text <- "rp"
  } else if (waifw_type == "p") {
    waifw_text <- "pm"
  }

  for (race in c("hisp", "white", "black", "other")) {
    if (race == "hisp") {
      data <- l_df_foi$hisp
    } else if (race == "white") {
      data <- l_df_foi$white
    } else if (race == "black") {
      data <- l_df_foi$black
    } else if (race == "other") {
      data <- l_df_foi$other
    }
    pdf(paste0("results/plot_", race, "_foi_", waifw_text, ".pdf"))
    plot_APChexamap(dat = data, y_var = "foi")
    dev.off()
  }
}

plot_hexamap_foi(results$df_all_long_assort, "a")
plot_hexamap_foi(results$df_all_long_rp, "r")
plot_hexamap_foi(results$df_all_long_pm, "p")

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
      lower_prev = quantile(prev, c(0.25)),
      upper_prev = quantile(prev, c(0.75))
      )
  df_notrt <- df_notrt %>%
    # filter(strain == "I") %>%
    mutate(strain = ifelse(strain == "I", "All Strains",
                           ifelse(strain == "Ir", "Treatment Resistant",
                                  "Treatment Sensitive")))
  plot_1 <- ggplot(data = df_notrt,
                   aes(x = time, color = as.factor(scenario))) +
    geom_ribbon(mapping = aes(ymin = lower_prev, ymax = upper_prev,
                              alpha = 0.1, fill = scenario)) +
    geom_line(mapping = aes(y = mean_prev)) +
    facet_grid(race ~ strain,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_fill_nejm() +
    scale_color_nejm() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.8) ,
                       breaks = c(0.0, 0.2, 0.4, 0.6, 0.8),
                       labels = scales::percent_format(accuracy = 1)) +
    # scale_y_continuous() +
    xlab("Year") + ylab("Prevalence") +
    guides(fill = guide_legend(title = "Scenario:"),
           color = "none",
           alpha = "none") +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))

  return(plot_1)
}

# Test-and-treat all
plot_prev_assort <- graph_results_grid(df_all_long_assort) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Assortative WAIFW")
plot_prev_rp     <- graph_results_grid(df_all_long_rp) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Random Partnership WAIFW")
plot_prev_pm     <- graph_results_grid(df_all_long_pm) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Proportional Mixing WAIFW")
ggsave("results/plot_ch3_prev_assort.png", plot_prev_assort,
       width = 10, height = 10)
ggsave("results/plot_ch3_prev_rp.png", plot_prev_rp,
       width = 10, height = 10)
ggsave("results/plot_ch3_prev_pm.png", plot_prev_pm,
       width = 10, height = 10)

# Targeted Policy
plot_prev_assort_race <- graph_results_grid(df_all_long_assort_race) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Assortative WAIFW")
plot_prev_rp_race     <- graph_results_grid(df_all_long_rp_race) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Random Partnership WAIFW")
plot_prev_pm_race     <- graph_results_grid(df_all_long_pm_race) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Proportional Mixing WAIFW")
ggsave("results/plot_ch3_prev_assort_race.png", plot_prev_assort_race,
       width = 10, height = 10)
ggsave("results/plot_ch3_prev_rp_race.png", plot_prev_rp_race,
       width = 10, height = 10)
ggsave("results/plot_ch3_prev_pm_race.png", plot_prev_pm_race,
       width = 10, height = 10)

# Combined
df_assort_race_mod <- df_all_long_assort_race %>%
  filter(scenario == "18+: Test and Treat") %>%
  mutate(scenario = ifelse(scenario == "18+: Test and Treat",
                           "Targeted Test and Treat", .))
df_rp_race_mod <- df_all_long_rp_race %>%
  filter(scenario == "18+: Test and Treat") %>%
  mutate(scenario = ifelse(scenario == "18+: Test and Treat",
                           "Targeted Test and Treat", .))
df_pm_race_mod <- df_all_long_pm_race %>%
  filter(scenario == "18+: Test and Treat") %>%
  mutate(scenario = ifelse(scenario == "18+: Test and Treat",
                           "Targeted Test and Treat", .))

df_assort_comb <- rbind(df_all_long_assort, df_assort_race_mod)
df_rp_comb     <- rbind(df_all_long_rp, df_rp_race_mod)
df_pm_comb     <- rbind(df_all_long_pm, df_pm_race_mod)

plot_prev_assort_comb <- graph_results_grid(df_assort_comb) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Assortative WAIFW")
plot_prev_rp_comb     <- graph_results_grid(df_rp_comb) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Random Partnership WAIFW")
plot_prev_pm_comb     <- graph_results_grid(df_pm_comb) +
  ggtitle("Prevalence by Strain and Race/Ethnicity: Proportional Mixing WAIFW")
ggsave("results/plot_ch3_prev_assort_comb.png", plot_prev_assort_comb,
       width = 10, height = 10)
ggsave("results/plot_ch3_prev_rp_comb.png", plot_prev_rp_comb,
       width = 10, height = 10)
ggsave("results/plot_ch3_prev_pm_comb.png", plot_prev_pm_comb,
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
  # df_wide <- df_wide_I
  colnames(df_wide) <- c("waifw", "race", "id", "t_0_sq",
                         "t_50_sq", "t_0_p1", "t_50_p1", "strain")

  df_wide <- df_wide %>%
    # filter(scenario != "No Treatment") %>%
    dplyr::mutate(pct_diff = (t_50_p1 - t_50_sq) / t_50_sq) %>%
    dplyr::mutate(net_diff = (t_50_p1 - t_50_sq)) # %>%
    # group_by(waifw, race) %>%
    # summarize(
    #   n = sum(!is.na(net_diff)),
    #   mu_net = mean(net_diff, na.rm = TRUE),
    #   se_net = sd(net_diff, na.rm = TRUE) / sqrt(n),
    #   lb_net = mu_net - (1.96 * se_net),
    #   ub_net = mu_net + (1.96 * se_net),
    #   mu_rel = mean(pct_diff, na.rm = TRUE),
    #   se_rel = sd(pct_diff, na.rm = TRUE) / sqrt(n),
    #   lb_rel = mu_rel - (1.96 * se_rel),
    #   ub_rel = mu_rel + (1.96 * se_rel)
    # )

  plot1 <- ggplot(data = df_wide, aes(x = strain, y = net_diff, fill = waifw)) +
    # geom_bar(stat = "identity") +
    # geom_errorbar(aes(ymin = lb_net,
    #                   ymax = ub_net),
    #               position = position_dodge(.9), width = 0.3) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(. ~ race,
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
          legend.position = "bottom")

  plot2 <- ggplot(data = df_wide, aes(x = strain, y = pct_diff, fill = waifw)) +
    # geom_bar(stat = "identity") +
    # geom_errorbar(aes(ymin = lb_rel,
    #                   ymax = ub_rel),
    #               position = position_dodge(.9), width = 0.3) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(. ~ race,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Strain") + ylab("Relative Change in Prevalence (%)") +
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
          legend.position = "bottom")

  return(list(plot1, plot2))
}

df_all_long_assort$waifw <- "Assortative"
df_all_long_rp$waifw     <- "Random Partnership"
df_all_long_pm$waifw     <- "Proportional Mixing"

df_all_long_assort_race$waifw <- "Assortative"
df_all_long_rp_race$waifw     <- "Random Partnership"
df_all_long_pm_race$waifw     <- "Proportional Mixing"

scenario_results_all <- as.data.frame(rbind(df_all_long_assort,
                                            df_all_long_rp,
                                            df_all_long_pm))
scenario_results_all_race <- as.data.frame(rbind(df_all_long_assort_race,
                                                 df_all_long_rp_race,
                                                 df_all_long_pm_race))


pct_plots      <- pct_impact(scenario_results_all)
pct_plots_race <- pct_impact(scenario_results_all_race)

pct_plots[[1]] <- pct_plots[[1]] +
  ggtitle("Impact of Test-and-Treat Policy: Year 50 (All Subgroups)")

pct_plots_race[[1]] <- pct_plots_race[[1]] +
  ggtitle("Impact of Test-and-Treat Policy: Year 50 (Targeted Policy)")

ggsave(plot = pct_plots[[1]], filename = "results/plot_ch3_net.png",
       height = 10, width = 10)
ggsave(plot = pct_plots_race[[1]], filename = "results/plot_ch3_net_race.png",
       height = 10, width = 10)

ggsave(plot = pct_plots[[2]], filename = "results/plot_ch3_rel.pdf",
       height = 8, width = 12)

# scatter comparing race and ethnicity groups, facet by scenario & alphas
plot_scatter <- function(df_all, df_alpha = df_q) {
  require(ggplot2)
  require(dplyr)
  require(ggsci)

  df_all <- df_all %>%
    filter(time == 50)

  df_wide <- dcast(setDT(df_all), race + scenario + run ~ scenario,
                   value.var = "I")
  colnames(df_wide) <- c("race", "scenario", "run", "no_trt",
                         "scen_1", "scen_2")

  df_wide <- df_wide %>%
    group_by(run, race) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(
      sq = max(no_trt),
      eff = ifelse(scenario == "18+: One-time Mass Treatment",
                   sq - scen_1, sq - scen_2),
      rel_eff = eff / sq
    ) %>%
    ungroup() %>%
    dplyr::filter(scenario != "No Treatment") %>%
    arrange(desc(run))

  df_plot <- data.table(df_wide, key="run")[
    data.table(df_alpha, key="run"),
    allow.cartesian=TRUE
  ]

  plot1 <- ggplot(data = df_plot, aes(x = Q, y = eff, color = race)) +
    geom_point() +
    geom_smooth(method = "lm", color = "black") +
    facet_grid(race ~ scenario,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
  xlab("Degree of Assortativity") +
  ylab("Net Reduction in Prevalence (%)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  # scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
  theme(strip.text = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.position = "none")

  return(list(plot1, df_plot))
}
colnames(m_alphas) <- c("a_11", "a_22", "a_33", "a_44", "a_12",
                     "a_13", "a_14", "a_23", "a_24", "a_34")
df_q <- as.data.frame(m_alphas) %>%
  dplyr::select(c("a_11", "a_22", "a_33", "a_44")) %>%
  dplyr::mutate(run = row_number()) %>%
  rowwise() %>%
  dplyr::mutate(Q = mean(c(a_11, a_22, a_33, a_44)))

plot_scatter_alphas <- plot_scatter(results$df_all_long)
ggsave("results/plot_ch3_alpha_scatter.png", plot_scatter_alphas[[1]])

## Scatter by strain
plot_scatter_strain <- function(df_all, df_alpha = df_q) {
  require(ggplot2)
  require(dplyr)
  require(ggsci)

  df_all <- df_all %>%
    filter(time == 50)

  df_wide_i  <- dcast(setDT(df_all), race + scenario + run ~ scenario,
                      value.var = "I")
  df_wide_iw <- dcast(setDT(df_all), race + scenario + run ~ scenario,
                      value.var = "Iw")
  df_wide_ir <- dcast(setDT(df_all), race + scenario + run ~ scenario,
                      value.var = "Ir")
  colnames(df_wide_i)  <- c("race", "scenario", "run", "no_trt_i",
                            "scen_1_i", "scen_2_i")
  colnames(df_wide_iw) <- c("race", "scenario", "run", "no_trt_iw",
                            "scen_1_iw", "scen_2_iw")
  colnames(df_wide_ir) <-  c("race", "scenario", "run", "no_trt_ir",
                             "scen_1_ir", "scen_2_ir")
  df_wide_all <- cbind(df_wide_i, df_wide_iw$no_trt_iw, df_wide_iw$scen_1_iw,
                       df_wide_iw$scen_2_iw, df_wide_ir$no_trt_ir,
                       df_wide_ir$scen_1_ir, df_wide_ir$scen_2_ir)
  colnames(df_wide_all) <-  c("race", "scenario", "run", "no_trt_i",
                              "scen_1_i", "scen_2_i", "no_trt_iw",
                              "scen_1_iw", "scen_2_iw", "no_trt_ir",
                             "scen_1_ir", "scen_2_ir")

  df_wide <- df_wide_all %>%
    group_by(run, race) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(
      sq_iw = max(no_trt_iw),
      eff_iw = ifelse(scenario == "18+: One-time Mass Treatment",
                      sq_iw - scen_1_iw, sq_iw - scen_2_iw),
      rel_eff_iw = eff_iw / sq_iw,
      sq_ir = max(no_trt_ir),
      eff_ir = ifelse(scenario == "18+: One-time Mass Treatment",
                      sq_ir - scen_1_ir, sq_ir - scen_2_ir),
      rel_eff_ir = eff_ir / sq_ir
    ) %>%
    ungroup() %>%
    dplyr::filter(scenario != "No Treatment") %>%
    arrange(desc(run))

  df_plot <- data.table(df_wide, key="run")[
    data.table(df_alpha, key="run"),
    allow.cartesian=TRUE
  ]

  df_plot <- df_plot %>%
    pivot_longer(
      cols = starts_with("eff"),
      names_to = "strain",
      names_prefix = "eff_",
      values_to = "eff"
    ) %>%
    mutate(strain = ifelse(strain == "ir", "Treatment Resistant",
                           "Treatment Sensitive"))

  plot1 <- ggplot(data = df_plot, aes(x = Q, y = eff, color = strain)) +
    geom_point() +
    facet_grid(race ~ scenario,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
    xlab("Degree of Assortativity") +
    ylab("Net Reduction in Prevalence (%)") +
    guides(color = guide_legend(title="Strain")) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "bottom")

  return(plot1)
}

plot_scatter_alphas_strain <- plot_scatter_strain(results$df_all_long)
plot_scatter_alphas_strain
ggsave("results/plot_ch3_alpha_scatter_strain.png", plot_scatter_alphas_strain)


# Prevalence difference from Status Quo (mean 95% Credible Interval)
graph_results_grid_diff <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)
  require(stringr)

  df_plot <- df_all %>%
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
    pivot_wider(names_from = c(strain, scenario), values_from = c(prev),
                id_cols = c(id, race, time)) %>%
    mutate(I_diff_all = .[[5]] - .[[4]]) %>%
    mutate(I_diff_tar = .[[6]] - .[[4]]) %>%
    mutate(Ir_diff_all = .[[8]] - .[[7]]) %>%
    mutate(Ir_diff_tar = .[[9]] - .[[7]]) %>%
    mutate(Iw_diff_all = .[[11]] - .[[10]]) %>%
    mutate(Iw_diff_tar = .[[12]] - .[[10]]) %>%
    gather(strain_scen, diff_prev, I_diff_all:Iw_diff_tar) %>%
    dplyr::mutate(strain = case_when(
      str_detect(strain_scen, "I_diff") ~ "All Strains",
      str_detect(strain_scen, "Ir_diff") ~ "Treatment Resistant",
      str_detect(strain_scen, "Iw_diff") ~ "Treatment Sensitive"
    )) %>%
    dplyr::mutate(scenario = case_when(
      str_detect(strain_scen, "all") ~ "18+ Test and Treat All",
      str_detect(strain_scen, "tar") ~ "18+ Targeted Test and Treat",
    )) %>%
    group_by(time, race, strain, scenario) %>%
    summarize(
      mean_prev = mean(diff_prev),
      lower_prev = quantile(diff_prev, c(0.25)),
      upper_prev = quantile(diff_prev, c(0.75))
    )

  plot_1 <- ggplot(data = df_plot,
                   aes(x = time, color = as.factor(scenario))) +
    geom_ribbon(mapping = aes(ymin = lower_prev, ymax = upper_prev,
                              alpha = 0.1, fill = scenario)) +
    geom_line(mapping = aes(y = mean_prev)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    facet_grid(race ~ strain,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_fill_nejm() +
    scale_color_nejm() +
    theme_bw() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    # scale_y_continuous(limits = c(0, 0.8) ,
    #                    breaks = c(0.0, 0.2, 0.4, 0.6, 0.8)) +
    xlab("Year") + ylab("Prevalence") +
    guides(fill = guide_legend(title = "Scenario:"),
           color = "none",
           alpha = "none") +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))

  return(plot_1)
}

plot_diff_a <- graph_results_grid_diff(df_assort_comb) +
  ggtitle("Prevalence by Strain and Race/Ethnicity vs. Status Quo",
          subtitle = "Assortative WAIFW ")
plot_diff_rp <- graph_results_grid_diff(df_rp_comb) +
  ggtitle("Prevalence by Strain and Race/Ethnicity vs. Status Quo",
          subtitle = "Random Partnership WAIFW ")
plot_diff_pm <- graph_results_grid_diff(df_pm_comb) +
  ggtitle("Prevalence by Strain and Race/Ethnicity vs. Status Quo",
          subtitle = "Proportional Mixing WAIFW ")
ggsave("results/plot_ch3_diff_a.png", plot_diff_a,
       width = 10, height = 10)
ggsave("results/plot_ch3_diff_rp.png", plot_diff_rp,
       width = 10, height = 10)
ggsave("results/plot_ch3_diff_pm.png", plot_diff_pm,
       width = 10, height = 10)



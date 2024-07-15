### Script to run scenarios for 3rd chapter
library(readr)
library(rootSolve)
library(tidyverse)
library(APCtools)
library(data.table)

source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

##### Functions #####
# Function to get prevalence by age
get_prev <- function(df, race) {
  v_totals <- paste0(c("i0w_", "i0r_", "i1r_"), race)
  prev_strain <- dtplyr::lazy_dt(df) %>%
    dplyr::select(., lubridate::intersect(dplyr::starts_with("i"),
                                          dplyr::ends_with(race))) %>%
    dplyr::select(-dplyr::one_of(v_totals)) %>%
    dplyr::as_tibble()

  l_prev <- purrr::map(purrr::set_names(c("I0w", "I0r", "I1r")),
                       ~dplyr::select(prev_strain, dplyr::starts_with(.x)))
  l_prev <- lapply(l_prev, as.matrix)

  m_prev <- l_prev[[1]] + l_prev[[2]] + l_prev[[3]]
  return(m_prev)
}

get_foi <- function(m_prev, m_waifw, scen_n) {
  m_foi <- matrix(nrow = nrow(m_prev), ncol = ncol(m_prev))
  for (j in seq_len(nrow(m_prev))) {
    m_foi[j, ] <- m_waifw %*% m_prev[j, ]
  }
  df_foi <- dtplyr::lazy_dt(m_foi) %>%
    dplyr::mutate(time = dplyr::row_number() - 1) %>%
    dplyr::mutate(scen = scen_n) %>%
    dplyr::as_tibble()

  return(df_foi)
}

# Assortative WAIFW parameters
waifw_assort <- as.matrix(read_csv("results/waifw/assort_Beta_hat_w2_ga.csv"))
waifw_assort <- waifw_assort[-321:-322, -1]

# RP WAIFW parameters
waifw_rp <- as.matrix(read_csv("results/waifw/rp_Beta_hat_w2_ga.csv",
                               col_names = FALSE))

### Assortative
v_race       <- c("Hispanic", "NH Black", "NH White", "Other")
v_race_names <- c("hisp", "black", "white", "other")

v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race,
  waifw = waifw_assort,
  trt_year = 50000000,
  end_t = 1000,
  ages = seq(1, 80),
  prob = FALSE
)

# Make copy of v_parameter and update for burn-in
v_parameter_burn         <- v_parameter
v_parameter_burn$sigma   <- 0
v_parameter_burn$v_alpha <- rep(0, length(v_parameter_burn$v_alpha))
v_parameter_burn$v_psi   <- 0
v_parameter_burn$v_eta   <- rep(0, length(v_parameter_burn$v_eta))
v_parameter_burn$burn    <- TRUE

# Run burn-in period and store starting states
burn_results <- runsteady(y = v_parameter_burn$v_state,
                          times = c(0, 1E5), func = sis_abr_model,
                          parms = v_parameter_burn,
                          verbose = TRUE)
v_burn_state <- burn_results$y
v_burn_state_assort <- v_burn_state

# Simulate forward to 2024 from 1991
v_parameter_2024            <- v_parameter
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
v_current_state_assort <- v_current_state
v_results <- list()
# Scenario 1: No test and treat
v_parameter_1         <- v_parameter
v_parameter_1$v_state <- v_current_state
v_parameter_1$v_time  <- seq(0, 50, by = 1)
v_parameter_1$v_eta   <- rep(0, length(v_parameter_1$v_eta))

scen_1_results <- get_sis_abr_model_results_all(
  v_params = v_parameter_1,
  v_race_names = v_race_names
)
v_results[[1]] <- scen_1_results
# Scenario 2: One-time treat everyone
v_parameter_2                 <- v_parameter
v_parameter_2$v_state         <- v_current_state
v_parameter_2$v_time          <- seq(0, 50, by = 1)
v_parameter_2$trt_year        <- 1
v_parameter_2$one_time_policy <- TRUE
v_parameter_2$trt_sus_policy  <- TRUE
v_parameter_2$v_eta           <- rep(0, length(v_parameter_2$v_eta))

scen_2_results <- get_sis_abr_model_results_all(
  v_params = v_parameter_2,
  v_race_names = v_race_names
)
v_results[[2]] <- scen_2_results

# Scenario 3: 18+, test & treat adults
v_parameter_3          <- v_parameter
v_parameter_3$v_state  <- v_current_state
v_parameter_3$v_time   <- seq(0, 50, by = 1)

scen_3_results <- get_sis_abr_model_results_all(
  v_params = v_parameter_3,
  v_race_names = v_race_names
)
v_results[[3]] <- scen_3_results

# Make simple results data frame
v_results <- mapply(function(x, y) "[<-"(x, "scenario", value = y),
                    v_results, seq(1, 3), SIMPLIFY = FALSE)
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
  dplyr::mutate(scenario = factor(scenario, levels = c(1, 2, 3),
                                  labels = c("No Treatment",
                                             "18+: One-time Mass Treatment",
                                             "18+: Test and Treat")))

# FOI by age, race, time, scenario
# Scenario 1
m_prev_hisp_1  <- get_prev(scen_1_results, "hisp")
m_prev_white_1 <- get_prev(scen_1_results, "white")
m_prev_black_1 <- get_prev(scen_1_results, "black")
m_prev_other_1 <- get_prev(scen_1_results, "other")
m_prev_all_1 <- cbind(m_prev_hisp_1, m_prev_white_1,
                      m_prev_black_1, m_prev_other_1)
df_foi_1 <- get_foi(m_prev = m_prev_all_1, m_waifw = waifw_assort, scen_n = 1)
# Scenario 2
m_prev_hisp_2  <- get_prev(scen_2_results, "hisp")
m_prev_white_2 <- get_prev(scen_2_results, "white")
m_prev_black_2 <- get_prev(scen_2_results, "black")
m_prev_other_2 <- get_prev(scen_2_results, "other")
m_prev_all_2 <- cbind(m_prev_hisp_2, m_prev_white_2,
                      m_prev_black_2, m_prev_other_2)
df_foi_2 <- get_foi(m_prev = m_prev_all_2, m_waifw = waifw_assort, scen_n = 2)
# Scenario 3
m_prev_hisp_3  <- get_prev(scen_3_results, "hisp")
m_prev_white_3 <- get_prev(scen_3_results, "white")
m_prev_black_3 <- get_prev(scen_3_results, "black")
m_prev_other_3 <- get_prev(scen_3_results, "other")
m_prev_all_3 <- cbind(m_prev_hisp_3, m_prev_white_3,
                      m_prev_black_3, m_prev_other_3)
df_foi_3 <- get_foi(m_prev = m_prev_all_3, m_waifw = waifw_assort, scen_n = 3)

# Combine
df_foi <- rbind(df_foi_1, df_foi_2, df_foi_3)
results_assort <- list(df_foi = df_foi, df_all_long = df_all_long,
                       v_parameter = v_parameter)


### Random Partnership
v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race,
  waifw = waifw_rp,
  trt_year = 50000000,
  end_t = 1000,
  ages = seq(1, 80),
  prob = FALSE
)

# Make copy of v_parameter and update for burn-in
v_parameter_burn         <- v_parameter
v_parameter_burn$sigma   <- 0
v_parameter_burn$v_alpha <- rep(0, length(v_parameter_burn$v_alpha))
v_parameter_burn$v_psi   <- 0
v_parameter_burn$v_eta   <- rep(0, length(v_parameter_burn$v_eta))
v_parameter_burn$burn    <- TRUE

# Run burn-in period and store starting states
burn_results <- runsteady(y = v_parameter_burn$v_state,
                          times = c(0, 1E5), func = sis_abr_model,
                          parms = v_parameter_burn,
                          verbose = TRUE)
v_burn_state <- burn_results$y
v_burn_state_rp <- v_burn_state
# Simulate forward to 2024 from 1991
v_parameter_2024            <- v_parameter
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
v_current_state_rp <- v_current_state
v_results <- list()
# Scenario 1: No test and treat
v_parameter_1         <- v_parameter
v_parameter_1$v_state <- v_current_state
v_parameter_1$v_time  <- seq(0, 50, by = 1)
v_parameter_1$v_eta   <- rep(0, length(v_parameter_1$v_eta))

scen_1_results <- get_sis_abr_model_results_all(
  v_params = v_parameter_1,
  v_race_names = v_race_names
)
v_results[[1]] <- scen_1_results
# Scenario 2: One-time treat everyone
v_parameter_2                 <- v_parameter
v_parameter_2$v_state         <- v_current_state
v_parameter_2$v_time          <- seq(0, 50, by = 1)
v_parameter_2$trt_year        <- 1
v_parameter_2$one_time_policy <- TRUE
v_parameter_2$trt_sus_policy  <- TRUE
v_parameter_2$v_eta           <- rep(0, length(v_parameter_2$v_eta))

scen_2_results <- get_sis_abr_model_results_all(
  v_params = v_parameter_2,
  v_race_names = v_race_names
)
v_results[[2]] <- scen_2_results

# Scenario 3: 18+, test & treat adults
v_parameter_3          <- v_parameter
v_parameter_3$v_state  <- v_current_state
v_parameter_3$v_time   <- seq(0, 50, by = 1)

scen_3_results <- get_sis_abr_model_results_all(
  v_params = v_parameter_3,
  v_race_names = v_race_names
)
v_results[[3]] <- scen_3_results

# Make simple results data frame
v_results <- mapply(function(x, y) "[<-"(x, "scenario", value = y),
                    v_results, seq(1, 3), SIMPLIFY = FALSE)
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
  dplyr::mutate(scenario = factor(scenario, levels = c(1, 2, 3),
                                  labels = c("No Treatment",
                                             "18+: One-time Mass Treatment",
                                             "18+: Test and Treat")))

# FOI by age, race, time, scenario
# Scenario 1
m_prev_hisp_1  <- get_prev(scen_1_results, "hisp")
m_prev_white_1 <- get_prev(scen_1_results, "white")
m_prev_black_1 <- get_prev(scen_1_results, "black")
m_prev_other_1 <- get_prev(scen_1_results, "other")
m_prev_all_1 <- cbind(m_prev_hisp_1, m_prev_white_1,
                      m_prev_black_1, m_prev_other_1)
df_foi_1 <- get_foi(m_prev = m_prev_all_1, m_waifw = waifw_rp, scen_n = 1)
# Scenario 2
m_prev_hisp_2  <- get_prev(scen_2_results, "hisp")
m_prev_white_2 <- get_prev(scen_2_results, "white")
m_prev_black_2 <- get_prev(scen_2_results, "black")
m_prev_other_2 <- get_prev(scen_2_results, "other")
m_prev_all_2 <- cbind(m_prev_hisp_2, m_prev_white_2,
                      m_prev_black_2, m_prev_other_2)
df_foi_2 <- get_foi(m_prev = m_prev_all_2, m_waifw = waifw_rp, scen_n = 2)
# Scenario 3
m_prev_hisp_3  <- get_prev(scen_3_results, "hisp")
m_prev_white_3 <- get_prev(scen_3_results, "white")
m_prev_black_3 <- get_prev(scen_3_results, "black")
m_prev_other_3 <- get_prev(scen_3_results, "other")
m_prev_all_3 <- cbind(m_prev_hisp_3, m_prev_white_3,
                      m_prev_black_3, m_prev_other_3)
df_foi_3 <- get_foi(m_prev = m_prev_all_3, m_waifw = waifw_rp, scen_n = 3)
# Combine
df_foi <- rbind(df_foi_1, df_foi_2, df_foi_3)
results_rp <- list(df_foi = df_foi, df_all_long = df_all_long,
                   v_parameter = v_parameter)

## Save results
write.csv(results_assort$df_foi,
          file = "results/df_method_policy_foi_a.csv", row.names = FALSE)
write.csv(results_rp$df_foi,
          file = "results/df_method_policy_foi_r.csv", row.names = FALSE)
write.csv(results_assort$df_all_long,
          file = "results/df_method_policy_long_a.csv", row.names = FALSE)
write.csv(results_rp$df_all_long,
          file = "results/df_method_policy_long_r.csv", row.names = FALSE)

# FOI Figures (Assort)
df_foi_long <- results_assort$df_foi
colnames(df_foi_long) <- c(paste0(c(rep("hisp", 80),
                                    rep("white", 80),
                                    rep("black", 80),
                                    rep("other", 80)), "_", rep(1:80, 4)),
                           "time", "scen")
df_foi_long <- df_foi_long %>%
  gather(age_race, foi, hisp_1:other_80)

df_foi_long[c("race", "age")] <- str_split_fixed(df_foi_long$age_race, "_", 2)
df_foi_long$period <- df_foi_long$time + 2024
df_foi_long$age <- as.numeric(df_foi_long$age)

df_foi_long_1 <- df_foi_long %>%
  filter(scen == 1)
df_foi_long_2 <- df_foi_long %>%
  filter(scen == 2)
df_foi_long_3 <- df_foi_long %>%
  filter(scen == 3)
l_df_foi_1 <- split(df_foi_long_1, df_foi_long_1$race)
l_df_foi_2 <- split(df_foi_long_2, df_foi_long_2$race)
l_df_foi_3 <- split(df_foi_long_3, df_foi_long_3$race)

pdf("results/plot_ch3_foi_assort_sq.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi_1$hisp, y_var = "foi",
                                  legend_title = "FOI - Hispanic")
plot_white_foi <- plot_APChexamap(dat = l_df_foi_1$white, y_var = "foi",
                                  legend_title = "FOI - NH White")
plot_black_foi <- plot_APChexamap(dat = l_df_foi_1$black, y_var = "foi",
                                  legend_title = "FOI - NH Black")
plot_other_foi <- plot_APChexamap(dat = l_df_foi_1$other, y_var = "foi",
                                  legend_title = "FOI - NH Other")
dev.off()

pdf("results/plot_ch3_foi_assort_p1.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi_2$hisp, y_var = "foi",
                                  legend_title = "FOI - Hispanic")
plot_white_foi <- plot_APChexamap(dat = l_df_foi_2$white, y_var = "foi",
                                  legend_title = "FOI - NH White")
plot_black_foi <- plot_APChexamap(dat = l_df_foi_2$black, y_var = "foi",
                                  legend_title = "FOI - NH Black")
plot_other_foi <- plot_APChexamap(dat = l_df_foi_2$other, y_var = "foi",
                                  legend_title = "FOI - NH Other")
dev.off()

pdf("results/plot_ch3_foi_assort_p2.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi_3$hisp, y_var = "foi",
                                  legend_title = "FOI - Hispanic")
plot_white_foi <- plot_APChexamap(dat = l_df_foi_3$white, y_var = "foi",
                                  legend_title = "FOI - NH White")
plot_black_foi <- plot_APChexamap(dat = l_df_foi_3$black, y_var = "foi",
                                  legend_title = "FOI - NH Black")
plot_other_foi <- plot_APChexamap(dat = l_df_foi_3$other, y_var = "foi",
                                  legend_title = "FOI - NH Other")
dev.off()

# FOI Figures (RP)
df_foi_long <- results_rp$df_foi
colnames(df_foi_long) <- c(paste0(c(rep("hisp", 80),
                                    rep("white", 80),
                                    rep("black", 80),
                                    rep("other", 80)), "_", rep(1:80, 4)),
                           "time", "scen")
df_foi_long <- df_foi_long %>%
  gather(age_race, foi, hisp_1:other_80)

df_foi_long[c("race", "age")] <- str_split_fixed(df_foi_long$age_race, "_", 2)
df_foi_long$period <- df_foi_long$time + 2024
df_foi_long$age <- as.numeric(df_foi_long$age)

df_foi_long_1 <- df_foi_long %>%
  filter(scen == 1)
df_foi_long_2 <- df_foi_long %>%
  filter(scen == 2)
df_foi_long_3 <- df_foi_long %>%
  filter(scen == 3)
l_df_foi_1 <- split(df_foi_long_1, df_foi_long_1$race)
l_df_foi_2 <- split(df_foi_long_2, df_foi_long_2$race)
l_df_foi_3 <- split(df_foi_long_3, df_foi_long_3$race)

pdf("results/plot_ch3_foi_rp_sq.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi_1$hisp, y_var = "foi",
                                  legend_title = "FOI - Hispanic")
plot_white_foi <- plot_APChexamap(dat = l_df_foi_1$white, y_var = "foi",
                                  legend_title = "FOI - NH White")
plot_black_foi <- plot_APChexamap(dat = l_df_foi_1$black, y_var = "foi",
                                  legend_title = "FOI - NH Black")
plot_other_foi <- plot_APChexamap(dat = l_df_foi_1$other, y_var = "foi",
                                  legend_title = "FOI - NH Other")
dev.off()

pdf("results/plot_ch3_foi_rp_p1.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi_2$hisp, y_var = "foi",
                                  legend_title = "FOI - Hispanic")
plot_white_foi <- plot_APChexamap(dat = l_df_foi_2$white, y_var = "foi",
                                  legend_title = "FOI - NH White")
plot_black_foi <- plot_APChexamap(dat = l_df_foi_2$black, y_var = "foi",
                                  legend_title = "FOI - NH Black")
plot_other_foi <- plot_APChexamap(dat = l_df_foi_2$other, y_var = "foi",
                                  legend_title = "FOI - NH Other")
dev.off()

pdf("results/plot_ch3_foi_rp_p2.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi_3$hisp, y_var = "foi",
                                  legend_title = "FOI - Hispanic")
plot_white_foi <- plot_APChexamap(dat = l_df_foi_3$white, y_var = "foi",
                                  legend_title = "FOI - NH White")
plot_black_foi <- plot_APChexamap(dat = l_df_foi_3$black, y_var = "foi",
                                  legend_title = "FOI - NH Black")
plot_other_foi <- plot_APChexamap(dat = l_df_foi_3$other, y_var = "foi",
                                  legend_title = "FOI - NH Other")
dev.off()

# Prevalence mean 95% Credible Interval
graph_results_grid <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  df_notrt <- df_all %>%
    tidyr::gather(strain, prev, I:Iw) %>%
    dplyr::mutate(strain = ifelse(strain == "I", "Total Prevalence",
                                  ifelse(strain == "Ir", "Resistant Strain",
                                         "Sensitive Strain"))) %>%
    dplyr::mutate(time = time + 2024)

  plot_1 <- ggplot(data = df_notrt,
                   aes(x = time, color = as.factor(scenario))) +
    geom_line(mapping = aes(y = prev)) +
    facet_grid(race ~ strain,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.6),
                       breaks = c(0, 0.2, 0.4, 0.6)) +
    scale_x_continuous(breaks = c(2024, 2044, 2064)) +
    xlab("Year") + ylab("Prevalence") +
    guides(color = guide_legend(title = "Scenario:")) +
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

# df_all_long <- read.csv("results/df_all_long.csv")
plot_assort_prev <- graph_results_grid(results_assort$df_all_long)
ggsave("results/plot_ch3_assort_prev.pdf", plot_assort_prev)

plot_rp_prev <- graph_results_grid(results_rp$df_all_long)
ggsave("results/plot_ch3_rp_prev.pdf", plot_rp_prev)

pct_impact <- function(df_all) {
  require(ggplot2)
  require(dplyr)
  require(ggsci)
  
  df_all <- df_all %>%
    filter(time == 0 | time == 50)
  
  df_wide <- dcast(setDT(df_all), scenario + waifw + race ~ time,
                   value.var = "I")
  colnames(df_wide) <- c("scenario", "waifw", "race", "t_0", "t_50")
  
  df_wide <- df_wide %>%
    filter(scenario != "No Treatment") %>%
    dplyr::mutate(pct_diff = abs(t_50 - t_0) / t_0) %>%
    dplyr::mutate(net_diff = abs(t_0 - t_50))
  
  plot1 <- ggplot(data = df_wide, aes(x = waifw, y = net_diff, fill = waifw)) +
    geom_bar(stat = "identity") +
    facet_grid(scenario ~ race,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Mixing Pattern") + ylab("Net Reduction in Prevalence (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")
  
  plot2 <- ggplot(data = df_wide, aes(x = waifw, y = pct_diff, fill = waifw)) +
    geom_bar(stat = "identity") +
    facet_grid(scenario ~ race,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Mixing Pattern") + ylab("Relative Reduction in Prevalence (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")
  
  return(list(plot1, plot2))
}

results_assort$df_all_long$waifw <- "Assortative"
results_rp$df_all_long$waifw <- "Random Partnership"
scenario_results_all <- rbind(results_assort$df_all_long, results_rp$df_all_long) 

pct_plots <- pct_impact(scenario_results_all)
ggsave(plot = pct_plots[[1]], filename = "results/plot_ch3_net.pdf",
       height = 8, width = 12)
ggsave(plot = pct_plots[[2]], filename = "results/plot_ch3_rel.pdf",
       height = 8, width = 12)

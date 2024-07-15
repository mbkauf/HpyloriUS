##### Load Libraries #####
library(doParallel)
library(doRNG)
library(foreach)
library(dtplyr)
library(dplyr)
library(readr)
library(rootSolve)
library(tidyverse)
library(APCtools)

source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

##### Functions #####
get_waifw_draws <- function(m_cov, v_betas) {
  chol_mat <- chol(m_cov)
  v_in <- qnorm(runif(nrow(m_cov), 0, 1), 0, 1)  # vector of random inverse normal
  v_tz <- rep(0, nrow(m_cov))

  for (i in 1:nrow(m_cov)) {
    v_tz[i] <- v_in %*% chol_mat[i,]
  }

  v_betas_rand <- v_betas + v_tz
  return(v_betas_rand)
}

get_new_waifw <- function(old_waifw, old_betas, new_betas) {
  new_waifw <- old_waifw
  for (i in 1:length(old_betas)) {
    new_waifw[new_waifw == old_betas[i]] <- new_betas[i]
  }
  return(new_waifw)
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


##### Get Parameter Sets #####
# assort_beta_hat_2 <- read_csv("results/waifw/assort_beta_hat_2.csv")
# v_betas_test <- as.numeric(assort_beta_hat_2[2, 2:33])

# assort_beta_hat_cov_2 <- read_csv("results/waifw/assort_beta_hat_cov_2.csv")
# m_cov_test <- assort_beta_hat_cov_2[, 34:65]
# m_cov_test <- as.matrix(m_cov_test)

# test_waifw_betas <- get_waifw_draws(m_cov = m_cov_test,
#                                     v_betas = v_betas_test)
# # test replace matrix values
# assort_Beta_hat_w2_2 <- as.matrix(read_csv("results/waifw/assort_Beta_hat_w2_2.csv"))
# assort_Beta_hat_w2_2 <- assort_Beta_hat_w2_2[, -1]

# rand_waifw <- get_new_waifw(old_waifw = assort_Beta_hat_w2_2,
#                             old_betas = v_betas_test,
#                             new_betas = test_waifw_betas)

# Assortative WAIFW parameters
assort_beta_hat_2 <- read_csv("results/waifw/assort_beta_hat_ga.csv")
v_betas_assort <- as.numeric(assort_beta_hat_2[2, 2:33])

assort_beta_hat_cov_2 <- read_csv("results/waifw/assort_beta_cov_w2_ga.csv")
m_cov_assort <- assort_beta_hat_cov_2[, -1]
m_cov_assort <- as.matrix(m_cov_assort)

waifw_assort <- as.matrix(read_csv("results/waifw/assort_Beta_hat_w2_ga.csv"))
waifw_assort <- waifw_assort[-321, -1]

# RP WAIFW parameters
rp_beta_hat_2 <- read_csv("results/waifw/rp_beta_hat_ga.csv")
v_betas_rp <- as.numeric(rp_beta_hat_2[2, 2:81])

rp_beta_hat_cov_2 <- read_csv("results/waifw/rp_beta_cov_w2_ga.csv")
m_cov_rp <- rp_beta_hat_cov_2[, -1]
m_cov_rp <- as.matrix(m_cov_rp)

waifw_rp <- as.matrix(read_csv("results/waifw/rp_Beta_hat_w2_ga.csv"))
waifw_rp <- waifw_rp[, -1]

##### Parallel Code #####
comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}

n <- 100
n_cores <- 50

my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my_cluster)
doRNG::registerDoRNG(seed = 12345)

results_assort <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
                          .packages = c("deSolve", "dplyr",
                                        "rootSolve", "dtplyr",
                                        "purrr", "lubridate")) %dopar% {
  # Outputs:
  # 1) FOI by age and time
  # 2)
  v_race       <- c("Hispanic", "NH Black", "NH White", "Other")
  v_race_names <- c("hisp", "black", "white", "other")

  v_rand_betas <- get_waifw_draws(m_cov = m_cov_assort,
                                  v_betas = v_betas_assort)
  rand_waifw <- get_new_waifw(old_waifw = waifw_assort,
                              old_betas = v_betas_assort,
                              new_betas = v_rand_betas)

  v_parameter <- load_sis_abr_model_params_all(
    v_race = v_race,
    waifw = rand_waifw,
    trt_year = 50000000,
    end_t = 1000,
    ages = seq(1, 80),
    prob = TRUE
  )

  # Make copy of v_parameter and update for burn-in
  v_parameter_burn         <- v_parameter
  v_parameter_burn$sigma   <- 0
  v_parameter_burn$v_alpha <- 0
  v_parameter_burn$v_psi   <- 0

  # Run burn-in period and store starting states
  burn_results <- runsteady(y = v_parameter_burn$v_state,
                            times = c(0,1E5), func = sis_abr_model,
                            parms = v_parameter_burn)
  v_burn_state <- burn_results$y

  # Simulate forward to 2024 from 1991
  v_parameter$v_time <- seq(0, 33, by = 1)
  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names
  )
  v_current_state <- as.numeric(tmp_results[nrow(tmp_results),
                                            2:1601])

  v_results <- list()
  # Scenario 1: No test and treat
  v_parameter_1         <- v_parameter
  v_parameter_1$v_state <- v_current_state
  v_parameter_1$v_time  <- seq(0, 50, by = 1)

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
                                               "18+: ?",
                                               "18+: ??")))
  df_all_long$id <- i

  # FOI by age, race, time, scenario
  # Scenario 1
  m_prev_hisp_1  <- get_prev(scen_1_results, "hisp")
  m_prev_white_1 <- get_prev(scen_1_results, "white")
  m_prev_black_1 <- get_prev(scen_1_results, "black")
  m_prev_other_1 <- get_prev(scen_1_results, "other")
  m_prev_all_1 <- cbind(m_prev_hisp_1, m_prev_white_1,
                        m_prev_black_1, m_prev_other_1)

  df_foi_1 <- get_foi(m_prev = m_prev_all_1, m_waifw = rand_waifw, scen_n = 1)

  # Scenario 2
  m_prev_hisp_2  <- get_prev(scen_2_results, "hisp")
  m_prev_white_2 <- get_prev(scen_2_results, "white")
  m_prev_black_2 <- get_prev(scen_2_results, "black")
  m_prev_other_2 <- get_prev(scen_2_results, "other")
  m_prev_all_2 <- cbind(m_prev_hisp_2, m_prev_white_2,
                        m_prev_black_2, m_prev_other_2)

  df_foi_2 <- get_foi(m_prev = m_prev_all_2, m_waifw = rand_waifw, scen_n = 2)

  # Scenario 3
  m_prev_hisp_3  <- get_prev(scen_3_results, "hisp")
  m_prev_white_3 <- get_prev(scen_3_results, "white")
  m_prev_black_3 <- get_prev(scen_3_results, "black")
  m_prev_other_3 <- get_prev(scen_3_results, "other")
  m_prev_all_3 <- cbind(m_prev_hisp_3, m_prev_white_3,
                        m_prev_black_3, m_prev_other_3)

  df_foi_3 <- get_foi(m_prev = m_prev_all_3, m_waifw = rand_waifw, scen_n = 3)

  # Combine
  df_foi <- rbind(df_foi_1, df_foi_2, df_foi_3)

  list(df_foi = df_foi, df_all_long = df_all_long, v_parameter = v_parameter)
}
parallel::stopCluster(cl = my_cluster)

# Random Partnership
my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my_cluster)
doRNG::registerDoRNG(seed = 12345)

results_rp <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
                      .packages = c("deSolve", "dplyr",
                                    "rootSolve", "dtplyr",
                                    "purrr", "lubridate")) %dopar% {
  # Outputs:
  # 1) FOI by age and time
  # 2)
  v_race       <- c("Hispanic", "NH Black", "NH White", "Other")
  v_race_names <- c("hisp", "black", "white", "other")

  v_rand_betas <- get_waifw_draws(m_cov = m_cov_rp,
                                  v_betas = v_betas_rp)
  rand_waifw <- get_new_waifw(old_waifw = waifw_rp,
                              old_betas = v_betas_rp,
                              new_betas = v_rand_betas)

  v_parameter <- load_sis_abr_model_params_all(
    v_race = v_race,
    waifw = rand_waifw,
    trt_year = 50000000,
    end_t = 1000,
    ages = seq(1, 80),
    prob = TRUE
  )

  # Make copy of v_parameter and update for burn-in
  v_parameter_burn         <- v_parameter
  v_parameter_burn$sigma   <- 0
  v_parameter_burn$v_alpha <- 0
  v_parameter_burn$v_psi   <- 0

  # Run burn-in period and store starting states
  burn_results <- runsteady(y = v_parameter_burn$v_state,
                            times = c(0,1E5), func = sis_abr_model,
                            parms = v_parameter_burn)
  v_burn_state <- burn_results$y

  # Simulate forward to 2024 from 1991
  v_parameter$v_time <- seq(0, 33, by = 1)
  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names
  )
  v_current_state <- as.numeric(tmp_results[nrow(tmp_results),
                                            2:1601])

  v_results <- list()
  # Scenario 1: No test and treat
  v_parameter_1         <- v_parameter
  v_parameter_1$v_state <- v_current_state
  v_parameter_1$v_time  <- seq(0, 50, by = 1)

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
                                               "18+: ?",
                                               "18+: ??")))
  df_all_long$id <- i

  # FOI by age, race, time, scenario
  # Scenario 1
  m_prev_hisp_1  <- get_prev(scen_1_results, "hisp")
  m_prev_white_1 <- get_prev(scen_1_results, "white")
  m_prev_black_1 <- get_prev(scen_1_results, "black")
  m_prev_other_1 <- get_prev(scen_1_results, "other")
  m_prev_all_1 <- cbind(m_prev_hisp_1, m_prev_white_1,
                        m_prev_black_1, m_prev_other_1)

  df_foi_1 <- get_foi(m_prev = m_prev_all_1, m_waifw = rand_waifw, scen_n = 1)

  # Scenario 2
  m_prev_hisp_2  <- get_prev(scen_2_results, "hisp")
  m_prev_white_2 <- get_prev(scen_2_results, "white")
  m_prev_black_2 <- get_prev(scen_2_results, "black")
  m_prev_other_2 <- get_prev(scen_2_results, "other")
  m_prev_all_2 <- cbind(m_prev_hisp_2, m_prev_white_2,
                        m_prev_black_2, m_prev_other_2)

  df_foi_2 <- get_foi(m_prev = m_prev_all_2, m_waifw = rand_waifw, scen_n = 2)

  # Scenario 3
  m_prev_hisp_3  <- get_prev(scen_3_results, "hisp")
  m_prev_white_3 <- get_prev(scen_3_results, "white")
  m_prev_black_3 <- get_prev(scen_3_results, "black")
  m_prev_other_3 <- get_prev(scen_3_results, "other")
  m_prev_all_3 <- cbind(m_prev_hisp_3, m_prev_white_3,
                        m_prev_black_3, m_prev_other_3)

  df_foi_3 <- get_foi(m_prev = m_prev_all_3, m_waifw = rand_waifw, scen_n = 3)

  # Combine
  df_foi <- rbind(df_foi_1, df_foi_2, df_foi_3)

  list(df_foi = df_foi, df_all_long = df_all_long, v_parameter = v_parameter)
}
parallel::stopCluster(cl = my_cluster)

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
                           "time", "scen", "id")
df_foi_long <- df_foi_long %>%
  gather(age_race, foi, hisp_1:other_80)

df_foi_long[c("race", "age")] <- str_split_fixed(df_foi_long$age_race, "_", 2)
df_foi_long$period <- df_foi_long$time + 2024
df_foi_long$age <- as.numeric(df_foi_long$age)

df_foi_long_1 <- df_foi_long %>%
  filter(scen == 1)
l_df_foi <- split(df_foi_long_1, df_foi_long_1$race)

pdf("results/plot_black_foi_assort.pdf")
plot_black_foi <- plot_APChexamap(dat = l_df_foi$black, y_var = "foi")
dev.off()

pdf("results/plot_hisp_foi_assort.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi$hisp, y_var = "foi")
dev.off()

pdf("results/plot_white_foi_assort.pdf")
plot_white_foi <- plot_APChexamap(dat = l_df_foi$white, y_var = "foi")
dev.off()

pdf("results/plot_other_foi_assort.pdf")
plot_other_foi <- plot_APChexamap(dat = l_df_foi$other, y_var = "foi")
dev.off()

# FOI Figures (RP)
df_foi_long <- results_rp$df_foi
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

pdf("results/plot_black_foi_rp.pdf")
plot_black_foi <- plot_APChexamap(dat = l_df_foi$black, y_var = "foi")
dev.off()

pdf("results/plot_hisp_foi_rp.pdf")
plot_hisp_foi  <- plot_APChexamap(dat = l_df_foi$hisp, y_var = "foi")
dev.off()

pdf("results/plot_white_foi_rp.pdf")
plot_white_foi <- plot_APChexamap(dat = l_df_foi$white, y_var = "foi")
dev.off()

pdf("results/plot_other_foi_rp.pdf")
plot_other_foi <- plot_APChexamap(dat = l_df_foi$other, y_var = "foi")
dev.off()

# Prevalence mean 95% Credible Interval
graph_results_grid <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  df_notrt <- df_all %>%
    # filter(scen == 0) %>%
    # mutate(scenario_name = "Status Quo") %>%
    gather(strain, prev, I:Iw) %>%
    group_by(time, scen, race, strain, scenario) %>%
    summarize(
      mean_prev = mean(prev),
      lower_prev = quantile(prev, c(0.025)),
      upper_prev = quantile(prev, c(0.975))
      )

  plot_1 <- ggplot(data = df_notrt, 
                   aes(x = time, color = as.factor(scenario))) +
    geom_ribbon(mapping = aes(ymin = lower_prev, ymax = upper_prev,
                              alpha = 0.2, fill = scenario)) +
    geom_line(mapping = aes(y = mean_prev)) +
    facet_grid(race ~ strain,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
    # scale_y_continuous(limits = c(0, 0.6) ,
    #                    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    xlab("Year") + ylab("Prevalence by Strain") +
    guides(fill = guide_legend(title = "Scenario:"),
           linetype = guide_legend(title = "Strain:"),
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

df_all_long <- read.csv("results/df_all_long.csv")
test_plot <- graph_results_grid(df_all_long)
test_plot

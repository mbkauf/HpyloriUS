##### This file is for calibrating alpha values #####
### Load libraries
library(readr)
library(rootSolve)
library(tidyverse)
library(APCtools)
library(doParallel)
library(foreach)
library(data.table)


### Functions
# Bounds function
bound <- function(v_var, min, max) {
  vec <- pmin(v_var, max)
  vec <- pmax(vec, min)
  return(vec)
}

sample_alphas <- function(n) {
  # Hispanic
  a_12 <- runif(n, min = 0.048, max = 0.686)
  a_13 <- runif(n, min = 0.016, max = 0.299) * (1 - a_12)
  a_13 <- bound(a_13, min = 0.016, max = 0.299)
  a_14 <- runif(n, min = 0.028, max = 0.349) * (1 - a_12 - a_13)
  a_14 <- bound(a_14, min = 0.028, max = 0.349)
  a_11 <- 1 - a_12 - a_13 - a_14
  a_11 <- bound(a_11, min = 0.044, max = 0.660)

  # NH White
  a_23 <- runif(n, min = 0.012, max = 0.388) * (1 - a_12)
  a_23 <- bound(a_23, min = 0.012, max = 0.388)
  a_24 <- runif(n, min = 0.049, max = 0.409) * (1 - a_12 - a_23)
  a_24 <- bound(a_24, min = 0.049, max = 0.409)
  a_22 <- 1 - a_12 - a_23 - a_24
  a_22 <- bound(a_22, min = 0.138, max = 0.648)

  # NH Black
  a_34 <- runif(n, min = 0.029, max = 0.157) * (1 - a_13 - a_23)
  a_34 <- bound(a_34, min = 0.029, max = 0.157)
  a_33 <- 1 - a_13 - a_23 - a_34
  a_33 <- bound(a_33, min = 0.005, max = 0.671)

  # NH Other
  a_44 <- 1 - a_14 - a_24 - a_34
  a_44 <- bound(a_44, min = 0.000, max = 0.661)

  m_alphas <- matrix(c(a_11, a_12, a_13, a_14, a_22,
                       a_23, a_24, a_33, a_34, a_44), ncol = 10)
  a_1j <- m_alphas[, 1] + m_alphas[, 2] + m_alphas[, 3] + m_alphas[, 4]
  a_2j <- m_alphas[, 2] + m_alphas[, 5] + m_alphas[, 6] + m_alphas[, 7]
  a_3j <- m_alphas[, 3] + m_alphas[, 6] + m_alphas[, 8] + m_alphas[, 9]
  a_4j <- m_alphas[, 4] + m_alphas[, 7] + m_alphas[, 9] + m_alphas[, 10]

  m_alphas <- cbind(m_alphas, a_1j, a_2j, a_3j, a_4j)
  m_alphas <- m_alphas[m_alphas[, 11] == 1, ]
  m_alphas <- m_alphas[m_alphas[, 12] == 1, ]
  m_alphas <- m_alphas[m_alphas[, 13] == 1, ]
  m_alphas <- m_alphas[m_alphas[, 14] == 1, ]
  m_alphas <- m_alphas[, 1:10]
  return(m_alphas)
}

# Create matrix of alphas to overlay WAIFWs
foi_transition_matrix_alpha <- function(v_alpha, w = w2,
                                        waifw_breaks = c(1, 5, 15, 25, 45,
                                                         55, 65, 70, 81)) {
  alpha_1  <- v_alpha[1]
  alpha_2  <- v_alpha[2]
  alpha_3  <- v_alpha[3]
  alpha_4  <- v_alpha[4]
  alpha_5  <- v_alpha[5]
  alpha_6  <- v_alpha[6]
  alpha_7  <- v_alpha[7]
  alpha_8  <- v_alpha[8]
  alpha_9  <- v_alpha[9]
  alpha_10 <- v_alpha[10]

  beta1  <- get_transmission_matrix(betas = rep(alpha_1, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta2  <- get_transmission_matrix(betas = rep(alpha_2, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta3  <- get_transmission_matrix(betas = rep(alpha_3, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta4  <- get_transmission_matrix(betas = rep(alpha_4, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta5  <- get_transmission_matrix(betas = rep(alpha_5, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta6  <- get_transmission_matrix(betas = rep(alpha_6, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta7  <- get_transmission_matrix(betas = rep(alpha_7, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta8  <- get_transmission_matrix(betas = rep(alpha_8, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta9  <- get_transmission_matrix(betas = rep(alpha_9, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta10 <- get_transmission_matrix(betas = rep(alpha_10, 8),  #nolint
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

# Function to get FOI by age
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

# Load other functions
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

### Get parameter sets
groups <- c(1:80)
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

v_race       <- c("Hispanic", "NH White", "NH Black", "Other")
v_race_names <- c("hisp", "white", "black", "other")
v_demo_vars <- get_demographic_vars_all(v_race = v_race,
                                        ages = groups)
v_prev_vars <- get_prevalence_vars_all(v_race = v_race,
                                       ages = groups)
v_inf <- v_prev_vars$prevalence * v_demo_vars$v_age_prop

# Assortative WAIFW parameters
waifw_assort <- as.matrix(read_csv("results/waifw/assort_Beta_hat_w2_ga.csv"))
waifw_assort <- waifw_assort[-321:-322, -1]

# RP WAIFW parameters
waifw_rp <- as.matrix(read_csv("results/waifw/rp_Beta_hat_w2_ga.csv",
                               col_names = FALSE))

# Sample alphas
m_alphas <- sample_alphas(10000)
nrow(m_alphas)

l_waifws <- list()
for (i in seq.int(nrow(m_alphas))) {
  v_alphas_within     <- m_alphas[i, 1:4]
  v_alphas_within_off <- 1 - m_alphas[i, 1:4]
  v_alphas_between    <- m_alphas[i, 5:10]
  m_within      <- alpha_matrix_within(v_alphas_within)
  m_within_off  <- alpha_matrix_within(v_alphas_within_off)
  m_between     <- alpha_matrix_between(v_alphas_between)
  waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
    (m_between * waifw_rp)
  l_waifws[[i]]  <- waifw_p
}

# nll <- c(NA, nrow(m_alphas))
# for (i in seq.int(nrow(m_alphas))) {
#   v_alphas <- unname(m_alphas[i, ])
#   v_one_minus_alphas <- 1 - v_alphas
#   m_waifw_alphas <- foi_transition_matrix_alpha(v_alphas)
#   m_one_minus_waifw <- foi_transition_matrix_alpha(v_one_minus_alphas)
#   waifw_p  <- (m_one_minus_waifw * waifw_rp) +
#     (m_waifw_alphas * waifw_assort)
#
#   foi_hat <- waifw_p %*% v_inf
#   nll[i]  <- -1 * sum(dnorm(x = foi_hat,
#                             mean = v_prev_vars$foi,
#                             sd = v_prev_vars$foi_sd,
#                             log = TRUE))
# }
#
# test_foi_hat <- waifw_assort %*% v_inf

### Set up parallel runs
comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}

n <- length(l_waifws)
n_cores <- 80

my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my_cluster)

results <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
                   .packages = c("deSolve", "dplyr",
                                  "rootSolve", "dtplyr",
                                  "purrr", "lubridate")) %dopar% {


  v_parameter <- load_sis_abr_model_params_all(
    v_race = v_race,
    waifw = l_waifws[[i]],
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
  df_foi$run <- i
  df_all_long$run <- i
  results <- list(df_foi = df_foi, df_all_long = df_all_long)
}
parallel::stopCluster(my_cluster)

graph_results_grid <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  df_plot <- df_all %>%
    tidyr::gather(strain, prev, I:Iw) %>%
    dplyr::mutate(strain = ifelse(strain == "I", "Total Prevalence",
                                  ifelse(strain == "Ir", "Resistant Strain",
                                         "Sensitive Strain"))) %>%
    dplyr::mutate(time = time + 2024)

  plot_1 <- ggplot(data = df_plot,
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




# Visualize the alpha posterior
library(psych)
colnames(m_alphas) <- c("a_11", "a_12", "a_13", "a_14", "a_22",
                        "a_23", "a_24", "a_33", "a_34", "a_44")
pdf("results/plot_ch3_alpha_posterior.pdf")
pairs.panels(m_alphas,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = TRUE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4)       # Histograms color
dev.off()

# Boxplot comparing race and ethnicity groups, facet by intervention
plot_box <- function(df_all) {
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
    dplyr::filter(scenario != "No Treatment") %>%
    group_by(race, run) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(
      sq = max(no_trt),
      eff = ifelse(scenario == "18+: One-time Mass Treatment",
                   scen_1 - sq, scen_2 - sq)
    ) %>%
    ungroup()

  plot1 <- ggplot(data = df_wide, aes(x = race, y = eff, fill = race)) +
    geom_boxplot() +
    facet_grid(scenario ~ .,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Race and Ethnicity") + ylab("Net Reduction in Prevalence (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")

  return(plot1)
}

plot_alpha_range_scen <- plot_box(results$df_all_long)
ggsave("results/plot_ch3_alpha_box.pdf", plot_alpha_range_scen)

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

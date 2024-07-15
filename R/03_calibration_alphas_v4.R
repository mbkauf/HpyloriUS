##### This file is for calibrating alpha values #####
### Load libraries
library(readr)
library(rootSolve)
library(tidyverse)
library(APCtools)
library(doParallel)
library(foreach)
library(data.table)
library(DEoptim)
library(dampack)
library(LaplacesDemon)
library(texreg)
library(lhs)


### Functions
# Create matrix of alphas to overlay WAIFWs
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

# Load other functions
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)


# Run calibration of alpha parameters
# Number of age groups
groups <- c(1:80)
n_ages <- length(groups)

# Population growth
q <- 0

# Load demographic and prevalence variables
v_race_calibrate <- c("Hispanic", "NH White", "NH Black", "Other")
v_break <- waifw_breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

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

v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race_calibrate,
  waifw = NULL,
  trt_year = 50000000,
  end_t = 1000,
  ages = groups,
  prob = FALSE,
  p_trans = 0.005
)

v_race_prop <- v_parameter$v_pop / sum(v_parameter$v_pop)
v_race_prop <- rep(v_race_prop, each = length(groups))
v_inf <- v_parameter$v_prev * v_race_prop * v_parameter$v_age_prop
v_foi <- v_parameter$v_foi

med_age <- c(2, 10, 20, 35, 50, 60, 67, 75)
v_med_ages <- c(med_age, med_age + 80, med_age + 160, med_age + 240)

v_p <- seq(from = 0.005, to = 0.05, by = 0.005)

# Define alpha bounds
a_11_bounds <- c(0.072, 0.928)
a_12_bounds <- c(0.046, 0.421)
a_13_bounds <- c(0.023, 0.195)
a_14_bounds <- c(0.019, 0.204)
a_21_bounds <- c(0.201, 0.612)
a_22_bounds <- c(0.119, 0.881)
a_23_bounds <- c(0.112, 0.353)
a_24_bounds <- c(0.121, 0.372)
a_31_bounds <- c(0.008, 0.117)
a_32_bounds <- c(0.008, 0.148)
a_33_bounds <- c(0.026, 0.974)
a_34_bounds <- c(0.002, 0.048)
a_41_bounds <- c(0.003, 0.108)
a_42_bounds <- c(0.002, 0.119)
a_43_bounds <- c(0.000, 0.032)
a_44_bounds <- c(0.069, 0.931)

l_a_bounds <- list(a_11_bounds, a_12_bounds, a_13_bounds, a_14_bounds,
                   a_21_bounds, a_22_bounds, a_23_bounds, a_24_bounds,
                   a_31_bounds, a_32_bounds, a_33_bounds, a_34_bounds,
                   a_41_bounds, a_42_bounds, a_43_bounds, a_44_bounds)


### Test IMIS
library(IMIS)
# sample.prior -- draws samples, and we have already created this
sample.prior <- function(n.samp) {
  # n.samp: the number of samples desired
  draws0 <- randomLHS(n = n.samp, k = 16)
  draws <- matrix(nrow = n.samp, ncol = 16)

  draws <- data.frame(
    a_11 <- qunif(draws0[, 1],  min = l_a_bounds[[1]][1],
                  max = l_a_bounds[[1]][2]),
    a_12 <- qunif(draws0[, 2],  min = l_a_bounds[[2]][1],
                  max = l_a_bounds[[2]][2]),
    a_13 <- qunif(draws0[, 3],  min = l_a_bounds[[3]][1],
                  max = l_a_bounds[[3]][2]),
    a_14 <- qunif(draws0[, 4],  min = l_a_bounds[[4]][1],
                  max = l_a_bounds[[4]][2]),
    a_21 <- qunif(draws0[, 5],  min = l_a_bounds[[5]][1],
                  max = l_a_bounds[[5]][2]),
    a_22 <- qunif(draws0[, 6],  min = l_a_bounds[[6]][1],
                  max = l_a_bounds[[6]][2]),
    a_23 <- qunif(draws0[, 7],  min = l_a_bounds[[7]][1],
                  max = l_a_bounds[[7]][2]),
    a_24 <- qunif(draws0[, 8],  min = l_a_bounds[[8]][1],
                  max = l_a_bounds[[8]][2]),
    a_31 <- qunif(draws0[, 9],  min = l_a_bounds[[9]][1],
                  max = l_a_bounds[[9]][2]),
    a_32 <- qunif(draws0[, 10], min = l_a_bounds[[10]][1],
                  max = l_a_bounds[[10]][2]),
    a_33 <- qunif(draws0[, 11], min = l_a_bounds[[11]][1],
                  max = l_a_bounds[[11]][2]),
    a_34 <- qunif(draws0[, 12], min = l_a_bounds[[12]][1],
                  max = l_a_bounds[[12]][2]),
    a_41 <- qunif(draws0[, 13], min = l_a_bounds[[13]][1],
                  max = l_a_bounds[[13]][2]),
    a_42 <- qunif(draws0[, 14], min = l_a_bounds[[14]][1],
                  max = l_a_bounds[[14]][2]),
    a_43 <- qunif(draws0[, 15], min = l_a_bounds[[15]][1],
                  max = l_a_bounds[[15]][2]),
    a_44 <- qunif(draws0[, 16], min = l_a_bounds[[16]][1],
                  max = l_a_bounds[[16]][2])
  )
  colnames(draws) <- c("a_11", "a_12", "a_13", "a_14",
                       "a_21", "a_22", "a_23", "a_24",
                       "a_31", "a_32", "a_33", "a_34",
                       "a_41", "a_42", "a_43", "a_44")
  return(as.matrix(draws))
}

# prior -- evaluates prior density of a parameter set or sets
l_prior <- function(params, v_bounds = l_a_bounds) {
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params)
  }

  v_ll <- rep(0, nrow(params))
  for (i in 1:ncol(params)) {
    v_ll <- v_ll + dunif(params[, i], min = l_a_bounds[[i]][1],
                         max = l_a_bounds[[i]][2], log = TRUE)
  }

  # ll_a_11 <- dunif(params[,1],  min = l_a_bounds[[1]][1],
  #                  max = l_a_bounds[[1]][2], log = TRUE)
  # ll_a_12 <- dunif(params[,2],  min = l_a_bounds[[2]][1],
  #                  max = l_a_bounds[[2]][2], log = TRUE)
  # ll_a_13 <- dunif(params[,3],  min = l_a_bounds[[3]][1],
  #                  max = l_a_bounds[[3]][2], log = TRUE)
  # ll_a_14 <- dunif(params[,4],  min = l_a_bounds[[4]][1],
  #                  max = l_a_bounds[[4]][2], log = TRUE)
  # ll_a_21 <- dunif(params[,5],  min = l_a_bounds[[5]][1],
  #                  max = l_a_bounds[[5]][2], log = TRUE)
  # ll_a_22 <- dunif(params[,6],  min = l_a_bounds[[6]][1],
  #                  max = l_a_bounds[[6]][2], log = TRUE)
  # ll_a_23 <- dunif(params[,7],  min = l_a_bounds[[7]][1],
  #                  max = l_a_bounds[[7]][2], log = TRUE)
  # ll_a_24 <- dunif(params[,8],  min = l_a_bounds[[8]][1],
  #                  max = l_a_bounds[[8]][2], log = TRUE)
  # ll_a_31 <- dunif(params[,9],  min = l_a_bounds[[9]][1],
  #                  max = l_a_bounds[[9]][2], log = TRUE)
  # ll_a_32 <- dunif(params[,10], min = l_a_bounds[[10]][1],
  #                  max = l_a_bounds[[10]][2], log = TRUE)
  # ll_a_33 <- dunif(params[,11], min = l_a_bounds[[11]][1],
  #                  max = l_a_bounds[[11]][2], log = TRUE)
  # ll_a_34 <- dunif(params[,12], min = l_a_bounds[[12]][1],
  #                  max = l_a_bounds[[12]][2], log  = TRUE)
  # ll_a_41 <- dunif(params[,13], min = l_a_bounds[[13]][1],
  #                  max = l_a_bounds[[13]][2], log = TRUE)
  # ll_a_42 <- dunif(params[,14], min = l_a_bounds[[14]][1],
  #                  max = l_a_bounds[[14]][2], log = TRUE)
  # ll_a_43 <- dunif(params[,15], min = l_a_bounds[[15]][1],
  #                  max = l_a_bounds[[15]][2], log = TRUE)
  # ll_a_44 <- dunif(params[,16], min = l_a_bounds[[16]][1],
  #                  max = l_a_bounds[[16]][2], log = TRUE)
  #
  # ll <- sum(ll_a_11, ll_a_12, ll_a_13, ll_a_14,
  #           ll_a_21, ll_a_22, ll_a_23, ll_a_24,
  #           ll_a_31, ll_a_32, ll_a_33, ll_a_34,
  #           ll_a_41, ll_a_42, ll_a_43, ll_a_44)

  return(v_ll)
}

prior <- function(par_vector) {
  v_prior <- exp(l_prior(par_vector))
  v_prior <- replace(v_prior, v_prior == 0, 1.1e-100)
  return(v_prior)
}


# likelihood -- evaluates likelihood of a parameter set or sets
l_likelihood <- function(params, p_trans = p_i, assort_resamp = m_samp_assort,
                         rp_resamp = m_samp_rp, inf = v_inf) {
  if (is.null(dim(params))) { # If vector, change to matrix
    params <- t(params)
  }
  llik <- rep(0, nrow(params))
  for (j in seq_len(nrow(params))) {
    jj <- tryCatch({
      waifw_assort <- foi_transition_matrix_assort(assort_resamp[j, ])
      waifw_rp     <- foi_transition_matrix_rp(rp_resamp[j, ])

      v_alphas_within     <- as.numeric(params[j, c(1, 6, 11, 16)])
      v_alphas_within_off <- as.numeric(1 - params[j, c(1, 6, 11, 16)])
      v_alphas_between    <- as.numeric(params[j, c(2:5, 7:10, 12:15)])
      m_within      <- alpha_matrix_within(v_alphas_within)
      m_within_off  <- alpha_matrix_within(v_alphas_within_off)
      m_between     <- alpha_matrix_between(v_alphas_between)
      waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
        (m_between * waifw_rp)
      foi_hat <- (p_trans * waifw_p) %*% inf
      foi_hat_med_age <- foi_hat[v_med_ages]
      llik[j]  <- sum(dnorm(x = v_parameter$v_foi[v_med_ages],
                            mean = foi_hat_med_age,
                            sd = v_parameter$v_foi_sd[v_med_ages],
                            log = TRUE))
    }, error = function(e) NA)
    if (is.na(jj)) (llik[j] <- -Inf)
  }
  # llik <- llik - 1750
  return(llik)
}

likelihood <- function(params) {
  exp(l_likelihood(params))
}

# Load first stage calibrated resamples
## Assortative
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

## Random Partnership
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

# Run IMIS
seed_num <- 135135
n_b  <- 1000
n_b0 <- n_b * 10
l_imis_alpha <- list()
for (i in 1:length(v_p)) {
  p_i <- v_p[i]
  imis_assort <- l_m_resamp_assort[[i]]
  imis_rp     <- l_m_resamp_rp[[i]]
  set.seed(seed_num)
  m_samp_assort <- imis_assort[sample(nrow(imis_assort),
                                      size = n_b0,
                                      replace = TRUE), ]
  m_samp_rp <- imis_rp[sample(nrow(imis_rp),
                              size = n_b0,
                              replace = TRUE), ]
  set.seed(seed_num + i)
  imis_res <- IMIS(B = n_b, B.re = 1e3, number_k = 500, D = 0)
  l_imis_alpha[[i]] <- imis_res
  print(i)
}

for(i in 1:length(v_p)) {
  p_i  <- v_p[i]
  # Unique parameter sets
  print(length(unique(l_imis_alpha[[i]]$resample[, 1])))

  # Effective sample size
  print(sum(table(l_imis_alpha[[i]]$resample[, 1]))^2 / sum(table(l_imis_alpha[[i]]$resample[, 1])^2))

  # Max weight
  print(max(table(l_imis_alpha[[i]]$resample[, 1])) / sum(table(l_imis_alpha[[i]]$resample[, 1])))

  # Max log-likelihood
  print(max(test_l_likelihood(l_imis_alpha[[i]]$resample, p_trans = p_i,
                              inf = v_inf_rp,
                              breaks = v_break, w = w2)))
}

### Save
for (i in 1:length(l_imis_alpha)) {
  fwrite(l_imis_alpha[[i]]$stat, paste0("results/waifw/imis_stat_", i,
                                         "_alphas_w2.csv"))
  fwrite(l_imis_alpha[[i]]$resample, paste0("results/waifw/imis_resample_", i,
                                             "_alphas_w2.csv"))
  fwrite(l_imis_alpha[[i]]$center, paste0("results/waifw/imis_center_", i,
                                           "_alphas_w2.csv"))
}

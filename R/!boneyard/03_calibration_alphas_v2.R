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

  m_alphas <- matrix(c(a_11, a_22, a_33, a_44, a_12,
                       a_13, a_14, a_23, a_24, a_34), ncol = 10)
  a_1j <- m_alphas[, 1] + m_alphas[, 5] + m_alphas[, 6] + m_alphas[, 7]
  a_2j <- m_alphas[, 2] + m_alphas[, 5] + m_alphas[, 8] + m_alphas[, 9]
  a_3j <- m_alphas[, 3] + m_alphas[, 6] + m_alphas[, 8] + m_alphas[, 10]
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
alpha_matrix_within <- function(v_alpha, w = w2,
                                waifw_breaks = c(1, 5, 15, 25, 45,
                                                 55, 65, 70, 81)) {
  alpha_11  <- v_alpha[1]
  alpha_22  <- v_alpha[2]
  alpha_33  <- v_alpha[3]
  alpha_44  <- v_alpha[4]

  beta1  <- get_transmission_matrix(betas = rep(alpha_11, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta2  <- get_transmission_matrix(betas = rep(alpha_22, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta3  <- get_transmission_matrix(betas = rep(alpha_33, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta4  <- get_transmission_matrix(betas = rep(alpha_44, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta0  <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)

  Beta <- rbind(cbind(beta1, beta0, beta0, beta0),
                cbind(beta0, beta2, beta0, beta0),
                cbind(beta0, beta0, beta3, beta0),
                cbind(beta0, beta0, beta0, beta4))
  return(Beta)
}

alpha_matrix_between <- function(v_alpha, w = w2,
                                 waifw_breaks = c(1, 5, 15, 25, 45,
                                                  55, 65, 70, 81)) {
  alpha_12  <- v_alpha[1]
  alpha_13  <- v_alpha[2]
  alpha_14  <- v_alpha[3]
  alpha_23  <- v_alpha[4]
  alpha_24  <- v_alpha[5]
  alpha_34  <- v_alpha[6]


  beta1  <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta2  <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta3  <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta4  <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta5  <- get_transmission_matrix(betas = rep(alpha_12, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta6  <- get_transmission_matrix(betas = rep(alpha_13, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta7  <- get_transmission_matrix(betas = rep(alpha_14, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta8  <- get_transmission_matrix(betas = rep(alpha_23, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta9  <- get_transmission_matrix(betas = rep(alpha_24, 8),  #nolint
                                    breaks = waifw_breaks,
                                    w = w,
                                    upper = TRUE)
  beta10 <- get_transmission_matrix(betas = rep(alpha_34, 8),  #nolint
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

# Dirichlet Parameters
mean_a_11 <- (0.660 - 0.044) / 2
mean_a_22 <- (0.648 - 0.138) / 2
mean_a_33 <- (0.671 - 0.005) / 2
mean_a_44 <- (0.661 - 0.000) / 2
mean_a_12 <- (0.686 - 0.048) / 2
mean_a_13 <- (0.299 - 0.016) / 2
mean_a_14 <- (0.349 - 0.028) / 2
mean_a_23 <- (0.388 - 0.012) / 2
mean_a_24 <- (0.409 - 0.049) / 2
mean_a_34 <- (0.157 - 0.029) / 2

sd_a_11 <- (mean_a_11 - 0.044) / 1.96
sd_a_22 <- (mean_a_22 - 0.138) / 1.96
sd_a_33 <- (mean_a_33 - 0.005) / 1.96
sd_a_44 <- (mean_a_44 - 0.000) / 1.96
sd_a_12 <- (mean_a_12 - 0.048) / 1.96
sd_a_13 <- (mean_a_13 - 0.016) / 1.96
sd_a_14 <- (mean_a_14 - 0.028) / 1.96
sd_a_23 <- (mean_a_23 - 0.012) / 1.96
sd_a_24 <- (mean_a_24 - 0.049) / 1.96
sd_a_34 <- (mean_a_34 - 0.029) / 1.96

v_means_1 <- c(mean_a_11, mean_a_12, mean_a_13, mean_a_14)
v_sds_1   <- c(sd_a_11, sd_a_12, sd_a_13, sd_a_14)

v_means_2 <- c(mean_a_22, mean_a_12, mean_a_23, mean_a_24)
v_sds_2   <- c(sd_a_22, sd_a_12, sd_a_23, sd_a_24)

v_means_3 <- c(mean_a_33, mean_a_13, mean_a_23, mean_a_34)
v_sds_3   <- c(sd_a_33, sd_a_13, sd_a_23, sd_a_34)

v_means_4 <- c(mean_a_44, mean_a_14, mean_a_24, mean_a_34)
v_sds_4   <- c(sd_a_44, sd_a_14, sd_a_24, sd_a_34)

v_params_1 <- dirichlet_params(v_means_1, v_sds_1)
v_params_2 <- dirichlet_params(v_means_2, v_sds_2)
v_params_3 <- dirichlet_params(v_means_3, v_sds_3)
v_params_4 <- dirichlet_params(v_means_4, v_sds_4)

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

nl_likelihood <- function(params) {
  v_alphas_within     <- params[1:4]
  v_alphas_within_off <- 1 - params[1:4]
  v_alphas_between    <- params[5:10]
  m_within      <- alpha_matrix_within(v_alphas_within)
  m_within_off  <- alpha_matrix_within(v_alphas_within_off)
  m_between     <- alpha_matrix_between(v_alphas_between)
  waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
    (m_between * waifw_rp)
  foi_hat <- waifw_p %*% v_inf
  nll  <- -1 * sum(dnorm(x = foi_hat,
                         mean = v_prev_vars$foi,
                         sd = v_prev_vars$foi_sd,
                         log = TRUE))
  return(nll)
}

# l_prior <- function(params) {
#   # Block diagonal alphas
#   ll_a_11 <- dunif(params[1], min = 0.044, max = 0.660, log = TRUE)
#   ll_a_22 <- dunif(params[2], min = 0.138, max = 0.648, log = TRUE)
#   ll_a_33 <- dunif(params[3], min = 0.005, max = 0.671, log = TRUE)
#   ll_a_44 <- dunif(params[4], min = 0.000, max = 0.661, log = TRUE)
#
#   # Off block diagonal
#   ll_a_12 <- dunif(params[5], min = 0.048, max = 0.686, log = TRUE)
#   ll_a_13 <- dunif(params[6], min = 0.016, max = 0.299, log = TRUE)
#   ll_a_14 <- dunif(params[7], min = 0.028, max = 0.349, log = TRUE)
#   ll_a_23 <- dunif(params[8], min = 0.012, max = 0.388, log = TRUE)
#   ll_a_24 <- dunif(params[9], min = 0.049, max = 0.409, log = TRUE)
#   ll_a_34 <- dunif(params[10], min = 0.029, max = 0.157, log  = TRUE)
#
#   nll <- -1 * sum(ll_a_11, ll_a_22, ll_a_33, ll_a_44, ll_a_12,
#                   ll_a_13, ll_a_14, ll_a_23, ll_a_24, ll_a_34)
#   return(nll)
# }

l_prior <- function(params) {
  ll_1 <- LaplacesDemon::ddirichlet(params[c(1, 5, 6, 7)],
                                    v_params_1, log = TRUE)
  ll_2 <- LaplacesDemon::ddirichlet(params[c(2, 5, 8, 9)],
                                    v_params_2, log = TRUE)
  ll_3 <- LaplacesDemon::ddirichlet(params[c(3, 6, 8, 10)],
                                    v_params_3, log = TRUE)
  ll_4 <- LaplacesDemon::ddirichlet(params[c(4, 7, 9, 10)],
                                    v_params_4, log = TRUE)

  nll <- -1 * sum(ll_1, ll_2, ll_3, ll_4)
  return(nll)
}

l_posterior <- function(params) {
  return(l_prior(params) + nl_likelihood(params))
}

calibrated_alphas <- DEoptim(
  fn = l_posterior,
  lower = c(0.044, 0.138, 0.005, 0.000, 0.048,
            0.016, 0.028, 0.012, 0.049, 0.029),
  upper = c(0.660, 0.648, 0.671, 0.661, 0.686,
            0.299, 0.349, 0.388, 0.409, 0.157),
  control = list(itermax = 1000,
                parallelType = "auto",
                parVar = list("v_inf", "w2", "l_prior", "l_posterior",
                              "nl_likelihood", "groups", "v_prev_vars",
                              "waifw_assort", "waifw_rp")))
alpha_hess <- numDeriv::hessian(l_posterior, calibrated_alphas$optim$bestmem)
alpha_cov <- solve(alpha_hess)
alpha_se <- sqrt(diag(alpha_cov))

alpha_names <- paste0(paste0("$", paste("\\alpha", c("{11}", "{22}", "{33}",
                                                     "{44}", "{12}", "{13}",
                                                     "{14}", "{23}",
                                                     "{24}", "{34}"),
                                        sep =  "_")), "$")
alpha_tex <- createTexreg(coef.names = alpha_names,
                          coef = calibrated_alphas$optim$bestmem,
                          se = alpha_se,
                          gof = calibrated_alphas$optim$bestval,
                          gof.names = "NLL")
alpha_tex_txt <- texreg(alpha_tex, digits = 4, stars = numeric(0),
                        booktabs = TRUE, single.row = TRUE,
                        custom.model.names = c("DEoptim"))

# Save DEoptim results
fwrite(alpha_cov, file = "results/waifw/alpha_cov.csv")
fwrite(as.data.table(alpha_se), file = "results/waifw/alpha_se.csv")
fwrite(as.data.table(calibrated_alphas$optim$bestmem),
       file = "results/waifw/alpha_se.csv")
write.table(alpha_tex_txt, file = "results/waifw/alpha_tex.txt")

## Try SIR
library(lhs)
library(magrittr)
sample_prior_lhs <- function(n) {
  # n: the number of samples desired
  draws0 <- randomLHS(n=n,k=10)
  draws <- data.frame(
    a_11 <- qunif(draws0[,1], min = 0.044, max = 0.660),
    a_22 <- qunif(draws0[,2], min = 0.138, max = 0.648),
    a_33 <- qunif(draws0[,3], min = 0.005, max = 0.671),
    a_44 <- qunif(draws0[,4], min = 0.000, max = 0.661),
    a_12 <- qunif(draws0[,5], min = 0.048, max = 0.686),
    a_13 <- qunif(draws0[,6], min = 0.016, max = 0.299),
    a_14 <- qunif(draws0[,7], min = 0.028, max = 0.349),
    a_23 <- qunif(draws0[,8], min = 0.012, max = 0.388),
    a_24 <- qunif(draws0[,9], min = 0.049, max = 0.409),
    a_34 <- qunif(draws0[,10], min = 0.029, max = 0.157)
  )
  return(as.matrix(draws))
}

n_samp <- 100000
m_samp <- sample_prior_lhs(n_samp)

v_llik <- rep(NA, n_samp)
for(i in 1:n_samp) {
  v_llik[i] <- -1 * nl_likelihood(as.numeric(m_samp[i,]))
  if(i/100==round(i/100,0)) {
    cat('\r',paste(i/nrow(m_samp)*100,"% done",sep=""))
  }
}

wt <- exp(v_llik - max(v_llik))
id_samp <- sample.int(n_samp, replace=T, prob=wt)
post_sir <- m_samp[id_samp, ]

# Unique parameter sets
length(unique(id_samp))

# Effective sample size
sum(table(id_samp))^2 / sum(table(id_samp)^2)

# Max weight
max(table(id_samp)) / sum(table(id_samp))

### Test IMIS
library(IMIS)
# sample.prior -- draws samples, and we have already created this
sample.prior <- function(n.samp) {
  # n.samp: the number of samples desired
  draws0 <- randomLHS(n=n.samp,k=10)
  draws <- matrix(nrow = n.samp, ncol = 10)
  # draws <- data.frame(
  #   a_11 <- qunif(draws0[,1], min = 0.044, max = 0.660),
  #   a_22 <- qunif(draws0[,2], min = 0.138, max = 0.648),
  #   a_33 <- qunif(draws0[,3], min = 0.005, max = 0.671),
  #   a_44 <- qunif(draws0[,4], min = 0.000, max = 0.661),
  #   a_12 <- qunif(draws0[,5], min = 0.048, max = 0.686),
  #   a_13 <- qunif(draws0[,6], min = 0.016, max = 0.299),
  #   a_14 <- qunif(draws0[,7], min = 0.028, max = 0.349),
  #   a_23 <- qunif(draws0[,8], min = 0.012, max = 0.388),
  #   a_24 <- qunif(draws0[,9], min = 0.049, max = 0.409),
  #   a_34 <- qunif(draws0[,10], min = 0.029, max = 0.157)
  # )
  # colnames(draws) <- c("a_11", "a_22", "a_33", "a_44", "a_12",
  #                      "a_13", "a_14", "a_23", "a_24", "a_34")

  draws[,1] <- qunif(draws0[,1], min = 0.044, max = 0.660)
  draws[,2] <- qunif(draws0[,2], min = 0.138, max = 0.648)
  draws[,3] <- qunif(draws0[,3], min = 0.005, max = 0.671)
  draws[,4] <- qunif(draws0[,4], min = 0.000, max = 0.661)
  draws[,5] <- qunif(draws0[,5], min = 0.048, max = 0.686)
  draws[,6] <- qunif(draws0[,6], min = 0.016, max = 0.299)
  draws[,7] <- qunif(draws0[,7], min = 0.028, max = 0.349)
  draws[,8] <- qunif(draws0[,8], min = 0.012, max = 0.388)
  draws[,9] <- qunif(draws0[,9], min = 0.049, max = 0.409)
  draws[,10] <- qunif(draws0[,10], min = 0.029, max = 0.157)

  return(draws)
}

# prior -- evaluates prior density of a parameter set or sets
l_prior <- function(par_vector) {
  # Block diagonal alphas
  ll_a_11 <- dunif(par_vector[,1], min = 0.044, max = 0.660, log = TRUE)
  ll_a_22 <- dunif(par_vector[,2], min = 0.138, max = 0.648, log = TRUE)
  ll_a_33 <- dunif(par_vector[,3], min = 0.005, max = 0.671, log = TRUE)
  ll_a_44 <- dunif(par_vector[,4], min = 0.000, max = 0.661, log = TRUE)

  # Off block diagonal
  ll_a_12 <- dunif(par_vector[,5], min = 0.048, max = 0.686, log = TRUE)
  ll_a_13 <- dunif(par_vector[,6], min = 0.016, max = 0.299, log = TRUE)
  ll_a_14 <- dunif(par_vector[,7], min = 0.028, max = 0.349, log = TRUE)
  ll_a_23 <- dunif(par_vector[,8], min = 0.012, max = 0.388, log = TRUE)
  ll_a_24 <- dunif(par_vector[,9], min = 0.049, max = 0.409, log = TRUE)
  ll_a_34 <- dunif(par_vector[,10], min = 0.029, max = 0.157, log  = TRUE)

  nll <- sum(ll_a_11, ll_a_22, ll_a_33, ll_a_44, ll_a_12,
             ll_a_13, ll_a_14, ll_a_23, ll_a_24, ll_a_34)
  return(nll)
}

prior <- function(par_vector) {
  exp(l_prior(par_vector))
}


# likelihood -- evaluates likelihood of a parameter set or sets
l_likelihood <- function(par_vector) {
  llik <- rep(0,nrow(par_vector))
  for(j in 1:nrow(par_vector)) {
    jj <- tryCatch( {
      v_alphas_within     <- par_vector[j,1:4]
      v_alphas_within_off <- 1 - par_vector[j,1:4]
      v_alphas_between    <- par_vector[j,5:10]
      m_within      <- alpha_matrix_within(v_alphas_within)
      m_within_off  <- alpha_matrix_within(v_alphas_within_off)
      m_between     <- alpha_matrix_between(v_alphas_between)
      waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
        (m_between * waifw_rp)
      foi_hat <- waifw_p %*% v_inf
      llik[j]  <- sum(dnorm(x = foi_hat,
                             mean = v_prev_vars$foi,
                             sd = v_prev_vars$foi_sd,
                             log = TRUE))
    }, error = function(e) NA)
    if(is.na(jj)) { llik[j] <- -Inf }
  }
  return(llik)
}

likelihood <- function(par_vector) {
  exp(l_likelihood(par_vector))
}

# Run IMIS
set.seed(1234)
imis_res <- IMIS(B = 100, B.re = 1e7, number_k = 400, D = 1)
imis_res <- IMIS_test(B = 10, B.re = 1e4, number_k = 400, D = 1)
# Unique parameter sets
length(unique(imis_res$resample[, 1]))

# Effective sample size
sum(table(imis_res$resample[, 1]))^2 / sum(table(imis_res$resample[, 1])^2)

# Max weight
max(table(imis_res$resample[, 1])) / sum(table(imis_res$resample[, 1]))

likelihood(test)

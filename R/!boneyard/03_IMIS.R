### IMIS Calibration ###
library(IMIS)
library(lhs)
library(psych)
library(rootSolve)
library(parallel)
library(doParallel)
library(doSNOW)

# Functions needed for IMIS
sample.prior <- function(n_samp, n_params = n_params_assort,
                         param_names = params_assort,
                         lb = lb_assort, ub = ub_assort) {
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_params)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_params)
  colnames(m_param_samp) <- param_names
  for (i in 1:n_params) {
    m_param_samp[, i] <- qunif(m_lhs_unit[, i],
                               min = lb[i],
                               max = ub[i])
  }
  return(m_param_samp)
}

f_log_prior <- function(v_params, param_names = params_assort,
                        lb = lb_assort, ub = ub_assort,
                        n_param = n_params_assort){
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_samp <- nrow(v_params)
  colnames(v_params) <- param_names
  lprior <- rep(0, n_samp)
  for (i in 1:n_param){
    lprior <- lprior + dunif(v_params[, i],
                             min = lb[i],
                             max = ub[i],
                             log = T)
  }
  return(lprior)
}

prior <- function(v_params, param_names = params_assort,
                  lb = lb_assort, ub = ub_assort,
                  n_param = n_params_assort) {
  exp(f_log_prior(v_params, param_names, lb, ub, n_param))
}



f_llik <- function(v_params){
  # par_vector: a vector (or matrix) of model parameters
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_samp <- nrow(v_params)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_targets)
  llik_overall <- numeric(n_samp)

  # Set-up parallel runs
  n_cores <- round(parallel::detectCores() * 0.8, 0)
  my_cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(my_cl)
  opts <- list(attachExportEnv = TRUE)
  foreach::foreach(j = 1:n_samp, .combine = "cbind",
                   .export = ls(globalenv()),
                   .packages = c("rootSolve", "MHadaptive"),
                   .options.snow = opts) %dopar% {
                     tmp_llik <- generate_loglik_assort(param = v_params[j, ])
                     for (i in 1:length(tmp_llik)) {
                       v_llik[j, i] <- tmp_llik[[i]]
                     }
                   }
  parallel::stopCluster(my_cl)
  # for(j in 1:n_samp) { # j=1
  # jj <- tryCatch( {
  # tmp_llik <- generate_loglik_assort(param = v_params[j, ])
  # for (i in 1:length(tmp_llik)) {
  #   v_llik[j, i] <- tmp_llik[[i]]
  # }
  # }, error = function(e) NA)
  # if (is.na(jj)) {llik_overall <- -Inf }
  # } # End loop over sampled parameter sets
  # return LLIK
  return(llik_overall)
}

generate_loglik_rp <- function(param, breaks = waifw_breaks, w = w2,
                               upper = TRUE, inf = v_inf) {
  # Create WAIFW matrix
  n_groups  <- length(param) / 8

  if (n_groups == 1) {
    v_betas_1 <- param[1:8]
    beta <- get_transmission_matrix(betas = v_betas_1,
                                    breaks = breaks,
                                    w = w,
                                    upper = upper)
  } else if (n_groups == 3) {
    v_betas_1 <- param[1:8]
    v_betas_2 <- param[9:16]
    v_betas_3 <- param[17:24]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta <- rbind(cbind(beta1, beta3),
                  cbind(beta3, beta2))
  } else if (n_groups == 6) {
    v_betas_1 <- param[1:8]
    v_betas_2 <- param[9:16]
    v_betas_3 <- param[17:24]
    v_betas_4 <- param[25:32]
    v_betas_5 <- param[33:40]
    v_betas_6 <- param[41:48]

    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta5 <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta6 <- get_transmission_matrix(betas = v_betas_6,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta <- rbind(cbind(beta1, beta4, beta5),
                  cbind(beta4, beta2, beta6),
                  cbind(beta5, beta6, beta3))
  } else if (n_groups == 10) {
    v_betas_1  <- param[1:8]
    v_betas_2  <- param[9:16]
    v_betas_3  <- param[17:24]
    v_betas_4  <- param[25:32]
    v_betas_5  <- param[33:40]
    v_betas_6  <- param[41:48]
    v_betas_7  <- param[49:56]
    v_betas_8  <- param[57:64]
    v_betas_9  <- param[65:72]
    v_betas_10 <- param[73:80]

    beta1  <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta2  <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta3  <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta4  <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta5  <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta6  <- get_transmission_matrix(betas = v_betas_6,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta7  <- get_transmission_matrix(betas = v_betas_7,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta8  <- get_transmission_matrix(betas = v_betas_8,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta9  <- get_transmission_matrix(betas = v_betas_9,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)
    beta10 <- get_transmission_matrix(betas = v_betas_10,  #nolint
                                      breaks = breaks,
                                      w = w,
                                      upper = upper)

    beta <- rbind(cbind(beta1, beta5, beta6,  beta7),
                  cbind(beta5, beta2, beta8,  beta9),
                  cbind(beta6, beta8, beta3,  beta10),
                  cbind(beta7, beta9, beta10, beta4))

  } else if (n_groups == 15) {
    v_betas_1 <- param[1:8]
    v_betas_2 <- param[9:16]
    v_betas_3 <- param[17:24]
    v_betas_4 <- param[25:32]
    v_betas_5 <- param[33:40]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta5 <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta <- rbind(cbind(beta1, beta0, beta0, beta0, beta0),
                  cbind(beta0, beta2, beta0, beta0, beta0),
                  cbind(beta0, beta0, beta3, beta0, beta0),
                  cbind(beta0, beta0, beta0, beta4, beta0),
                  cbind(beta0, beta0, beta0, beta0, beta5))
  }

  v_parameter <- load_sis_abr_model_params_all(
    v_race = v_race,
    waifw = beta,
    trt_year = 50000000,
    end_t = 1000,
    ages = seq(1, 80))

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

  v_burn_inf   <- v_burn_state[321:640]
  foi_hat      <- beta %*% v_burn_inf

  ll_prev <- sum(dnorm(x = v_parameter$prevalence,
                       mean = v_burn_inf,
                       sd = v_prev_vars$prevalence_sd,
                       log = TRUE))

  ll_foi <- sum(dnorm(x = v_prev_vars$foi,
                      mean = foi_hat,
                      sd = v_prev_vars$foi_sd,
                      log = TRUE))
  return(list(ll_prev = ll_prev, ll_foi = ll_foi))
}

generate_loglik_assort <- function(param, breaks = waifw_breaks, w = w2,
                                   upper = TRUE, inf = v_inf) {
  # Create WAIFW matrix
  n_groups  <- length(param) / 8

  if (n_groups == 1) {
    v_betas_1 <- param[1:8]
    beta <- get_transmission_matrix(betas = v_betas_1,
                                    breaks = breaks,
                                    w = w,
                                    upper = upper)
  } else if (n_groups == 2) {
    v_betas_1 <- param[1:8]
    v_betas_2 <- param[9:16]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta <- rbind(cbind(beta1, beta0),
                  cbind(beta0, beta2))
  } else if (n_groups == 3) {
    v_betas_1 <- param[1:8]
    v_betas_2 <- param[9:16]
    v_betas_3 <- param[17:24]

    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta <- rbind(cbind(beta1, beta0, beta0),
                  cbind(beta0, beta2, beta0),
                  cbind(beta0, beta0, beta3))
  } else if (n_groups == 4) {
    v_betas_1 <- param[1:8]
    v_betas_2 <- param[9:16]
    v_betas_3 <- param[17:24]
    v_betas_4 <- param[25:32]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta <- rbind(cbind(beta1, beta0, beta0, beta0),
                  cbind(beta0, beta2, beta0, beta0),
                  cbind(beta0, beta0, beta3, beta0),
                  cbind(beta0, beta0, beta0, beta4))
  } else if (n_groups == 5) {
    v_betas_1 <- param[1:8]
    v_betas_2 <- param[9:16]
    v_betas_3 <- param[17:24]
    v_betas_4 <- param[25:32]
    v_betas_5 <- param[33:40]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta5 <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = breaks,
                                     w = w,
                                     upper = upper)
    beta <- rbind(cbind(beta1, beta0, beta0, beta0, beta0),
                  cbind(beta0, beta2, beta0, beta0, beta0),
                  cbind(beta0, beta0, beta3, beta0, beta0),
                  cbind(beta0, beta0, beta0, beta4, beta0),
                  cbind(beta0, beta0, beta0, beta0, beta5))
  }

  v_parameter <- load_sis_abr_model_params_all(
    v_race = v_race,
    waifw = beta,
    trt_year = 50000000,
    end_t = 1000,
    ages = seq(1, 80))

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

  v_burn_inf   <- v_burn_state[321:640] / v_parameter$v_age_prop
  foi_hat      <- beta %*% v_burn_inf

  v_ll_prev <- dnorm(x = v_parameter$v_prev,
                     mean = v_burn_inf,
                     sd = v_parameter$v_prev_sd,
                     log = TRUE)

  v_ll_foi <- dnorm(x = v_parameter$v_foi,
                    mean = foi_hat,
                    sd = v_parameter$v_foi_sd,
                    log = TRUE)

  v_ll_prev <- v_ll_prev[-c(1, 81, 161, 241)]
  v_ll_foi  <- v_ll_foi[-c(1, 81, 161, 241)]

  ll_prev <- sum(v_ll_prev)
  ll_foi  <- sum(v_ll_foi)

  return(list(ll_prev = ll_prev, ll_foi = ll_foi))
}

likelihood <- function(v_params) {
  exp(f_llik(v_params))
}

f_log_post <- function(v_params, param_names = params_assort,
                       lb = lb_assort, ub = ub_assort,
                       n_param = n_params_assort) {
  lpost <- f_log_prior(v_params, param_names, lb, ub, n_param) +
    f_llik(v_params)
  return(lpost)
}

## Variables to run IMIS
# WAIFW matrices
# w1 <- matrix(c(1, 9, 9, 9, 9, 9, 9, 9,
#                0, 2, 9, 9, 9, 9, 9, 9,
#                0, 0, 3, 9, 9, 9, 9, 9,
#                0, 0, 0, 4, 9, 9, 9, 9,
#                0, 0, 0, 0, 5, 9, 9, 9,
#                0, 0, 0, 0, 0, 6, 9, 9,
#                0, 0, 0, 0, 0, 0, 7, 9,
#                0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)

w2 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)

w3 <- matrix(c(1, 1, 1, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)

# v_waifw_structure <- list(w1, w2, w3)
v_waifw_structure <- list(w2, w3)
waifw_breaks      <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
# Age groups
groups <- c(1:80)

# Number of age groups
n_ages <- length(groups)

# Population growth
q <- 0

# IMIS parameters
n_resamp          <- 100

params_assort <- as.character(sapply(list("hisp_hisp", "white_white",
                                          "black_black", "other_other"),
                                     FUN = function(x) paste0(x, 1:8)))

params_rp <- as.character(sapply(list("hisp_hisp", "white_white",
                                      "black_black", "other_other",
                                      "hisp_white", "hisp_black",
                                      "hisp_other", "white_black",
                                      "white_other", "black_other"),
                                 FUN = function(x) paste0(x, 1:8)))

n_params_assort <- length(params_assort)
n_params_rp     <- length(params_rp)

lb_assort <- rep(c(0, 0, 0, 0, 0, 0, 0, 0), 4)
ub_assort <- rep(c(1, 1, 1, 1, 1, 1, 1, 1), 4)

lb_rp <- rep(c(0, 0, 0, 0, 0, 0, 0, 0), 10)
ub_rp <- rep(c(1, 1, 1, 1, 1, 1, 1, 1), 10)

target_names <- c("FOI", "Prevalence")
n_targets    <- length(target_names)

v_race <- c("Hispanic", "NH White", "NH Black", "Other")
## Run IMIS
m_param_samp_assort <- sample.prior(n_samp = 10, n_params = n_params_assort,
                                    param_names = params_assort, lb = lb_assort,
                                    ub = ub_assort)
m_param_samp_rp     <- sample.prior(n_samp = 10, n_params = n_params_rp,
                                    param_names = params_rp, lb = lb_rp,
                                    ub = ub_rp)

fit_imis <- IMIS(B = 100, # the incremental sample size at each iteration of IMIS.
                 B.re = n_resamp, # the desired posterior sample size
                 number_k = 10, # the maximum number of iterations in IMIS.
                 D = 0)

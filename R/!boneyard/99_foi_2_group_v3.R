library(scales)
library(ggplot2)
library(Matrix)
library(parallel)
library(optimParallel)
library(openxlsx)

generate_nloglik_assort <- function(param, breaks, w,
                                    upper = TRUE, inf = v_inf,
                                    v_model_param = v_parameter) {
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

  foi_hat <- beta %*% inf
  nll_foi <- -1 * sum(dnorm(x = foi_hat,
                            mean = v_prev_vars$foi,
                            sd = v_prev_vars$foi_sd,
                            log = TRUE))

  v_model_param$waifw <- beta
  burn_results <- runsteady(y = v_model_param$v_state,
                            times = c(0, 1E5), func = sis_abr_model,
                            parms = v_parameter_burn)
  nll_prev <- -1 * sum(dnorm(x = burn_results$y,
                             mean = v_model_param$prev,
                             sd = v_model_param$prev_sd,
                             log = TRUE))

  out <- nll_foi + nll_prev
  return(out)
}

estimate_beta_assort  <- function(v_beta_0,
                                  w, i,
                                  upper = TRUE,
                                  beta_names = v_beta_names,
                                  v_breaks = waifw_breaks) {
  require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)
  require(optimParallel)

  waifw_ll   <- optimParallel(par = v_beta_0,
                              fn = generate_nloglik_assort,
                              breaks = v_breaks, w = w,
                              inf = v_inf,
                              hessian = TRUE,  method = "L-BFGS-B",
                              lower = rep(0, length(v_beta_0)),
                              upper = rep(50, length(v_beta_0)),
                              control = list(maxit = 1e4))
  beta_llk   <- waifw_ll$value
  v_beta_hat <- as.vector(waifw_ll$par)

  ## Check if HESSIAN is Positive Definite
  ## If not, make covariance Positive Definite
  ## Is Positive Definite?
  if (MHadaptive::isPositiveDefinite(waifw_ll$hessian) == FALSE) {
    print("Hessian is NOT Positive Definite")
    m_hess <- Matrix::nearPD(waifw_ll$hessian)$mat
    beta_hat_cov <- solve(m_hess)

  } else {
    print("Hessian IS Positive Definite")
    print("No additional adjustment to COV matrix")
    beta_hat_cov <- solve(waifw_ll$hessian)
  }

  beta_hat_cov <- as.matrix(beta_hat_cov)
  ### Plot correlation matrix
  beta_hat_cor <- cov2cor(beta_hat_cov)

  ### Compute SE
  beta_hat_se <- sqrt(diag(beta_hat_cov))

  ### Generate big WAIFw matrices
  n_groups  <- length(v_beta_hat) / 8

  if (n_groups == 1) {
    v_betas_1 <- v_beta_hat[1:8]
    beta_hat <- get_transmission_matrix(betas = v_betas_1,
                                        breaks = v_breaks,
                                        w = w,
                                        upper = upper)
  } else if (n_groups == 2) {
    v_betas_1 <- v_beta_hat[1:8]
    v_betas_2 <- v_beta_hat[9:16]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta_hat <- rbind(cbind(beta1, beta0),
                      cbind(beta0, beta2))
  } else if (n_groups == 3) {
    v_betas_1 <- v_beta_hat[1:8]
    v_betas_2 <- v_beta_hat[9:16]
    v_betas_3 <- v_beta_hat[17:24]

    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta_hat <- rbind(cbind(beta1, beta0, beta0),
                      cbind(beta0, beta2, beta0),
                      cbind(beta0, beta0, beta3))
  } else if (n_groups == 4) {
    v_betas_1 <- v_beta_hat[1:8]
    v_betas_2 <- v_beta_hat[9:16]
    v_betas_3 <- v_beta_hat[17:24]
    v_betas_4 <- v_beta_hat[25:32]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta_hat <- rbind(cbind(beta1, beta0, beta0, beta0),
                      cbind(beta0, beta2, beta0, beta0),
                      cbind(beta0, beta0, beta3, beta0),
                      cbind(beta0, beta0, beta0, beta4))
  } else if (n_groups == 5) {
    v_betas_1 <- v_beta_hat[1:8]
    v_betas_2 <- v_beta_hat[9:16]
    v_betas_3 <- v_beta_hat[17:24]
    v_betas_4 <- v_beta_hat[25:32]
    v_betas_5 <- v_beta_hat[33:40]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta5 <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta_hat <- rbind(cbind(beta1, beta0, beta0, beta0, beta0),
                      cbind(beta0, beta2, beta0, beta0, beta0),
                      cbind(beta0, beta0, beta3, beta0, beta0),
                      cbind(beta0, beta0, beta0, beta4, beta0),
                      cbind(beta0, beta0, beta0, beta0, beta5))
  }

  # Create latex code for coefficient table
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat, se = beta_hat_se,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
              beta_hat_se, beta_hat, beta_tex))
}

calibrate_betas_assort <- function(n_age_groups = 8,
                                   v_waifw_structure,
                                   v_race_names,
                                   v_beta0) {
  n_waifw <- length(v_waifw_structure)
  n_betas <- n_age_groups * length(v_race_names)
  ## Beta and WAIFW latex names
  v_race_beta_names <- as.vector(outer(1:8, v_race_names, paste, sep = "_"))
  v_beta_names  <- paste0(paste0("$", paste("\\beta", v_race_beta_names,
                                            sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("w", 1:n_waifw, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_waifw)
  v_beta_llk     <- numeric(n_waifw)
  m_beta_hat     <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  m_beta_hat_se  <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cov <- array(0, dim = c(n_betas, n_betas, n_waifw),
                          dimnames = list(paste0("b", seq(1:n_betas)),
                                          paste0("b", seq(1:n_betas)),
                                          paste0("w", seq(1:n_waifw))))
  a_beta_hat_cor <- a_beta_hat_cov
  v_beta_hat     <- vector("list", n_waifw)
  v_beta_tex     <- vector("list", n_waifw)
  names(v_beta_tex) <- v_waifw_names

  for (j in 1:n_waifw) {
    results_beta <-  estimate_beta_assort(w = v_waifw_structure[[j]],
                                          i = j,
                                          beta_names = v_beta_names,
                                          v_beta_0 = v_beta0)
    v_waifw_ll[[j]]       <- results_beta[[1]]
    v_beta_llk[j]         <- results_beta[[2]]
    m_beta_hat[j, ]       <- results_beta[[3]]
    a_beta_hat_cov[, , j] <- results_beta[[4]]
    a_beta_hat_cor[, , j] <- results_beta[[5]]
    m_beta_hat_se[j, ]    <- results_beta[[6]]
    v_beta_hat[[j]]       <- results_beta[[7]]
    v_beta_tex[[j]]       <- results_beta[8]
  }

  return(list(v_waifw_ll = v_waifw_ll,
              v_beta_llk = v_beta_llk,
              m_beta_hat = m_beta_hat,
              a_beta_hat_cov = a_beta_hat_cov,
              a_beta_hat_cor = a_beta_hat_cor,
              m_beta_hat_se = m_beta_hat_se,
              v_Beta_hat = v_beta_hat,
              v_beta_tex = v_beta_tex))
}

generate_nloglik_rp <- function(param, breaks, w,
                                upper = TRUE, inf = v_inf,
                                v_model_param = v_parameter) {
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

  foi_hat <- beta %*% inf

  out <- -1 * sum(dnorm(x = foi_hat,
                        mean = v_prev_vars$foi,
                        sd = v_prev_vars$foi_sd,
                        log = TRUE))
  return(out)
}

estimate_beta_rp  <- function(v_beta_0,
                              w, i,
                              beta_names = v_beta_names,
                              upper = TRUE,
                              v_breaks = waifw_breaks) {
  require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)
  require(optimParallel)

  waifw_ll   <- optimParallel(par = v_beta_0,
                              fn = generate_nloglik_rp,
                              breaks = v_breaks, w = w,
                              inf = v_inf,
                              hessian = TRUE,  method = "L-BFGS-B",
                              lower = rep(0, length(v_beta_0)),
                              upper = rep(50, length(v_beta_0)),
                              control = list(maxit = 1e4))
  beta_llk   <- waifw_ll$value
  v_beta_hat <- as.vector(waifw_ll$par)

  ## Check if HESSIAN is Positive Definite
  ## If not, make covariance Positive Definite
  ## Is Positive Definite?
  if (MHadaptive::isPositiveDefinite(waifw_ll$hessian) == FALSE) {
    print("Hessian is NOT Positive Definite")
    m_hess <- Matrix::nearPD(waifw_ll$hessian)$mat
    beta_hat_cov <- solve(m_hess)

  } else {
    print("Hessian IS Positive Definite")
    print("No additional adjustment to COV matrix")
    beta_hat_cov <- solve(waifw_ll$hessian)
  }

  beta_hat_cov <- as.matrix(beta_hat_cov)
  ### Plot correlation matrix
  beta_hat_cor <- cov2cor(beta_hat_cov)

  ### Compute SE
  beta_hat_se <- sqrt(diag(beta_hat_cov))

  ### Generate big WAIFW matrices
  n_groups  <- length(v_beta_hat) / 8

  if (n_groups == 1) {
    v_betas_1 <- v_beta_hat[1:8]
    beta_hat <- get_transmission_matrix(betas = v_betas_1,
                                        breaks = v_breaks,
                                        w = w,
                                        upper = upper)
  } else if (n_groups == 3) {
    v_betas_1 <- v_beta_hat[1:8]
    v_betas_2 <- v_beta_hat[9:16]
    v_betas_3 <- v_beta_hat[17:24]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta_hat <- rbind(cbind(beta1, beta3),
                      cbind(beta3, beta2))
  } else if (n_groups == 6) {
    v_betas_1 <- v_beta_hat[1:8]
    v_betas_2 <- v_beta_hat[9:16]
    v_betas_3 <- v_beta_hat[17:24]
    v_betas_4 <- v_beta_hat[25:32]
    v_betas_5 <- v_beta_hat[33:40]
    v_betas_6 <- v_beta_hat[41:48]

    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta5 <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta6 <- get_transmission_matrix(betas = v_betas_6,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta_hat <- rbind(cbind(beta1, beta4, beta5),
                      cbind(beta4, beta2, beta6),
                      cbind(beta5, beta6, beta3))
  } else if (n_groups == 10) {
    v_betas_1  <- v_beta_hat[1:8]
    v_betas_2  <- v_beta_hat[9:16]
    v_betas_3  <- v_beta_hat[17:24]
    v_betas_4  <- v_beta_hat[25:32]
    v_betas_5  <- v_beta_hat[33:40]
    v_betas_6  <- v_beta_hat[41:48]
    v_betas_7  <- v_beta_hat[49:56]
    v_betas_8  <- v_beta_hat[57:64]
    v_betas_9  <- v_beta_hat[65:72]
    v_betas_10 <- v_beta_hat[73:80]

    beta1  <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta2  <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta3  <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta4  <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta5  <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta6  <- get_transmission_matrix(betas = v_betas_6,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta7  <- get_transmission_matrix(betas = v_betas_7,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta8  <- get_transmission_matrix(betas = v_betas_8,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta9  <- get_transmission_matrix(betas = v_betas_9,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)
    beta10 <- get_transmission_matrix(betas = v_betas_10,  #nolint
                                      breaks = v_breaks,
                                      w = w,
                                      upper = upper)

    beta_hat <- rbind(cbind(beta1, beta5, beta6,  beta7),
                      cbind(beta5, beta2, beta8,  beta9),
                      cbind(beta6, beta8, beta3,  beta10),
                      cbind(beta7, beta9, beta10, beta4))

  } else if (n_groups == 15) {
    v_betas_1 <- v_beta_hat[1:8]
    v_betas_2 <- v_beta_hat[9:16]
    v_betas_3 <- v_beta_hat[17:24]
    v_betas_4 <- v_beta_hat[25:32]
    v_betas_5 <- v_beta_hat[33:40]
    beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta4 <- get_transmission_matrix(betas = v_betas_4,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta5 <- get_transmission_matrix(betas = v_betas_5,  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta0 <- get_transmission_matrix(betas = rep(0, 8),  #nolint
                                     breaks = v_breaks,
                                     w = w,
                                     upper = upper)
    beta_hat <- rbind(cbind(beta1, beta0, beta0, beta0, beta0),
                      cbind(beta0, beta2, beta0, beta0, beta0),
                      cbind(beta0, beta0, beta3, beta0, beta0),
                      cbind(beta0, beta0, beta0, beta4, beta0),
                      cbind(beta0, beta0, beta0, beta0, beta5))
  }

  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat, se = beta_hat_se,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
              beta_hat_se, beta_hat, beta_tex))
}

calibrate_betas_rp <- function(n_age_groups = 8,
                               v_waifw_structure,
                               v_race_names,
                               v_beta0) {
  n_waifw <- length(v_waifw_structure)
  n_betas <- length(v_beta0)
  ## Beta and WAIFW latex names
  v_race_beta_names <- as.vector(outer(1:8, v_race_names, paste, sep = "_"))
  v_beta_names  <- paste0(paste0("$", paste("\\beta", v_race_beta_names,
                                            sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("w", 1:n_waifw, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_waifw)
  v_beta_llk     <- numeric(n_waifw)
  m_beta_hat     <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  m_beta_hat_se  <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cov <- array(0, dim = c(n_betas, n_betas, n_waifw),
                          dimnames = list(paste0("b", seq(1:n_betas)),
                                          paste0("b", seq(1:n_betas)),
                                          paste0("w", seq(1:n_waifw))))
  a_beta_hat_cor <- a_beta_hat_cov
  v_beta_hat     <- vector("list", n_waifw)
  v_beta_tex     <- vector("list", n_waifw)
  names(v_beta_tex) <- v_waifw_names

  for (j in 1:n_waifw) {
    results_beta <-  estimate_beta_rp(w = v_waifw_structure[[j]],
                                      i = j,
                                      beta_names = v_beta_names,
                                      v_beta_0 = v_beta0)
    v_waifw_ll[[j]]       <- results_beta[[1]]
    v_beta_llk[j]         <- results_beta[[2]]
    m_beta_hat[j, ]       <- results_beta[[3]]
    a_beta_hat_cov[, , j] <- results_beta[[4]]
    a_beta_hat_cor[, , j] <- results_beta[[5]]
    m_beta_hat_se[j, ]    <- results_beta[[6]]
    v_beta_hat[[j]]       <- results_beta[[7]]
    v_beta_tex[[j]]       <- results_beta[8]
  }

  return(list(v_waifw_ll = v_waifw_ll,
              v_beta_llk = v_beta_llk,
              m_beta_hat = m_beta_hat,
              a_beta_hat_cov = a_beta_hat_cov,
              a_beta_hat_cor = a_beta_hat_cor,
              m_beta_hat_se = m_beta_hat_se,
              v_Beta_hat = v_beta_hat,
              v_beta_tex = v_beta_tex))
}

## Variables needed for calibration
# WAIFW matrices
w1 <- matrix(c(1, 9, 9, 9, 9, 9, 9, 9,
               0, 2, 9, 9, 9, 9, 9, 9,
               0, 0, 3, 9, 9, 9, 9, 9,
               0, 0, 0, 4, 9, 9, 9, 9,
               0, 0, 0, 0, 5, 9, 9, 9,
               0, 0, 0, 0, 0, 6, 9, 9,
               0, 0, 0, 0, 0, 0, 7, 9,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)

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

v_waifw_structure <- list(w1, w2, w3)

## Age groups
groups <- c(1:80)

## Number of age groups
n_ages <- length(groups)

## Population growth
q <- 0

## Load demographic and prevalence variables
v_race_calibrate <- c("Hispanic", "NH White", "NH Black", "Other")

v_demo_vars <- get_demographic_vars_all(v_race = v_race_calibrate,
                                        ages = groups)
v_prev_vars <- get_prevalence_vars_all(v_race = v_race_calibrate,
                                       ages = groups)
v_inf <- v_prev_vars$prevalence * v_demo_vars$v_age_prop
v_foi <- v_prev_vars$foi

## betas and breaks to feed calibration
l_beta0 <- list(rep(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                    length(v_race_calibrate)),
                rep(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                    length(v_race_calibrate)),
                rep(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
                    length(v_race_calibrate)),
                rep(c(0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.05, 0.01),
                    length(v_race_calibrate)))
waifw_breaks      <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

start_time <- Sys.time()
parallel::detectCores()
n_cores <- 80

my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

setDefaultCluster(cl=my_cluster)
clusterExport(my_cluster,
              varlist = list("get_transmission_matrix", "groups",
                             "v_prev_vars"))

x <- 0
for (betas in l_beta0) {
  x <- x + 1
  tmp <- calibrate_betas_assort(v_waifw_structure = v_waifw_structure,
                                v_beta0 = betas,
                                v_race_names = v_race_calibrate)
  assign(paste0("calibrated_assort_", x), tmp)
}
parallel::stopCluster(cl = my_cluster)

end_time <- Sys.time()
end_time - start_time

# Save Assortative Calibration Results
# Latex
all_tex_1 <- list(calibrated_assort_1[["v_beta_tex"]][["$w_1$"]][[1]],
                  calibrated_assort_1[["v_beta_tex"]][["$w_2$"]][[1]],
                  calibrated_assort_1[["v_beta_tex"]][["$w_3$"]][[1]])
all_tex_2 <- list(calibrated_assort_2[["v_beta_tex"]][["$w_1$"]][[1]],
                  calibrated_assort_2[["v_beta_tex"]][["$w_2$"]][[1]],
                  calibrated_assort_2[["v_beta_tex"]][["$w_3$"]][[1]])
all_tex_3 <- list(calibrated_assort_3[["v_beta_tex"]][["$w_1$"]][[1]],
                  calibrated_assort_3[["v_beta_tex"]][["$w_2$"]][[1]],
                  calibrated_assort_3[["v_beta_tex"]][["$w_3$"]][[1]])
all_latex_1 <- texreg(all_tex_1, digits = 4, stars = numeric(0), booktabs = TRUE,
                      single.row = TRUE, custom.model.names = c("W1", "W2", "W3"))
all_latex_2 <- texreg(all_tex_2, digits = 4, stars = numeric(0), booktabs = TRUE,
                      single.row = TRUE, custom.model.names = c("W1", "W2", "W3"))
all_latex_3 <- texreg(all_tex_3, digits = 4, stars = numeric(0), booktabs = TRUE,
                      single.row = TRUE, custom.model.names = c("W1", "W2", "W3"))
write.table(all_latex_1, file = "results/assort_betas_1.txt")
write.table(all_latex_2, file = "results/assort_betas_2.txt")
write.table(all_latex_3, file = "results/assort_betas_3.txt")

# Tables with beta hat
write.csv(calibrated_assort_1$m_beta_hat, file = "results/assort_beta_hat_1.csv")
write.csv(calibrated_assort_2$m_beta_hat, file = "results/assort_beta_hat_2.csv")
write.csv(calibrated_assort_3$m_beta_hat, file = "results/assort_beta_hat_3.csv")

# Tables with beta hat se
write.csv(calibrated_assort_1$m_beta_hat_se, file = "results/assort_beta_hat_se_1.csv")
write.csv(calibrated_assort_2$m_beta_hat_se, file = "results/assort_beta_hat_se_2.csv")
write.csv(calibrated_assort_3$m_beta_hat_se, file = "results/assort_beta_hat_se_3.csv")

# Tables with cov
beta_hat_cov_1 <- as.data.frame(calibrated_assort_1$a_beta_hat_cov)
beta_hat_cov_2 <- as.data.frame(calibrated_assort_2$a_beta_hat_cov)
beta_hat_cov_3 <- as.data.frame(calibrated_assort_3$a_beta_hat_cov)
write.csv(beta_hat_cov_1, file = "results/assort_beta_hat_cov_1.csv")
write.csv(beta_hat_cov_2, file = "results/assort_beta_hat_cov_2.csv")
write.csv(beta_hat_cov_3, file = "results/assort_beta_hat_cov_3.csv")

# Tables with cor
beta_hat_cor_1 <- as.data.frame(calibrated_assort_1$a_beta_hat_cor)
beta_hat_cor_2 <- as.data.frame(calibrated_assort_2$a_beta_hat_cor)
beta_hat_cor_3 <- as.data.frame(calibrated_assort_3$a_beta_hat_cor)
write.csv(beta_hat_cor_1, file = "results/assort_beta_hat_cor_1.csv")
write.csv(beta_hat_cor_2, file = "results/assort_beta_hat_cor_2.csv")
write.csv(beta_hat_cor_3, file = "results/assort_beta_hat_cor_3.csv")

# Calibrated full WAIFWs
write.csv(calibrated_assort_1$v_Beta_hat[[1]],
          file = "results/assort_Beta_hat_w1_1.csv")
write.csv(calibrated_assort_1$v_Beta_hat[[2]],
          file = "results/assort_Beta_hat_w2_1.csv")
write.csv(calibrated_assort_1$v_Beta_hat[[3]],
          file = "results/assort_Beta_hat_w3_1.csv")

write.csv(calibrated_assort_2$v_Beta_hat[[1]],
          file = "results/assort_Beta_hat_w1_2.csv")
write.csv(calibrated_assort_2$v_Beta_hat[[2]],
          file = "results/assort_Beta_hat_w2_2.csv")
write.csv(calibrated_assort_2$v_Beta_hat[[3]],
          file = "results/assort_Beta_hat_w3_2.csv")

write.csv(calibrated_assort_3$v_Beta_hat[[1]],
          file = "results/assort_Beta_hat_w1_3.csv")
write.csv(calibrated_assort_3$v_Beta_hat[[2]],
          file = "results/assort_Beta_hat_w2_3.csv")
write.csv(calibrated_assort_3$v_Beta_hat[[3]],
          file = "results/assort_Beta_hat_w3_3.csv")

## Random partnership WAIFWs
# Load demographic and prevalence variables
v_race_calibrate <- c("Hispanic", "NH White", "NH Black", "Other")

v_demo_vars <- get_demographic_vars_all(v_race = v_race_calibrate,
                                        ages = groups)
v_prev_vars <- get_prevalence_vars_all(v_race = v_race_calibrate,
                                       ages = groups)
v_inf <- v_prev_vars$prevalence * v_demo_vars$v_age_prop
v_foi <- v_prev_vars$foi

l_beta0 <- list(rep(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1), 10))
waifw_breaks      <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

v_race_comb_names <- c("hisp_hisp", "white_white", "black_black", "other_other",
                       "hisp_white", "hisp_black", "hisp_other", "white_black",
                       "white_other", "black_other")
start_time <- Sys.time()
parallel::detectCores()
n_cores <- 80

my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

setDefaultCluster(cl=my_cluster)
clusterExport(my_cluster,
              varlist = list("get_transmission_matrix", "groups",
                             "v_prev_vars"))

x <- 0
for (betas in l_beta0) {
  x <- x + 1
  tmp <- calibrate_betas_rp(v_waifw_structure = v_waifw_structure,
                            v_beta0 = betas,
                            v_race_names = v_race_comb_names)
  assign(paste0("calibrated_rp_", x), tmp)
}
parallel::stopCluster(cl = my_cluster)

end_time <- Sys.time()
end_time - start_time

# Save Random Partnership Calibration Results
# Latex
all_tex_1 <- list(calibrated_rp_1[["v_beta_tex"]][["$w_1$"]][[1]],
                  calibrated_rp_1[["v_beta_tex"]][["$w_2$"]][[1]],
                  calibrated_rp_1[["v_beta_tex"]][["$w_3$"]][[1]])
all_latex_1 <- texreg(all_tex_1, digits = 4, stars = numeric(0), booktabs = TRUE,
                      single.row = TRUE, custom.model.names = c("W1", "W2", "W3"))
write.table(all_latex_1, file = "results/rp_betas_1.txt")

# Tables with beta hat
write.csv(calibrated_rp_1$m_beta_hat, file = "results/rp_beta_hat_1.csv")

# Tables with beta hat se
write.csv(calibrated_rp_1$m_beta_hat_se, file = "results/rp_beta_hat_se_1.csv")

# Tables with cov
beta_hat_cov_1 <- as.data.frame(calibrated_rp_1$a_beta_hat_cov)
write.csv(beta_hat_cov_1, file = "results/rp_beta_hat_cov_1.csv")

# Tables with cor
beta_hat_cor_1 <- as.data.frame(calibrated_rp_1$a_beta_hat_cor)
write.csv(beta_hat_cor_1, file = "results/rp_beta_hat_cor_1.csv")


# Calibrated full WAIFWs
write.csv(calibrated_rp_1$v_Beta_hat[[1]],
          file = "results/rp_Beta_hat_w1_1.csv")
write.csv(calibrated_rp_1$v_Beta_hat[[2]],
          file = "results/rp_Beta_hat_w2_1.csv")
write.csv(calibrated_rp_1$v_Beta_hat[[3]],
          file = "results/rp_Beta_hat_w3_1.csv")

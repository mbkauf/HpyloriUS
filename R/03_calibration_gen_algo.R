library(scales)
library(ggplot2)
library(ggsci)
library(Matrix)
library(parallel)
library(optimParallel)
library(openxlsx)
library(DEoptim)
library(parallelly)
library(numDeriv)

str_stack <- function(x) {
  x %>% str_split("") %>% map(~ .x %>% paste(collapse = "\n")) #nolint
}

foi_transition_matrix_assort <- function(v_betas, w) {
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

foi_transition_matrix_rp <- function(v_betas, w) {
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

generate_nloglik_assort <- function(param, breaks = v_break, w = w2,
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

  foi_hat <- beta %*% inf

  out <- -1 * sum(dnorm(x = foi_hat,
                        mean = v_parameter$v_foi,
                        sd = v_parameter$v_foi_sd,
                        log = TRUE))
  return(out)
}

estimate_beta_assort  <- function(w, i,
                                  v_upper,
                                  beta_names = v_beta_names,
                                  v_breaks = waifw_breaks) {
  require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)
  require(DEoptim)

  if (i == 1) {
    v_upper <- c(16, 2, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 16, 2, 1, 0.5, 0.5, 0.5, 0.5, 0.2,
                 16, 2, 1, 0.2, 0.2, 0.2, 0.2, 0.1,
                 16, 2, 1, 0.2, 0.2, 0.2, 0.2, 0.1)
  }
  waifw_ll   <- DEoptim(fn = generate_nloglik_assort,
                        w = w,
                        lower = rep(0, length(v_upper)),
                        upper = v_upper,
                        control = list(itermax = 1000,
                                       parallelType = "auto",
                                       parVar = list("v_break", "v_inf",
                                                     "get_transmission_matrix",
                                                     "groups", "v_parameter")))

  beta_llk   <- waifw_ll$optim$bestval
  v_beta_hat <- as.vector(waifw_ll$optim$bestmem)

  # Trace of betas
  m_betas    <- waifw_ll$member$bestmemit

  # Trace of nll
  v_nll      <- waifw_ll$member$bestvalit

  # Calculate Hessian to get covariance and SE
  m_hess <- numDeriv::hessian(generate_nloglik_assort, v_beta_hat)
  if (MHadaptive::isPositiveDefinite(m_hess) == FALSE) {
    print("Hessian is NOT Positive Definite")
    m_hess <- as.matrix(Matrix::nearPD(m_hess)$mat)
  }
  m_cov  <- solve(m_hess)
  v_se   <- sqrt(diag(m_cov))

  # Generate big WAIFW matrices
  n_groups  <- length(v_beta_hat) / 8

  Beta_hat <- foi_transition_matrix_assort(v_beta_hat, w)

  # Create latex code for coefficient table
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, m_betas, Beta_hat, beta_tex,
              v_se, m_cov, v_nll))
}

calibrate_betas_assort <- function(v_waifw_structure,
                                   v_waifw_breaks = waifw_breaks,
                                   v_race_names,
                                   v_upper) {

  n_age_groups <- length(v_waifw_breaks) - 1
  n_waifw <- length(v_waifw_structure)
  n_betas <- n_age_groups * length(v_race_names)
  ## Beta and WAIFW latex names
  v_race_beta_names <- as.vector(outer(1:n_age_groups, v_race_names,
                                       paste, sep = "_"))
  v_beta_names  <- paste0(paste0("$", paste("\\beta", v_race_beta_names,
                                            sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("w", 1:n_waifw, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_waifw)
  v_beta_llk     <- numeric(n_waifw)
  m_beta_hat     <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  a_betas <- array(0, dim = c(n_betas, 1000, n_waifw),
                   dimnames = list(paste0("b", seq(1:n_betas)),
                                   paste0("r", seq(1:1000)),
                                   paste0("w", seq(1:n_waifw))))
  m_beta_hat_se  <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cor <- array(0, dim = c(n_betas, n_betas, n_waifw),
                          dimnames = list(paste0("b", seq(1:n_betas)),
                                          paste0("b", seq(1:n_betas)),
                                          paste0("w", seq(1:n_waifw))))

  v_Beta_hat     <- vector("list", n_waifw)
  v_beta_tex     <- vector("list", n_waifw)
  names(v_beta_tex) <- v_waifw_names

  m_nll     <- matrix(0, nrow = 1000, ncol = n_waifw,
                      dimnames = list(paste0("r", seq(1:1000)),
                                      paste0("b", seq(1:n_waifw))))

  for (j in 1:n_waifw) {
    results_beta <-  estimate_beta_assort(w = v_waifw_structure[[j]],
                                          i = j,
                                          beta_names = v_beta_names,
                                          v_upper = v_upper)
    v_waifw_ll[[j]]       <- results_beta[[1]]
    v_beta_llk[j]         <- results_beta[[2]]
    m_beta_hat[j, ]       <- results_beta[[3]]
    a_betas[, , j]        <- results_beta[[4]]
    v_Beta_hat[[j]]       <- results_beta[[5]]
    v_beta_tex[[j]]       <- results_beta[[6]]
    m_beta_hat_se[j, ]    <- results_beta[[7]]
    a_beta_hat_cor[, , j] <- results_beta[[8]]
    m_nll[, j]            <- results_beta[[9]]
  }

  return(list(v_waifw_ll = v_waifw_ll,
              v_beta_llk = v_beta_llk,
              m_beta_hat = m_beta_hat,
              a_betas = a_betas,
              v_Beta_hat = v_Beta_hat,
              v_beta_tex = v_beta_tex,
              m_beta_hat_se = m_beta_hat_se,
              a_beta_hat_cor = a_beta_hat_cor,
              m_nll = m_nll))
}

generate_nloglik_rp <- function(param, breaks = v_break, w = w2,
                                upper = TRUE, inf = v_inf) {
  ### USE THIS IF CALIBRATING LOG OF PARAMETERS ###
  # param <- exp(param)

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
                        mean = v_parameter$v_foi,
                        sd = v_parameter$v_foi_sd,
                        log = TRUE))
  return(out)
}

estimate_beta_rp  <- function(w, i,
                              v_upper,
                              beta_names = v_beta_names,
                              v_breaks = waifw_breaks) {
  require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)
  require(DEoptim)
  if (i == 1) {
    v_upper <- c(4, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 4, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 5, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 5, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 3, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 7, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 6, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 1, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 1, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 7, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1)
    v_lower <- rep(0, length(v_upper))
  } else if (i == 2) {
    v_upper <- c(0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1)
    v_lower <- rep(0, length(v_upper))
  } else {
    v_upper <- c(0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1)
    v_lower <- rep(0, length(v_upper))
  }


  ### USE THIS IF CALIBRATING LOG OF PARAMETERS ###
  # if (i == 1) {
  #   v_lower <- c(-14, -5, -1, -298, -4, -3, -127, -4, -239, -126, -89, -320, -2,
  #                -168, -137, -4, -176, -1, -287, -172, -25, -342, -8, -5, 2,
  #                -345, -387, -5, -4, -298,   -4, -274,   -1, -349, -119,  -67,
  #                -9, -100, -2, -448, 2, -36, -16, -3, -3, -221, -5, -284, -219,
  #                -1, -329, -4, -4, -4, -277, -199, 0, -2, -2, -149, -375, -3,
  #                -225, -17, -7, -348,  -37, -3, -324, -4, -36, -4, 0, -171, -1,
  #                -310, -423, -242, -3, -315)
  #   v_upper <- c(-3, -3, 0, -50, -3, -2, -39, -3, -17, -15, -21, -108, -1, -13,
  #                -64, -2, -43, 0, -51, -30, -17, -106, -5, -4, 3, -115, -53, -3,
  #                -2, -133, -3, -153, 0, -124, -39, -38, -8, -18, -1, -210, 3,
  #                -7, -6, -2, -2, -66, -4, -40, -61, 0, -86, -3, -3, -3, -150,
  #                -37, 1, -1, -1, -50, -150, -2, -61, -7, -2, -90, -8, -2, -203,
  #                -3, -17, -3, 1, -70, 0, -148, -117, -102, -2, -262)
  # } else if (i == 2) {
  #   v_lower <- c(-3, -2, -400, -10, -10, -20, -8, -300, -500, -10, -60, -3.4,
  #                -140, -4, -3.6, -400, -4, -2, -500, -200, -10, -200, -50, -300,
  #                -10, -2, -500, -200, -20, -15, -500, -40, -500, -200, -7, -300,
  #                -400, -350, -450, -4.4, -300, -400, -8, -10, -6, -380, -500,
  #                -8, -450, -500, -4, -400, -430, -450, -10, -300, -400, -300,
  #                -9.5, -500, -350, -280, -500, -500, -4.8, -300, -70, -150,
  #                -4.5, -500, -500, -100, -400, -500, -5, -15, -100, -20, -15,
  #                -500)

  #   v_upper <- c(2, 0, -30, 0, 0, 2, 0, -150, 0, 0, 2, -3, 0, -3, -3.3, 0, 1, 0,
  #                -10, 0, 0, 1, 0, -10, 2, 0, -10, 0, 1, 1, 0, 1, 0, 3, 0, 0, 1,
  #                1, 0, -4, 2, 2, 1, -2, -3, -10, 0, -5.3, 2, 0, -1.5, 0, -50, 2,
  #                -1, 1, -50, 2, -2, 2, 1, 1, 1, -10, -3.5, 1, 1, 1, -4.2, -50,
  #                -20, 0, 1, -20, 0, 1, 1, 1, 0, 0)
  # } else if (i == 3) {
  #   v_lower <- rep(-500, length(v_upper))
  #   v_upper <- rep(3, length(v_upper))
  # }

  waifw_ll   <- DEoptim(fn = generate_nloglik_rp,
                        w = w,
                        lower = v_lower,
                        upper = v_upper,
                        control = list(itermax = 2000,
                                       parallelType = "auto",
                                       parVar = list("v_break", "v_inf",
                                                     "get_transmission_matrix",
                                                     "groups", "v_parameter")))

  beta_llk   <- waifw_ll$optim$bestval
  v_beta_hat <- as.vector(waifw_ll$optim$bestmem)

  # Trace of betas
  m_betas    <- waifw_ll$member$bestmemit

  # Trace of nll
  v_nll      <- waifw_ll$member$bestvalit

  # Calculate Hessian to get covariance and SE
  m_hess <- numDeriv::hessian(generate_nloglik_rp, v_beta_hat)
  if (MHadaptive::isPositiveDefinite(m_hess) == FALSE) {
    print("Hessian is NOT Positive Definite")
    m_hess <- as.matrix(Matrix::nearPD(m_hess)$mat)
  }
  m_cov  <- solve(m_hess)
  v_se   <- sqrt(diag(m_cov))

  ### Generate big WAIFw matrices
  n_groups  <- length(v_beta_hat) / 8

  Beta_hat <- foi_transition_matrix_rp(v_beta_hat, w)

  # Create latex code for coefficient table
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, m_betas, Beta_hat, beta_tex,
              v_se, m_cov, v_nll))
}

calibrate_betas_rp <- function(v_waifw_structure,
                               v_waifw_breaks = waifw_breaks,
                               v_race_names,
                               v_upper) {

  n_age_groups <- length(v_waifw_breaks) - 1
  n_waifw <- length(v_waifw_structure)
  n_betas <- length(v_upper)

  ## Beta and WAIFW latex names
  v_race_beta_names <- as.vector(outer(1:n_age_groups, v_race_names,
                                       paste, sep = "_"))
  v_beta_names  <- paste0(paste0("$", paste("\\beta", v_race_beta_names,
                                            sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("w", 1:n_waifw, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_waifw)
  v_beta_llk     <- numeric(n_waifw)
  m_beta_hat     <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))

  a_betas <- array(0, dim = c(n_betas, 2000, n_waifw),
                   dimnames = list(paste0("b", seq(1:n_betas)),
                                   paste0("r", seq(1:2000)),
                                   paste0("w", seq(1:n_waifw))))
  m_beta_hat_se  <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cov <- array(0, dim = c(n_betas, n_betas, n_waifw),
                          dimnames = list(paste0("b", seq(1:n_betas)),
                                          paste0("b", seq(1:n_betas)),
                                          paste0("w", seq(1:n_waifw))))
  v_Beta_hat     <- vector("list", n_waifw)
  v_beta_tex     <- vector("list", n_waifw)
  names(v_beta_tex) <- v_waifw_names

  m_nll     <- matrix(0, nrow = 2000, ncol = n_waifw,
                      dimnames = list(paste0("r", seq(1:2000)),
                                      paste0("b", seq(1:n_waifw))))

  for (j in 1:n_waifw) {
    results_beta <-  estimate_beta_rp(w = v_waifw_structure[[j]],
                                      i = j,
                                      beta_names = v_beta_names,
                                      v_upper = v_upper)
    v_waifw_ll[[j]]       <- results_beta[[1]]
    v_beta_llk[j]         <- results_beta[[2]]
    m_beta_hat[j, ]       <- results_beta[[3]]
    a_betas[, , j]        <- results_beta[[4]]
    v_Beta_hat[[j]]       <- results_beta[[5]]
    v_beta_tex[[j]]       <- results_beta[[6]]
    m_beta_hat_se[j, ]    <- results_beta[[7]]
    a_beta_hat_cov[, , j] <- results_beta[[8]]
    m_nll[, j]            <- results_beta[[9]]
  }

  return(list(v_waifw_ll = v_waifw_ll,
              v_beta_llk = v_beta_llk,
              m_beta_hat = m_beta_hat,
              a_betas = a_betas,
              v_Beta_hat = v_Beta_hat,
              v_beta_tex = v_beta_tex,
              m_beta_hat_se = m_beta_hat_se,
              a_beta_hat_cov = a_beta_hat_cov,
              m_nll = m_nll))
}

estimate_beta_rp_converge  <- function(w, i,
                                       v_start,
                                       beta_names = v_beta_names,
                                       v_breaks = waifw_breaks) {
  require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)
  require(DEoptim)

  if (i == 1) {
    v_upper <- c(4, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 4, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 5, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 5, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1,
                 3, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 7, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 6, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 1, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 1, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1,
                 7, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.1)
    v_lower <- rep(0, length(v_upper))
    v_init  <- v_start
  } else if (i == 2) {
    v_upper <- c(0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1)
    v_lower <- rep(0, length(v_upper))
    v_init  <- v_start
  } else {
    v_upper <- c(0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                 0.4, 0.4, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1)
    v_lower <- rep(0, length(v_upper))
    v_init  <- v_start
  }

  ### USE THIS IF CALIBRATING LOG OF PARAMETERS ###
  # if (i == 1) {
  #   v_lower <- c(-14, -5, -1, -298, -4, -3, -127, -4, -239, -126, -89, -320, -2,
  #                -168, -137, -4, -176, -1, -287, -172, -25, -342, -8, -5, 2,
  #                -345, -387, -5, -4, -298,   -4, -274,   -1, -349, -119,  -67,
  #                -9, -100, -2, -448, 2, -36, -16, -3, -3, -221, -5, -284, -219,
  #                -1, -329, -4, -4, -4, -277, -199, 0, -2, -2, -149, -375, -3,
  #                -225, -17, -7, -348,  -37, -3, -324, -4, -36, -4, 0, -171, -1,
  #                -310, -423, -242, -3, -315)
  #   v_upper <- c(-3, -3, 0, -50, -3, -2, -39, -3, -17, -15, -21, -108, -1, -13,
  #                -64, -2, -43, 0, -51, -30, -17, -106, -5, -4, 3, -115, -53, -3,
  #                -2, -133, -3, -153, 0, -124, -39, -38, -8, -18, -1, -210, 3,
  #                -7, -6, -2, -2, -66, -4, -40, -61, 0, -86, -3, -3, -3, -150,
  #                -37, 1, -1, -1, -50, -150, -2, -61, -7, -2, -90, -8, -2, -203,
  #                -3, -17, -3, 1, -70, 0, -148, -117, -102, -2, -262)
  # } else if (i == 2) {
  #   v_lower <- c(-3, -2, -400, -10, -10, -20, -8, -300, -500, -10, -60, -3.4,
  #                -140, -4, -3.6, -400, -4, -2, -500, -200, -10, -200, -50, -300,
  #                -10, -2, -500, -200, -20, -15, -500, -40, -500, -200, -7, -300,
  #                -400, -350, -450, -4.4, -300, -400, -8, -10, -6, -380, -500,
  #                -8, -450, -500, -4, -400, -430, -450, -10, -300, -400, -300,
  #                -9.5, -500, -350, -280, -500, -500, -4.8, -300, -70, -150,
  #                -4.5, -500, -500, -100, -400, -500, -5, -15, -100, -20, -15,
  #                -500)

  #   v_upper <- c(2, 0, -30, 0, 0, 2, 0, -150, 0, 0, 2, -3, 0, -3, -3.3, 0, 1, 0,
  #                -10, 0, 0, 1, 0, -10, 2, 0, -10, 0, 1, 1, 0, 1, 0, 3, 0, 0, 1,
  #                1, 0, -4, 2, 2, 1, -2, -3, -10, 0, -5.3, 2, 0, -1.5, 0, -50, 2,
  #                -1, 1, -50, 2, -2, 2, 1, 1, 1, -10, -3.5, 1, 1, 1, -4.2, -50,
  #                -20, 0, 1, -20, 0, 1, 1, 1, 0, 0)
  # } else if (i == 3) {
  #   v_lower <- rep(-500, length(v_upper))
  #   v_upper <- rep(3, length(v_upper))
  # }

  m_pop <- as.matrix(rbind(v_init, v_init, v_init, v_init, v_init, 
                           v_init, v_init, v_init, v_init, v_init))

  waifw_ll   <- DEoptim(fn = generate_nloglik_rp,
                        w = w,
                        lower = v_lower,
                        upper = v_upper,
                        control = list(itermax = 5000,
                                       trace = 50,
                                       # steptol = 200,
                                       # reltol = 0.00001,
                                       parallelType = "auto",
                                       parVar = list("v_break", "v_inf",
                                                     "get_transmission_matrix",
                                                     "groups", "v_parameter")))
  ### USE THIS OTHERWISE
  # waifw_ll   <- DEoptim(fn = generate_nloglik_rp,
  #                       w = w,
  #                       lower = rep(0, length(v_upper)),
  #                       upper = v_upper,
  #                       control = list(itermax = 2000,
  #                                      parallelType = "auto",
  #                                      parVar = list("v_break", "v_inf",
  #                                                    "get_transmission_matrix",
  #                                                    "groups", "v_prev_vars")))

  beta_llk   <- waifw_ll$optim$bestval
  v_beta_hat <- as.vector(waifw_ll$optim$bestmem)

  # Trace of betas
  m_betas    <- waifw_ll$member$bestmemit

  # Trace of nll
  v_nll      <- waifw_ll$member$bestvalit

  # Calculate Hessian to get covariance and SE
  m_hess <- numDeriv::hessian(generate_nloglik_rp, v_beta_hat)
  if (MHadaptive::isPositiveDefinite(m_hess) == FALSE) {
    print("Hessian is NOT Positive Definite")
    m_hess <- as.matrix(Matrix::nearPD(m_hess)$mat)
  }
  m_cov  <- solve(m_hess)
  v_se   <- sqrt(diag(m_cov))

  ### Generate big WAIFw matrices
  n_groups  <- length(v_beta_hat) / 8

  Beta_hat <- foi_transition_matrix_rp(v_beta_hat, w)

  # Create latex code for coefficient table
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, m_betas, Beta_hat, beta_tex,
              v_se, m_cov, v_nll))
}

calibrate_betas_rp_converge <- function(v_waifw_structure,
                                        v_waifw_breaks = waifw_breaks,
                                        v_race_names,
                                        m_start) {

  n_age_groups <- length(v_waifw_breaks) - 1
  n_waifw <- length(v_waifw_structure)
  n_betas <- 80

  ## Beta and WAIFW latex names
  v_race_beta_names <- as.vector(outer(1:n_age_groups, v_race_names,
                                       paste, sep = "_"))
  v_beta_names  <- paste0(paste0("$", paste("\\beta", v_race_beta_names,
                                            sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("w", 1:n_waifw, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_waifw)
  v_beta_llk     <- numeric(n_waifw)
  m_beta_hat     <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))

  a_betas <- array(0, dim = c(n_betas, 5000, n_waifw),
                   dimnames = list(paste0("b", seq(1:n_betas)),
                                   paste0("r", seq(1:5000)),
                                   paste0("w", seq(1:n_waifw))))
  m_beta_hat_se  <- matrix(0, ncol = n_betas, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_betas))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cov <- array(0, dim = c(n_betas, n_betas, n_waifw),
                          dimnames = list(paste0("b", seq(1:n_betas)),
                                          paste0("b", seq(1:n_betas)),
                                          paste0("w", seq(1:n_waifw))))
  v_Beta_hat     <- vector("list", n_waifw)
  v_beta_tex     <- vector("list", n_waifw)
  names(v_beta_tex) <- v_waifw_names

  m_nll     <- matrix(0, nrow = 5000, ncol = n_waifw,
                      dimnames = list(paste0("r", seq(1:5000)),
                                      paste0("b", seq(1:n_waifw))))

  for (j in 1:n_waifw) {
    results_beta <-  estimate_beta_rp_converge(w = v_waifw_structure[[j]],
                                               i = j,
                                               beta_names = v_beta_names,
                                               v_start = m_start[j, ncol(m_beta_hat_rp)-1:ncol(m_beta_hat_rp)])
    
    # if(length(results_beta[[9]]) < 3000) {
    #   results_beta[[9]] <- c(results_beta[[9]], 
    #                          rep(NA, length(results_beta[[9]]) - 3000))
    # }
    v_waifw_ll[[j]]       <- results_beta[[1]]
    v_beta_llk[j]         <- results_beta[[2]]
    m_beta_hat[j, ]       <- results_beta[[3]]
    # a_betas[, , j]        <- results_beta[[4]]
    v_Beta_hat[[j]]       <- results_beta[[5]]
    v_beta_tex[[j]]       <- results_beta[[6]]
    m_beta_hat_se[j, ]    <- results_beta[[7]]
    a_beta_hat_cov[, , j] <- results_beta[[8]]
    # m_nll[, j]            <- results_beta[[9]]
  }

  return(list(v_waifw_ll = v_waifw_ll,
              v_beta_llk = v_beta_llk,
              m_beta_hat = m_beta_hat,
              a_betas = a_betas,
              v_Beta_hat = v_Beta_hat,
              v_beta_tex = v_beta_tex,
              m_beta_hat_se = m_beta_hat_se,
              a_beta_hat_cov = a_beta_hat_cov,
              m_nll = m_nll))
}

get_foi_plot <- function(v_betas, m_cov, waifw, n_samp = 1000) {
  m_foi_hat <- matrix(nrow = n_samp, ncol = length(v_parameter$v_foi))

  for (i in seq.int(n_samp)) {
    rand_betas <- get_waifw_draws(m_cov = m_cov, v_betas = v_betas)
    rand_waifw <- get_new_waifw(old_waifw = waifw, old_betas = v_betas,
                                new_betas = rand_betas)
    tmp_foi <- rand_waifw %*% v_inf
    m_foi_hat[i, ] <- tmp_foi
  }
  df_foi_rand <- as.data.frame(m_foi_hat)
  foi_hat_lb <- apply(df_foi_rand, 2, quantile, probs = 0.025)
  foi_hat_ub <- apply(df_foi_rand, 2, quantile, probs = 0.975)
  foi_hat <- waifw %*% v_inf

  df_foi <- as.data.frame(cbind(v_parameter$v_foi,
                                v_parameter$v_foi_lb,
                                v_parameter$v_foi_ub,
                                c(rep("Hispanic", length(groups)),
                                  rep("NH White", length(groups)),
                                  rep("NH Black", length(groups)),
                                  rep("Other", length(groups))),
                                rep(groups, length(v_race)),
                                foi_hat, foi_hat_lb, foi_hat_ub))

  colnames(df_foi) <- c("foi_est", "foi_est_lb", "foi_est_ub", "race",
                        "age", "foi_model", "foi_model_lb", "foi_model_ub")

  df_foi$foi_est <- as.numeric(df_foi$foi_est)
  df_foi$foi_est_lb <- as.numeric(df_foi$foi_est_lb)
  df_foi$foi_est_ub <- as.numeric(df_foi$foi_est_ub)
  df_foi$foi_model <- as.numeric(df_foi$foi_model)
  df_foi$foi_model_lb <- as.numeric(df_foi$foi_model_lb)
  df_foi$foi_model_ub <- as.numeric(df_foi$foi_model_ub)
  df_foi$age <- as.numeric(df_foi$age)
  df_foi$race2 <- df_foi$race

  plot <- ggplot(data = df_foi) +
    facet_wrap(. ~ race, scales = "free") +
    ggsci::scale_color_nejm(guide = "none") +
    ggsci::scale_fill_nejm(guide = "none") +
    theme_bw() +
    geom_line(aes(x = age, y = foi_est, color = race, alpha = "Observed")) +
    geom_ribbon(aes(x = age, y = foi_est, ymax = foi_est_ub,
                    ymin = foi_est_lb, fill = race),
                alpha = 0.3, color = NA) +
    geom_point(aes(x = age, y = foi_model, color = race,
                   alpha = "Model Predicted")) +
    geom_errorbar(aes(x = age, ymin = foi_hat_lb, ymax = foi_hat_ub,
                      color = race)) +
    theme(legend.position = "right",
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 11),
          strip.text = element_text(size = 14)) +
    xlab("Age") + ylab("Force of Infection (FOI)") +
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
    scale_alpha_manual(name = NULL,
                       values = c(1, 1),
                       breaks = c("Model Predicted", "Observed"),
                       guide = guide_legend(override.aes =
                                              list(linetype = c(0, 1),
                                                   shape = c(16, NA),
                                                   color = "black")))

  return(plot)
}

get_waifw_draws <- function(m_cov, v_betas) {
  chol_mat <- chol(m_cov)
  # vector of random inverse normal
  v_in <- qnorm(runif(nrow(m_cov), 0, 1), 0, 1)
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

get_cov_se <- function(v_betas, assort = TRUE) {
  if (assort == TRUE) {
    m_hess <- numDeriv::hessian(generate_nloglik_assort, v_betas)
    if (MHadaptive::isPositiveDefinite(m_hess) == FALSE) {
      print("Hessian is NOT Positive Definite")
      m_hess <- Matrix::nearPD(m_hess)$mat
    }
    m_cov  <- solve(m_hess)
    v_se   <- sqrt(diag(m_cov))
  } else {
    m_hess <- numDeriv::hessian(generate_nloglik_rp, v_betas)
    if (MHadaptive::isPositiveDefinite(m_hess) == FALSE) {
      print("Hessian is NOT Positive Definite")
      m_hess <- Matrix::nearPD(m_hess)$mat
    }
    m_cov  <- solve(m_hess)
    v_se   <- sqrt(diag(m_cov))
  }
  return(list(m_cov = m_cov, v_se = v_se))
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

v_waifw <- list(w1, w2, w3)

# Age groups
groups <- c(1:80)

# Number of age groups
n_ages <- length(groups)

# Population growth
q <- 0

# Load demographic and prevalence variables
v_race_calibrate <- c("Hispanic", "NH White", "NH Black", "Other")

v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race_calibrate,
  waifw = NULL,
  trt_year = 50000000,
  end_t = 1000,
  ages = groups,
  prob = FALSE
)

v_race_prop <- rep(v_parameter$v_pop, n_ages)
v_inf <- v_parameter$v_prev * v_parameter$v_age_prop * v_race_prop
v_foi <- v_parameter$v_foi

v_break <- waifw_breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

v_race_comb_names <- c("hisp_hisp", "white_white", "black_black", "other_other",
                       "hisp_white", "hisp_black", "hisp_other", "white_black",
                       "white_other", "black_other")

## Assortative WAIFWs
v_upper <- c(0.7, 0.6, 0.3, 0.15, 0.1, 0.1, 0.1, 0.1,
             0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
             0.7, 0.6, 0.3, 0.15, 0.1, 0.1, 0.1, 0.1,
             0.7, 0.6, 0.3, 0.15, 0.1, 0.1, 0.1, 0.1)

start_time <- Sys.time()

calibrated_assort <- calibrate_betas_assort(v_waifw_structure = v_waifw,
                                            v_race_names = v_race_calibrate,
                                            v_upper = v_upper)

end_time <- Sys.time()
end_time - start_time

### Test RP WAIFW ###
v_upper <- c(0.2, 0.2, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05,
             0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
             0.3, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05,
             0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05,
             0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
             0.3, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05,
             0.3, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05,
             0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
             0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
             0.3, 0.2, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05)

start_time <- Sys.time()

calibrated_rp <- calibrate_betas_rp(v_waifw_structure = v_waifw,
                                    v_race_names = v_race_comb_names,
                                    v_upper = v_upper)

end_time <- Sys.time()
end_time - start_time

## Save Assortative Calibration Results
# Latex
all_tex_assort <- list(calibrated_assort$v_beta_tex$`$w_1$`,
                       calibrated_assort$v_beta_tex$`$w_2$`,
                       calibrated_assort$v_beta_tex$`$w_3$`)

all_latex_assort <- texreg(all_tex_assort, digits = 4, stars = numeric(0),
                           booktabs = TRUE, single.row = TRUE,
                           custom.model.names = c("W1", "W2", "W3"))

write.table(all_latex_assort, file = "results/waifw/assort_betas_ga.txt")

# Table with beta hat
write.csv(calibrated_assort$m_beta_hat,
          file = "results/waifw/assort_beta_hat_ga.csv")

# Table with SEs
write.csv(calibrated_assort$m_beta_hat_se,
          file = "results/waifw/assort_beta_se_ga.csv")

# Correlation matrices
beta_hat_cov_w1 <- as.data.frame(calibrated_assort$a_beta_hat_cov[, , 1])
beta_hat_cov_w2 <- as.data.frame(calibrated_assort$a_beta_hat_cov[, , 2])
beta_hat_cov_w3 <- as.data.frame(calibrated_assort$a_beta_hat_cov[, , 3])

write.csv(beta_hat_cov_w1, file = "results/waifw/assort_beta_cov_w1_ga.csv")
write.csv(beta_hat_cov_w2, file = "results/waifw/assort_beta_cov_w2_ga.csv")
write.csv(beta_hat_cov_w3, file = "results/waifw/assort_beta_cov_w3_ga.csv")

# Calibrated full WAIFWs
write.csv(calibrated_assort$v_Beta_hat[[1]],
          file = "results/waifw/assort_Beta_hat_w1_ga.csv")

write.csv(calibrated_assort$v_Beta_hat[[2]],
          file = "results/waifw/assort_Beta_hat_w2_ga.csv")

write.csv(calibrated_assort$v_Beta_hat[[3]],
          file = "results/waifw/assort_Beta_hat_w3_ga.csv")

# Negative Log Likelihood Trace
write.csv(calibrated_assort$m_nll[[1]],
          file = "results/waifw/assort_nll_w1_ga.csv")

write.csv(calibrated_assort$m_nll[[2]],
          file = "results/waifw/assort_nll_w2_ga.csv")

write.csv(calibrated_assort$m_nll[[3]],
          file = "results/waifw/assort_nll_w3_ ga.csv")

## Save RP Calibration Results
# Latex
all_tex_rp <- list(calibrated_rp$v_beta_tex$`$w_1$`,
                   calibrated_rp$v_beta_tex$`$w_2$`,
                   calibrated_rp$v_beta_tex$`$w_3$`)

all_latex_rp <- texreg(all_tex_rp, digits = 4, stars = numeric(0),
                       booktabs = TRUE, single.row = TRUE,
                       custom.model.names = c("W1", "W2", "W3"))

write.table(all_latex_rp, file = "results/waifw/rp_betas_ga.txt")

# Tables with beta hat
write.csv(calibrated_rp$m_beta_hat, file = "results/waifw/rp_beta_hat_ga.csv")

# Table with SEs
write.csv(calibrated_rp$m_beta_hat_se,
          file = "results/waifw/rp_beta_se_ga.csv")

# Covariance matrices
beta_hat_cov_w1 <- as.data.frame(calibrated_rp$a_beta_hat_cov[, , 1])
beta_hat_cov_w2 <- as.data.frame(calibrated_rp$a_beta_hat_cov[, , 2])
beta_hat_cov_w3 <- as.data.frame(calibrated_rp$a_beta_hat_cov[, , 3])

write.csv(beta_hat_cov_w1, file = "results/waifw/rp_beta_cov_w1_ga.csv")
write.csv(beta_hat_cov_w2, file = "results/waifw/rp_beta_cov_w2_ga.csv")
write.csv(beta_hat_cov_w3, file = "results/waifw/rp_beta_cov_w3_ga.csv")

# Calibrated full WAIFWs
write.csv(calibrated_rp$v_Beta_hat[[1]],
          file = "results/waifw/rp_Beta_hat_w1_ga.csv")

write.csv(calibrated_rp$v_Beta_hat[[2]],
          file = "results/waifw/rp_Beta_hat_w2_ga.csv")

write.csv(calibrated_rp$v_Beta_hat[[3]],
          file = "results/waifw/rp_Beta_hat_w3_ga.csv")

# Negative Log Likelihood Trace
write.csv(calibrated_rp$m_nll[[1]],
          file = "results/waifw/rp_nll_w1_ga.csv")

write.csv(calibrated_rp$m_nll[[2]],
          file = "results/waifw/rp_nll_w2_ga.csv")

write.csv(calibrated_rp$m_nll[[3]],
          file = "results/waifw/rp_nll_w3_ ga.csv")


### Refine calibration
m_beta_hat_rp <- fread(file = "results/waifw/rp_beta_hat_ga.csv")
v_beta_hat_rp_w1 <- m_beta_hat_rp[1, ncol(m_beta_hat_rp)-1:ncol(m_beta_hat_rp)]
v_beta_hat_rp_w2 <- m_beta_hat_rp[2, ncol(m_beta_hat_rp)-1:ncol(m_beta_hat_rp)]
v_beta_hat_rp_w3 <- m_beta_hat_rp[3, ncol(m_beta_hat_rp)-1:ncol(m_beta_hat_rp)]

calibrated_rp_converge <- calibrate_betas_rp_converge(
  v_waifw_structure = v_waifw,
  v_race_names = v_race_comb_names,
  m_start = m_beta_hat_rp
)

## Check steady state prevalence
library(rootSolve)

# Assortative
v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race_calibrate,
  waifw = calibrated_assort$v_Beta_hat[[2]],
  trt_year = 50000000,
  end_t = 1000,
  ages = seq(1, 80),
  prob = FALSE
)

# Make copy of v_parameter and update for burn-in
v_parameter_burn         <- v_parameter
v_parameter_burn$sigma   <- 0
v_parameter_burn$v_alpha <- 0
v_parameter_burn$v_psi   <- 0

# Run burn-in period and store starting states
burn_results <- runsteady(y = v_parameter_burn$v_state,
                          times = c(0, 1E5), func = sis_abr_model,
                          parms = v_parameter_burn)
v_burn_state_assort <- burn_results$y

## RP
v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race_calibrate,
  waifw = calibrated_rp$v_Beta_hat[[2]],
  trt_year = 50000000,
  end_t = 1000,
  ages = seq(1, 80),
  prob = FALSE
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
v_burn_state_rp <- burn_results$y

## Check differences from predicted prevalence
v_diff_assort <- v_parameter$v_state - v_burn_state_assort
v_diff_rp     <- v_parameter$v_state - v_burn_state_rp
v_diff_steady <- v_burn_state_assort - v_burn_state_rp



## Create FOI Plots
v_race <- v_race_calibrate
assort_w2_hess <- numDeriv::hessian(generate_nloglik_assort,
                                    calibrated_assort$m_beta_hat[2, ])
m_cov_asort_2 <- solve(assort_w2_hess)



### RP ###
# W1
cov_w1      <- as.matrix(read.csv("results/waifw/rp_beta_cov_w1_ga.csv"))[-1, -1]
betas_w1    <- read.csv("results/waifw/rp_beta_hat_ga.csv")
betas_w1    <- as.numeric(betas_w1[1, 2:81])
waifw_w1_rp <- foi_transition_matrix_rp(exp(betas_w1), w=w1)

# W2
cov_w2      <- as.matrix(read.csv("results/waifw/rp_beta_cov_w2_ga.csv"))[-1, -1]
betas_w2    <- read.csv("results/waifw/rp_beta_hat_ga.csv")
betas_w2    <- as.numeric(betas_w2[2, 2:81])
waifw_w2_rp <- foi_transition_matrix_rp(exp(betas_w2), w=w2)

# W3
cov_w3      <- as.matrix(read.csv("results/waifw/rp_beta_cov_w3_ga.csv"))[-1, -1]
betas_w3    <- read.csv("results/waifw/rp_beta_hat_ga.csv")
betas_w3    <- as.numeric(betas_w3[3, 2:81])
waifw_w3_rp <- foi_transition_matrix_rp(exp(betas_w3), w=w3)

plot_foi_rp_w1 <- get_foi_plot(waifw = waifw_w1_rp) + ggtitle("RP - W1")
plot_foi_rp_w1

plot_foi_rp_w2 <- get_foi_plot(waifw = waifw_w2_rp) + ggtitle("RP - W2")
plot_foi_rp_w2

plot_foi_rp_w3 <- get_foi_plot(waifw = waifw_w3_rp) + ggtitle("RP - W3")
plot_foi_rp_w3

### Assort ###
waifw_w1_a <- as.matrix(read.csv("results/waifw/assort_Beta_hat_w1_ga.csv"))[-321, -1]
plot_foi_assort_w1 <- get_foi_plot(waifw = waifw_w1_a) + ggtitle("Assort - W1")
plot_foi_assort_w1



plot_foi_assort_w1 <- get_foi_plot(calibrated_assort$v_Beta_hat[[1]]) +
  ggtitle("Assortative - W1")
plot_foi_assort_w2 <- get_foi_plot(v_betas = calibrated_assort$m_beta_hat[2, ],
                                   m_cov = m_cov_asort_2,
                                   waifw = calibrated_assort$v_Beta_hat[[2]]) +
  ggtitle("Assortative - W2")

ggsave(filename = "results/FOI_validation_assort_w1_ga.pdf",
       plot = plot_foi_assort_w1, width = 10, height = 4)

ggsave(filename = "results/FOI_validation_rp_ga.pdf",
       plot = plot_foi_valid, width = 10, height = 4)


m_assort_w1 <- calibrated_assort$v_waifw_ll[[1]][["member"]][["bestmemit"]]
m_boot_assort_w1 <- m_assort_w1[sample(nrow(m_assort_w1),
                                       size = 10000, replace = TRUE), ]
v_boot_assort_w1_sd <- apply(m_boot_assort_w1, 2, sd)
m_boot_assort_w1_cor <- cor(m_boot_assort_w1)




### Generate visualizations of WAIFW matrices ###
simple_waifw <- function(w, v_betas) {
  if (w == "w1") {
    waifw <- matrix(c(v_betas[1], 0, 0, 0, 0, 0, 0, 0,
                      0, v_betas[2], 0, 0, 0, 0, 0, 0,
                      0, 0, v_betas[3], 0, 0, 0, 0, 0,
                      0, 0, 0, v_betas[4], 0, 0, 0, 0,
                      0, 0, 0, 0, v_betas[5], 0, 0, 0,
                      0, 0, 0, 0, 0, v_betas[6], 0, 0,
                      0, 0, 0, 0, 0, 0, v_betas[7], 0,
                      0, 0, 0, 0, 0, 0, 0, v_betas[8]),
                    ncol = 8, byrow = TRUE)
  } else if (w == "w2") {
    waifw <- matrix(c(v_betas[1], v_betas[1], v_betas[3], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[1], v_betas[2], v_betas[3], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[3], v_betas[3], v_betas[3], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[4], v_betas[4], v_betas[4], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[5], v_betas[5], v_betas[5], v_betas[5],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[6], v_betas[6], v_betas[6], v_betas[6],
                      v_betas[6], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[7], v_betas[7], v_betas[7], v_betas[7],
                      v_betas[7], v_betas[7], v_betas[7], v_betas[8],
                      v_betas[8], v_betas[8], v_betas[8], v_betas[8],
                      v_betas[8], v_betas[8], v_betas[8], v_betas[8]),
                    ncol = 8, byrow = TRUE)
  } else if (w == "w3") {
    waifw <- matrix(c(v_betas[1], v_betas[1], v_betas[1], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[1], v_betas[2], v_betas[3], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[1], v_betas[3], v_betas[3], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[4], v_betas[4], v_betas[4], v_betas[4],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[5], v_betas[5], v_betas[5], v_betas[5],
                      v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[6], v_betas[6], v_betas[6], v_betas[6],
                      v_betas[6], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[7], v_betas[7], v_betas[7], v_betas[7],
                      v_betas[7], v_betas[7], v_betas[7], v_betas[8],
                      v_betas[8], v_betas[8], v_betas[8], v_betas[8],
                      v_betas[8], v_betas[8], v_betas[8], v_betas[8]),
                    ncol = 8, byrow = TRUE)
  }
  rownames(waifw) <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")
  colnames(waifw) <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

  return(waifw)
}
## Read in WAIFW matrices
# Assortative
m_beta_hat_assort <- read.csv("results/waifw/assort_beta_hat_ga.csv")[-4,-1]
m_waifw_assort_w1 <- as.matrix(read.csv("results/waifw/assort_Beta_hat_w1_ga.csv"))[-321,-1]
m_waifw_assort_w2 <- as.matrix(read.csv("results/waifw/assort_Beta_hat_w2_ga.csv"))[-321:-322,-1]
m_waifw_assort_w3 <- as.matrix(read.csv("results/waifw/assort_Beta_hat_w3_ga.csv"))[-321:-322,-1]

colnames(m_waifw_assort_w1) <-
  colnames(m_waifw_assort_w2) <-
  colnames(m_waifw_assort_w1) <- as.vector(outer(1:80, v_race_calibrate,
                                                 paste, sep = "_"))
rownames(m_waifw_assort_w1) <-
  rownames(m_waifw_assort_w2) <-
  rownames(m_waifw_assort_w1) <- as.vector(outer(1:80, v_race_calibrate,
                                                 paste, sep = "_"))
labels <- as.vector(outer(1:80, v_race_calibrate, paste, sep = "-"))
# RP
m_beta_hat_rp <- read.csv("results/waifw/rp_beta_hat_ga.csv")[-4,]
m_waifw_rp_w1 <- as.matrix(read.csv("results/waifw/rp_Beta_hat_w1_ga.csv"))[-321:-323,]
m_waifw_rp_w2 <- as.matrix(read.csv("results/waifw/rp_Beta_hat_w2_ga.csv"))[-321,]
m_waifw_rp_w3 <- as.matrix(read.csv("results/waifw/rp_Beta_hat_w3_ga.csv"))[-321,]

# Heatmap of WAIFW
library(plotly)
fig <- plot_ly(x = labels, y = labels, z = m_waifw_assort_w1, type = "heatmap") # %>% add_surface()
fig

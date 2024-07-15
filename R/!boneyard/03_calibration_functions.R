generate_nloglik <- function(betas, breaks, w, upper = TRUE, inf = v_inf) {
  # Create WAIFW matrix
  beta <- get_transmission_matrix(betas = betas, breaks = breaks, w = w, #nolint
                                  upper = upper)
  foi_hat <- beta %*% inf
  out <- -1 * sum(dnorm(x = foi_hat, mean = v_prev_vars$foi,
                        sd = v_prev_vars$foi_sd, log = TRUE))

  return(out)
}


estimate_beta  <- function(n_age_groups, w, i, beta_names = v_beta_names,  #nolint
                           v_breaks = waifw_breaks, group_names = groups,
                           beta0) {
  require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)

  ### Run optimization
  waifw_ll   <- optim(par = beta0, fn = generate_nloglik,
                      breaks = v_breaks, w = w,
                      inf = v_inf,
                      hessian = TRUE,  method = "L-BFGS-B",
                      lower = rep(0, n_age_groups),
                      upper = rep(100, n_age_groups))
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
  beta_hat <- get_transmission_matrix(betas = v_beta_hat,  #nolint
                                      breaks = v_breaks,
                                      w = w)
  # print(length(v_beta_hat))
  # print(length(beta_names))
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat, se = beta_hat_se,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
              beta_hat_se, beta_hat, beta_tex))
}


calibrate_betas <- function(n_age_groups, n_waifw, v_waifw_structure, beta0) {
  ## Beta and WAIFW latex names
  v_beta_names  <- paste0(paste0("$", paste("\\beta", 1:8, sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("w", 1:n_waifw, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_waifw)
  v_beta_llk     <- numeric(n_waifw)
  m_beta_hat     <- matrix(0, ncol = n_age_groups, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_age_groups))))
  m_beta_hat_se  <- matrix(0, ncol = n_age_groups, nrow = n_waifw,
                           dimnames = list(paste0("w", seq(1:n_waifw)),
                                           paste0("b", seq(1:n_age_groups))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cov <- array(0, dim = c(n_age_groups, n_age_groups, n_waifw),
                          dimnames = list(paste0("b", seq(1:n_age_groups)),
                                          paste0("b", seq(1:n_age_groups)),
                                          paste0("w", seq(1:n_waifw))))
  a_beta_hat_cor <- a_beta_hat_cov
  v_beta_hat     <- vector("list", n_waifw)
  v_beta_tex     <- vector("list", n_waifw)
  names(v_beta_tex) <- v_waifw_names

  for (i in 1:n_waifw) {
    results_beta <-  estimate_beta(n_age_groups = n_age_groups,
                                   w = v_waifw_structure[[i]],
                                   i = i,
                                   beta_names = v_beta_names,
                                   beta0 = beta0)
    v_waifw_ll[[i]]     <- results_beta[[1]]
    v_beta_llk[i]       <- results_beta[[2]]
    m_beta_hat[i, ]       <- results_beta[[3]]
    a_beta_hat_cov[, , i] <- results_beta[[4]]
    a_beta_hat_cor[, , i] <- results_beta[[5]]
    m_beta_hat_se[i, ]    <- results_beta[[6]]
    v_beta_hat[[i]]       <- results_beta[[7]]
    v_beta_tex[[i]]       <- results_beta[8]
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


### Plot calibration
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



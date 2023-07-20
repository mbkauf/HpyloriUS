get_transmission_matrix_mixing <- function(betas, breaks, W, upper = TRUE, group.names = groups) {

  ### Number of age groups = number of data to estimate
  n.a <- (length(breaks)-1)
  n.params <- (length(breaks)-1)*4

  ### Age.max
  age.max <- breaks[length(breaks)]-1

  ### Age group labels
  age.names <- paste(breaks[-length(breaks)], breaks[-1]-1, sep = "-")
  age.names <- rep(age.names, 4)

  if (max(W) > n.params ) {betas <- c(betas, 0)}
  ## Complete index matrix from upper triangular index matrix
  complete.index <- function(index.upper){
    index <- index.upper + t(index.upper)
    diag(index) <- diag(index.upper)
    return(index)
  }

  if (upper == TRUE) {
    index <- complete.index(index.upper = W)

    ### Generate age-group waifw
    bij <- matrix(betas[index], ncol = n.params) #, dimnames = list(age.names, age.names))

    ### Generate 1-year age-specific waifw
    index.full <- index[rep(rep(1:n.a, diff(breaks)), 4),
                        rep(rep(1:n.a, diff(breaks)), 4)]
    colnames(index.full) <- rep(age.names, rep(diff(breaks), 4))
    m_waifw <- matrix(betas[index.full], ncol = (length(group.names))*4)

  } else {
    index <- W
    ### Generate age-group waifw
    bij <- matrix(betas[index], ncol = n.params, dimnames = list(age.names, age.names))
    ### Generate 1-year age-specific waifw
    index.full <- index[rep(1:n.a, diff(breaks)),
                        rep(1:n.a, diff(breaks))]
    colnames(index.full) <- rep(age.names, diff(breaks))
    m_waifw <- matrix(betas[index.full], ncol = (length(group.names)))
  }

  colnames(m_waifw) <- rownames(m_waifw) <- rep(group.names, 4)
  return(m_waifw)
}

generate_nloglik_mix <- function(betas, breaks, W, upper = TRUE, inf = v_inf){
  # Create WAIFW matrix
  Beta <- get_transmission_matrix_mixing(betas = betas, breaks = breaks, W = W,
                                  upper = upper)
  foi_hat <- Beta %*% inf
  out <- -1 * sum(dnorm(x = foi_hat, mean = foi_mixing,
                        sd = foi_sd_mixing, log = T))
  return(out)
}


estimate_Beta  <- function(n_age_groups, W, i, beta_names = v_beta_names,
                           v_breaks = waifw.breaks, group.names = groups) {
  require(MHadaptive)
  require(texreg)

  ### Run optimization
  waifw_ll   <- optim(par = beta0, fn = generate_nloglik_mix,
                      breaks = v_breaks, W = W,
                      hessian = T,  method = "L-BFGS-B",
                      lower = rep(0, n_age_groups), upper = rep(15, n_age_groups))
  beta_llk   <- waifw_ll$value
  v_beta_hat <- as.vector(waifw_ll$par)
  ### Check if Hessian is positive definite
  print(paste0("Is W = ", i, " positive definite? ",
               MHadaptive::isPositiveDefinite(waifw_ll$hessian)))
  beta_hat_cov <- solve(waifw_ll$hessian)
  ### Plot correlation matrix
  beta_hat_cor <- cov2cor(beta_hat_cov)

  ### Compute SE
  beta_hat_se <- sqrt(diag(beta_hat_cov))

  ### Generate big WAIFW matrices
  Beta_hat <- get_transmission_matrix_mixing(betas = v_beta_hat,
                                      breaks = v_breaks,
                                      W = W)
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat, se = beta_hat_se,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
              beta_hat_se, Beta_hat, beta_tex))
}

calibrate_Betas_mixing <- function(n_age_groups, n_Betas, v_waifw_structure) {
  ## Beta and WAIFW latex names
  v_beta_names <- paste0(paste0("$", paste("\\beta", 1:8, sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("W", 1:n_Betas, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_Betas)
  v_beta_llk     <- numeric(n_Betas)
  m_beta_hat     <- matrix(0, ncol = n_age_groups, nrow = n_Betas,
                           dimnames = list(paste0("W", seq(1:n_Betas)),
                                           paste0("B", seq(1:n_age_groups))))
  m_beta_hat_se  <- matrix(0, ncol = n_age_groups, nrow = n_Betas,
                           dimnames = list(paste0("W", seq(1:n_Betas)),
                                           paste0("B", seq(1:n_age_groups))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cov <- array(0, dim = c(n_age_groups, n_age_groups, n_Betas),
                          dimnames = list(paste0("B", seq(1:n_age_groups)),
                                          paste0("B", seq(1:n_age_groups)),
                                          paste0("W", seq(1:n_Betas))))
  a_beta_hat_cor <- a_beta_hat_cov
  v_Beta_hat     <- vector("list", n_Betas)
  v_beta_tex     <- vector("list", n_Betas)
  names(v_beta_tex) <- v_waifw_names

  for (i in 1:n_Betas) {
    results_beta <-  estimate_Beta(n_age_groups = n_age_groups,
                                   W = v_waifw_structure[[i]],
                                   i = i,
                                   beta_names = v_beta_names)
    v_waifw_ll[[i]] <- results_beta[[1]]
    v_beta_llk[i]   <- results_beta[[2]]
    m_beta_hat[i,]  <- results_beta[[3]]

    a_beta_hat_cov[,,i] <- results_beta[[4]]
    a_beta_hat_cor[,,i] <- results_beta[[5]]
    m_beta_hat_se[i,]   <- results_beta[[6]]
    v_Beta_hat[[i]]     <- results_beta[[7]]
    v_beta_tex[[i]]     <- results_beta[8]
  }

  return(list(v_waifw_ll = v_waifw_ll,
              v_beta_llk = v_beta_llk,
              m_beta_hat = m_beta_hat,
              a_beta_hat_cov = a_beta_hat_cov,
              a_beta_hat_cor = a_beta_hat_cor,
              m_beta_hat_se = m_beta_hat_se,
              v_Beta_hat = v_Beta_hat,
              v_beta_tex = v_beta_tex))
}

W1_mixing <- matrix(c(seq(1,32),
                      c(0, seq(2,32)),
                      c(rep(0, 2),  seq(3,32)),
                      c(rep(0, 3),  seq(4,32)),
                      c(rep(0, 4),  seq(5,32)),
                      c(rep(0, 5),  seq(6,32)),
                      c(rep(0, 6),  seq(7,32)),
                      c(rep(0, 7),  seq(8,32)),
                      c(rep(0, 8),  seq(9,32)),
                      c(rep(0, 9),  seq(10,32)),
                      c(rep(0, 10), seq(11,32)),
                      c(rep(0, 11), seq(12,32)),
                      c(rep(0, 12), seq(13,32)),
                      c(rep(0, 13), seq(14,32)),
                      c(rep(0, 14), seq(15,32)),
                      c(rep(0, 15), seq(16,32)),
                      c(rep(0, 16), seq(17,32)),
                      c(rep(0, 17), seq(18,32)),
                      c(rep(0, 18), seq(19,32)),
                      c(rep(0, 19), seq(20,32)),
                      c(rep(0, 20), seq(21,32)),
                      c(rep(0, 21), seq(22,32)),
                      c(rep(0, 22), seq(23,32)),
                      c(rep(0, 23), seq(24,32)),
                      c(rep(0, 24), seq(25,32)),
                      c(rep(0, 25), seq(26,32)),
                      c(rep(0, 26), seq(27,32)),
                      c(rep(0, 27), seq(28,32)),
                      c(rep(0, 28), seq(29,32)),
                      c(rep(0, 29), seq(30,32)),
                      c(rep(0, 30), seq(31,32)),
                      c(rep(0, 31), 32)), ncol = 32, byrow = T)

beta0 <- rep(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1), 4)
waifw.breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_age_group_names <- rep(c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64", "65-69", "70-80"), 4)

Beta <- get_transmission_matrix_mixing(betas = beta0,
                                       breaks = waifw.breaks,
                                       W = W1_mixing,
                                       group.names = groups)

## Combine prevalence
v_prev_vars_white <- get_prevalence_vars(spec_groups = groups,
                                         Race = "Non-Hispanic White")
v_prev_vars_black <- get_prevalence_vars(spec_groups = groups,
                                         Race = "Non-Hispanic Black")
v_prev_vars_hisp  <- get_prevalence_vars(spec_groups = groups,
                                         Race = "Hispanic")
v_prev_vars_other <- get_prevalence_vars(spec_groups = groups,
                                         Race = "Other")
v_demo_white <- get_demographic_vars(spec_groups = groups, Race = "Non-Hispanic White")
v_demo_black <- get_demographic_vars(spec_groups = groups, Race = "Non-Hispanic Black")
v_demo_hisp  <- get_demographic_vars(spec_groups = groups, Race = "Hispanic")
v_demo_other <- get_demographic_vars(spec_groups = groups, Race = "Other")

prev_white <- v_prev_vars_white$prevalence * v_demo_white$v_age_prop
prev_black <- v_prev_vars_black$prevalence * v_demo_black$v_age_prop
prev_hisp  <- v_prev_vars_hisp$prevalence * v_demo_hisp$v_age_prop
prev_other <- v_prev_vars_other$prevalence * v_demo_other$v_age_prop

v_inf <- c(prev_white*0.25, prev_black*0.25, prev_hisp*0.25, prev_other*0.25)

foi_mixing    <- c(v_prev_vars_white$foi, v_prev_vars_black$foi, v_prev_vars_hisp$foi, v_prev_vars_other$foi)
foi_sd_mixing <- c(v_prev_vars_white$foi_sd, v_prev_vars_black$foi_sd, v_prev_vars_hisp$foi_sd, v_prev_vars_other$foi_sd)

n_age_groups <- (length(waifw.breaks) - 1)*4

Beta_mixing <- calibrate_Betas_mixing(n_age_groups = n_age_groups,
                                      n_Betas = 1,
                                      v_waifw_structure = list(W1_mixing))

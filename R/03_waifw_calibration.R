################################################################################
# This file will run the genetic algorithm Differential Evolution to
# calibrate the transmission parameters of the model.
# This takes place by calibrating the contact matrices with a fixed
# probability of transmission. Our likelihood function for each fixed
# probability of transmissions is the joint likelihood of the assortative and
# random partnership contact matrices. This procedure will be repeated for
# various values of the probability of transmission.
################################################################################
### Load libraries
library(scales)
library(ggplot2)
library(ggsci)
library(Matrix)
library(parallel)
library(openxlsx)
library(DEoptim)
library(parallelly)
library(numDeriv)

### Load model functions
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

### Calibration functions
foi_transition_matrix_assort <- function(v_betas, w = m_w) {
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

foi_transition_matrix_rp <- function(v_betas, w = m_w) {
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

generate_nloglik <- function(param,
                             p_trans,
                             inf,
                             breaks = v_break,
                             w = m_w) {
  # param <- exp(param)  ## Use if calibrating log scale parameters
  if (length(param) == 32) {
    waifw <- foi_transition_matrix_assort(param)
  } else if (length(param) == 80) {
    waifw <- foi_transition_matrix_rp(param)
  }

  foi_hat <- (p_trans * waifw) %*% inf
  out     <- -1 * sum(dnorm(x = foi_hat,
                            mean = v_parameter$v_foi,
                            sd = v_parameter$v_foi_sd,
                            log = TRUE))
  return(out)
}

estimate_waifw  <- function(waifw_names = v_waifw_names,
                            v_breaks = waifw_breaks,
                            v_inf, iter, p_transmission,
                            v_upper, i) {
  require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)
  require(DEoptim)

  v_upper <- v_upper[i, ]

  waifw_ll   <- DEoptim(fn = generate_nloglik,
                        p_trans = p_transmission,
                        inf = v_inf,
                        lower = rep(-40, length(v_upper)),
                        upper = v_upper,
                        control = list(itermax = iter,
                                       parallelType = "parallel",
                                       parVar = list("v_break", "waifw_breaks",
                                                     "m_w",
                                                     "get_transmission_matrix",
                                                     "groups", "v_parameter",
                                                     "foi_transition_matrix_assort",
                                                     "foi_transition_matrix_rp")))

  waifw_llk   <- waifw_ll$optim$bestval
  v_waifw_hat <- as.vector(waifw_ll$optim$bestmem)

  # Trace of betas
  m_waifw    <- waifw_ll$member$bestmemit

  # Trace of nll
  v_nll      <- waifw_ll$member$bestvalit

  # Calculate Hessian to get covariance and SE
  m_hess <- numDeriv::hessian(generate_nloglik, v_waifw_hat,
                              p_trans = p_transmission,
                              inf = v_inf)
  if (MHadaptive::isPositiveDefinite(m_hess) == FALSE) {
    print("Hessian is NOT Positive Definite")
    m_hess <- as.matrix(Matrix::nearPD(m_hess)$mat)
  }
  m_cov  <- solve(m_hess)
  v_se   <- sqrt(diag(m_cov))

  # Generate big WAIFW matrices
  if (length(v_waifw_hat) == 32) {
    m_waifw_hat <- foi_transition_matrix_assort(v_waifw_hat)
  } else if (length(v_waifw_hat) == 80) {
    m_waifw_hat <- foi_transition_matrix_rp(v_waifw_hat)
  }

  # Percent text
  prob_text <- paste0(p_transmission * 100, "%")

  # Create latex code for coefficient table
  waifw_tex <- createTexreg(coef.names = waifw_names,
                            coef = v_waifw_hat,
                            se = v_se,
                            gof = waifw_llk,
                            gof.names = as.character(prob_text))

  return(list(waifw_ll, waifw_llk, v_waifw_hat, m_waifw,
              m_waifw_hat, waifw_tex, v_se, m_cov, v_nll))
}

calibrate_contacts <- function(v_waifw_breaks = waifw_breaks,
                               v_race_names,
                               v_p_trans,
                               v_ub,
                               v_inf,
                               n_iter = 1500) {

  n_age_groups <- length(v_waifw_breaks) - 1
  n_p_trans <- length(v_p_trans)
  n_betas <- ncol(v_ub)

  ## Beta and WAIFW latex names
  v_race_beta_names <- as.vector(outer(1:n_age_groups, v_race_names,
                                       paste, sep = "_"))
  v_waifw_names  <- paste0(paste0("$", paste("\\beta", v_race_beta_names,
                                             sep = "_")), "$")
  v_p_names <- paste0(paste0("$", paste("p", 1:n_p_trans, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_p_trans)
  v_beta_llk     <- numeric(n_p_trans)
  m_beta_hat     <- matrix(0, ncol = n_betas, nrow = n_p_trans,
                           dimnames = list(paste0("p", seq(1:n_p_trans)),
                                           paste0("b", seq(1:n_betas))))
  a_betas <- array(0, dim = c(n_betas, n_iter, n_p_trans),
                   dimnames = list(paste0("b", seq(1:n_betas)),
                                   paste0("r", seq(1:n_iter)),
                                   paste0("p", seq(1:n_p_trans))))
  m_beta_hat_se  <- matrix(0, ncol = n_betas, nrow = n_p_trans,
                           dimnames = list(paste0("p", seq(1:n_p_trans)),
                                           paste0("b", seq(1:n_betas))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_p_names

  a_beta_hat_cor <- array(0, dim = c(n_betas, n_betas, n_p_trans),
                          dimnames = list(paste0("b", seq(1:n_betas)),
                                          paste0("b", seq(1:n_betas)),
                                          paste0("p", seq(1:n_p_trans))))

  v_Beta_hat     <- vector("list", n_p_trans)
  v_beta_tex     <- vector("list", n_p_trans)
  names(v_beta_tex) <- v_p_names

  m_nll     <- matrix(0, nrow = n_iter, ncol = n_p_trans,
                      dimnames = list(paste0("r", seq(1:n_iter)),
                                      paste0("p", seq(1:n_p_trans))))

  for (j in seq.int(n_p_trans)) {
    results_beta <-  estimate_waifw(v_inf = v_inf,
                                    iter = n_iter,
                                    p_transmission = v_p_trans[j],
                                    waifw_names = v_waifw_names,
                                    v_upper = v_ub,
                                    i = j)
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

### Variables needed for calibration
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
w3 <- matrix(c(1, 1, 2, 3, 3, 5, 5, 5,
               0, 1, 2, 3, 3, 5, 5, 5,
               0, 0, 2, 4, 4, 6, 7, 8,
               0, 0, 0, 4, 4, 6, 7, 8,
               0, 0, 0, 0, 4, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)

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
  prob = FALSE,
  p_trans = NULL
)

v_race_prop <- v_parameter$v_pop / sum(v_parameter$v_pop)
v_race_prop <- rep(v_race_prop, each = length(groups))
v_inf <- v_parameter$v_prev * v_race_prop * v_parameter$v_age_prop
v_foi <- v_parameter$v_foi

v_break <- waifw_breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

v_race_comb_names <- c("hisp_hisp", "white_white", "black_black", "other_other",
                       "hisp_white", "hisp_black", "hisp_other", "white_black",
                       "white_other", "black_other")


### Run calibration for multiple fixed probabilities of transmission
v_p <- seq(from = 0.005, to = 0.05, by = 0.005)

m_ub <- rbind(rep(1500, 32),
              rep(1000, 32),
              rep(800, 32),
              rep(600, 32),
              rep(600, 32),
              rep(400, 32),
              rep(400, 32),
              rep(400, 32),
              rep(200, 32),
              rep(200, 32))

# m_w <- w1
# calibrate_contacts_assort_w1 <- calibrate_contacts(
#   v_race_names = v_race_calibrate,
#   v_p_trans = v_p,
#   v_ub = m_ub,
#   v_inf = v_inf
# )

m_w <- w2
calibrate_contacts_assort_w2 <- calibrate_contacts(
  v_race_names = v_race_calibrate,
  v_p_trans = v_p,
  v_ub = m_ub,
  v_inf = v_inf
)

m_w <- w3
calibrate_contacts_assort_w3 <- calibrate_contacts(
  v_race_names = v_race_calibrate,
  v_p_trans = v_p,
  v_ub = m_ub,
  v_inf = v_inf,
  n_iter = 3500
)

m_ub <- rbind(rep(1500, 80),
              rep(1000, 80),
              rep(800, 80),
              rep(600, 80),
              rep(600, 80),
              rep(400, 80),
              rep(400, 80),
              rep(400, 80),
              rep(200, 80),
              rep(200, 80))

# m_w <- w1
# calibrate_contacts_rp_w1 <- calibrate_contacts(
#   v_race_names = v_race_comb_names,
#   v_p_trans = v_p,
#   v_ub = m_ub,
#   v_inf = v_inf,
#   n_iter = 3500,
# )

m_w <- w2
calibrate_contacts_rp_w2 <- calibrate_contacts(
  v_race_names = v_race_comb_names,
  v_p_trans = v_p,
  v_ub = m_ub,
  v_inf = v_inf,
  n_iter = 3500,
)

# m_ub <- rbind(rep(7, 80),
#               rep(6, 80),
#               rep(6, 80),
#               rep(5, 80),
#               rep(5, 80),
#               rep(4, 80),
#               rep(4, 80),
#               rep(4, 80),
#               rep(4, 80),
#               rep(4, 80))
# v_p <- c(0.005)
m_w <- w3
calibrate_contacts_rp_w3 <- calibrate_contacts(
  v_race_names = v_race_comb_names,
  v_p_trans = v_p,
  v_ub = m_ub,
  v_inf = v_inf,
  n_iter = 5000,
)

???### Save calibration results
save_results <- function(l_results, type = "assort", w = c("W1", "w2", "w3")) {
  require(data.table)

  if (type == "assort") {
    file_path <- "results/waifw/assort_"
  } else {
    file_path <- "results/waifw/rp_"
  }

  # Latex
  all_latex <- texreg(l_results$v_beta_tex, digits = 4, stars = numeric(0),
                      booktabs = TRUE, single.row = TRUE,
                      custom.model.names = as.character(v_p))
  write.table(all_latex, file = paste0(file_path, "contacts_ga_", w, ".txt"))

  # Table with beta hat
  fwrite(l_results$m_beta_hat,
         file = paste0(file_path, "contacts_hat_ga_", w, ".csv"))

  # Table with SEs
  fwrite(l_results$m_beta_hat_se,
         file = paste0(file_path, "contacts_se_ga_", w, ".csv"))

  # Covariance matrices, full WAIFWs, contacts trace
  for (i in 1:length(v_p)) {
    fwrite(as.data.frame(l_results$a_beta_hat_cor[, , i]),
           file = paste0(file_path, "contacts_cov_", i, "_ga_", w, ".csv"))
    fwrite(as.data.frame(l_results$a_betas[, , i]),
           file = paste0(file_path, "contacts_trace_", i, "_ga_", w, ".csv"))
    fwrite(as.data.frame(l_results$v_Beta_hat[[i]]),
           file = paste0(file_path, "Contacts_hat_", i, "_ga_", w, ".csv"))
  }

  # Negative log-likelihood trace
  fwrite(l_results$m_nll,
         file = paste0(file_path, "m_nll_ga_", w, ".csv"))

}

save_results(calibrate_contacts_assort_w2, w = "w2")
save_results(calibrate_contacts_rp_w2, type = "rp", w = "w2")

save_results(calibrate_contacts_assort_w3, w = "w3")
save_results(calibrate_contacts_rp_w3, type = "rp", w = "w3")

##### Try IMIS #####
library(IMIS)
library(lhs)
library(truncnorm)
# sample.prior -- draws samples, and we have already created this
sample.prior <- function(n.samp, n.params, v_mean = v_mean_i, v_se = v_se_i) {
  # n.samp: the number of samples desired
  draws0 <- randomLHS(n = n.samp, k = n.params)
  draws <- matrix(nrow = n.samp, ncol = n.params)

  lb <- rep(0, n.params)
  for (i in seq_len(ncol(draws))) {
  draws[, i] <- qtruncnorm(draws0[, i], a = lb[i],
                           mean = v_mean[i], sd = v_se[i])
  }

  return(draws)
}

# prior -- evaluates prior density of a parameter set or sets
l_prior <- function(params, v_mean = v_mean_i, v_se = v_se_i) {
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params)
  }
  v_ll <- rep(0, nrow(params))
  for (i in 1:ncol(params)) {
    v_ll <- v_ll + log(dtruncnorm(params[, i], a = 0, mean = v_mean[i],
                              sd = v_se[i]))
  }

  return(v_ll)
}

prior <- function(par_vector) {
  v_prior <- exp(l_prior(par_vector))
  v_prior <- replace(v_prior, v_prior == 0, 1.1e-100)
  return(v_prior)
}


# likelihood -- evaluates likelihood of a parameter set or sets
l_likelihood <- function(params,
                         p_trans,
                         inf,
                         breaks = v_break,
                         w = m_w) {
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params)
  }

  llik <- rep(0, nrow(params))

  for (j in 1:nrow(params)) {
    jj <- tryCatch({
      if (ncol(params) == 32) {
        waifw <- foi_transition_matrix_assort(params[j, ])
      } else if (ncol(params) == 80) {
        waifw <- foi_transition_matrix_rp(params[j, ])
      }
      foi_hat         <- (p_trans * waifw) %*% inf
      foi_hat_med_age <- foi_hat[v_med_ages]
      llik[j]  <- sum(dnorm(x = foi_hat_med_age,
                            mean = v_parameter$v_foi[v_med_ages],
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

# Run IMIS
set.seed(12345)
v_p <- seq(from = 0.005, to = 0.05, by = 0.005)
med_age <- c(2, 10, 20, 35, 50, 60, 67, 75)
v_med_ages <- c(med_age, med_age + 80, med_age + 160, med_age + 240)


# W1
# m_w <- w1
# m_assort    <- calibrate_contacts_assort_w1$m_beta_hat
# m_assort_se <- calibrate_contacts_assort_w1$m_beta_hat_se

# m_rp    <- calibrate_contacts_rp_w1$m_beta_hat
# m_rp_se <- calibrate_contacts_rp_w1$m_beta_hat_se

# l_imis_assort_w1 <- list()
# l_imis_rp_w1     <- list()
# for (i in 1:length(v_p)) {
#   v_ub_i <- as.numeric(m_assort[i, ]) * 2
#   p_i  <- v_p[i]
#   v_assort <- as.numeric(m_assort[i, ])
#   v_assort_se <- as.numeric(m_assort_se[i, ]) * 2
#   # Set default values for both WAIFW structures
#   formals(l_likelihood) <- alist(params = , p_trans = p_i, inf = v_inf,
#                                  breaks = v_break, w = m_w)

#   # Assortative
#   formals(sample.prior) <- alist(n.samp = , n.params = 32, v_mean = v_assort,
#                                  v_se = v_assort_se)
#   formals(l_prior) <- alist(params = , v_mean = v_assort, v_se = v_assort_se)
#   imis_assort <- IMIS(B = 1000, B.re = 1e3, number_k = 400, D = 1)
#   l_imis_assort_w1[[i]] <- imis_assort

#   # Random Partnership
#   v_ub_i <- as.numeric(m_rp[i, ]) * 3
#   v_rp <- as.numeric(m_rp[i, ])
#   v_rp_se <- v_rp * 0.5
#   formals(sample.prior) <- alist(n.samp = , n.params = 80, v_mean = v_rp,
#                                  v_se = v_rp_se)
#   formals(l_prior) <- alist(params = , v_mean = v_rp, v_se = v_rp_se)

#   imis_rp     <- IMIS(B = 1000, B.re = 1e3, number_k = 300, D = 3)
#   l_imis_rp_w1[[i]] <- imis_rp
#   print(paste0("W1_", i))
# }

# W2
m_w <- w2
m_assort    <- calibrate_contacts_assort_w2$m_beta_hat
m_assort_se <- calibrate_contacts_assort_w2$m_beta_hat_se

m_rp    <- calibrate_contacts_rp_w2$m_beta_hat
m_rp_se <- calibrate_contacts_rp_w2$m_beta_hat_se

l_imis_assort_w2 <- list()
l_imis_rp_w2     <- list()
for (i in 1:length(v_p)) {
  v_ub_i <- as.numeric(m_assort[i, ]) * 2
  p_i  <- v_p[i]
  v_assort <- as.numeric(m_assort[i, ])
  v_assort_se <- as.numeric(m_assort_se[i, ]) * 2
  # Set default values for both WAIFW structures
  formals(l_likelihood) <- alist(params = , p_trans = p_i, inf = v_inf,
                                 breaks = v_break, w = m_w)

  # Assortative
  formals(sample.prior) <- alist(n.samp = , n.params = 32, v_mean = v_assort,
                                 v_se = v_assort_se)
  formals(l_prior) <- alist(params = , v_mean = v_assort, v_se = v_assort_se)
  imis_assort <- IMIS(B = 1000, B.re = 1e3, number_k = 400, D = 1)
  l_imis_assort_w2[[i]] <- imis_assort

  # Random Partnership
  v_ub_i <- as.numeric(m_rp[i, ]) * 3
  v_rp <- as.numeric(m_rp[i, ])
  v_rp_se <- as.numeric(m_rp_se[i, ]) * 0.2
  formals(sample.prior) <- alist(n.samp = , n.params = 80, v_mean = v_rp,
                                 v_se = v_rp_se)
  formals(l_prior) <- alist(params = , v_mean = v_rp, v_se = v_rp_se)

  imis_rp     <- IMIS(B = 1000, B.re = 1e3, number_k = 300, D = 3)
  l_imis_rp_w2[[i]] <- imis_rp
  print(paste0("W2_", i))
}


# W3
m_w <- w3
m_assort    <- calibrate_contacts_assort_w3$m_beta_hat
m_assort_se <- calibrate_contacts_assort_w3$m_beta_hat_se

m_rp    <- calibrate_contacts_rp_w3$m_beta_hat
m_rp_se <- calibrate_contacts_rp_w3$m_beta_hat_se

l_imis_assort_w3 <- list()
l_imis_rp_w3     <- list()
for (i in 1:length(v_p)) {
  # v_ub_i <- as.numeric(m_assort[i, ]) * 2
  # p_i  <- v_p[i]
  # v_assort <- as.numeric(m_assort[i, ])
  # v_assort_se <- as.numeric(m_assort_se[i, ]) * 0.1
  # # Set default values for both WAIFW structures
  # formals(l_likelihood) <- alist(params = , p_trans = p_i, inf = v_inf,
  #                                breaks = v_break, w = m_w)
  #
  # # Assortative
  # formals(sample.prior) <- alist(n.samp = , n.params = 32, v_mean = v_assort,
  #                                v_se = v_assort_se)
  # formals(l_prior) <- alist(params = , v_mean = v_assort, v_se = v_assort_se)
  # imis_assort <- IMIS(B = 1000, B.re = 1e3, number_k = 400, D = 1)
  # l_imis_assort_w3[[i]] <- imis_assort

  # Random Partnership
  v_ub_i <- as.numeric(m_rp[i, ]) * 3
  v_rp <- as.numeric(m_rp[i, ])
  # v_rp_se <-as.numeric(m_rp_se[i, ]) * 0.005
  v_rp_se <- as.numeric(m_rp[i, ]) * (1/5)
  formals(sample.prior) <- alist(n.samp = , n.params = 80, v_mean = v_rp,
                                 v_se = v_rp_se)
  formals(l_prior) <- alist(params = , v_mean = v_rp, v_se = v_rp_se)

  imis_rp     <- IMIS(B = 1000, B.re = 1e3, number_k = 300, D = 3)
  l_imis_rp_w3[[i]] <- imis_rp
  print(paste0("W3_", i))
}

### Save
for (i in 1:length(l_imis_assort_w3)) {
  fwrite(l_imis_assort_w3[[i]]$stat, paste0("results/waifw/imis_stat_", i,
                                         "_assort_w3.csv"))
  fwrite(l_imis_assort_w3[[i]]$resample, paste0("results/waifw/imis_resample_", i,
                                             "_assort_w3.csv"))
  fwrite(l_imis_assort_w3[[i]]$center, paste0("results/waifw/imis_center_", i,
                                           "_assort_w3.csv"))
}

for (i in 1:length(l_imis_rp_w3)) {
  fwrite(l_imis_rp_w2[[i]]$stat, paste0("results/waifw/imis_stat_", i,
                                         "_rp_w2.csv"))
  fwrite(l_imis_rp_w2[[i]]$resample, paste0("results/waifw/imis_resample_", i,
                                             "_rp_w2.csv"))
  fwrite(l_imis_rp_w2[[i]]$center, paste0("results/waifw/imis_center_", i,
                                           "_rp_w2.csv"))
}

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

foi_transition_matrix_assort <- function(v_betas, w = m_w) {
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

foi_transition_matrix_rp <- function(v_betas, w = m_w) {
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

# W1
## Assortative
# m_w <- w1
# l_m_resamp_assort <- list(l_imis_assort_w1[[1]]$resample,
#                           l_imis_assort_w1[[2]]$resample,
#                           l_imis_assort_w1[[3]]$resample,
#                           l_imis_assort_w1[[4]]$resample,
#                           l_imis_assort_w1[[5]]$resample,
#                           l_imis_assort_w1[[6]]$resample,
#                           l_imis_assort_w1[[7]]$resample,
#                           l_imis_assort_w1[[8]]$resample,
#                           l_imis_assort_w1[[9]]$resample,
#                           l_imis_assort_w1[[10]]$resample)

# ## Random Partnership
# l_m_resamp_rp <- list(l_imis_rp_w1[[1]]$resample, l_imis_rp_w1[[2]]$resample,
#                       l_imis_rp_w1[[3]]$resample, l_imis_rp_w1[[4]]$resample,
#                       l_imis_rp_w1[[5]]$resample, l_imis_rp_w1[[6]]$resample,
#                       l_imis_rp_w1[[7]]$resample, l_imis_rp_w1[[8]]$resample,
#                       l_imis_rp_w1[[9]]$resample, l_imis_rp_w1[[10]]$resample)

# seed_num <- 135135
# n_b  <- 1000
# n_b0 <- n_b * 10
# l_imis_alpha_w1 <- list()
# for (i in 1:length(v_p)) {
#   p_i <- v_p[i]
#   imis_assort <- l_m_resamp_assort[[i]]
#   imis_rp     <- l_m_resamp_rp[[i]]
#   set.seed(seed_num)
#   m_samp_assort <- imis_assort[sample(nrow(imis_assort),
#                                       size = n_b0,
#                                       replace = TRUE), ]
#   m_samp_rp <- imis_rp[sample(nrow(imis_rp),
#                               size = n_b0,
#                               replace = TRUE), ]
#   set.seed(seed_num + i)
#   imis_res <- IMIS(B = n_b, B.re = 1e3, number_k = 500, D = 0)
#   l_imis_alpha_w1[[i]] <- imis_res
#   print(paste0("W1_", i))
# }

# W2
## Assortative
m_w <- w2
l_m_resamp_assort <- list(l_imis_assort_w2[[1]]$resample,
                          l_imis_assort_w2[[2]]$resample,
                          l_imis_assort_w2[[3]]$resample,
                          l_imis_assort_w2[[4]]$resample,
                          l_imis_assort_w2[[5]]$resample,
                          l_imis_assort_w2[[6]]$resample,
                          l_imis_assort_w2[[7]]$resample,
                          l_imis_assort_w2[[8]]$resample,
                          l_imis_assort_w2[[9]]$resample,
                          l_imis_assort_w2[[10]]$resample)

## Random Partnership
l_m_resamp_rp <- list(l_imis_rp_w2[[1]]$resample, l_imis_rp_w2[[2]]$resample,
                      l_imis_rp_w2[[3]]$resample, l_imis_rp_w2[[4]]$resample,
                      l_imis_rp_w2[[5]]$resample, l_imis_rp_w2[[6]]$resample,
                      l_imis_rp_w2[[7]]$resample, l_imis_rp_w2[[8]]$resample,
                      l_imis_rp_w2[[9]]$resample, l_imis_rp_w2[[10]]$resample)

seed_num <- 458456
n_b  <- 1000
n_b0 <- n_b * 10
l_imis_alpha_w2 <- list()
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
  l_imis_alpha_w2[[i]] <- imis_res
  print(paste0("W2_", i))
}

# W3
## Assortative
m_w <- w3
l_m_resamp_assort <- list(l_imis_assort_w3[[1]]$resample,
                          l_imis_assort_w3[[2]]$resample,
                          l_imis_assort_w3[[3]]$resample,
                          l_imis_assort_w3[[4]]$resample,
                          l_imis_assort_w3[[5]]$resample,
                          l_imis_assort_w3[[6]]$resample,
                          l_imis_assort_w3[[7]]$resample,
                          l_imis_assort_w3[[8]]$resample,
                          l_imis_assort_w3[[9]]$resample,
                          l_imis_assort_w3[[10]]$resample)

## Random Partnership
l_m_resamp_rp <- list(l_imis_rp_w3[[1]]$resample, l_imis_rp_w3[[2]]$resample,
                      l_imis_rp_w3[[3]]$resample, l_imis_rp_w3[[4]]$resample,
                      l_imis_rp_w3[[5]]$resample, l_imis_rp_w3[[6]]$resample,
                      l_imis_rp_w3[[7]]$resample, l_imis_rp_w3[[8]]$resample,
                      l_imis_rp_w3[[9]]$resample, l_imis_rp_w3[[10]]$resample)

seed_num <- 456789
n_b  <- 1000
n_b0 <- n_b * 10
l_imis_alpha_w3 <- list()
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
  l_imis_alpha_w3[[i]] <- imis_res
  print(paste0("W3_", i))
}

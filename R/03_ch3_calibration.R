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

generate_nloglik <- function(param,
                             p_trans,
                             inf,
                             breaks = v_break,
                             w = w2) {
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
                        lower = rep(0, length(v_upper)),
                        upper = v_upper,
                        control = list(itermax = iter,
                                       parallelType = "parallel",
                                       parVar = list("v_break", "waifw_breaks",
                                                     "w2",
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
w2 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
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
v_inf_rp <- v_parameter$v_prev * v_race_prop * v_parameter$v_age_prop
v_inf_assort <- v_parameter$v_prev * v_parameter$v_age_prop
v_foi <- v_parameter$v_foi

v_break <- waifw_breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

v_race_comb_names <- c("hisp_hisp", "white_white", "black_black", "other_other",
                       "hisp_white", "hisp_black", "hisp_other", "white_black",
                       "white_other", "black_other")

##### Run calibration for a single fixed probability of transmission #####
# v_p         <- c(0.05)
# v_ub_assort <- rep(100, 32)
# v_ub_rp     <- rep(0.8, 80)
# 
# calibrate_contacts_assort <- calibrate_contacts(
#   v_race_names = v_race_calibrate,
#   v_p_trans = v_p,
#   v_ub = v_ub_assort,
#   v_inf = v_inf_assort
# )
# 
# calibrate_contacts_rp2 <- calibrate_contacts(
#   v_race_names = v_race_comb_names,
#   v_p_trans = v_p,
#   v_ub = v_ub_rp,
#   v_inf = v_inf_rp,
#   n_iter = 3500,
# )
# 
# ### Run calibration for multiple fixed probabilities of transmission
# v_p <- seq(from = 0.005, to = 0.05, by = 0.005)
# 
# m_ub <- rbind(rep(1500, 32),
#               rep(1000, 32),
#               rep(800, 32),
#               rep(600, 32),
#               rep(600, 32),
#               rep(400, 32),
#               rep(400, 32),
#               rep(400, 32),
#               rep(200, 32),
#               rep(200, 32))
# 
# calibrate_contacts_assort <- calibrate_contacts(
#   v_race_names = v_race_calibrate,
#   v_p_trans = v_p,
#   v_ub = m_ub,
#   v_inf = v_inf_rp
# )
# 
# m_ub <- rbind(rep(1500, 80),
#               rep(1000, 80),
#               rep(800, 80),
#               rep(600, 80),
#               rep(600, 80),
#               rep(400, 80),
#               rep(400, 80),
#               rep(400, 80),
#               rep(200, 80),
#               rep(200, 80))
# 
# calibrate_contacts_rp <- calibrate_contacts(
#   v_race_names = v_race_comb_names,
#   v_p_trans = v_p,
#   v_ub = m_ub,
#   v_inf = v_inf_rp,
#   n_iter = 3500,
# )

# ### Save calibration results
# save_results <- function(l_results, type = "assort") {
#   require(data.table)
# 
#   if (type == "assort") {
#     file_path <- "results/waifw/assort_"
#   } else {
#     file_path <- "results/waifw/rp_"
#   }
# 
#   # Latex
#   all_latex <- texreg(l_results$v_beta_tex, digits = 4, stars = numeric(0),
#                       booktabs = TRUE, single.row = TRUE,
#                       custom.model.names = as.character(v_p))
#   write.table(all_latex, file = paste0(file_path, "contacts_ga.txt"))
# 
#   # Table with beta hat
#   fwrite(l_results$m_beta_hat,
#          file = paste0(file_path, "contacts_hat_ga.csv"))
# 
#   # Table with SEs
#   fwrite(l_results$m_beta_hat_se,
#          file = paste0(file_path, "contacts_se_ga.csv"))
# 
#   # Covariance matrices, full WAIFWs, contacts trace
#   for (i in 1:length(v_p)) {
#     fwrite(as.data.frame(l_results$a_beta_hat_cor[, , i]),
#            file = paste0(file_path, "contacts_cov_", i, "_ga.csv"))
#     fwrite(as.data.frame(l_results$a_betas[, , i]),
#            file = paste0(file_path, "contacts_trace_", i, "_ga.csv"))
#     fwrite(as.data.frame(l_results$v_Beta_hat[[i]]),
#            file = paste0(file_path, "Contacts_hat_", i, "_ga.csv"))
#   }
# 
#   # Negative log-likelihood trace
#   fwrite(l_results$m_nll,
#          file = paste0(file_path, "m_nll_ga.csv"))
# 
# }
# 
# save_results(calibrate_contacts_assort)
# save_results(calibrate_contacts_rp, type = "rp")

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
l_prior <- function(params, v_ub = v_ub_i) {
  if(is.null(dim(params))) { # If vector, change to matrix
    params <- t(params)
  }
  v_ll <- rep(0, nrow(params))
  for (i in 1:ncol(params)) {
    v_ll <- v_ll + dunif(params[, i], min = 0, max = v_ub[i], log = TRUE)
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
                         w = w2) {
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

# l_likelihood_par <- function(params,
#                              ...) {
#   if (is.null(dim(params))) { # If vector, change to matrix
#     params <- t(params)
#   }
#
#   n_samp <- nrow(params)
#   no_cores <- parallel::detectCores() - 1
#
#   # Initialize cluster object
#   cl <- parallel::makeCluster(no_cores)
#   doParallel::registerDoParallel(cl)
#   opts <- list(attachExportEnv = TRUE)
#   v_llk <- foreach::foreach(i = 1:n_samp, .combine = c,
#                             .export = ls(globalenv()),
#                             .packages = c(),
#                             .options.snow = opts) %dopar% {
#                               l_likelihood(params[i, ])
#                             }
#   parallel::stopCluster(cl)
#   return(v_llk)
# }

likelihood <- function(params) {
  exp(l_likelihood(params))
}

# Run IMIS
set.seed(12345)
v_p <- seq(from = 0.005, to = 0.05, by = 0.005)
m_assort <- read.csv("results/waifw/assort_contacts_hat_ga.csv")
m_rp <- read.csv("results/waifw/rp_contacts_hat_ga.csv")

m_assort_se <- read.csv("results/waifw/assort_contacts_se_ga.csv")
m_rp_se <- read.csv("results/waifw/rp_contacts_se_ga.csv")

# l_imis_assort <- list()
# l_imis_rp     <- list()
# for (i in 1:length(v_p)) {
#   v_ub_i <- as.numeric(m_assort[i, ]) * 2
#   p_i  <- v_p[i]
#   v_assort <- as.numeric(m_assort[i, ])
#   v_assort_se <- as.numeric(m_assort_se[i, ])*1.5
#   # m_assort_cov <- as.matrix(read.csv("results/waifw/assort_contacts_cov_1_ga.csv"))
#   # Set default values for both WAIFW structures
#   formals(l_prior)      <- alist(params = , v_ub = v_ub_i)
#   formals(l_likelihood) <- alist(params = , p_trans = p_i, inf = v_inf_rp,
#                                  breaks = v_break, w = w2)
#
#   # Assortative
#   # # formals(sample.prior) <- alist(n.samp = , n.params = 32, v_ub = v_ub_i)
#   # formals(sample.prior) <- alist(n.samp = , n.params = 32, v_mean = v_assort,
#   #                                v_se = v_assort_se)
#   # # formals(sample.prior) <- alist(n.samp = , n.params = 32, v_mean = v_assort,
#   # #                                v_cov = m_assort_cov)
#   # imis_assort <- IMIS(B = 1000, B.re = 1e3, number_k = 200, D = 3)
#   # l_imis_assort[[i]] <- imis_assort
#
#   # Random Partnership
#   v_ub_i <- as.numeric(m_rp[i, ]) * 3
#   v_rp <- as.numeric(m_rp[i, ])
#   # v_rp_se <- as.numeric(m_rp_se[i, ]) * 0.005
#   v_rp_se <- v_rp * 0.3
#   formals(l_prior)      <- alist(params = , v_ub = v_ub_i)
#   formals(sample.prior) <- alist(n.samp = , n.params = 80, v_mean = v_rp,
#                                 v_se = v_rp_se)
#
#   imis_rp     <- IMIS(B = 1000, B.re = 1e3, number_k = 300, D = 3)
#   l_imis_rp[[i]] <- imis_rp
#   print(i)
# }

### Run with fewer calibration targets
med_age <- c(2, 10, 20, 35, 50, 60, 67, 75)
v_med_ages <- c(med_age, med_age + 80, med_age + 160, med_age + 240)

l_imis_assort <- list()
l_imis_rp     <- list()
for (i in 1:length(v_p)) {
  v_ub_i <- as.numeric(m_assort[i, ]) * 2
  p_i  <- v_p[i]
  v_assort <- as.numeric(m_assort[i, ])
  v_assort_se <- as.numeric(m_assort_se[i, ]) * 2
  # Set default values for both WAIFW structures
  formals(l_likelihood) <- alist(params = , p_trans = p_i, inf = v_inf_rp,
                                 breaks = v_break, w = w2)

  # Assortative
  formals(sample.prior) <- alist(n.samp = , n.params = 32, v_mean = v_assort,
                                 v_se = v_assort_se)
  imis_assort <- IMIS(B = 1000, B.re = 1e3, number_k = 400, D = 1)
  l_imis_assort[[i]] <- imis_assort

  # Random Partnership
  # v_ub_i <- as.numeric(m_rp[i, ]) * 3
  # v_rp <- as.numeric(m_rp[i, ])
  # v_rp_se <- v_rp * 0.5
  # formals(sample.prior) <- alist(n.samp = , n.params = 80, v_mean = v_rp,
  #                                v_se = v_rp_se)
  #
  # imis_rp     <- IMIS(B = 1000, B.re = 1e3, number_k = 300, D = 3)
  # l_imis_rp[[i]] <- imis_rp
  print(i)
}



# Assortative results
for(i in 1:length(v_p)) {
  p_i  <- v_p[i]
  # Unique parameter sets
  print(length(unique(l_imis_assort[[i]]$resample[, 1])))

  # Effective sample size
  print(sum(table(l_imis_assort[[i]]$resample[, 1]))^2 / sum(table(l_imis_assort[[i]]$resample[, 1])^2))

  # Max weight
  print(max(table(l_imis_assort[[i]]$resample[, 1])) / sum(table(l_imis_assort[[i]]$resample[, 1])))

  # Max log-likelihood
  print(max(test_l_likelihood(l_imis_assort[[i]]$resample, p_trans = p_i,
                              inf = v_inf_rp,
                              breaks = v_break, w = w2)))
}

# Random Partnership results
for(i in 1:length(v_p)) {
  p_i  <- v_p[i]
  # Unique parameter sets
  print(length(unique(l_imis_rp[[i]]$resample[, 1])))

  # Effective sample size
  print(sum(table(l_imis_rp[[i]]$resample[, 1]))^2 / sum(table(l_imis_rp[[i]]$resample[, 1])^2))

  # Max weight
  print(max(table(l_imis_rp[[i]]$resample[, 1])) / sum(table(l_imis_rp[[i]]$resample[, 1])))

  # Max log-likelihood
  print(max(test_l_likelihood(l_imis_rp[[i]]$resample, p_trans = p_i,
                              inf = v_inf_rp,
                              breaks = v_break, w = w2)))
}

# Save IMIS results
test_l_likelihood <- function(params,
                         p_trans,
                         inf,
                         breaks = v_break,
                         w = w2) {
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
      foi_hat <- (p_trans * waifw) %*% inf
      llik[j]  <- sum(dnorm(x = foi_hat,
                            mean = v_parameter$v_foi,
                            sd = v_parameter$v_foi_sd,
                            log = TRUE))
    }, error = function(e) NA)
    if (is.na(jj)) (llik[j] <- -Inf)
  }

  return(llik)
}

i <- 1
max(test_l_likelihood(imis_assort$resample, p_trans = p_i, inf = v_inf_rp,
                      breaks = v_break, w = w2))
test_l_likelihood(imis_rp$resample, p_trans = p_i, inf = v_inf_rp,
                  breaks = v_break, w = w2)

### Save
for (i in 1:length(l_imis_assort)) {
  fwrite(l_imis_assort[[i]]$stat, paste0("results/waifw/imis_stat_", i,
                                         "_assort.csv"))
  fwrite(l_imis_assort[[i]]$resample, paste0("results/waifw/imis_resample_", i,
                                             "_assort.csv"))
  fwrite(l_imis_assort[[i]]$center, paste0("results/waifw/imis_center_", i,
                                           "_assort.csv"))
}

for (i in 1:length(l_imis_rp)) {
  fwrite(l_imis_rp[[i]]$stat, paste0("results/waifw/imis_stat_", i,
                                         "_rp.csv"))
  fwrite(l_imis_rp[[i]]$resample, paste0("results/waifw/imis_resample_", i,
                                             "_rp.csv"))
  fwrite(l_imis_rp[[i]]$center, paste0("results/waifw/imis_center_", i,
                                           "_rp.csv"))
}

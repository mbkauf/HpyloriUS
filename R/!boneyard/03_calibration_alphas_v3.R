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

# Load other functions
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

# Define alpha bounds
a_11_bounds <- c(0.072, 0.928)
a_12_bounds <- c(0.086, 0.476)
a_13_bounds <- c(0.086, 0.326)
a_14_bounds <- c(0.086, 0.334)
a_21_bounds <- c(0.388, 0.835)
a_22_bounds <- c(0.119, 0.881)
a_23_bounds <- c(0.388, 0.770)
a_24_bounds <- c(0.388, 0.776)
a_31_bounds <- c(0.017, 0.138)
a_32_bounds <- c(0.017, 0.162)
a_33_bounds <- c(0.026, 0.974)
a_34_bounds <- c(0.017, 0.090)
a_41_bounds <- c(0.003, 0.110)
a_42_bounds <- c(0.003, 0.120)
a_43_bounds <- c(0.003, 0.055)
a_44_bounds <- c(0.069, 0.931)

l_a_bounds <- list(a_11_bounds, a_12_bounds, a_13_bounds, a_14_bounds,
                   a_21_bounds, a_22_bounds, a_23_bounds, a_24_bounds,
                   a_31_bounds, a_32_bounds, a_33_bounds, a_34_bounds,
                   a_41_bounds, a_42_bounds, a_43_bounds, a_44_bounds)

nl_likelihood <- function(params, p_trans, waifw_assort, waifw_rp, inf) {
  v_alphas_within     <- params[c(1, 6, 11, 16)]
  v_alphas_within_off <- 1 - params[c(1, 6, 11, 16)]
  v_alphas_between    <- params[c(2:5, 7:10, 12:15)]
  m_within      <- alpha_matrix_within(v_alphas_within)
  m_within_off  <- alpha_matrix_within(v_alphas_within_off)
  m_between     <- alpha_matrix_between(v_alphas_between)

  waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
    (m_between * waifw_rp)
  foi_hat <- (p_trans * waifw_p) %*% inf

  nll  <- -1 * sum(dnorm(x = v_foi,
                         mean = foi_hat,
                         sd = v_parameter$v_foi_sd,
                         log = TRUE))

  return(nll)
}

l_prior <- function(params, v_bounds = l_a_bounds) {
  ll_a_11 <- dunif(params[1],  min = l_a_bounds[[1]][1],
                   max = l_a_bounds[[1]][2], log = TRUE)
  ll_a_12 <- dunif(params[2],  min = l_a_bounds[[2]][1],
                   max = l_a_bounds[[2]][2], log = TRUE)
  ll_a_13 <- dunif(params[3],  min = l_a_bounds[[3]][1],
                   max = l_a_bounds[[3]][2], log = TRUE)
  ll_a_14 <- dunif(params[4],  min = l_a_bounds[[4]][1],
                   max = l_a_bounds[[4]][2], log = TRUE)
  ll_a_21 <- dunif(params[5],  min = l_a_bounds[[5]][1],
                   max = l_a_bounds[[5]][2], log = TRUE)
  ll_a_22 <- dunif(params[6],  min = l_a_bounds[[6]][1],
                   max = l_a_bounds[[6]][2], log = TRUE)
  ll_a_23 <- dunif(params[7],  min = l_a_bounds[[7]][1],
                   max = l_a_bounds[[7]][2], log = TRUE)
  ll_a_24 <- dunif(params[8],  min = l_a_bounds[[8]][1],
                   max = l_a_bounds[[8]][2], log = TRUE)
  ll_a_31 <- dunif(params[9],  min = l_a_bounds[[9]][1],
                   max = l_a_bounds[[9]][2], log = TRUE)
  ll_a_32 <- dunif(params[10], min = l_a_bounds[[10]][1],
                   max = l_a_bounds[[10]][2], log = TRUE)
  ll_a_33 <- dunif(params[11], min = l_a_bounds[[11]][1],
                   max = l_a_bounds[[11]][2], log = TRUE)
  ll_a_34 <- dunif(params[12], min = l_a_bounds[[12]][1],
                   max = l_a_bounds[[12]][2], log  = TRUE)
  ll_a_41 <- dunif(params[13], min = l_a_bounds[[13]][1],
                   max = l_a_bounds[[13]][2], log = TRUE)
  ll_a_42 <- dunif(params[14], min = l_a_bounds[[14]][1],
                   max = l_a_bounds[[14]][2], log = TRUE)
  ll_a_43 <- dunif(params[15], min = l_a_bounds[[15]][1],
                   max = l_a_bounds[[15]][2], log = TRUE)
  ll_a_44 <- dunif(params[16], min = l_a_bounds[[16]][1],
                   max = l_a_bounds[[16]][2], log = TRUE)

  nll <- -1 * sum(ll_a_11, ll_a_12, ll_a_13, ll_a_14,
                  ll_a_21, ll_a_22, ll_a_23, ll_a_24,
                  ll_a_31, ll_a_32, ll_a_33, ll_a_34,
                  ll_a_41, ll_a_42, ll_a_43, ll_a_44)
  return(nll)
}

l_posterior <- function(params, p_trans, waifw_assort, waifw_rp, inf) {
  return(l_prior(params) +
           nl_likelihood(params, p_trans, waifw_assort, waifw_rp, inf))
}

estimate_alphas <- function(v_inf, iter, p_transmission, i,
                            waifw_assort, waifw_rp) {
  # require(devtools)
  require(MHadaptive)  # devtools::install_github("cjbayesian/MHadaptive")
  require(texreg)
  require(Matrix)
  require(DEoptim)

  calibrated_alphas <- DEoptim(
    fn = l_posterior,
    p_trans = p_transmission,
    waifw_assort = waifw_assort,
    waifw_rp = waifw_rp,
    inf = v_inf,
    lower = c(a_11_bounds[1], a_12_bounds[1], a_13_bounds[1], a_14_bounds[1],
              a_21_bounds[1], a_22_bounds[1], a_23_bounds[1], a_24_bounds[1],
              a_31_bounds[1], a_32_bounds[1], a_33_bounds[1], a_34_bounds[1],
              a_41_bounds[1], a_42_bounds[1], a_43_bounds[1], a_44_bounds[1]),
    upper = c(a_11_bounds[2], a_12_bounds[2], a_13_bounds[2], a_14_bounds[2],
              a_21_bounds[2], a_22_bounds[2], a_23_bounds[2], a_24_bounds[2],
              a_31_bounds[2], a_32_bounds[2], a_33_bounds[2], a_34_bounds[2],
              a_41_bounds[2], a_42_bounds[2], a_43_bounds[2], a_44_bounds[2]),
    control = list(itermax = iter,
                   parallelType = "parallel",
                   parVar = list("v_inf_rp", "l_prior", "l_posterior",
                                 "nl_likelihood", "groups", "v_parameter",
                                 "l_a_bounds", "alpha_matrix_within",
                                 "alpha_matrix_between", "v_foi"))
  )

  alpha_llk   <- calibrated_alphas$optim$bestval
  v_alpha_hat <- as.vector(calibrated_alphas$optim$bestmem)

  # Trace of betas
  m_alpha    <- calibrated_alphas$member$bestmemit

  # Trace of nll
  v_nll      <- calibrated_alphas$member$bestvalit

  alpha_hess <- numDeriv::hessian(nl_likelihood, calibrated_alphas$optim$bestmem,
                                  p_trans = p_transmission,
                                  waifw_assort = waifw_assort,
                                  waifw_rp = waifw_rp,
                                  inf = v_inf)

  if (MHadaptive::isPositiveDefinite(alpha_hess) == FALSE) {
    print("Hessian is NOT Positive Definite")
    alpha_hess <- as.matrix(Matrix::nearPD(alpha_hess)$mat)
  }

  alpha_cov  <- solve(alpha_hess)
  alpha_se   <- sqrt(diag(alpha_cov))

  # Percent text
  prob_text <- paste0(p_transmission * 100, "%")

  alpha_names <- paste0(paste0("$", paste("\\alpha",
                                          c("{11}", "{12}", "{13}", "{14}",
                                            "{21}", "{22}", "{23}", "{24}",
                                            "{31}", "{32}", "{33}", "{34}",
                                            "{41}", "{42}", "{43}", "{44}"),
                                          sep =  "_")), "$")
  alpha_tex <- createTexreg(coef.names = alpha_names,
                            coef = calibrated_alphas$optim$bestmem,
                            se = alpha_se,
                            gof = calibrated_alphas$optim$bestval,
                            gof.names = as.character(prob_text))

  return(list(calibrated_alphas, alpha_llk, v_alpha_hat, m_alpha,
              alpha_tex, alpha_se, alpha_cov, v_nll))
}

calibrate_alphas <- function(v_p_trans,
                             v_inf,
                             l_waifw_assort,
                             l_waifw_rp,
                             n_iter = 1500) {

  n_p_trans    <- length(v_p_trans)
  n_alphas     <- 16

  ## Beta and WAIFW latex names

  v_p_names <- paste0(paste0("$", paste("p", 1:n_p_trans, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_alphas_ll     <- vector("list", n_p_trans)
  v_alphas_llk    <- numeric(n_p_trans)
  m_alphas_hat     <- matrix(0, ncol = n_alphas, nrow = n_p_trans,
                             dimnames = list(paste0("p", seq(1:n_p_trans)),
                                             paste0("b", seq(1:n_alphas))))
  a_alphas <- array(0, dim = c(n_alphas, n_iter, n_p_trans),
                    dimnames = list(paste0("b", seq(1:n_alphas)),
                                    paste0("r", seq(1:n_iter)),
                                    paste0("p", seq(1:n_p_trans))))
  m_alphas_hat_se  <- matrix(0, ncol = n_alphas, nrow = n_p_trans,
                             dimnames = list(paste0("p", seq(1:n_p_trans)),
                                             paste0("b", seq(1:n_alphas))))
  rownames(m_alphas_hat) <- rownames(m_alphas_hat_se) <- v_p_names

  a_alphas_hat_cov <- array(0, dim = c(n_alphas, n_alphas, n_p_trans),
                            dimnames = list(paste0("b", seq(1:n_alphas)),
                                            paste0("b", seq(1:n_alphas)),
                                            paste0("p", seq(1:n_p_trans))))

  v_alphas_tex     <- vector("list", n_p_trans)
  names(v_alphas_tex) <- v_p_names

  m_nll     <- matrix(0, nrow = n_iter, ncol = n_p_trans,
                      dimnames = list(paste0("r", seq(1:n_iter)),
                                      paste0("p", seq(1:n_p_trans))))

  for (j in seq.int(n_p_trans)) {
    results_alpha <-  estimate_alphas(v_inf = v_inf,
                                      iter = n_iter,
                                      p_transmission = v_p_trans[j],
                                      waifw_assort = l_waifw_assort[[j]],
                                      waifw_rp = l_waifw_rp[[j]],
                                      i = j)
    v_alphas_ll[[j]]        <- results_alpha[[1]]
    v_alphas_llk[j]         <- results_alpha[[2]]
    m_alphas_hat[j, ]       <- results_alpha[[3]]
    a_alphas[, , j]         <- results_alpha[[4]]
    v_alphas_tex[[j]]       <- results_alpha[[5]]
    m_alphas_hat_se[j, ]    <- results_alpha[[6]]
    a_alphas_hat_cov[, , j] <- results_alpha[[7]]
    m_nll[, j]              <- results_alpha[[8]]
  }

  return(list(v_alphas_ll = v_alphas_ll,
              v_alphas_llk = v_alphas_llk,
              m_alphas_hat = m_alphas_hat,
              a_alphas = a_alphas,
              v_alphas_tex = v_alphas_tex,
              m_alphas_hat_se = m_alphas_hat_se,
              a_alphas_hat_cov = a_alphas_hat_cov,
              m_nll = m_nll))
}

# Run calibration of alpha parameters
# Number of age groups
groups <- c(1:80)
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

v_race_prop <- v_parameter$v_pop / sum(v_parameter$v_pop)
v_race_prop <- rep(v_race_prop, each = length(groups))
v_inf_rp <- v_parameter$v_prev * v_race_prop
v_inf_assort <- v_parameter$v_prev
v_foi <- v_parameter$v_foi

v_p <- seq(from = 0.005, to = 0.05, by = 0.005)

# Load contact matrices
v_m_contacts_assort <- c("assort_Contacts_hat_1_ga", "assort_Contacts_hat_2_ga",
                         "assort_Contacts_hat_3_ga", "assort_Contacts_hat_4_ga",
                         "assort_Contacts_hat_5_ga", "assort_Contacts_hat_6_ga",
                         "assort_Contacts_hat_7_ga", "assort_Contacts_hat_8_ga",
                         "assort_Contacts_hat_9_ga",
                         "assort_Contacts_hat_10_ga")
v_m_contacts_rp <- c("rp_Contacts_hat_1_ga", "rp_Contacts_hat_2_ga",
                     "rp_Contacts_hat_3_ga", "rp_Contacts_hat_4_ga",
                     "rp_Contacts_hat_5_ga", "rp_Contacts_hat_6_ga",
                     "rp_Contacts_hat_7_ga", "rp_Contacts_hat_8_ga",
                     "rp_Contacts_hat_9_ga", "rp_Contacts_hat_10_ga")


v_m_contacts_all <- c(v_m_contacts_assort, v_m_contacts_rp)

for (i in seq_along(v_m_contacts_all)) {
  assign(v_m_contacts_all[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_contacts_all[i], ".csv"))))
}

l_m_contacts_assort <- list(assort_Contacts_hat_1_ga, assort_Contacts_hat_2_ga,
                            assort_Contacts_hat_3_ga, assort_Contacts_hat_4_ga,
                            assort_Contacts_hat_5_ga, assort_Contacts_hat_6_ga,
                            assort_Contacts_hat_7_ga, assort_Contacts_hat_8_ga,
                            assort_Contacts_hat_9_ga,
                            assort_Contacts_hat_10_ga)
l_m_contacts_rp <- list(rp_Contacts_hat_1_ga, rp_Contacts_hat_2_ga,
                        rp_Contacts_hat_3_ga, rp_Contacts_hat_4_ga,
                        rp_Contacts_hat_5_ga, rp_Contacts_hat_6_ga,
                        rp_Contacts_hat_7_ga, rp_Contacts_hat_8_ga,
                        rp_Contacts_hat_9_ga, rp_Contacts_hat_10_ga)

results_alphas <- calibrate_alphas(v_p_trans = v_p,
                                   v_inf = v_inf_rp,
                                   l_waifw_assort = l_m_contacts_assort,
                                   l_waifw_rp = l_m_contacts_rp)

### Edits names to make results
alpha_names <- paste0(paste("a", c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4,4)),
                            sep = "_"), rep(seq(1:4), 4))
colnames(results_alphas$m_alphas_hat) <-
  colnames(results_alphas$m_alphas_hat_se) <- alpha_names

file_path <- "results/waifw/"

# Latex
l_tex_alphas <- list(results_alphas[["v_alphas_tex"]][["$p_1$"]],
                     results_alphas[["v_alphas_tex"]][["$p_2$"]],
                     results_alphas[["v_alphas_tex"]][["$p_3$"]],
                     results_alphas[["v_alphas_tex"]][["$p_4$"]],
                     results_alphas[["v_alphas_tex"]][["$p_5$"]],
                     results_alphas[["v_alphas_tex"]][["$p_6$"]],
                     results_alphas[["v_alphas_tex"]][["$p_7$"]],
                     results_alphas[["v_alphas_tex"]][["$p_8$"]],
                     results_alphas[["v_alphas_tex"]][["$p_9$"]],
                     results_alphas[["v_alphas_tex"]][["$p_10$"]])
all_latex <- texreg(l_tex_alphas, digits = 3, stars = numeric(0),
                    booktabs = TRUE, single.row = TRUE,
                    custom.model.names = as.character(v_p))
write.table(all_latex, file = paste0(file_path, "alphas_tex_ga.txt"))

# Table with alphas hat
fwrite(results_alphas$m_alphas_hat,
       file = paste0(file_path, "alphas_hat_ga.csv"))

# Table with SEs
fwrite(results_alphas$m_alphas_hat_se,
       file = paste0(file_path, "alphas_se_ga.csv"))

# Covariance matrices, full WAIFWs, contacts trace
for (i in 1:length(v_p)) {
  fwrite(as.data.frame(results_alphas$a_alphas_hat_cov[, , i]),
         file = paste0(file_path, "alphas_cov_", i, "_ga.csv"))
  fwrite(as.data.frame(results_alphas$a_alphas[, , i]),
         file = paste0(file_path, "alphas_trace_", i, "_ga.csv"))
}

# Negative log-likelihood trace
fwrite(results_alphas$m_nll,
       file = paste0(file_path, "m_nll_alphas_ga.csv"))

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
  draws0 <- randomLHS(n = n.samp, k = 16)
  draws <- matrix(nrow = n.samp, ncol = 16)

  draws <- data.frame(
    a_11 <- qunif(draws0[, 1],  min = l_a_bounds[[1]][1],
                  max = l_a_bounds[[1]][2], log = TRUE),
    a_12 <- qunif(draws0[, 2],  min = l_a_bounds[[2]][1],
                  max = l_a_bounds[[2]][2], log = TRUE),
    a_13 <- qunif(draws0[, 3],  min = l_a_bounds[[3]][1],
                  max = l_a_bounds[[3]][2], log = TRUE),
    a_14 <- qunif(draws0[, 4],  min = l_a_bounds[[4]][1],
                  max = l_a_bounds[[4]][2], log = TRUE),
    a_21 <- qunif(draws0[, 5],  min = l_a_bounds[[5]][1],
                  max = l_a_bounds[[5]][2], log = TRUE),
    a_22 <- qunif(draws0[, 6],  min = l_a_bounds[[6]][1],
                  max = l_a_bounds[[6]][2], log = TRUE),
    a_23 <- qunif(draws0[, 7],  min = l_a_bounds[[7]][1],
                  max = l_a_bounds[[7]][2], log = TRUE),
    a_24 <- qunif(draws0[, 8],  min = l_a_bounds[[8]][1],
                  max = l_a_bounds[[8]][2], log = TRUE),
    a_31 <- qunif(draws0[, 9],  min = l_a_bounds[[9]][1],
                  max = l_a_bounds[[9]][2], log = TRUE),
    a_32 <- qunif(draws0[, 10], min = l_a_bounds[[10]][1],
                  max = l_a_bounds[[10]][2], log = TRUE),
    a_33 <- qunif(draws0[, 11], min = l_a_bounds[[11]][1],
                  max = l_a_bounds[[11]][2], log = TRUE),
    a_34 <- qunif(draws0[, 12], min = l_a_bounds[[12]][1],
                  max = l_a_bounds[[12]][2], log  = TRUE),
    a_41 <- qunif(draws0[, 13], min = l_a_bounds[[13]][1],
                  max = l_a_bounds[[13]][2], log = TRUE),
    a_42 <- qunif(draws0[, 14], min = l_a_bounds[[14]][1],
                  max = l_a_bounds[[14]][2], log = TRUE),
    a_43 <- qunif(draws0[, 15], min = l_a_bounds[[15]][1],
                  max = l_a_bounds[[15]][2], log = TRUE),
    a_44 <- qunif(draws0[, 16], min = l_a_bounds[[16]][1],
                  max = l_a_bounds[[16]][2], log = TRUE)
  )
  colnames(draws) <- c("a_11", "a_12", "a_13", "a_14",
                       "a_21", "a_22", "a_23", "a_24",
                       "a_31", "a_32", "a_33", "a_34",
                       "a_41", "a_42", "a_43", "a_44")
  return(draws)
}

# prior -- evaluates prior density of a parameter set or sets
l_prior <- function(params, v_bounds = l_a_bounds) {
  ll_a_11 <- dunif(params[1],  min = l_a_bounds[[1]][1],
                   max = l_a_bounds[[1]][2], log = TRUE)
  ll_a_12 <- dunif(params[2],  min = l_a_bounds[[2]][1],
                   max = l_a_bounds[[2]][2], log = TRUE)
  ll_a_13 <- dunif(params[3],  min = l_a_bounds[[3]][1],
                   max = l_a_bounds[[3]][2], log = TRUE)
  ll_a_14 <- dunif(params[4],  min = l_a_bounds[[4]][1],
                   max = l_a_bounds[[4]][2], log = TRUE)
  ll_a_21 <- dunif(params[5],  min = l_a_bounds[[5]][1],
                   max = l_a_bounds[[5]][2], log = TRUE)
  ll_a_22 <- dunif(params[6],  min = l_a_bounds[[6]][1],
                   max = l_a_bounds[[6]][2], log = TRUE)
  ll_a_23 <- dunif(params[7],  min = l_a_bounds[[7]][1],
                   max = l_a_bounds[[7]][2], log = TRUE)
  ll_a_24 <- dunif(params[8],  min = l_a_bounds[[8]][1],
                   max = l_a_bounds[[8]][2], log = TRUE)
  ll_a_31 <- dunif(params[9],  min = l_a_bounds[[9]][1],
                   max = l_a_bounds[[9]][2], log = TRUE)
  ll_a_32 <- dunif(params[10], min = l_a_bounds[[10]][1],
                   max = l_a_bounds[[10]][2], log = TRUE)
  ll_a_33 <- dunif(params[11], min = l_a_bounds[[11]][1],
                   max = l_a_bounds[[11]][2], log = TRUE)
  ll_a_34 <- dunif(params[12], min = l_a_bounds[[12]][1],
                   max = l_a_bounds[[12]][2], log  = TRUE)
  ll_a_41 <- dunif(params[13], min = l_a_bounds[[13]][1],
                   max = l_a_bounds[[13]][2], log = TRUE)
  ll_a_42 <- dunif(params[14], min = l_a_bounds[[14]][1],
                   max = l_a_bounds[[14]][2], log = TRUE)
  ll_a_43 <- dunif(params[15], min = l_a_bounds[[15]][1],
                   max = l_a_bounds[[15]][2], log = TRUE)
  ll_a_44 <- dunif(params[16], min = l_a_bounds[[16]][1],
                   max = l_a_bounds[[16]][2], log = TRUE)

  ll <- sum(ll_a_11, ll_a_12, ll_a_13, ll_a_14,
            ll_a_21, ll_a_22, ll_a_23, ll_a_24,
            ll_a_31, ll_a_32, ll_a_33, ll_a_34,
            ll_a_41, ll_a_42, ll_a_43, ll_a_44)
  return(ll)
}

prior <- function(par_vector) {
  exp(l_prior(par_vector))
}


# likelihood -- evaluates likelihood of a parameter set or sets
l_likelihood <- function(params, p_trans, waifw_assort, waifw_rp, inf) {
  llik <- rep(0, nrow(params))
  for (j in seq_len(nrow(params))) {
    jj <- tryCatch({
      v_alphas_within     <- params[c(1, 6, 11, 16)]
      v_alphas_within_off <- 1 - params[c(1, 6, 11, 16)]
      v_alphas_between    <- params[c(2:5, 7:10, 12:15)]
      m_within      <- alpha_matrix_within(v_alphas_within)
      m_within_off  <- alpha_matrix_within(v_alphas_within_off)
      m_between     <- alpha_matrix_between(v_alphas_between)

      waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
        (m_between * waifw_rp)
      foi_hat <- (p_trans * waifw_p) %*% inf
      llik[j]  <- sum(dnorm(x = foi_hat,
                            mean = v_prev_vars$foi,
                            sd = v_prev_vars$foi_sd,
                            log = TRUE))
    }, error = function(e) NA)
    if (is.na(jj)) (llik[j] <- -Inf)
  }
  return(llik)
}

likelihood <- function(params, p_trans, waifw_assort, waifw_rp, inf) {
  exp(l_likelihood(params, p_trans, waifw_assort, waifw_rp, inf))
}

# Run IMIS
set.seed(1234)
for (i in seq_len(v_p)) {
  formals(l_likelihood) <- formals(likelihood) <-
    alist(params = , p_trans = v_p[i], waifw_assort = l_m_contacts_assort[[i]],
          waifw_rp = l_m_contacts_rp[[i]], inf = v_inf)
  imis_res <- IMIS(B = 1000, B.re = 1e4, number_k = 10, D = 1)
}

# Unique parameter sets
length(unique(imis_res$resample[, 1]))

# Effective sample size
sum(table(imis_res$resample[, 1]))^2 / sum(table(imis_res$resample[, 1])^2)

# Max weight
max(table(imis_res$resample[, 1])) / sum(table(imis_res$resample[, 1]))

likelihood(test)

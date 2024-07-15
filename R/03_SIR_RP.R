### Functions to run SIR for to get sample for RP WAIFW ###
foi_transition_matrix_rp <- function(v_betas, w,
                                     waifw_breaks = c(1, 5, 15, 25, 45,
                                                      55, 65, 70, 81)) {
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

get_waifw_draws <- function(m_cov, v_betas) {
  chol_mat <- chol(m_cov)
  v_in <- qnorm(runif(nrow(m_cov), 0, 1), 0, 1)  # vector of random inverse normal
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

v_race       <- c("Hispanic", "NH Black", "NH White", "Other")
v_race_names <- c("hisp", "black", "white", "other")

m_beta_hat_rp <- as.matrix(read.csv("results/waifw/rp_beta_hat_ga.csv"))
v_betas_rp_w1 <- as.numeric(m_beta_hat_rp[1, ])
v_betas_rp_w2 <- as.numeric(m_beta_hat_rp[2, ])
v_betas_rp_w3 <- as.numeric(m_beta_hat_rp[3, ])

m_cov_rp_w1 <- as.matrix(read.csv("results/waifw/rp_beta_cov_w1_ga.csv"))
m_cov_rp_w2 <- as.matrix(read.csv("results/waifw/rp_beta_cov_w2_ga.csv"))
m_cov_rp_w3 <- as.matrix(read.csv("results/waifw/rp_beta_cov_w3_ga.csv"))

m_waifw_rp_w1 <- as.matrix(read.csv("results/waifw/rp_Beta_hat_w1_ga.csv"))
m_waifw_rp_w2 <- as.matrix(read.csv("results/waifw/rp_Beta_hat_w2_ga.csv"))
m_waifw_rp_w3 <- as.matrix(read.csv("results/waifw/rp_Beta_hat_w3_ga.csv"))

##### Parallel Code #####
comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}

n <- 100000
n_cores <- 80

my_cluster <- parallel::makeCluster(
  n_cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my_cluster)
doRNG::registerDoRNG(seed = 12345)

out <- foreach(i = 1:n, .combine = "comb", .multicombine = TRUE,
               .packages = c("deSolve", "dplyr", "rootSolve", "dtplyr",
                             "purrr", "lubridate")) %dopar% {


v_rand_betas <- get_waifw_draws(m_cov = m_cov_rp,
                                v_betas = v_betas_rp)
rand_waifw <- get_new_waifw(old_waifw = waifw_rp,
                            old_betas = v_betas_rp,
                            new_betas = v_rand_betas)

  list(v_rand_betas_w1 = v_rand_betas_w1,
       v_rand_betas_w2 = v_rand_betas_w2,
       v_rand_betas_w3 = v_rand_betas_w3,
       nll_w1 
  )
}

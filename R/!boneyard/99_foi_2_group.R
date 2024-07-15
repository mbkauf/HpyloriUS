library(CVXR)

## Age groups
groups <- c(1:80)

## WAIFW parameters
waifw_breaks      <- c(1, 20, 55, 81)
beta0             <- c(0.8, 0.5, 0.3, 0.8, 0.5, 0.3, 0.6, 0.2, 0.1)
v_age_group_names <- c("1-19", "20-54", "55-80")

# w2 <- matrix(c(1, 2, 3, 7, 8, 9,
#                0, 2, 3, 8, 8, 9,
#                0, 0, 3, 9, 9, 9,
#                0, 0, 0, 4, 5, 6,
#                0, 0, 0, 0, 5, 6,
#                0, 0, 0, 0, 0, 6), ncol = 6, byrow = TRUE)

w2 <- matrix(c(1, 2, 3,
               0, 2, 3,
               0, 0, 3), ncol = 3, byrow = TRUE)


generate_nloglik <- function(betas, breaks, w, alphas,
                             upper = TRUE, inf = v_inf) {
  # Create WAIFW matrix
  betas  <- betas@value
  alphas <- alphas@value

  v_betas_1 <- betas[1:3] * alphas[1]
  v_betas_2 <- betas[4:6] * alphas[2]
  v_betas_3 <- betas[7:9] * alphas[3]

  print(v_betas_1)

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

  foi_hat <- beta %*% inf
  print(foi_hat)
  out <- -1 * sum(dnorm(x = foi_hat,
                        mean = c(v_prev_vars_white$foi, v_prev_vars_black$foi),
                        sd = c(v_prev_vars_white$foi_sd, v_prev_vars_black$foi_sd),
                        log = TRUE))
  return(out)
}

# get_nll <- function(beta, inf = v_inf, v_foi, v_foi_sd) {
#   foi_hat <- beta %*% inf
#   nll <- -1 * sum(dnorm(x = foi_hat, mean = v_foi,
#                   sd = v_foi_sd, log = TRUE))
#   return(nll)
# }

estimate_beta  <- function(w, i = NULL,
                           # beta_names = v_beta_names,  #nolint
                           v_breaks = waifw_breaks) {
  # require(MHadaptive)
  # require(texreg)
  # require(Matrix)
  require(CVXR)

  ### Set up variables for convex optimization
  v_waifw         <- Variable(9)
  value(v_waifw)  <- beta0
  v_alphas        <- Variable(3)
  value(v_alphas) <- c(0.2, 0.2, 0.2)

  l_constraints <- list(v_waifw >= 0,
                        v_waifw <= 50,
                        v_alphas <= 1 / 2,
                        v_alphas >= 0,
                        v_alphas[1] + v_alphas[2] + v_alphas[3] <= 1)
  nll <- generate_nloglik(betas = v_waifw, breaks = waifw_breaks,
                          w = w2, alphas = v_alphas)
  prob <- Problem(Minimize(nll), l_constraints)

  waifw_ll    <- CVXR::solve(prob)
  beta_llk    <- waifw_ll$value
  v_beta_hat  <- waifw_ll$getValue(v_waifw)
  v_alpha_hat <- waifw_ll$getValue(v_alphas)

  return(list(waifw_ll, beta_llk, v_beta_hat, v_alpha_hat))
}


v_demo_vars_white <- get_demographic_vars(spec_groups = groups,
                               race = "NH White")
v_prev_vars_white <- get_prevalence_vars(spec_groups = groups,
                                   race = "NH White")
v_inf_white <- v_prev_vars_white$prevalence * v_demo_vars_white$v_age_prop

v_demo_vars_black <- get_demographic_vars(spec_groups = groups,
                                          race = "NH Black")
v_prev_vars_black <- get_prevalence_vars(spec_groups = groups,
                                         race = "NH Black")
v_inf_black <- v_prev_vars_black$prevalence * v_demo_vars_black$v_age_prop

v_inf <- c(v_inf_white, v_inf_black)
v_foi <- c(v_prev_vars_white$foi, v_prev_vars_black$foi)
n_age_groups <- length(waifw_breaks) - 1

test <- estimate_beta(w = w2)


## Test
v_alphas <- c(0.2, 0.2, 0.2)
v_betas_1 <- beta0[1:3] * v_alphas[1]
v_betas_2 <- beta0[4:6] * v_alphas[2]
v_betas_3 <- beta0[7:9] * v_alphas[3]

beta1 <- get_transmission_matrix(betas = v_betas_1,  #nolint
                                 breaks = waifw_breaks,
                                 w = w2,
                                 upper = T)
beta2 <- get_transmission_matrix(betas = v_betas_2,  #nolint
                                 breaks = breaks,
                                 w = w,
                                 upper = upper)
beta3 <- get_transmission_matrix(betas = v_betas_3,  #nolint
                                 breaks = breaks,
                                 w = w,
                                 upper = upper)

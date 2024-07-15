library(scales)
library(ggplot2)
library(Matrix)

## Age groups
groups <- c(1:80)

## WAIFW parameters (3 age groups)
# waifw_breaks      <- c(1, 20, 55, 81)
# v_beta_0             <- c(0.8, 0.5, 0.3, 0.8, 0.5, 0.3, 0.6, 0.2, 0.1)
# # v_alphas_0 <- c(0.2, 0.2, 0.2)
# v_age_group_names <- c("1-19", "20-54", "55-80")
# w2 <- matrix(c(1, 2, 3,
#                0, 2, 3,
#                0, 0, 3), ncol = 3, byrow = TRUE)



## WAIFW parameters (8 age groups)
waifw_breaks      <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_beta_0          <- rep(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1), 3)

v_beta_0          <- c(0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1,
                       0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
                       0,0,0,0,0,0,0,0)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")
w2 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)


generate_nloglik <- function(param, breaks, w,
                             upper = TRUE, inf = v_inf) {
  # Create WAIFW matrix
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

  foi_hat <- beta %*% inf
  # print(foi_hat)
  out <- -1 * sum(dnorm(x = foi_hat,
                        mean = c(v_prev_vars_white$foi, v_prev_vars_black$foi),
                        sd = c(v_prev_vars_white$foi_sd, v_prev_vars_black$foi_sd),
                        log = TRUE))
  return(out)
}

estimate_beta  <- function(n_age_groups, w, i, beta_names = v_beta_names,  #nolint
                           v_breaks = waifw_breaks, group_names = groups) {
  require(MHadaptive)
  require(texreg)
  require(Matrix)

  ### Run optimization
  # waifw_ll   <- optim(par = c(v_beta_0, v_alphas_0),
  #                     fn = generate_nloglik,
  #                     breaks = v_breaks, w = w,
  #                     inf = v_inf,
  #                     hessian = TRUE,  method = "L-BFGS-B",
  #                     lower = c(rep(0, n_age_groups), c(0, 0, 0)),
  #                     upper = c(rep(100, n_age_groups), c(0.4, 0.4, 0.2)))

  waifw_ll   <- optim(par = v_beta_0,
                      fn = generate_nloglik,
                      breaks = v_breaks, w = w,
                      inf = v_inf,
                      hessian = TRUE,  method = "L-BFGS-B",
                      lower = rep(0, n_age_groups),
                      upper = rep(100, n_age_groups),
                      control = list(maxit = 1e4))
  # beta_llk   <- waifw_ll$value
  # v_beta_hat <- as.vector(waifw_ll$par)
  #
  # ## Check if HESSIAN is Positive Definite
  # ## If not, make covariance Positive Definite
  # ## Is Positive Definite?
  # if (MHadaptive::isPositiveDefinite(waifw_ll$hessian) == FALSE) {
  #   print("Hessian is NOT Positive Definite")
  #   m_hess <- Matrix::nearPD(waifw_ll$hessian)$mat
  #   beta_hat_cov <- solve(m_hess)
  #
  # } else {
  #   print("Hessian IS Positive Definite")
  #   print("No additional adjustment to COV matrix")
  #   beta_hat_cov <- solve(waifw_ll$hessian)
  # }
  #
  # beta_hat_cov <- as.matrix(beta_hat_cov)
  # ### Plot correlation matrix
  # beta_hat_cor <- cov2cor(beta_hat_cov)
  #
  # ### Compute SE
  # beta_hat_se <- sqrt(diag(beta_hat_cov))
  #
  # ### Generate big wAIFw matrices
  # beta_hat <- get_transmission_matrix(betas = v_beta_hat,  #nolint
  #                                     breaks = v_breaks,
  #                                     w = w)
  # beta_tex <- createTexreg(coef.names = beta_names,
  #                          coef = v_beta_hat, se = beta_hat_se,
  #                          gof = beta_llk, gof.names = as.character(i))
  #
  # return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
  #             beta_hat_se, beta_hat, beta_tex))
  return(waifw_ll)
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

test3 <- estimate_beta(w = w2)

### Calibrate alphas
## Assortative WAIFWs
beta_white_a <- get_transmission_matrix(betas = test4[["par"]][1:8],
                                        breaks = waifw_breaks,
                                        w = w2,
                                        group_names = groups)
beta_black_a <- get_transmission_matrix(betas = test4[["par"]][9:16],
                                        breaks = waifw_breaks,
                                        w = w2,
                                        group_names = groups)
beta_mix_a   <- get_transmission_matrix(betas = rep(0, 8),
                                        breaks = waifw_breaks,
                                        w = w2,
                                        group_names = groups)

## Random partnership WAIFWs
beta_white_r <- get_transmission_matrix(betas = test3[["par"]][1:8],
                                        breaks = waifw_breaks,
                                        w = w2,
                                        group_names = groups)
beta_black_r <- get_transmission_matrix(betas = test3[["par"]][9:16],
                                        breaks = waifw_breaks,
                                        w = w2,
                                        group_names = groups)
beta_mix_r   <- get_transmission_matrix(betas = test3[["par"]][17:24],
                                        breaks = waifw_breaks,
                                        w = w2,
                                        group_names = groups)

# Functions to calibrate alphas
generate_nloglik_alpha <- function(param, inf = v_inf,
                                   beta_1_a = beta_white_a,
                                   beta_2_a = beta_black_a,
                                   beta_3_a = beta_mix_a,
                                   beta_1_r = beta_white_r,
                                   beta_2_r = beta_black_r,
                                   beta_3_r = beta_mix_r) {

  beta1 <- ((1 - param[1]) *  beta_1_r) + (param[1] * beta_1_a)
  beta2 <- ((1 - param[1]) *  beta_2_r) + (param[1] * beta_2_a)
  beta3 <- ((1 - param[1]) *  beta_3_r)

  beta <- rbind(cbind(beta1, beta3),
                cbind(beta3, beta2))

  foi_hat <- beta %*% inf

  out <- -1 * sum(dnorm(x = foi_hat,
                        mean = c(v_prev_vars_white$foi, v_prev_vars_black$foi),
                        sd = c(v_prev_vars_white$foi_sd, v_prev_vars_black$foi_sd),
                        log = TRUE))
  return(out)
}

estimate_beta  <- function(beta_names = v_beta_names) {
  require(MHadaptive)
  require(texreg)
  require(Matrix)

  ### Run optimization
  # waifw_ll   <- optim(par = c(v_beta_0, v_alphas_0),
  #                     fn = generate_nloglik,
  #                     breaks = v_breaks, w = w,
  #                     inf = v_inf,
  #                     hessian = TRUE,  method = "L-BFGS-B",
  #                     lower = c(rep(0, n_age_groups), c(0, 0, 0)),
  #                     upper = c(rep(100, n_age_groups), c(0.4, 0.4, 0.2)))

  waifw_ll   <- optim(par = 0,
                      fn = generate_nloglik_alpha,
                      inf = v_inf,
                      hessian = TRUE,  method = "L-BFGS-B",
                      lower = 0,
                      upper = 0.1)
  # beta_llk   <- waifw_ll$value
  # v_beta_hat <- as.vector(waifw_ll$par)
  #
  # ## Check if HESSIAN is Positive Definite
  # ## If not, make covariance Positive Definite
  # ## Is Positive Definite?
  # if (MHadaptive::isPositiveDefinite(waifw_ll$hessian) == FALSE) {
  #   print("Hessian is NOT Positive Definite")
  #   m_hess <- Matrix::nearPD(waifw_ll$hessian)$mat
  #   beta_hat_cov <- solve(m_hess)
  #
  # } else {
  #   print("Hessian IS Positive Definite")
  #   print("No additional adjustment to COV matrix")
  #   beta_hat_cov <- solve(waifw_ll$hessian)
  # }
  #
  # beta_hat_cov <- as.matrix(beta_hat_cov)
  # ### Plot correlation matrix
  # beta_hat_cor <- cov2cor(beta_hat_cov)
  #
  # ### Compute SE
  # beta_hat_se <- sqrt(diag(beta_hat_cov))
  #
  # ### Generate big wAIFw matrices
  # beta_hat <- get_transmission_matrix(betas = v_beta_hat,  #nolint
  #                                     breaks = v_breaks,
  #                                     w = w)
  # beta_tex <- createTexreg(coef.names = beta_names,
  #                          coef = v_beta_hat, se = beta_hat_se,
  #                          gof = beta_llk, gof.names = as.character(i))
  #
  # return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
  #             beta_hat_se, beta_hat, beta_tex))
  return(waifw_ll)
}

alpha_calib <- estimate_beta()

## Plot FOI
alpha <- alpha_calib$par

alpha <- 0
beta1 <- ((1 - alpha) *  beta_white_r) + (alpha * beta_white_a)
beta2 <- ((1 - alpha) *  beta_black_r) + (alpha * beta_black_a)
beta3 <- ((1 - alpha) *  beta_mix_r)

beta <- rbind(cbind(beta1, beta3),
              cbind(beta3, beta2))

foi_hat <- beta %*% v_inf
-1 * sum(dnorm(x = foi_hat,
               mean = c(v_prev_vars_white$foi, v_prev_vars_black$foi),
               sd = c(v_prev_vars_white$foi_sd, v_prev_vars_black$foi_sd),
               log = TRUE))

# v_model_foi_white <- calibrate_white[["v_Beta_hat"]][[1]] %*% prev_white
# v_model_foi_black <- calibrate_black[["v_Beta_hat"]][[1]] %*% prev_black



df_foi <- as.data.frame(rbind(cbind(v_prev_vars_white$foi,
                                    v_prev_vars_white$foi_lb,
                                    v_prev_vars_white$foi_ub,
                                    rep("NH White", length(groups)),
                                    groups, foi_hat[1:80,]),
                              cbind(v_prev_vars_black$foi,
                                    v_prev_vars_black$foi_lb,
                                    v_prev_vars_black$foi_ub,
                                    rep("NH Black", length(groups)),
                                    groups,
                                    foi_hat[81:160,])))

colnames(df_foi) <- c("foi_est", "foi_est_lb", "foi_est_ub", "race",
                      "age", "foi_model")

df_foi$foi_est <- as.numeric(df_foi$foi_est)
df_foi$foi_est_lb <- as.numeric(df_foi$foi_est_lb)
df_foi$foi_est_ub <- as.numeric(df_foi$foi_est_ub)
df_foi$foi_model <- as.numeric(df_foi$foi_model)
df_foi$age <- as.numeric(df_foi$age)
df_foi$race2 <- df_foi$race

library(ggsci)
library(stringr)
library(dplyr)

str_stack <- function(x) {
  x %>% str_split("") %>% map(~ .x %>% paste(collapse = "\n")) #nolint
}

plot_foi_mix <- ggplot(data = df_foi) +
  facet_wrap(. ~ race, scales = "free") +
  scale_color_nejm(guide = "none") +
  scale_fill_nejm(guide = "none") +
  theme_bw() +
  geom_line(aes(x = age, y = foi_est, color = race, alpha = "Observed")) +
  geom_ribbon(aes(x = age, y = foi_est, ymax = foi_est_ub,
                  ymin = foi_est_lb, fill = race),
              alpha = 0.3, color = NA) +
  geom_point(aes(x = age, y = foi_model, color = race, alpha = "Model Predicted")) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  xlab("Age") + ylab("FOI") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  # ylim(0, 0.06) +
  scale_alpha_manual(name = NULL,
                     values = c(1, 1),
                     breaks = c("Model Predicted", "Observed"),
                     guide = guide_legend(override.aes =
                                            list(linetype = c(0, 1),
                                                 shape = c(16, NA),
                                                 color = "black")))
plot_foi_mix
ggsave(filename = paste0("results/FOI_mix_alpha_", alpha, ".png"),
       plot = plot_foi_mix)

# library(lattice)
library(pheatmap)
library(RColorBrewer)
library(viridis)

colnames(beta) <- seq(1:160)
rownames(beta) <- seq(1:160)
# levelplot(beta, ylim = c(160, 1), xlim = c(1, 160), tri.upper = T)
# List with colors for each annotation.
quantile_breaks <- function(xs, n = 80) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(beta, n = 81)


pheatmap(mat = beta, cluster_rows=F, cluster_cols=F, color = inferno(23),
         breaks = mat_breaks, labels_row = rep("", 180),
         labels_col = rep("", 180))

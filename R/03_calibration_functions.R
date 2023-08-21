library(scales)
library(ggplot2)
library(Matrix)

generate_nloglik <- function(betas, breaks, w, upper = TRUE, inf = v_inf) {
  # Create wAIFw matrix
  beta <- get_transmission_matrix(betas = betas, breaks = breaks, w = w, #nolint
                                  upper = upper)
  foi_hat <- beta %*% inf
  out <- -1 * sum(dnorm(x = foi_hat, mean = v_prev_vars$foi,
                        sd = v_prev_vars$foi_sd, log = TRUE))
  return(out)
}


estimate_beta  <- function(n_age_groups, w, i, beta_names = v_beta_names,  #nolint
                           v_breaks = waifw_breaks, group_names = groups) {
  require(MHadaptive)
  require(texreg)
  require(Matrix)

  ### Run optimization
  waifw_ll   <- optim(par = beta0, fn = generate_nloglik,
                      breaks = v_breaks, w = w,
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

  ### Generate big wAIFw matrices
  beta_hat <- get_transmission_matrix(betas = v_beta_hat,  #nolint
                                      breaks = v_breaks,
                                      w = w)
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat, se = beta_hat_se,
                           gof = beta_llk, gof.names = as.character(i))
  print(typeof(beta_hat_cov))
  return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
              beta_hat_se, beta_hat, beta_tex))
}


calibrate_betas <- function(n_age_groups, n_betas, v_waifw_structure) {
  ## Beta and wAIFw latex names
  v_beta_names <- paste0(paste0("$", paste("\\beta", 1:8, sep = "_")), "$")
  v_waifw_names <- paste0(paste0("$", paste("w", 1:n_betas, sep = "_")), "$")

  ## Initialize vectors, matrices and arrays
  v_waifw_ll     <- vector("list", n_betas)
  v_beta_llk     <- numeric(n_betas)
  m_beta_hat     <- matrix(0, ncol = n_age_groups, nrow = n_betas,
                           dimnames = list(paste0("w", seq(1:n_betas)),
                                           paste0("b", seq(1:n_age_groups))))
  m_beta_hat_se  <- matrix(0, ncol = n_age_groups, nrow = n_betas,
                           dimnames = list(paste0("w", seq(1:n_betas)),
                                           paste0("b", seq(1:n_age_groups))))
  rownames(m_beta_hat) <- rownames(m_beta_hat_se) <- v_waifw_names

  a_beta_hat_cov <- array(0, dim = c(n_age_groups, n_age_groups, n_betas),
                          dimnames = list(paste0("b", seq(1:n_age_groups)),
                                          paste0("b", seq(1:n_age_groups)),
                                          paste0("w", seq(1:n_betas))))
  a_beta_hat_cor <- a_beta_hat_cov
  v_beta_hat     <- vector("list", n_betas)
  v_beta_tex     <- vector("list", n_betas)
  names(v_beta_tex) <- v_waifw_names

  for (i in 1:n_betas) {
    results_beta <-  estimate_beta(n_age_groups = n_age_groups,
                                   w = v_waifw_structure[[i]],
                                   i = i,
                                   beta_names = v_beta_names)
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



## Variables needed for calibration
# wAIFw matricies
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

# w4 <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1,
#                2, 2, 2, 2, 2, 2, 2, 2,
#                3, 3, 3, 3, 3, 3, 3, 3,
#                4, 4, 4, 4, 4, 4, 4, 4,
#                5, 5, 5, 5, 5, 5, 5, 5,
#                6, 6, 6, 6, 6, 6, 6, 6,
#                7, 7, 7, 7, 7, 7, 7, 7,
#                8, 8, 8, 8, 8, 8, 8, 8), ncol = 8, byrow = T)

v_waifw_structure <- list(w1, w2, w3) # , w4

## WAIFW parameters
waifw_breaks      <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
beta0             <- c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

## Age groups
groups <- c(1:80)

## Number of age groups
n_ages <- length(groups)

## Population growth
q <- 0

n_betas <- length(v_waifw_structure)
n_age_groups <- length(waifw_breaks) - 1

v_race <- c("Hispanic", "NH White", "NH Black", "Other")

# v_race <- c("Hispanic", "Other")
# v_race <- c("Non-Hispanic Black")
# v_waifw_structure = list(w3)

for (race in v_race) {
  if (race == "NH White") {
    race_text <- "white"
  } else if (race == "NH Black") {
    race_text <- "black"
  } else if (race == "Hispanic") {
    race_text <- "hispanic"
  } else if (race == "Other") {
    race_text <- "other"
  }
  print(race)
  v_demo <- get_demographic_vars(spec_groups = groups,
                                 race = race)
  v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                     race = race,
                                     birth_cohort = 1960)
  v_inf <- v_prev_vars$prevalence * v_demo$v_age_prop
  tmp <- calibrate_betas(n_age_groups = n_age_groups,
                         n_betas = n_betas,
                         v_waifw_structure = v_waifw_structure)
  assign(paste0("calibrate_", race_text), tmp)
}

### Plot calibration
simple_waifw <- function(w, v_betas) {
  if(w == "w1"){
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

white_waifw_w1 <- simple_waifw(w = "w1",
                               v_betas = calibrate_white[["m_beta_hat"]][1, ])
white_waifw_w2 <- simple_waifw(w = "w2",
                               v_betas = calibrate_white[["m_beta_hat"]][2, ])
white_waifw_w3 <- simple_waifw(w = "w3",
                               v_betas = calibrate_white[["m_beta_hat"]][3, ])

df_white <- rbind(as.data.frame(as.table(white_waifw_w1)),
                  as.data.frame(as.table(white_waifw_w2)),
                  as.data.frame(as.table(white_waifw_w3))) %>%
  mutate(w = c(rep("w1", 64), rep("w2", 64), rep("w3", 64))) %>%
  mutate(race = "NH white")


black_waifw_w1 <- simple_waifw(w = "w1",
                               v_betas = calibrate_black[["m_beta_hat"]][1, ])
black_waifw_w2 <- simple_waifw(w = "w2",
                               v_betas = calibrate_black[["m_beta_hat"]][2, ])
black_waifw_w3 <- simple_waifw(w = "w3",
                               v_betas = calibrate_black[["m_beta_hat"]][3, ])

df_black <- rbind(as.data.frame(as.table(black_waifw_w1)),
                  as.data.frame(as.table(black_waifw_w2)),
                  as.data.frame(as.table(black_waifw_w3))) %>%
  mutate(w = c(rep("w1", 64), rep("w2", 64), rep("w3", 64))) %>%
  mutate(race = "NH Black")

hispanic_waifw_w1 <- simple_waifw(w = "w1",
                                  v_betas = calibrate_hispanic
                                  [["m_beta_hat"]][1, ])
hispanic_waifw_w2 <- simple_waifw(w = "w2",
                                  v_betas = calibrate_hispanic
                                  [["m_beta_hat"]][2, ])
hispanic_waifw_w3 <- simple_waifw(w = "w3",
                                  v_betas = calibrate_hispanic
                                  [["m_beta_hat"]][3, ])

df_hispanic <- rbind(as.data.frame(as.table(hispanic_waifw_w1)),
                     as.data.frame(as.table(hispanic_waifw_w2)),
                     as.data.frame(as.table(hispanic_waifw_w3))) %>%
  mutate(w = c(rep("w1", 64), rep("w2", 64), rep("w3", 64))) %>%
  mutate(race = "Hispanic")

other_waifw_w1 <- simple_waifw(w = "w1",
                               v_betas = calibrate_other[["m_beta_hat"]][1, ])
other_waifw_w2 <- simple_waifw(w = "w2",
                               v_betas = calibrate_other[["m_beta_hat"]][2, ])
other_waifw_w3 <- simple_waifw(w = "w3",
                               v_betas = calibrate_other[["m_beta_hat"]][3, ])

df_other <- rbind(as.data.frame(as.table(other_waifw_w1)),
                  as.data.frame(as.table(other_waifw_w2)),
                  as.data.frame(as.table(other_waifw_w3))) %>%
  mutate(w = c(rep("w1", 64), rep("w2", 64), rep("w3", 64))) %>%
  mutate(race = "Other")


df_all <- rbind(df_white, df_black, df_hispanic, df_other)

plot_all <- ggplot(data = df_all, aes(x = Var1,
                                      y = Var2,
                                      fill = Freq)) +
  geom_tile() + facet_grid(w ~ race) +
  scale_fill_gradientn(name = "Beta Value", trans = "log",
                       # breaks = c(0.0025, 0.0498, 1.000),
                       colors = c("white", "#FFF5F0", "#FEE0D2", "#FCBBA1",
                                  "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D",
                                  "#A50F15", "#67000D")) +
  scale_y_discrete(limits = rev) + scale_x_discrete() + theme_bw() +
  xlab("Age Group") + ylab("Age Group") +
  theme(axis.text.x = element_text(size = 8, angle = -45,
                                   vjust = 0.2, hjust = 0.4),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.height = unit(1.5, "cm"))
plot_all

ggsave(filename = "results/wAIFw_calibration.png",
       plot = plot_all, width = 10, height = 4)

## Plot FOI
v_prev_vars_white <- get_prevalence_vars(spec_groups = groups,
                                         race = "NH white")
v_prev_vars_black <- get_prevalence_vars(spec_groups = groups,
                                         race = "NH Black")
v_prev_vars_hisp  <- get_prevalence_vars(spec_groups = groups,
                                         race = "Hispanic")
v_prev_vars_other <- get_prevalence_vars(spec_groups = groups,
                                         race = "Other")
v_demo_white <- get_demographic_vars(spec_groups = groups, race = "NH white")
v_demo_black <- get_demographic_vars(spec_groups = groups, race = "NH Black")
v_demo_hisp  <- get_demographic_vars(spec_groups = groups, race = "Hispanic")
v_demo_other <- get_demographic_vars(spec_groups = groups, race = "Other")


prev_white <- v_prev_vars_white$prevalence * v_demo_white$v_age_prop
prev_black <- v_prev_vars_black$prevalence * v_demo_black$v_age_prop
prev_hisp  <- v_prev_vars_hisp$prevalence * v_demo_hisp$v_age_prop
prev_other <- v_prev_vars_other$prevalence * v_demo_other$v_age_prop

v_model_foi_white <- calibrate_white[["v_Beta_hat"]][[1]] %*% prev_white
v_model_foi_black <- calibrate_black[["v_Beta_hat"]][[1]] %*% prev_black
v_model_foi_hisp  <- calibrate_hispanic[["v_Beta_hat"]][[2]] %*% prev_hisp
v_model_foi_other <- calibrate_other[["v_Beta_hat"]][[1]] %*% prev_other

df_foi <- as.data.frame(rbind(cbind(v_prev_vars_white$foi,
                                    v_prev_vars_white$foi_lb,
                                    v_prev_vars_white$foi_ub,
                                    rep("NH white", length(groups)),
                                    groups, v_model_foi_white),
                              cbind(v_prev_vars_black$foi,
                                    v_prev_vars_black$foi_lb,
                                    v_prev_vars_black$foi_ub,
                                    rep("NH Black", length(groups)),
                                    groups, v_model_foi_black),
                              cbind(v_prev_vars_hisp$foi,
                                    v_prev_vars_hisp$foi_lb,
                                    v_prev_vars_hisp$foi_ub,
                                    rep("Hispanic", length(groups)),
                                    groups, v_model_foi_hisp),
                              cbind(v_prev_vars_other$foi,
                                    v_prev_vars_other$foi_lb,
                                    v_prev_vars_other$foi_ub,
                                    rep("Other", length(groups)),
                                    groups, v_model_foi_other)))

colnames(df_foi) <- c("foi_est", "foi_est_lb", "foi_est_ub","race",
                      "age", "foi_model")

df_foi$foi_est <- as.numeric(df_foi$foi_est)
df_foi$foi_est_lb <- as.numeric(df_foi$foi_est_lb)
df_foi$foi_est_ub <- as.numeric(df_foi$foi_est_ub)
df_foi$foi_model <- as.numeric(df_foi$foi_model)
df_foi$age <- as.numeric(df_foi$age)
df_foi$race2 <- df_foi$race

library(ggsci)
library(tidyverse)
library(dplyr)

str_stack <- function(x) {
  x %>% str_split("") %>% map(~ .x %>% paste(collapse = "\n")) #nolint
}

plot_foi_valid <- ggplot(data = df_foi) +
  facet_wrap(. ~ race, scales = "free") +
  scale_color_nejm(guide = "none") +
  scale_fill_nejm(guide = "none") +
  theme_bw() +
  geom_line(aes(x = age, y = foi_est, color = race, alpha = "Observed")) +
  geom_ribbon(aes(x = age, y = foi_est, ymax = foi_est_ub,
                  ymin = foi_est_lb, fill = race),
              alpha=0.3, color = NA) +
  geom_point(aes(x = age, y = foi_model, color = race, alpha = "Fitted")) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  xlab("Age") + ylab("FOI") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_alpha_manual(name = NULL,
                     values = c(1, 1),
                     breaks = c("Fitted", "Observed"),
                     guide = guide_legend(override.aes =
                                            list(linetype = c(0, 1),
                                                 shape = c(16, NA),
                                                 color = "black")))
plot_foi_valid

ggsave(filename = "results/FOI_validation.png",
       plot = plot_foi_valid, width = 10, height = 4)

### Generate Latex Table of Results
all_tex <- list(calibrate_hispanic[["v_beta_tex"]][["$w_1$"]][[1]],
                calibrate_hispanic[["v_beta_tex"]][["$w_2$"]][[1]],
                calibrate_hispanic[["v_beta_tex"]][["$w_3$"]][[1]],
                calibrate_black[["v_beta_tex"]][["$w_1$"]][[1]],
                calibrate_black[["v_beta_tex"]][["$w_2$"]][[1]],
                calibrate_black[["v_beta_tex"]][["$w_3$"]][[1]],
                calibrate_white[["v_beta_tex"]][["$w_1$"]][[1]],
                calibrate_white[["v_beta_tex"]][["$w_2$"]][[1]],
                calibrate_white[["v_beta_tex"]][["$w_3$"]][[1]],
                calibrate_other[["v_beta_tex"]][["$w_1$"]][[1]],
                calibrate_other[["v_beta_tex"]][["$w_2$"]][[1]],
                calibrate_other[["v_beta_tex"]][["$w_3$"]][[1]])

waifw_names_hisp  <- paste0("Hispanic - ",
                            paste0("$", paste("w", 1:3, sep = "_")), "$")
waifw_names_black <- paste0("NH Black - ",
                            paste0("$", paste("w", 1:3, sep = "_")), "$")
waifw_names_white <- paste0("NH white - ",
                            paste0("$", paste("w", 1:3, sep = "_")), "$")
waifw_names_other <- paste0("Other - ",
                            paste0("$", paste("w", 1:3, sep = "_")), "$")
waifw_names_all   <- c(waifw_names_hisp,
                       waifw_names_black,
                       waifw_names_white,
                       waifw_names_other)

texreg(all_tex, digits = 3, stars = numeric(0), booktabs = TRUE,
       custom.model.names = waifw_names_all)

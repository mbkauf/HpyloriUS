library(scales)
library(ggplot2)
library(Matrix)

generate_nloglik <- function(betas, breaks, W, upper = TRUE, inf = v_inf){
  # Create WAIFW matrix
  Beta <- get_transmission_matrix(betas = betas, breaks = breaks, W = W,
                                  upper = upper)
  foi_hat <- Beta %*% inf
  out <- -1 * sum(dnorm(x = foi_hat, mean = v_prev_vars$foi,
                        sd = v_prev_vars$foi_sd, log = T))
  # out <- -1 * sum(dgamma(x = foi_hat,
  #                        shape = (v_prev_vars$foi)^2/(v_prev_vars$foi_sd^2),
  #                        rate = (v_prev_vars$foi_sd^2)/v_prev_vars$foi, log = T))
  return(out)
}


estimate_Beta  <- function(n_age_groups, W, i, beta_names = v_beta_names,
                           v_breaks = waifw.breaks, group.names = groups) {
  require(MHadaptive)
  require(texreg)
  require(Matrix)

  ### Run optimization
  waifw_ll   <- optim(par = beta0, fn = generate_nloglik,
                      breaks = v_breaks, W = W,
                      hessian = T,  method = "L-BFGS-B",
                      lower = rep(0, n_age_groups), upper = rep(100, n_age_groups))
  beta_llk   <- waifw_ll$value
  v_beta_hat <- as.vector(waifw_ll$par)

  ## Check if HESSIAN is Positive Definite; If not, make covariance Positive Definite
  ## Is Positive Definite?
  if(MHadaptive::isPositiveDefinite(waifw_ll$hessian)==FALSE){
    print("Hessian is NOT Positive Definite")
    m.hess <- Matrix::nearPD(waifw_ll$hessian)$mat
    beta_hat_cov <- solve(m.hess)

  } else{
    print("Hessian IS Positive Definite")
    print("No additional adjustment to COV matrix")
    beta_hat_cov <- solve(waifw_ll$hessian)
  }

  print(beta_hat_cov)
  ### Plot correlation matrix
  beta_hat_cor <- cov2cor(beta_hat_cov)

  ### Compute SE
  beta_hat_se <- sqrt(diag(beta_hat_cov))

  ### Generate big WAIFW matrices
  Beta_hat <- get_transmission_matrix(betas = v_beta_hat,
                                      breaks = v_breaks,
                                      W = W)
  beta_tex <- createTexreg(coef.names = beta_names,
                           coef = v_beta_hat, se = beta_hat_se,
                           gof = beta_llk, gof.names = as.character(i))

  return(list(waifw_ll, beta_llk, v_beta_hat, beta_hat_cov, beta_hat_cor,
              beta_hat_se, Beta_hat, beta_tex))
}


calibrate_Betas <- function(n_age_groups, n_Betas, v_waifw_structure) {
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



## Variables needed for calibration
# WAIFW matricies
W1 <- matrix(c(1, 9, 9, 9, 9, 9, 9, 9,
               0, 2, 9, 9, 9, 9, 9, 9,
               0, 0, 3, 9, 9, 9, 9, 9,
               0, 0, 0, 4, 9, 9, 9, 9,
               0, 0, 0, 0, 5, 9, 9, 9,
               0, 0, 0, 0, 0, 6, 9, 9,
               0, 0, 0, 0, 0, 0, 7, 9,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = T)

W2 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = T)

W3 <- matrix(c(1, 1, 1, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = T)

# W4 <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1,
#                2, 2, 2, 2, 2, 2, 2, 2,
#                3, 3, 3, 3, 3, 3, 3, 3,
#                4, 4, 4, 4, 4, 4, 4, 4,
#                5, 5, 5, 5, 5, 5, 5, 5,
#                6, 6, 6, 6, 6, 6, 6, 6,
#                7, 7, 7, 7, 7, 7, 7, 7,
#                8, 8, 8, 8, 8, 8, 8, 8), ncol = 8, byrow = T)

v_waifw_structure = list(W1, W2, W3) # , W4

n_Betas <- length(v_waifw_structure)
n_age_groups <- length(waifw.breaks) - 1

v_race <- c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other")

v_race <- c("Hispanic", "Other")
v_race <- c("Non-Hispanic Black")
v_waifw_structure = list(W3)

for(race in v_race) {
  if(race=="Non-Hispanic White") {
    race_text <- "White"
  } else if(race=="Non-Hispanic Black") {
    race_text <- "Black"
  } else if(race=="Hispanic") {
    race_text <- "Hispanic"
  } else if(race=="Other") {
    race_text <- "Other"
  }
    v_demo <- get_demographic_vars(spec_groups = groups,
                                   Race = race)
    v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                       Race = race,
                                       birth_cohort = 1960)
    v_inf <- v_prev_vars$prevalence * v_demo$v_age_prop
    tmp <- calibrate_Betas(n_age_groups = n_age_groups,
                           n_Betas = n_Betas,
                           v_waifw_structure = v_waifw_structure)
    assign(paste0("calibrate", race_text), tmp)
}

### Plot calibration
simple_waifw <- function(w, v_betas) {
  if(w=="W1"){
    waifw <- matrix(c(v_betas[1], 0, 0, 0, 0, 0, 0, 0,
                      0, v_betas[2], 0, 0, 0, 0, 0, 0,
                      0, 0, v_betas[3], 0, 0, 0, 0, 0,
                      0, 0, 0, v_betas[4], 0, 0, 0, 0,
                      0, 0, 0, 0, v_betas[5], 0, 0, 0,
                      0, 0, 0, 0, 0, v_betas[6], 0, 0,
                      0, 0, 0, 0, 0, 0, v_betas[7], 0,
                      0, 0, 0, 0, 0, 0, 0, v_betas[8]),
                    ncol = 8, byrow = T)
  } else if (w=="W2") {
    waifw <- matrix(c(v_betas[1], v_betas[1], v_betas[3], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[1], v_betas[2], v_betas[3], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[3], v_betas[3], v_betas[3], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[4], v_betas[4], v_betas[4], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[5], v_betas[5], v_betas[5], v_betas[5], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[6], v_betas[6], v_betas[6], v_betas[6], v_betas[6], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[8],
                      v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8]),
                    ncol = 8, byrow = T)
  } else if (w=="W3") {
    waifw <- matrix(c(v_betas[1], v_betas[1], v_betas[1], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[1], v_betas[2], v_betas[3], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[1], v_betas[3], v_betas[3], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[4], v_betas[4], v_betas[4], v_betas[4], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[5], v_betas[5], v_betas[5], v_betas[5], v_betas[5], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[6], v_betas[6], v_betas[6], v_betas[6], v_betas[6], v_betas[6], v_betas[7], v_betas[8],
                      v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[7], v_betas[8],
                      v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8], v_betas[8]),
                    ncol = 8, byrow = T)
  }
  rownames(waifw) <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64", "65-69", "70-80")
  colnames(waifw) <-c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64", "65-69", "70-80")

  return(waifw)
}

white_waifw_w1 <- simple_waifw(w="W1", v_betas = calibrateWhite[["m_beta_hat"]][1,])
white_waifw_w2 <- simple_waifw(w="W2", v_betas = calibrateWhite[["m_beta_hat"]][2,])
white_waifw_w3 <- simple_waifw(w="W3", v_betas = calibrateWhite[["m_beta_hat"]][3,])

df_white <- rbind(as.data.frame(as.table(white_waifw_w1)),
                       as.data.frame(as.table(white_waifw_w2)),
                       as.data.frame(as.table(white_waifw_w3))) %>%
  mutate(w = c(rep("W1", 64), rep("W2", 64), rep("W3", 64))) %>%
  mutate(race = "NH White")


black_waifw_w1 <- simple_waifw(w="W1", v_betas = calibrateBlack[["m_beta_hat"]][1,])
black_waifw_w2 <- simple_waifw(w="W2", v_betas = calibrateBlack[["m_beta_hat"]][2,])
black_waifw_w3 <- simple_waifw(w="W3", v_betas = calibrateBlack[["m_beta_hat"]][3,])

df_black <- rbind(as.data.frame(as.table(black_waifw_w1)),
                       as.data.frame(as.table(black_waifw_w2)),
                       as.data.frame(as.table(black_waifw_w3))) %>%
  mutate(w = c(rep("W1", 64), rep("W2", 64), rep("W3", 64))) %>%
  mutate(race = "NH Black")

hispanic_waifw_w1 <- simple_waifw(w="W1", v_betas = calibrateHispanic[["m_beta_hat"]][1,])
hispanic_waifw_w2 <- simple_waifw(w="W2", v_betas = calibrateHispanic[["m_beta_hat"]][2,])
hispanic_waifw_w3 <- simple_waifw(w="W3", v_betas = calibrateHispanic[["m_beta_hat"]][3,])

df_hispanic <- rbind(as.data.frame(as.table(hispanic_waifw_w1)),
                     as.data.frame(as.table(hispanic_waifw_w2)),
                     as.data.frame(as.table(hispanic_waifw_w3))) %>%
  mutate(w = c(rep("W1", 64), rep("W2", 64), rep("W3", 64))) %>%
  mutate(race = "Hispanic")

other_waifw_w1 <- simple_waifw(w="W1", v_betas = calibrateOther[["m_beta_hat"]][1,])
other_waifw_w2 <- simple_waifw(w="W2", v_betas = calibrateOther[["m_beta_hat"]][2,])
other_waifw_w3 <- simple_waifw(w="W3", v_betas = calibrateOther[["m_beta_hat"]][3,])

df_other <- rbind(as.data.frame(as.table(other_waifw_w1)),
                  as.data.frame(as.table(other_waifw_w2)),
                  as.data.frame(as.table(other_waifw_w3))) %>%
  mutate(w = c(rep("W1", 64), rep("W2", 64), rep("W3", 64))) %>%
  mutate(race = "Other")


df_all <- rbind(df_white, df_black, df_hispanic, df_other)

plot_all <- ggplot(data = df_all, aes(x = Var1,
                                      y = Var2,
                                      fill = Freq)) +
  geom_tile() + facet_grid(w~race) +
  scale_fill_gradientn(name = "Beta Value", trans = "log", # oob=squish_infinite,
                       # breaks = c(0.0025, 0.0498, 1.000),
                       colors = c("white", "#FFF5F0", "#FEE0D2", "#FCBBA1",
                                  "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D",
                                  "#A50F15", "#67000D")) +
  scale_y_discrete(limits=rev) + scale_x_discrete() + theme_bw() +
  xlab("Age Group") + ylab("Age Group") +
  theme(axis.text.x = element_text(size = 8, angle = -45, vjust = 0.2, hjust = 0.4),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.height = unit(1.5, "cm"))
plot_all

ggsave(filename = "results/WAIFW_calibration.png",
       plot = plot_all , width = 10, height = 4)

## Plot FOI
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

v_model_foi_white <- calibrateWhite[["v_Beta_hat"]][[1]] %*% prev_white
v_model_foi_black <- calibrateBlack[["v_Beta_hat"]][[1]] %*% prev_black
v_model_foi_hisp  <- calibrateHispanic[["v_Beta_hat"]][[2]] %*% prev_hisp
v_model_foi_other <- calibrateOther[["v_Beta_hat"]][[1]] %*% prev_other

df_foi <- as.data.frame(rbind(cbind(v_prev_vars_white$foi,
                                    v_prev_vars_white$foi_lb,
                                    v_prev_vars_white$foi_ub,
                                    rep("NH White", length(groups)),
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
colnames(df_foi) <- c("foi_est", "foi_est_lb", "foi_est_ub","race", "age", "foi_model")
df_foi$foi_est <- as.numeric(df_foi$foi_est)
df_foi$foi_est_lb <- as.numeric(df_foi$foi_est_lb)
df_foi$foi_est_ub <- as.numeric(df_foi$foi_est_ub)
df_foi$foi_model <- as.numeric(df_foi$foi_model)
df_foi$age <- as.numeric(df_foi$age)
df_foi$race2 <- df_foi$race

library(ggsci)
library(tidyverse)
str_stack <- function(x) {
  x %>% str_split("") %>% map(~ .x %>% paste(collapse = "\n"))
}
plot_foi_valid <- ggplot(data = df_foi) +
  facet_wrap(.~race, scales="free") +
  scale_color_nejm(guide = "none") +
  scale_fill_nejm(guide = "none") +
  theme_bw() +
  geom_line(aes(x = age, y = foi_est, color = race, alpha = "Observed")) +
  geom_ribbon(aes(x=age, y=foi_est, ymax=foi_est_ub, ymin=foi_est_lb, fill = race), alpha=0.3, color = NA) +
  geom_point(aes(x=age, y=foi_model, color = race, alpha = "Fitted")) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  xlab("Age") + ylab("FOI") +
  scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_alpha_manual(name = NULL,
                     values = c(1, 1),
                     breaks = c("Fitted", "Observed"),
                     guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                              shape = c(16, NA),
                                                              color = "black")))
plot_foi_valid

ggsave(filename = "results/FOI_validation.png",
       plot = plot_foi_valid , width = 10, height = 4)

### Generate Latex Table of Results
all_tex <- list(calibrateHispanic[["v_beta_tex"]][["$W_1$"]][[1]],
                calibrateHispanic[["v_beta_tex"]][["$W_2$"]][[1]],
                calibrateHispanic[["v_beta_tex"]][["$W_3$"]][[1]],
                calibrateBlack[["v_beta_tex"]][["$W_1$"]][[1]],
                calibrateBlack[["v_beta_tex"]][["$W_2$"]][[1]],
                calibrateBlack[["v_beta_tex"]][["$W_3$"]][[1]],
                calibrateWhite[["v_beta_tex"]][["$W_1$"]][[1]],
                calibrateWhite[["v_beta_tex"]][["$W_2$"]][[1]],
                calibrateWhite[["v_beta_tex"]][["$W_3$"]][[1]],
                calibrateOther[["v_beta_tex"]][["$W_1$"]][[1]],
                calibrateOther[["v_beta_tex"]][["$W_2$"]][[1]],
                calibrateOther[["v_beta_tex"]][["$W_3$"]][[1]])

waifw.names.hisp <- paste0("Hispanic - ", paste0("$",paste("W", 1:3, sep = "_")), "$")
waifw.names.black <- paste0("NH Black - ", paste0("$",paste("W", 1:3, sep = "_")), "$")
waifw.names.white <- paste0("NH White - ", paste0("$",paste("W", 1:3, sep = "_")), "$")
waifw.names.other <- paste0("Other - ", paste0("$",paste("W", 1:3, sep = "_")), "$")
waifw.names.all <- c(waifw.names.hisp,
                        waifw.names.black,
                        waifw.names.white,
                        waifw.names.other)

texreg(all_tex, digits = 3, stars = numeric(0), booktabs = T,
       custom.model.names = waifw.names.all)



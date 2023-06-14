
generate_nloglik <- function(betas, breaks, W, upper = TRUE, inf = v_inf){
  # Create WAIFW matrix
  Beta <- get_transmission_matrix(betas = betas, breaks = breaks, W = W,
                                  upper = upper)
  foi_hat <- Beta %*% inf
  out <- -sum(dnorm(x = v_prev_vars$foi, mean = foi_hat,
                    sd = v_prev_vars$foi_sd, log = T))
  return(out)
}


estimate_Beta  <- function(n_age_groups, W, i, beta_names = v_beta_names,
                           v_breaks = waifw.breaks, group.names = groups) {
  require(MHadaptive)
  require(texreg)

  ### Run optimization
  waifw_ll   <- optim(par = beta0, fn = generate_nloglik,
                      breaks = v_breaks, W = W,
                      hessian = T,  method = "L-BFGS-B",
                      lower = rep(0, n_age_groups), upper = rep(20, n_age_groups))
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

v_race <- c("Non-Hispanic White", "Non-Hispanic Black", "Other Hispanic",
            "Mexican-American", "Other")
v_sex  <- c("Male", "Female")

for(race in v_race) {
  if(race=="Non-Hispanic White") {
    race_text <- "White"
  } else if(race=="Non-Hispanic Black") {
    race_text <- "Black"
  } else if(race=="Other Hispanic") {
    race_text <- "Hispanic"
  } else if(race=="Mexican-American") {
    race_text <- "Mexican"
  } else if(race=="Other") {
    race_text <- "Other"
  }
  for(sex in v_sex) {
    v_demo <- get_demographic_vars(spec_groups = groups,
                                   Race = race, Sex = sex)
    v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                       Race = race,
                                       Sex = sex)
    v_inf <- v_prev_vars$prevalence * v_demo$v_age_prop
    tmp <- calibrate_Betas(n_age_groups = n_age_groups,
                           n_Betas = n_Betas,
                           v_waifw_structure = v_waifw_structure)
    assign(paste0("calibrate", race_text, sex), tmp)
  }
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
  rownames(waifw) <- c("0-4", "5-14", "15-24", "25-44", "45-54", "55-64", "65-69", "70-80")
  colnames(waifw) <-c("0-4", "5-14", "15-24", "25-44", "45-54", "55-64", "65-69", "70-80")

  return(waifw)
}

white_male_waifw_w1 <- simple_waifw(w="W1", v_betas = calibrateWhiteMale[["m_beta_hat"]][1,])
white_male_waifw_w2 <- simple_waifw(w="W2", v_betas = calibrateWhiteMale[["m_beta_hat"]][2,])
white_male_waifw_w3 <- simple_waifw(w="W3", v_betas = calibrateWhiteMale[["m_beta_hat"]][3,])

df_white_male <- rbind(as.data.frame(as.table(white_male_waifw_w1)),
                       as.data.frame(as.table(white_male_waifw_w2)),
                       as.data.frame(as.table(white_male_waifw_w3))) %>%
  mutate(w = c(rep("W1", 64), rep("W2", 64), rep("W3", 64))) %>%
  mutate(race_sex = "White - Males")


black_male_waifw_w1 <- simple_waifw(w="W1", v_betas = calibrateBlackMale[["m_beta_hat"]][1,])
black_male_waifw_w2 <- simple_waifw(w="W2", v_betas = calibrateBlackMale[["m_beta_hat"]][2,])
black_male_waifw_w3 <- simple_waifw(w="W3", v_betas = calibrateBlackMale[["m_beta_hat"]][3,])

df_black_male <- rbind(as.data.frame(as.table(black_male_waifw_w1)),
                       as.data.frame(as.table(black_male_waifw_w2)),
                       as.data.frame(as.table(black_male_waifw_w3))) %>%
  mutate(w = c(rep("W1", 64), rep("W2", 64), rep("W3", 64))) %>%
  mutate(race_sex = "Black - Males")

hispanic_male_waifw_w1 <- simple_waifw(w="W1", v_betas = calibrateHispanicMale[["m_beta_hat"]][1,])
hispanic_male_waifw_w2 <- simple_waifw(w="W2", v_betas = calibrateHispanicMale[["m_beta_hat"]][2,])
hispanic_male_waifw_w3 <- simple_waifw(w="W3", v_betas = calibrateHispanicMale[["m_beta_hat"]][3,])

df_hispanic_male <- rbind(as.data.frame(as.table(hispanic_male_waifw_w1)),
                       as.data.frame(as.table(hispanic_male_waifw_w2)),
                       as.data.frame(as.table(hispanic_male_waifw_w3))) %>%
  mutate(w = c(rep("W1", 64), rep("W2", 64), rep("W3", 64))) %>%
  mutate(race_sex = "Hispanic - Males")

df_all_male <- rbind(df_white_male, df_black_male, df_hispanic_male)

plot_all_male <- ggplot(data = df_all_male, aes(x = Var1,
                                                y = Var2,
                                                fill = Freq)) +
  geom_tile() + facet_grid(w~race_sex) +
  scale_fill_gradientn(name = "Beta Value", trans = "log", oob=squish_infinite,
                       breaks = c(0.007, 0.018, 0.050, 0.135, 0.368),
                       colors = c("white", "#FFF5F0", "#FEE0D2", "#FCBBA1",
                                  "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D",
                                  "#A50F15", "#67000D")) +
  scale_y_discrete(limits=rev) + scale_x_discrete() + theme_bw() +
  xlab("Age Group") + ylab("Age Group")
plot_all_male


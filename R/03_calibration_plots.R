library(dplyr)
library(ggplot2)
library(scales)
library(corrplot)
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

foi_transition_matrix_assort <- function(v_betas, w = w1) {
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

foi_transition_matrix_rp <- function(v_betas, w = w1) {
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

  beta <- rbind(cbind(beta1, beta5, beta6,  beta7),
                cbind(beta5, beta2, beta8,  beta9),
                cbind(beta6, beta8, beta3,  beta10),
                cbind(beta7, beta9, beta10, beta4))
  return(beta)
}

get_foi_plot <- function(l_m_samp, v_p_trans) {
  v_race_labels <- c(rep("Hispanic", length(groups)),
                     rep("NH White", length(groups)),
                     rep("NH Black", length(groups)),
                     rep("Other", length(groups)))
  v_age_labels <-  rep(groups, length(v_race))
  m_foi_hat <- matrix(nrow = 1000, ncol = 320)
  l_m_foi_hat <- list()
  for (j in 1:length(v_p_trans)) {
    m_samp <- l_m_samp[[j]]
    p_trans <- v_p_trans[j]

    for (i in seq.int(nrow(m_samp))) {
      if (ncol(m_samp) == 32) {
        tmp_waifw <- foi_transition_matrix_assort(m_samp[i, ])
      } else {
        tmp_waifw <- foi_transition_matrix_rp(m_samp[i, ])
      }
      tmp_foi <- (p_trans * tmp_waifw) %*% v_inf
      m_foi_hat[i, ] <- tmp_foi
    }
    # Observed
    df_foi_obs <- cbind(v_race_labels, v_age_labels, v_parameter$v_foi,
                        v_parameter$v_foi_lb, v_parameter$v_foi_ub,
                        rep("Observed", 320), rep(p_trans, 320))
    if (j == 1) {
      # Model Estimated
      df_foi_hat <- as.data.frame(m_foi_hat)
      foi_hat_lb     <- apply(df_foi_hat, 2, quantile, probs = 0.025)
      foi_hat_ub     <- apply(df_foi_hat, 2, quantile, probs = 0.975)
      foi_hat_median <- apply(df_foi_hat, 2, quantile, probs = 0.5)
      df_foi_hat <- cbind(v_race_labels, v_age_labels, foi_hat_median,
                          foi_hat_lb, foi_hat_ub, rep("Model Predicted", 320),
                          rep(p_trans, 320))

      df_foi_all <- rbind(df_foi_hat, df_foi_obs)
    }
    tmp_df_foi_hat <- as.data.frame(m_foi_hat)
    foi_hat_lb     <- apply(tmp_df_foi_hat, 2, quantile, probs = 0.025)
    foi_hat_ub     <- apply(tmp_df_foi_hat, 2, quantile, probs = 0.975)
    foi_hat_median <- apply(tmp_df_foi_hat, 2, quantile, probs = 0.5)
    tmp_df_foi_hat <- cbind(v_race_labels, v_age_labels, foi_hat_median,
                            foi_hat_lb, foi_hat_ub, rep("Model Predicted", 320),
                            rep(p_trans, 320))
    tmp_df_foi_all <- rbind(tmp_df_foi_hat, df_foi_obs)
    df_foi_all <- rbind(df_foi_all, tmp_df_foi_all)
  }

  # df_foi <- as.data.frame(rbind(cbind(c(rep("Hispanic", length(groups)),
  #                                       rep("NH White", length(groups)),
  #                                       rep("NH Black", length(groups)),
  #                                       rep("Other", length(groups))),
  #                                     rep(groups, length(v_race)),
  #                                     v_parameter$v_foi,
  #                                     v_parameter$v_foi_lb,
  #                                     v_parameter$v_foi_ub,
  #                                     rep("Observed", length(v_race)*length(groups))),
  #                               cbind(c(rep("Hispanic", length(groups)),
  #                                       rep("NH White", length(groups)),
  #                                       rep("NH Black", length(groups)),
  #                                       rep("Other", length(groups))),
  #                                     rep(groups, length(v_race)),
  #                                     foi_hat_median,
  #                                     foi_hat_lb,
  #                                     foi_hat_ub,
  #                                     rep("Model Predicted", length(v_race)*length(groups)))))
  colnames(df_foi_all) <- c("race", "age", "foi", "foi_lb",
                        "foi_ub", "model", "pr_trans")
  df_foi_all <- as.data.frame(df_foi_all)
  df_foi_all$foi       <- as.numeric(df_foi_all$foi)
  df_foi_all$foi_lb    <- as.numeric(df_foi_all$foi_lb)
  df_foi_all$foi_ub    <- as.numeric(df_foi_all$foi_ub)
  df_foi_all$age       <- as.numeric(df_foi_all$age)
  df_foi_all$pr_trans  <- factor(df_foi_all$pr_trans,
                                 levels = v_p_trans,
                                 labels = label_percent()(v_p_trans))
  df_foi_all$model     <- factor(df_foi_all$model,
                             levels = c("Observed", "Model Predicted"))
  print(head(df_foi_all))

  plot <- ggplot(data = df_foi_all) +
    facet_grid(race ~ pr_trans, scales = "free") +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    theme_bw() +
    geom_line(aes(x = age, y = foi, color = model), linewidth = 1.2) +
    geom_ribbon(aes(x = age, y = foi, ymax = foi_ub,
                    ymin = foi_lb, fill = model),
                alpha = 0.3, color = NA) +
    theme(legend.position = "bottom",
          title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 11),
          strip.text = element_text(size = 14),
          legend.title=element_blank()) +
    xlab("Age") + ylab("Force of Infection (FOI)") +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       sec.axis = sec_axis(~ . , name = "Probability of Transmission",
                                           breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = "Race/Ethnicity",
                                           breaks = NULL, labels = NULL))

  return(plot)
}

# Load other functions
source("R/01_model_inputs.R", echo = FALSE)
source("R/02_model_functions.R", echo = FALSE)

### Get parameter sets
groups <- c(1:80)
w1 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
               0, 2, 3, 4, 5, 6, 7, 8,
               0, 0, 3, 4, 5, 6, 7, 8,
               0, 0, 0, 4, 5, 6, 7, 8,
               0, 0, 0, 0, 5, 6, 7, 8,
               0, 0, 0, 0, 0, 6, 7, 8,
               0, 0, 0, 0, 0, 0, 7, 8,
               0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = TRUE)

v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")
v_break <- waifw_breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
v_race       <- c("Hispanic", "NH White", "NH Black", "Other")
v_race_names <- c("hisp", "white", "black", "other")

# Load demographic and prevalence variables=
v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race,
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

v_p <- seq(from = 0.005, to = 0.05, by = 0.005)

## Assortative
v_m_resamp_assort <- c("imis_resample_1_assort", "imis_resample_2_assort", "imis_resample_3_assort",
                       "imis_resample_4_assort", "imis_resample_5_assort", "imis_resample_6_assort",
                       "imis_resample_7_assort", "imis_resample_8_assort", "imis_resample_9_assort",
                       "imis_resample_10_assort")

for (i in seq_along(v_m_resamp_assort)) {
  assign(v_m_resamp_assort[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_assort[i], ".csv"))))
}

l_m_resamp_assort <- list(imis_resample_1_assort, imis_resample_2_assort, imis_resample_3_assort,
                       imis_resample_4_assort, imis_resample_5_assort, imis_resample_6_assort,
                       imis_resample_7_assort, imis_resample_8_assort, imis_resample_9_assort,
                       imis_resample_10_assort)

assort_foi <- get_foi_plot(l_m_resamp_assort, v_p) + ggtitle("Assortative")
assort_foi
ggsave(filename = "results/waifw/assort_foi_calibration.png",
       plot = assort_foi, width = 12, height = 8)

## Random Partnership
v_m_resamp_rp <- c("imis_resample_1_rp", "imis_resample_2_rp", "imis_resample_3_rp",
                   "imis_resample_4_rp", "imis_resample_5_rp", "imis_resample_6_rp",
                   "imis_resample_7_rp", "imis_resample_8_rp", "imis_resample_9_rp",
                   "imis_resample_10_rp")

for (i in seq_along(v_m_resamp_rp)) {
  assign(v_m_resamp_rp[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_rp[i], ".csv"))))
}

l_m_resamp_rp <- list(imis_resample_1_rp, imis_resample_2_rp, imis_resample_3_rp,
                          imis_resample_4_rp, imis_resample_5_rp, imis_resample_6_rp,
                          imis_resample_7_rp, imis_resample_8_rp, imis_resample_9_rp,
                          imis_resample_10_rp)

rp_foi <- get_foi_plot(l_m_resamp_rp, v_p) + ggtitle("Random Partnership")
rp_foi
ggsave(filename = "results/waifw/rp_foi_calibration.png",
       plot = rp_foi, width = 12, height = 8)

## Proportionate mixing
get_foi_plot_prop <- function(v_m_samp_assort, v_m_samp_rp,
                              v_m_samp_alphas, v_p_trans) {
  v_race_labels <- c(rep("Hispanic", length(groups)),
                     rep("NH White", length(groups)),
                     rep("NH Black", length(groups)),
                     rep("Other", length(groups)))
  v_age_labels <-  rep(groups, length(v_race))
  m_foi_hat <- matrix(nrow = 1000, ncol = 320)
  l_m_foi_hat <- list()

  for (j in 1:length(v_p_trans)) {
    m_samp_assort <- v_m_samp_assort[[j]]
    m_samp_rp     <- v_m_samp_rp[[j]]
    m_samp_alphas <- v_m_samp_alphas[[j]]
    p_trans       <- v_p_trans[j]

    for (i in seq.int(nrow(m_samp_assort))) {
      waifw_assort <- foi_transition_matrix_assort(m_samp_assort[i, ])
      waifw_rp     <- foi_transition_matrix_rp(m_samp_rp[i, ])

      v_alphas_within     <- as.numeric(m_samp_alphas[i, c(1, 6, 11, 16)])
      v_alphas_within_off <- as.numeric(1 - m_samp_alphas[i, c(1, 6, 11, 16)])
      v_alphas_between    <- as.numeric(m_samp_alphas[i, c(2:5, 7:10, 12:15)])
      m_within      <- alpha_matrix_within(v_alphas_within)
      m_within_off  <- alpha_matrix_within(v_alphas_within_off)
      m_between     <- alpha_matrix_between(v_alphas_between)
      waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
        (m_between * waifw_rp)
      tmp_foi <- (p_trans * waifw_p) %*% v_inf
      m_foi_hat[i, ] <- tmp_foi
    }

    # Observed
    df_foi_obs <- cbind(v_race_labels, v_age_labels, v_parameter$v_foi,
                        v_parameter$v_foi_lb, v_parameter$v_foi_ub,
                        rep("Observed", 320), rep(p_trans, 320))
    if (j == 1) {
      # Model Estimated
      df_foi_hat <- as.data.frame(m_foi_hat)
      foi_hat_lb     <- apply(df_foi_hat, 2, quantile, probs = 0.025)
      foi_hat_ub     <- apply(df_foi_hat, 2, quantile, probs = 0.975)
      foi_hat_median <- apply(df_foi_hat, 2, quantile, probs = 0.5)
      df_foi_hat <- cbind(v_race_labels, v_age_labels, foi_hat_median,
                          foi_hat_lb, foi_hat_ub, rep("Model Predicted", 320),
                          rep(p_trans, 320))

      df_foi_all <- rbind(df_foi_hat, df_foi_obs)
    }
    tmp_df_foi_hat <- as.data.frame(m_foi_hat)
    foi_hat_lb     <- apply(tmp_df_foi_hat, 2, quantile, probs = 0.025)
    foi_hat_ub     <- apply(tmp_df_foi_hat, 2, quantile, probs = 0.975)
    foi_hat_median <- apply(tmp_df_foi_hat, 2, quantile, probs = 0.5)
    tmp_df_foi_hat <- cbind(v_race_labels, v_age_labels, foi_hat_median,
                            foi_hat_lb, foi_hat_ub, rep("Model Predicted", 320),
                            rep(p_trans, 320))
    tmp_df_foi_all <- rbind(tmp_df_foi_hat, df_foi_obs)
    df_foi_all <- rbind(df_foi_all, tmp_df_foi_all)
  }

  colnames(df_foi_all) <- c("race", "age", "foi", "foi_lb",
                            "foi_ub", "model", "pr_trans")
  df_foi_all <- as.data.frame(df_foi_all)
  df_foi_all$foi       <- as.numeric(df_foi_all$foi)
  df_foi_all$foi_lb    <- as.numeric(df_foi_all$foi_lb)
  df_foi_all$foi_ub    <- as.numeric(df_foi_all$foi_ub)
  df_foi_all$age       <- as.numeric(df_foi_all$age)
  df_foi_all$pr_trans  <- factor(df_foi_all$pr_trans,
                                 levels = v_p_trans,
                                 labels = label_percent()(v_p_trans))
  df_foi_all$model     <- factor(df_foi_all$model,
                                 levels = c("Observed", "Model Predicted"))
  print(head(df_foi_all))

  plot <- ggplot(data = df_foi_all) +
    facet_grid(race ~ pr_trans, scales = "free") +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    theme_bw() +
    geom_line(aes(x = age, y = foi, color = model), linewidth = 1.2) +
    geom_ribbon(aes(x = age, y = foi, ymax = foi_ub,
                    ymin = foi_lb, fill = model),
                alpha = 0.3, color = NA) +
    theme(legend.position = "bottom",
          title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 11),
          strip.text = element_text(size = 14),
          legend.title=element_blank()) +
    xlab("Age") + ylab("Force of Infection (FOI)") +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80),
                       sec.axis = sec_axis(~ . , name = "Probability of Transmission",
                                           breaks = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = sec_axis(~ . , name = "Race/Ethnicity",
                                           breaks = NULL, labels = NULL))

  return(plot)
}


## Alphas
v_m_resamp_alphas <- c("imis_resample_1_alphas", "imis_resample_2_alphas",
                       "imis_resample_3_alphas", "imis_resample_4_alphas",
                       "imis_resample_5_alphas", "imis_resample_6_alphas",
                       "imis_resample_7_alphas", "imis_resample_8_alphas",
                       "imis_resample_9_alphas", "imis_resample_10_alphas")

for (i in seq_along(v_m_resamp_alphas)) {
  assign(v_m_resamp_alphas[i],
         as.matrix(read.csv(paste0("results/waifw/",
                                   v_m_resamp_alphas[i], ".csv"))))
}

l_m_resamp_alphas <- list(imis_resample_1_alphas, imis_resample_2_alphas,
                          imis_resample_3_alphas, imis_resample_4_alphas,
                          imis_resample_5_alphas, imis_resample_6_alphas,
                          imis_resample_7_alphas, imis_resample_8_alphas,
                          imis_resample_9_alphas, imis_resample_10_alphas)

alphas_foi <- get_foi_plot_prop(l_m_resamp_assort, l_m_resamp_rp,
                                l_m_resamp_alphas, v_p) +
  ggtitle("Proportional Mixing")
alphas_foi
ggsave(filename = "results/waifw/pm_foi_calibration.png",
       plot = alphas_foi, width = 12, height = 8)

## Alpha prior and posterior plots
plot_prior_posterior <- function(prior, posterior, num_rows) {
  df_prior <- as.data.frame(prior) %>%
    mutate(post = 0)
  df_posterior <- as.data.frame(posterior) %>%
    mutate(post = 1)
  df_all <- rbind(df_prior, df_posterior) %>%
    pivot_longer(!post, names_to = "parameter", values_to = "value")
  df_all$post <- factor(df_all$post, levels = c(0, 1), labels = c("Prior", "Posterior"))

  plot <- ggplot(df_all, aes(x = value, fill = post, group = post)) +
    geom_density(alpha = 0.4) +
    facet_wrap(vars(parameter), nrow = num_rows, scales = "free") +
    theme_bw() +
    xlab("Parameter Value") +
    ylab("Density") +
    scale_fill_viridis(discrete = TRUE, direction = -1) +
    scale_x_continuous(n.breaks = 4) +
    theme(legend.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6))
  return(plot)
}


v_p_text <- c("005", "01", "015", "02", "025", "03", "035", "04", "045", "05")
for (i in 1:10) {
  m_prior <- sample.prior(1000)
  tmp_plot <- plot_prior_posterior(m_prior, l_m_resamp_alphas[[i]], 4)
  ggsave(filename = paste0("results/waifw/plot_alphas_posterior_", v_p_text[i],".png"),
         plot = tmp_plot, width = 8, height = 8)
}

## Assortative prior and posterior plots
m_assort <- read.csv("results/waifw/assort_contacts_hat_ga.csv")
m_assort_se <- read.csv("results/waifw/assort_contacts_se_ga.csv")

v_names <- paste(rep(v_race, each = 8), v_age_group_names, sep = " ")
l_m_resamp_assort <- lapply(l_m_resamp_assort, function(x) {`colnames<-`(x, v_names)})

for (i in 1:10) {

  v_assort <- as.numeric(m_assort[i, ])
  v_assort_se <- as.numeric(m_assort_se[i, ])*1.5
  formals(sample.prior) <- alist(n.samp = , n.params = 32, v_mean = v_assort,
                                 v_se = v_assort_se)
  m_prior <- sample.prior(1000)
  colnames(m_prior) <- v_names

  tmp_plot <- plot_prior_posterior(m_prior, l_m_resamp_assort[[i]], num_rows = 4)
  ggsave(filename = paste0("results/waifw/plot_assort_posterior_", v_p_text[i],".png"),
         plot = tmp_plot, width = 12, height = 12)
}

plot_prior_posterior(m_prior, l_m_resamp_assort[[1]], num_rows = 4)

## Random Partnership
m_rp <- read.csv("results/waifw/rp_contacts_hat_ga.csv")
m_rp_se <- read.csv("results/waifw/rp_contacts_se_ga.csv")
v_race_comb_names <- c("Hisp. Hisp.", "White White", "Black Black", "Other Other",
                       "Hisp. White", "Hisp. Black", "Hisp. Other", "White Black",
                       "White Other", "Black Other")
v_names <- paste(rep(v_race_comb_names, each = 8), v_age_group_names, sep = " ")
l_m_resamp_rp <- lapply(l_m_resamp_rp, function(x) {`colnames<-`(x, v_names)})
for (i in 1:10) {
  v_rp <- as.numeric(m_rp[i, ])
  v_rp_se <- v_rp * 0.3
  formals(sample.prior) <- alist(n.samp = , n.params = 80, v_mean = v_rp,
                                 v_se = v_rp_se)
  m_prior <- sample.prior(1000)
  colnames(m_prior) <- v_names

  tmp_plot <- plot_prior_posterior(m_prior, l_m_resamp_rp[[i]], 10)
  ggsave(filename = paste0("results/waifw/plot_rp_posterior_", v_p_text[i],".png"),
         plot = tmp_plot, width = 12, height = 12)
}

## Plot overlay of all posteriors across all probs of transmission
library(viridis)
library(scales)
plot_posteriors <- function(l_posterior, num_rows, v_probs = v_p) {
  for (i in 1:length(v_probs)) {
    if (i == 1) {
      df_all <- as.data.frame(l_posterior[[i]]) %>%
        mutate(prob = v_probs[i])
    } else {
      tmp_df <- as.data.frame(l_posterior[[i]]) %>%
        mutate(prob = v_probs[i])
      df_all <- rbind(df_all, tmp_df)
    }
  }

  df_all <- df_all %>%
    pivot_longer(!prob, names_to = "parameter", values_to = "value")
  df_all$prob <- factor(df_all$prob, levels = v_probs,
                        labels = c("0.5%", "1%", "1.5%", "2%", "2.5%",
                                   "3%", "3.5%", "4%", "4.5%", "5%"))
  df_all$parameter <- factor(df_all$parameter,
                             levels = c("a_11", "a_12", "a_13", "a_14",
                                        "a_21", "a_22", "a_23", "a_24",
                                        "a_31", "a_32", "a_33", "a_34",
                                        "a_41", "a_42", "a_43", "a_44"),
                             labels = c("alpha[HH]", "alpha[HW]", "alpha[HB]", "alpha[HO]",
                                        "alpha[WH]", "alpha[WW]", "alpha[WB]", "alpha[WO]",
                                        "alpha[BH]", "alpha[BW]", "alpha[BB]", "alpha[BO]",
                                        "alpha[OH]", "alpha[OW]", "alpha[OB]", "alpha[OO]"))

  plot <- ggplot(df_all, aes(x = value, fill = prob, group = prob)) +
    geom_density(alpha = 0.4) +
    facet_wrap(vars(parameter), nrow = num_rows, scales = "free",
               labeller = label_parsed) +
    theme_bw() +
    xlab("Parameter Value") +
    ylab("Density") +
    scale_fill_viridis(discrete = TRUE, direction = -1) +
    scale_x_continuous(n.breaks = 4) +
    guides(fill = guide_legend(title = "Probability of Transmission",
                               position = "bottom", nrow = 1)) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8),
          strip.text.x = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10))
  return(plot)
}

plot_post_alphas <- plot_posteriors(l_m_resamp_alphas, num_rows = 4)
ggsave(filename = "results/waifw/alpha_posterior_all.png",
       plot = plot_post_alphas, width = 10, height = 10)

plot_post_assort <- plot_posteriors(l_m_resamp_assort, num_rows = 4)
ggsave(filename = "results/waifw/assort_posterior_all.png",
       plot = plot_post_assort, width = 12, height = 10)

plot_post_rp     <- plot_posteriors(l_m_resamp_rp, num_rows = 10)
ggsave(filename = "results/waifw/rp_posterior_all.png",
       plot = plot_post_rp, width = 15, height = 10)


### Correlation matrix by alpha parameter across probs
plot_alpha_cor <- function(df_alpha, alpha_index) {
  df_all <- as.data.frame(df_alpha) %>%
    mutate(run = row_number()) %>%
    pivot_longer(!run, names_to = "parameter", values_to = "value")
  df_all$parameter <- factor(df_all$parameter,
                           levels = c("a_11", "a_12", "a_13", "a_14",
                                      "a_21", "a_22", "a_23", "a_24",
                                      "a_31", "a_32", "a_33", "a_34",
                                      "a_41", "a_42", "a_43", "a_44"),
                           labels = c("alpha[HH]", "alpha[HW]", "alpha[HB]", "alpha[HO]",
                                      "alpha[WH]", "alpha[WW]", "alpha[WB]", "alpha[WO]",
                                      "alpha[BH]", "alpha[BW]", "alpha[BB]", "alpha[BO]",
                                      "alpha[OH]", "alpha[OW]", "alpha[OB]", "alpha[OO]"))
  df_all <- df_all %>%
    pivot_wider(names_from = "parameter", values_from = "value") %>%
    select(-run)
  # p <- GGally::ggcorr(df_all, parse = T, label = T) +
  #   theme_bw()
  # p[["geom_params"]]$parse <- TRUE
  # print(p)
  # corrplot(cor(df_all), xlab = c(expression(alpha[HH]), expression(alpha[HW]), expression(alpha[HB]), expression(alpha[HO]),
  #                                expression(alpha[WH]), expression(alpha[WW]), expression(alpha[WB]), expression(alpha[WO]),
  #                                expression(alpha[BH]), expression(alpha[BW]), expression(alpha[BB]), expression(alpha[BO]),
  #                                expression(alpha[OH]), expression(alpha[OW]), expression(alpha[OB]), expression(alpha[OO])))
  M <- cor(df_all)
  colnames(M) <- c("alpha[HH]", "alpha[HW]", "alpha[HB]", "alpha[HO]",
                   "alpha[WH]", "alpha[WW]", "alpha[WB]", "alpha[WO]",
                   "alpha[BH]", "alpha[BW]", "alpha[BB]", "alpha[BO]",
                   "alpha[OH]", "alpha[OW]", "alpha[OB]", "alpha[OO]")
  rownames(M) <- c("alpha[HH]", "alpha[HW]", "alpha[HB]", "alpha[HO]",
                   "alpha[WH]", "alpha[WW]", "alpha[WB]", "alpha[WO]",
                   "alpha[BH]", "alpha[BW]", "alpha[BB]", "alpha[BO]",
                   "alpha[OH]", "alpha[OW]", "alpha[OB]", "alpha[OO]")
  p <- corrplot(M)
  print(p)
  list(text = text, plot = plot)
}

for (i in 1:10) {
  tmp_plot <- plot_alpha_cor(l_m_resamp_alphas[[i]], i)
  ggsave(filename = paste0("results/waifw/alpha_cor_", i, ".png"),
         plot = tmp_plot, width = 12, height = 12)
}

test_p <- plot_alpha_cor(l_m_resamp_alphas[[1]], 1)

M <- cor(l_m_resamp_alphas[[1]])
colnames(M) <- c(":alpha[HH]", ":alpha[HW]", ":alpha[HB]", ":alpha[HO]",
                 ":alpha[WH]", ":alpha[WW]", ":alpha[WB]", ":alpha[WO]",
                 ":alpha[BH]", ":alpha[BW]", ":alpha[BB]", ":alpha[BO]",
                 ":alpha[OH]", ":alpha[OW]", ":alpha[OB]", ":alpha[OO]")
rownames(M) <- c(":alpha[HH]", ":alpha[HW]", ":alpha[HB]", ":alpha[HO]",
                 ":alpha[WH]", ":alpha[WW]", ":alpha[WB]", ":alpha[WO]",
                 ":alpha[BH]", ":alpha[BW]", ":alpha[BB]", ":alpha[BO]",
                 ":alpha[OH]", ":alpha[OW]", ":alpha[OB]", ":alpha[OO]")
pdf("results/waifw/alpha_cor_1.pdf", height = 10, width = 10)
corrplot(M, method = "color", type = "upper", tl.srt = 0, diag = F,
         addCoef.col = 'black', tl.col = "red", mar = c(1, 1, 1, 1),
         tl.offset = 0.8, tl.cex = 1.3)
dev.off()

png("results/waifw/alpha_cor_1.png", height = 800, width = 800)
corrplot(M, method = "color", type = "upper", tl.srt = 0, diag = F,
         addCoef.col = 'black', tl.col = "red", mar = c(1, 1, 1, 1),
         tl.offset = 0.8, tl.cex = 1.3)
dev.off()

### Correlation across all parameters
library(corrplot)

df_assort <- do.call(cbind, l_m_resamp_assort)
df_rp     <- do.call(cbind, l_m_resamp_rp)
df_alphas <- do.call(cbind, l_m_resamp_alphas)
df_all_cor <- cor(cbind(df_assort, df_rp, df_alphas))

png("results/waifw/plot_cor_posterior_all.png")
corrplot(df_all_cor, tl.pos = 'n', diag = FALSE,
         addgrid.col = NA, is.corr = FALSE)
dev.off()

png("results/waifw/plot_cor_posterior_assort.png")
corrplot(cor(df_assort), tl.pos = 'n', diag = FALSE,
         addgrid.col = NA, is.corr = FALSE)
dev.off()

png("results/waifw/plot_cor_posterior_rp.png")
corrplot(cor(df_rp), tl.pos = 'n', diag = FALSE,
         addgrid.col = NA, is.corr = FALSE)
dev.off()

for (i in 1:10) {
  tmp_df_all_cor <- cor(cbind(l_m_resamp_assort[[i]], l_m_resamp_rp[[i]],
                              l_m_resamp_alphas[[i]]))
  png(paste0("results/waifw/plot_cor_posterior_all_", v_p_text[i],".png"))
  corrplot(tmp_df_all_cor, tl.pos = 'n', diag = FALSE,
           addgrid.col = NA, is.corr = FALSE)
  dev.off()

  png(paste0("results/waifw/plot_cor_posterior_assort_", v_p_text[i],".png"))
  corrplot(cor(l_m_resamp_assort[[i]]), tl.pos = 'n', diag = FALSE,
           addgrid.col = NA, is.corr = FALSE)
  dev.off()

  png(paste0("results/waifw/plot_cor_posterior_rp_", v_p_text[i],".png"))
  corrplot(cor(l_m_resamp_rp[[i]]), tl.pos = 'n', diag = FALSE,
           addgrid.col = NA, is.corr = FALSE)
  dev.off()
}

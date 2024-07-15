library(scales)
library(ggplot2)
library(Matrix)
library(ggsci)
library(stringr)
library(dplyr)

## Variables needed for calibration
# WAIFW matrices
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

v_waifw_structure <- list(w1, w2, w3)

## WAIFW parameters
waifw_breaks      <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
beta0             <- c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
l_beta0           <- list(c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                          c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                          c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
                          c(0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.05, 0.01))
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
  v_demo <- get_demographic_vars(spec_groups = groups,
                                 race = race)
  v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                     race = race)
  v_inf <- v_prev_vars$prevalence * v_demo$v_age_prop

  x <- 0
  for (betas in l_beta0) {
    x <- x + 1
    tmp <- calibrate_betas(n_age_groups = n_age_groups,
                           n_betas = n_betas,
                           v_waifw_structure = v_waifw_structure,
                           beta0 = betas)
    assign(paste0("calibrate_", race_text, "_", x), tmp)
  }
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
  mutate(race = "NH White")


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
                       # breaks = c(0.0025, 0.0498, 1.000),  # nolint
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

ggsave(filename = "results/WAIFW_calibration.png",
       plot = plot_all, width = 10, height = 4)

## Plot FOI
v_prev_vars_white <- get_prevalence_vars(spec_groups = groups,
                                         race = "NH White")
v_prev_vars_black <- get_prevalence_vars(spec_groups = groups,
                                         race = "NH Black")
v_prev_vars_hisp  <- get_prevalence_vars(spec_groups = groups,
                                         race = "Hispanic")
v_prev_vars_other <- get_prevalence_vars(spec_groups = groups,
                                         race = "Other")
v_demo_white <- get_demographic_vars(spec_groups = groups, race = "NH White")
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

colnames(df_foi) <- c("foi_est", "foi_est_lb", "foi_est_ub", "race",
                      "age", "foi_model")

df_foi$foi_est <- as.numeric(df_foi$foi_est)
df_foi$foi_est_lb <- as.numeric(df_foi$foi_est_lb)
df_foi$foi_est_ub <- as.numeric(df_foi$foi_est_ub)
df_foi$foi_model <- as.numeric(df_foi$foi_model)
df_foi$age <- as.numeric(df_foi$age)
df_foi$race2 <- df_foi$race


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
              alpha = 0.3, color = NA) +
  geom_point(aes(x = age, y = foi_model, color = race, alpha = "Model Predicted")) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  xlab("Age") + ylab("Force of Infection (FOI)") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  # ylim(0, 0.065) +
  scale_alpha_manual(name = NULL,
                     values = c(1, 1),
                     breaks = c("Model Predicted", "Observed"),
                     guide = guide_legend(override.aes =
                                            list(linetype = c(0, 1),
                                                 shape = c(16, NA),
                                                 color = "black")))
plot_foi_valid

ggsave(filename = "results/FOI_validation_free_scale.pdf",
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
waifw_names_white <- paste0("NH White - ",
                            paste0("$", paste("w", 1:3, sep = "_")), "$")
waifw_names_other <- paste0("Other - ",
                            paste0("$", paste("w", 1:3, sep = "_")), "$")
waifw_names_all   <- c(waifw_names_hisp,
                       waifw_names_black,
                       waifw_names_white,
                       waifw_names_other)

texreg(all_tex, digits = 3, stars = numeric(0), booktabs = TRUE,
       custom.model.names = waifw_names_all)

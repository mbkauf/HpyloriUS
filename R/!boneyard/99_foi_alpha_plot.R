v_foi_white <- c()
v_foi_black <- c()
n <- 0
for (i in seq(0, 1, by = 0.1)) {
  n <- n + 1
  alpha <- i
  beta1 <- ((1 - alpha) *  beta_white_r) + (alpha * beta_white_a)
  beta2 <- ((1 - alpha) *  beta_black_r) + (alpha * beta_black_a)
  beta3 <- ((1 - alpha) *  beta_mix_r)

  beta <- rbind(cbind(beta1, beta3),
                cbind(beta3, beta2))

  foi_hat <- beta %*% v_inf

  v_foi_white <- c(v_foi_white, foi_hat[1:80,])
  v_foi_black <- c(v_foi_black, foi_hat[81:160,])
}

v_alpha <- c(rep(0, 80), rep(0.1, 80), rep(0.2, 80), rep(0.3, 80), rep(0.4, 80),
             rep(0.5, 80), rep(0.6, 80), rep(0.7, 80), rep(0.8, 80),
             rep(0.9, 80), rep(1, 80))

df_foi <- as.data.frame(rbind(cbind(rep(v_prev_vars_white$foi, 11),
                                    rep(v_prev_vars_white$foi_lb, 11),
                                    rep(v_prev_vars_white$foi_ub, 11),
                                    rep("NH White", length(v_foi_white)),
                                    rep(groups, 11),
                                    v_foi_white,
                                    v_alpha),
                              cbind(rep(v_prev_vars_black$foi, 11),
                                    rep(v_prev_vars_black$foi_lb, 11),
                                    rep(v_prev_vars_black$foi_ub, 11),
                                    rep("NH Black", length(v_foi_black)),
                                    rep(groups, 11),
                                    v_foi_black,
                                    v_alpha)))

colnames(df_foi) <- c("foi_est", "foi_est_lb", "foi_est_ub", "race",
                      "age", "foi_model", "alpha")


df_foi$foi_est <- as.numeric(df_foi$foi_est)
df_foi$foi_est_lb <- as.numeric(df_foi$foi_est_lb)
df_foi$foi_est_ub <- as.numeric(df_foi$foi_est_ub)
df_foi$foi_model <- as.numeric(df_foi$foi_model)
df_foi$age <- as.numeric(df_foi$age)
df_foi$race2 <- df_foi$race
df_foi$alpha <- factor(df_foi$alpha)


plot_foi_mix <- ggplot(data = df_foi) +
  facet_wrap(. ~ race, scales = "free") +
  scale_color_nejm(guide = "none") +
  scale_fill_nejm(guide = "none") +
  theme_bw() +
  geom_line(aes(x = age, y = foi_est, color = race, alpha = "Observed")) +
  geom_ribbon(aes(x = age, y = foi_est, ymax = foi_est_ub,
                  ymin = foi_est_lb, fill = race),
              alpha = 0.3, color = NA) +
  geom_point(aes(x = age, y = foi_model, color = race, shape = alpha,
                 alpha = "Model Predicted")) +
  theme(legend.position = "right",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 14)) +
  xlab("Age") + ylab("FOI") +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_shape_manual(values=1:nlevels(df_foi$alpha)) +
  scale_alpha_manual(name = NULL,
                     values = c(1, 1),
                     breaks = c("Model Predicted", "Observed"),
                     guide = guide_legend(override.aes =
                                            list(linetype = c(0, 1),
                                                 shape = c(16, NA),
                                                 color = "black")))
plot_foi_mix
ggsave(filename = paste0("results/FOI_black_white_alpha_range.png"),
       plot = plot_foi_mix)

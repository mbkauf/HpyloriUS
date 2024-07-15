get_foi_plot <- function(waifw_assort, waifw_rp, v_alphas, v_inf, p_trans) {
  v_alphas_within     <- v_alphas[c(1, 6, 11, 16)]
  v_alphas_within_off <- 1 - v_alphas[c(1, 6, 11, 16)]
  v_alphas_between    <- v_alphas[c(2:5, 7:10, 12:15)]
  m_within      <- alpha_matrix_within(v_alphas_within)
  m_within_off  <- alpha_matrix_within(v_alphas_within_off)
  m_between     <- alpha_matrix_between(v_alphas_between)
  
  waifw_p  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
    (m_between * waifw_rp)
  foi_hat <- (p_trans * waifw_p) %*% v_inf
  
  df_foi <- as.data.frame(cbind(v_parameter$v_foi,
                                v_parameter$v_foi_lb,
                                v_parameter$v_foi_ub,
                                c(rep("Hispanic", length(groups)),
                                  rep("NH White", length(groups)),
                                  rep("NH Black", length(groups)),
                                  rep("Other", length(groups))),
                                rep(groups, length(4)),
                                foi_hat))
  
  colnames(df_foi) <- c("foi_est", "foi_est_lb", "foi_est_ub", "race",
                        "age", "foi_model")
  
  df_foi$foi_est <- as.numeric(df_foi$foi_est)
  df_foi$foi_est_lb <- as.numeric(df_foi$foi_est_lb)
  df_foi$foi_est_ub <- as.numeric(df_foi$foi_est_ub)
  df_foi$foi_model <- as.numeric(df_foi$foi_model)
  df_foi$age <- as.numeric(df_foi$age)
  df_foi$race2 <- df_foi$race
  
  plot <- ggplot(data = df_foi) +
    facet_wrap(. ~ race, scales = "free") +
    ggsci::scale_color_nejm(guide = "none") +
    ggsci::scale_fill_nejm(guide = "none") +
    theme_bw() +
    geom_line(aes(x = age, y = foi_est, color = race, alpha = "Observed")) +
    geom_ribbon(aes(x = age, y = foi_est, ymax = foi_est_ub,
                    ymin = foi_est_lb, fill = race),
                alpha = 0.3, color = NA) +
    geom_point(aes(x = age, y = foi_model, color = race,
                   alpha = "Model Predicted")) +
    theme(legend.position = "right",
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.text = element_text(size = 11),
          strip.text = element_text(size = 14)) +
    xlab("Age") + ylab("Force of Infection (FOI)") +
    scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
    scale_alpha_manual(name = NULL,
                       values = c(1, 1),
                       breaks = c("Model Predicted", "Observed"),
                       guide = guide_legend(override.aes =
                                              list(linetype = c(0, 1),
                                                   shape = c(16, NA),
                                                   color = "black")))
  
  return(plot)
}

for (i in 1:length(v_p)) {
  tmp_plot <- get_foi_plot(waifw_assort = l_m_contacts_assort[[i]],
                           waifw_rp = l_m_contacts_rp[[i]],
                           v_alphas = results_alphas$m_alphas_hat[i,],
                           v_inf = v_inf_rp,
                           p_trans = v_p[i])
  assign(paste0("foi_plot_", i), tmp_plot)
}


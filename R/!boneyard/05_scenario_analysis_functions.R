get_sis_treat_scenario <- function(p_trt_eff = NULL,
                                   p_test_sens = NULL,
                                   p_test_spec = NULL,
                                   v_betas,
                                   v_race = c("Hispanic",
                                              "NH White",
                                              "NH Black",
                                              "Other")) {
  require(dplyr)
  require(rootSolve)
  v_results <- list()
  i <- 1
  for (x in v_race) {
    if (x == "Hispanic") {
      beta <- v_betas$m_beta_hispanic
    } else if (x == "NH White") {
      beta <- v_betas$m_beta_white
    } else if (x == "NH Black") {
      beta <- v_betas$m_beta_black
    } else if (x == "Other") {
      beta <- v_betas$m_beta_other
    }

    # Load prevalence and demographic variables
    v_prev <- get_prevalence_vars(spec_groups = groups,
                                  race = x)
    v_demo <- get_demographic_vars(spec_groups = groups,
                                   race = x)
    v_parameter <- load_sis_model_params(waifw = beta,
                                         demography_vars = v_demo,
                                         prevalence_vars = v_prev,
                                         ages = groups,
                                         trt_year = 50000000,
                                         end_t = 50000000)

    # Run burn-in period and store starting states
    burn_results <- runsteady(y = v_parameter$v_state,
                              times = c(0,1E5), func = sis_model,
                              parms = v_parameter)
    v_burn_state <- burn_results$y

    # Scenario 1: No treat
    v_parameter <- load_sis_model_params(waifw = beta,
                                         demography_vars = v_demo,
                                         prevalence_vars = v_prev,
                                         ages = groups,
                                         trt_year = 1000,
                                         end_t = 50,
                                         burn_state = TRUE,
                                         v_init_state = v_burn_state)

    tmp_results <- get_sis_model_results(v_params = v_parameter)
    v_results[[i]] <- tmp_results
    i <- i + 1

    # Scenario 2: One-time treatment, 99% effective
    v_parameter <- load_sis_model_params(waifw = beta,
                                         demography_vars = v_demo,
                                         prevalence_vars = v_prev,
                                         ages = groups,
                                         trt_year = 1,
                                         trt_eff = 0.99,
                                         end_t = 50,
                                         burn_state = TRUE,
                                         v_init_state = v_burn_state)
    tmp_results <- get_sis_model_results(v_params = v_parameter)
    v_results[[i]] <- tmp_results
    i <- i + 1

    # Scenario 3: One-time treatment, 80% effective
    v_parameter <- load_sis_model_params(waifw = beta,
                                         demography_vars = v_demo,
                                         prevalence_vars = v_prev,
                                         ages = groups,
                                         trt_year = 1,
                                         trt_eff = 0.8,
                                         end_t = 50,
                                         burn_state = TRUE,
                                         v_init_state = v_burn_state)
    tmp_results <- get_sis_model_results(v_params = v_parameter)
    v_results[[i]] <- tmp_results
    i <- i + 1
  }

  df_all <- rbind(
    v_results[[1]],
    v_results[[2]],
    v_results[[3]],
    v_results[[4]],
    v_results[[5]],
    v_results[[6]],
    v_results[[7]],
    v_results[[8]],
    v_results[[9]],
    v_results[[10]],
    v_results[[11]],
    v_results[[12]]
  )

  v_race_names <- c(rep("Hispanic", nrow(v_results[[1]])),
                    rep("Hispanic", nrow(v_results[[2]])),
                    rep("Hispanic", nrow(v_results[[3]])),
                    rep("NH White", nrow(v_results[[4]])),
                    rep("NH White", nrow(v_results[[5]])),
                    rep("NH White", nrow(v_results[[6]])),
                    rep("NH Black", nrow(v_results[[7]])),
                    rep("NH Black", nrow(v_results[[8]])),
                    rep("NH Black", nrow(v_results[[9]])),
                    rep("Other", nrow(v_results[[10]])),
                    rep("Other", nrow(v_results[[11]])),
                    rep("Other", nrow(v_results[[12]])))
  v_scenario <- c(rep(1, nrow(v_results[[1]])),
                  rep(2, nrow(v_results[[2]])),
                  rep(3, nrow(v_results[[3]])),
                  rep(1, nrow(v_results[[4]])),
                  rep(2, nrow(v_results[[5]])),
                  rep(3, nrow(v_results[[6]])),
                  rep(1, nrow(v_results[[7]])),
                  rep(2, nrow(v_results[[8]])),
                  rep(3, nrow(v_results[[9]])),
                  rep(1, nrow(v_results[[10]])),
                  rep(2, nrow(v_results[[11]])),
                  rep(3, nrow(v_results[[12]])))

  df_all <- cbind(df_all, v_race_names, v_scenario) %>%
    dplyr::select(c(time, s0, i0, s1, i1, v_race_names, v_scenario)) %>%
    dplyr::mutate(I = i0 + i1) %>%
    dplyr::mutate(scen = ifelse(v_scenario > 1, 1, 0)) %>%
    dplyr::mutate(scen = factor(scen, levels = c(0, 1))) %>%
    dplyr::mutate(v_scenario = factor(v_scenario, levels = c(1, 2, 3),
                                      labels = c("No Treatment",
                                                 "One-Time Treatment: 99% Effectiveness",
                                                 "One-Time Treatment: 80% Effectiveness")))


  return(df_all)
}


graph_results <- function(df_all) {
  require(ggplot2)

  l_scenario <- split(df_all, df_all$v_scenario)

  plot_1 <- ggplot(data = l_scenario[[1]], aes(color = v_race_names)) +
    geom_line(mapping = aes(x = time, y = I)) +
    ylim(0, 0.6) +
    xlab("Year") + ylab("Proportion Infected") +
    guides(color = guide_legend(title = "Race/Ethnicity")) +
    ggtitle("No Treatment") + theme(legend.position = "none")

  plot_2 <- ggplot(data = l_scenario[[2]], aes(color = v_race_names)) +
    geom_line(mapping = aes(x = time, y = I)) +
    ylim(0, 0.6) +
    xlab("Year") + ylab("Proportion Infected") +
    guides(color = guide_legend(title = "Race/Ethnicity")) +
    ggtitle("One-Time Treatment - 99% Effectiveness") + theme(legend.position = "none")

  plot_3 <- ggplot(data = l_scenario[[3]], aes(color = v_race_names)) +
    geom_line(mapping = aes(x = time, y = I)) +
    ylim(0, 0.6) +
    xlab("Year") + ylab("Proportion Infected") +
    guides(color = guide_legend(title = "Race/Ethnicity")) +
    ggtitle("One-Time Treatment - 80% Effectiveness")

  return(list(plot_1, plot_2, plot_3))
}

graph_results_grid <- function(df_all) {
  require(ggplot2)
  require(ggsci)

  df_notrt <- df_all %>%
    filter(scen == 0) %>%
    mutate(I_notrt = I) %>%
    dplyr::select(time, I_notrt, v_race_names) %>%
    mutate(scenario_name = "Status Quo")

  df_scen  <- df_all %>%
    filter(scen == 1) %>%
    mutate(scenario_name = "Treatment")

  df_plot  <- merge(df_scen, df_notrt, by = "v_race_names")

  plot_1 <- ggplot(data = df_scen, aes(color = v_race_names)) +
    geom_line(data = df_notrt,
              mapping = aes(x = time, y = I_notrt, size = scenario_name, linetype = scenario_name)) +
    geom_line(mapping = aes(x = time, y = I, size = scenario_name, linetype = scenario_name)) +
    scale_size_manual(values = c(0.5, 1), guide = 'none') +
    facet_wrap(. ~ v_scenario,
               labeller = labeller(v_scenario = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
    scale_linetype_manual(values = c(2, 1),
                          guide = guide_legend(override.aes=list(linetype=c(2,1),
                                                                 lwd=c(0.5,1)))) +
    scale_y_continuous(limits = c(0, 0.6) ,
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    xlab("Year") + ylab("Proportion Infected") +
    guides(color = guide_legend(title = "Race/Ethnicity"),
           linetype = guide_legend(title = element_blank())) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))

  return(plot_1)
}

beta_hisp <- get_transmission_matrix(betas = calibrate_hispanic_1$m_beta_hat[2, ],
                                     breaks = waifw_breaks,
                                     w = w2,
                                     group_names = groups)
beta_white <- get_transmission_matrix(betas = calibrate_white_1$m_beta_hat[2, ],
                                      breaks = waifw_breaks,
                                      w = w2,
                                      group_names = groups)
beta_black <- get_transmission_matrix(betas = calibrate_black_1$m_beta_hat[2, ],
                                      breaks = waifw_breaks,
                                      w = w2,
                                      group_names = groups)
beta_other <- get_transmission_matrix(betas = calibrate_other_1$m_beta_hat[2, ],
                                      breaks = waifw_breaks,
                                      w = w2,
                                      group_names = groups)

v_beta_race <- list(m_beta_hispanic = beta_hisp, m_beta_white = beta_white,
                    m_beta_black = beta_black, m_beta_other = beta_other)

scenario_results <- get_sis_treat_scenario(v_betas = v_beta_race)
# plot_results <- graph_results(scenario_results)
# ggsave(filename = "results/test_no_trt.pdf", plot = plot_results[[1]],
#        width = 4, height = 4)
# ggsave(filename = "results/test_trt_99.pdf", plot = plot_results[[2]],
#        width = 4, height = 4)
# ggsave(filename = "results/test_trt_80.pdf", plot = plot_results[[3]],
#        width = 6, height = 4)

plot_results_grid <- graph_results_grid(scenario_results)
plot_results_grid
ggsave(filename = "results/test_all_scenarios.pdf", plot = plot_results_grid,
       width = 10, height = 4)


## Create bar graph for relative impact of interventions
pct_impact <- function(df_all) {
  require(ggplot2)
  require(dplyr)
  require(ggsci)

  df_all <- df_all %>%
    filter(time == 0 | time == 50) %>%
    dplyr::mutate(I = i0 + i1)

  df_wide <- dcast(setDT(df_all), v_scenario +  v_race_names ~ time,
                   value.var = "I")
  colnames(df_wide) <- c("v_scenario", "v_race", "t_0", "t_50")

  df_wide <- df_wide %>%
    filter(v_scenario != "No Treatment") %>%
    dplyr::mutate(pct_diff = abs(t_50 - t_0) / t_0) %>%
    dplyr::mutate(net_diff = t_0 - t_50 )

  plot <- ggplot(data = df_wide, aes(x = v_race, y = net_diff, fill = v_race)) +
    geom_bar(stat = "identity") +
    facet_wrap(.~ v_scenario,
               labeller = labeller(v_scenario = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Race/Ethnicity") + ylab("Net Reduction in Prevalence (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")

  return(plot)
}

pct_plot <- pct_impact(scenario_results)

library(grid)
library(gridExtra)

poster <- grid.arrange(arrangeGrob(plot_results_grid, pct_plot, ncol = 2))
ggsave(filename = "results/scenario_fig.pdf", plot = poster,
       width = 12, height = 4)



get_sis_treat_scenario_all <- function(p_trt_eff = NULL,
                                       p_test_sens = NULL,
                                       p_test_spec = NULL,
                                       ages = groups,
                                       m_waifw,
                                       m_waifw_burn,
                                       v_race = c("Hispanic",
                                                  "NH White",
                                                  "NH Black",
                                                  "Other"),
                                       v_race_names = c("hisp", "black",
                                                        "white", "other")) {
  require(dplyr)
  require(rootSolve)
  v_results <- list()

  v_parameter <- load_sis_model_params_all(
    v_race = v_race,
    waifw = m_waifw_burn,
    trt_year = 50000000,
    end_t = 1000,
    ages = ages)

  model_results_sis_all <- get_sis_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  # Run burn-in period and store starting states
  burn_results <- runsteady(y = v_parameter$v_state,
                            times = c(0,1E5), func = sis_model_all,
                            parms = v_parameter)
  v_burn_state <- burn_results$y

  # Scenario 1: No treat
  v_parameter <- load_sis_model_params_all(
    v_race = v_race,
    waifw = m_waifw,
    trt_year = 100001,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = ages)

  tmp_results <- get_sis_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[1]] <- tmp_results

  # Scenario 2: One-time treatment, 99% effective
  v_parameter <- load_sis_model_params_all(
    v_race = v_race,
    waifw = m_waifw,
    trt_year = 1,
    trt_eff = 0.99,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = ages)

  tmp_results <- get_sis_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[2]] <- tmp_results

  # Scenario 3: One-time treatment, 80% effective
  v_parameter <- load_sis_model_params_all(
    v_race = v_race,
    waifw = m_waifw,
    trt_year = 1,
    trt_eff = 0.8,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = ages)

  tmp_results <- get_sis_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[3]] <- tmp_results

  # Combine Results
  df_all <- rbind(
    v_results[[1]],
    v_results[[2]],
    v_results[[3]]
  )

  m_hisp <- df_all %>%
    dplyr::select(c(time, s0_hisp, i0_hisp, s1_hisp, i1_hisp)) %>%
    dplyr::mutate(
      race = "Hispanic",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51))
    ) %>%
    dplyr::rename(
      s0 = s0_hisp,
      i0 = i0_hisp,
      s1 = s1_hisp,
      i1 = i1_hisp
    )

  m_black <- df_all %>%
    dplyr::select(c(time, s0_black, i0_black, s1_black, i1_black)) %>%
    dplyr::mutate(
      race = "NH Black",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51))
    ) %>%
    dplyr::rename(
      s0 = s0_black,
      i0 = i0_black,
      s1 = s1_black,
      i1 = i1_black
    )
  m_white <- df_all %>%
    dplyr::select(c(time, s0_white, i0_white, s1_white, i1_white)) %>%
    dplyr::mutate(
      race = "NH White",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51))
    ) %>%
    dplyr::rename(
      s0 = s0_white,
      i0 = i0_white,
      s1 = s1_white,
      i1 = i1_white
    )
  m_other <- df_all %>%
    dplyr::select(c(time, s0_other, i0_other, s1_other, i1_other)) %>%
    dplyr::mutate(
      race = "Other",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51))
    ) %>%
    dplyr::rename(
      s0 = s0_other,
      i0 = i0_other,
      s1 = s1_other,
      i1 = i1_other
    )

  df_all_long <- rbind(m_hisp, m_black, m_white, m_other) %>%
    # df_all_long <- rbind(m_black, m_white) %>%
    dplyr::mutate(I = i0 + i1) %>%
    dplyr::mutate(scen = ifelse(scenario > 1, 1, 0)) %>%
    dplyr::mutate(scen = factor(scen, levels = c(0, 1))) %>%
    dplyr::mutate(scenario = factor(scenario, levels = c(1, 2, 3),
                                    labels = c("No Treatment",
                                               "One-Time Treatment: 99% Effectiveness",
                                               "One-Time Treatment: 80% Effectiveness")))

  return(df_all_long)
}

graph_results <- function(df_all) {
  require(ggplot2)

  l_scenario <- split(df_all, df_all$scenario)

  plot_1 <- ggplot(data = l_scenario[[1]], aes(color = race)) +
    geom_line(mapping = aes(x = time, y = I)) +
    ylim(0, 0.6) +
    xlab("Year") + ylab("Proportion Infected") +
    guides(color = guide_legend(title = "Race/Ethnicity")) +
    ggtitle("No Treatment") + theme(legend.position = "none")

  plot_2 <- ggplot(data = l_scenario[[2]], aes(color = race)) +
    geom_line(mapping = aes(x = time, y = I)) +
    ylim(0, 0.6) +
    xlab("Year") + ylab("Proportion Infected") +
    guides(color = guide_legend(title = "Race/Ethnicity")) +
    ggtitle("One-Time Treatment - 99% Effectiveness") + theme(legend.position = "none")

  plot_3 <- ggplot(data = l_scenario[[3]], aes(color = race)) +
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
  require(ggpubr)

  df_notrt <- df_all %>%
    filter(scen == 0) %>%
    mutate(I_notrt = I) %>%
    dplyr::select(time, I_notrt, race) %>%
    mutate(scenario_name = "Status Quo")

  df_scen  <- df_all %>%
    filter(scen == 1) %>%
    mutate(scenario_name = "Treatment")

  df_plot  <- merge(df_scen, df_notrt, by = "race")

  plot_1 <- ggplot(data = df_scen, aes(color = race)) +
    geom_line(data = df_notrt,
              mapping = aes(x = time, y = I_notrt, size = scenario_name, linetype = scenario_name)) +
    geom_line(mapping = aes(x = time, y = I, size = scenario_name, linetype = scenario_name)) +
    scale_size_manual(values = c(0.5, 1), guide = 'none') +
    facet_wrap(. ~ scenario,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
    scale_linetype_manual(values = c(2, 1),
                          guide = guide_legend(override.aes=list(linetype=c(2,1),
                                                                 lwd=c(0.5,1)))) +
    scale_y_continuous(limits = c(0, 0.6) ,
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    xlab("Year") + ylab("Prevalence") +
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
graph_results_grid_all <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  df_notrt <- df_all %>%
    filter(scen == 0) %>%
    mutate(I_notrt = I) %>%
    # dplyr::select(time, I_notrt, race, scen_waifw) %>%
    mutate(scenario_name = "Status Quo")

  df_scen  <- df_all %>%
    filter(scen == 1) %>%
    mutate(scenario_name = "Treatment")

  # df_plot  <- merge(df_scen, df_notrt, by = "race")

  plot_1 <- ggplot(data = df_scen, aes(color = race)) +
    geom_line(data = df_notrt,
              mapping = aes(x = time, y = I_notrt, size = scenario, linetype = scenario)) +
    geom_line(mapping = aes(x = time, y = I, size = scenario, linetype = scenario)) +
    scale_size_manual(values = c(0.5, 0.5, 1.2, 1.2), guide = "none") +
    facet_grid(waifw ~ race,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    # scale_color_nejm() +
    scale_color_manual(values =
                         c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#20854EFF"),
                       guide = "none") +
    theme_bw() +
    scale_linetype_manual(values = c(2,3,1,5),
                          guide = guide_legend(override.aes=list(linetype=c(2,3,1,5),
                                                                 lwd=c(0.5,0.5,1.2,1.2)))) +
    scale_y_continuous(limits = c(0, 0.6) ,
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    xlab("Year") + ylab("Prevalence") +
    guides(# color = guide_legend(title = "Race/Ethnicity"),
           linetype = guide_legend(title = element_blank(),# "Scenario/WAIFW Structure",
                                   ncol = 4)) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))

  return(plot_1)
}

## Create WAIFW matrix
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
beta_mix_0 <- get_transmission_matrix(betas = rep(0, 8),
                                      breaks = waifw_breaks,
                                      w = w2,
                                      group_names = groups)
# beta_all <- rbind(cbind(beta_hisp,  beta_mix_0, beta_mix_0, beta_mix_0),
#                   cbind(beta_mix_0, beta_black, beta_mix_0, beta_mix_0),
#                   cbind(beta_mix_0, beta_mix_0, beta_white, beta_mix_0),
#                   cbind(beta_mix_0, beta_mix_0, beta_mix_0, beta_other))

beta_all <- as.matrix(read.csv(file = "results/assort_Beta_hat_w2_1.csv"))[,-1]

# alpha = 1 (No mixing)
beta_white_black_100 <- rbind(cbind(beta_white, beta_mix_0),
                              cbind(beta_mix_0, beta_black))

# alpha = 0.5 (some mixing)
alpha <- 0.5
beta1 <- ((1 - alpha) *  beta_white_r) + (alpha * beta_white_a)
beta2 <- ((1 - alpha) *  beta_black_r) + (alpha * beta_black_a)
beta3 <- ((1 - alpha) *  beta_mix_r)
beta_white_black_50 <- rbind(cbind(beta1, beta3),
                             cbind(beta3, beta2))


## Run Scenarios
# Age groups
groups <- c(1:80)
# Number of age groups
n_ages <- length(groups)
# Population growth
q <- 0
# WAIFWs
beta_rp <- as.matrix(read.csv(file = "results/waifw/rp_Beta_hat_w2_1.csv"))[,-1]
beta_assort <- as.matrix(read.csv(file = "results/waifw/assort_Beta_hat_w3_3.csv"))[,-1]

# Assortative
scenario_results <- get_sis_treat_scenario_all(
  m_waifw = beta_assort,
  m_waifw_burn = beta_rp,
  v_race = c("Hispanic", "NH White", "NH Black", "Other"),
  v_race_names = c("hisp", "white", "black", "other"))
plot_results_grid_assort <- graph_results_grid(scenario_results)
plot_results_grid_assort

# RP
scenario_results_rp <- get_sis_treat_scenario_all(
  m_waifw = beta_rp,
  m_waifw_burn = beta_rp,
  v_race = c("Hispanic", "NH White", "NH Black", "Other"),
  v_race_names = c("hisp", "white", "black", "other"))
plot_results_grid_rp <- graph_results_grid(scenario_results_rp)
plot_results_grid_rp

# Combine and Plot
# Change variable values to combine
scenario_results_assort <- scenario_results %>%
  mutate(waifw = "Assortative")

scenario_results_rp <- scenario_results_rp %>%
  mutate(waifw = "Random Partnership")

scenario_results_all <- rbind(scenario_results_assort, scenario_results_rp) %>%
  mutate(scen_waifw = ifelse(scen == 0 & waifw == "Assortative",
                             "Status Quo - Assortative",
                             ifelse(scen == 1 & waifw == "Assortative",
                                    "Treatment - Assortative",
                                    ifelse(scen == 0 & waifw == "Random Partnership",
                                           "Status Quo - Random Partnership",
                                           "Treatment - Random Partnership"))))

write.csv(scenario_results_all, "results/sis_scenario_all.csv")
scenario_results_all <- read.csv("results/sis_scenario_all.csv")
plot_sis_all <- graph_results_grid_all(scenario_results_all)
plot_sis_all
ggsave(plot = plot_sis_all, filename = "results/sis_scen_assort_rp_v3.pdf",
       height = 6, width = 12)

pct_impact <- function(df_all) {
  require(ggplot2)
  require(dplyr)
  require(ggsci)

  df_all <- df_all %>%
    filter(time == 0 | time == 50)

  df_wide <- dcast(setDT(df_test), scenario + waifw + race ~ time,
                   value.var = "I")
  colnames(df_wide) <- c("scenario", "waifw", "race", "t_0", "t_50")

  df_wide <- df_wide %>%
    filter(scenario != "No Treatment") %>%
    dplyr::mutate(pct_diff = abs(t_50 - t_0) / t_0) %>%
    dplyr::mutate(net_diff = t_0 - t_50 )

  plot1 <- ggplot(data = df_wide, aes(x = waifw, y = net_diff, fill = waifw)) +
    geom_bar(stat = "identity") +
    facet_grid(scenario ~ race,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Mixing Pattern") + ylab("Net Reduction in Prevalence (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")

  plot2 <- ggplot(data = df_wide, aes(x = waifw, y = pct_diff, fill = waifw)) +
    geom_bar(stat = "identity") +
    facet_grid(scenario ~ race,
               labeller = labeller(race = label_wrap_gen(21))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Mixing Pattern") + ylab("Relative Reduction in Prevalence (%)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")

  return(list(plot1, plot2))
}

pct_plots <- pct_impact(scenario_results_all)
ggsave(plot = pct_plots[[1]], filename = "results/sis_assort_rp_net.pdf",
       height = 8, width = 12)
ggsave(plot = pct_plots[[2]], filename = "results/sis_assort_rp_rel.pdf",
       height = 8, width = 12)

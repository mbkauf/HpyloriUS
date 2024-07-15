get_sis_abr_test_treat <- function(p_trt_eff = NULL,
                                   p_test_sens = NULL,
                                   p_test_spec = NULL,
                                   abr = TRUE,
                                   ages = groups,
                                   m_waifw = beta_all,
                                   m_waifw_burn,
                                   v_race = c("Hispanic",
                                              "NH Black",
                                              "NH White",
                                              "Other"),
                                   v_race_names = c("hisp", "black",
                                                    "white", "other")) {
  require(dplyr)
  require(rootSolve)
  v_results <- list()

  v_parameter <- load_sis_abr_model_params_all(
    v_race = c("Hispanic", "NH Black", "NH White", "Other"),
    waifw = m_waifw_burn,
    trt_year = 50000000,
    end_t = 1000,
    ages = groups)
  v_parameter$sigma <- 0
  v_parameter$v_alpha <- 0
  # Run burn-in period and store starting states
  burn_results <- runsteady(y = v_parameter$v_state,
                            times = c(0,1E5), func = sis_abr_model,
                            parms = v_parameter)
  v_burn_state <- burn_results$y

  # Scenario 1: No test and treat
  v_parameter <- load_sis_abr_model_params_all(
    v_race = c("Hispanic", "NH Black", "NH White", "Other"),
    waifw = m_waifw,
    trt_year = 50000000,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = groups)
  # if (abr == FALSE) (v_parameter$sigma <- 0)
  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[1]] <- tmp_results

  # Scenario 2: 90% test sensitivity
  v_parameter <- load_sis_abr_model_params_all(
    v_race = c("Hispanic", "NH Black", "NH White", "Other"),
    waifw = m_waifw,
    trt_year = 50000000,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = groups)
  v_parameter$test_sens <- 0.90

  # if (abr == FALSE) (v_parameter$sigma <- 0)
  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[2]] <- tmp_results

  # Scenario 3: 60% test sensitivity
  v_parameter <- load_sis_abr_model_params_all(
    v_race = c("Hispanic", "NH Black", "NH White", "Other"),
    waifw = m_waifw,
    trt_year = 50000000,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = groups)
  v_parameter$test_sens <- 0.60
  # if (abr == FALSE) (v_parameter$sigma <- 0)
  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[3]] <- tmp_results

  # Scenario 4: 90% test sensitivity (no testing before 18)
  v_parameter <- load_sis_abr_model_params_all(
    v_race = c("Hispanic", "NH Black", "NH White", "Other"),
    waifw = m_waifw,
    trt_year = 50000000,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = groups)
  v_parameter$test_sens <- 0.90
  v_parameter$age_based_test <- TRUE
  # if (abr == FALSE) (v_parameter$sigma <- 0)
  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[4]] <- tmp_results

  # Scenario 5: 60% test sensitivity (no testing before 18)
  v_parameter <- load_sis_abr_model_params_all(
    v_race = c("Hispanic", "NH Black", "NH White", "Other"),
    waifw = m_waifw,
    trt_year = 50000000,
    end_t = 50,
    burn_state = TRUE,
    v_init_state = v_burn_state,
    ages = groups)
  v_parameter$test_sens <- 0.60
  v_parameter$age_based_test <- TRUE
  # if (abr == FALSE) (v_parameter$sigma <- 0)
  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter,
    v_race_names = v_race_names)

  v_results[[5]] <- tmp_results

  # Combine Results
  # df_all <- rbind(
  #   v_results[[1]],
  #   v_results[[2]],
  #   v_results[[3]]
  # )
  df_all <- dplyr::bind_rows(v_results)


  m_hisp <- df_all %>%
    dplyr::select(c(time, s0_hisp, i0w_hisp, i0r_hisp, s1_hisp, i1r_hisp)) %>%
    dplyr::mutate(
      race = "Hispanic",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51), rep(4, 51), rep(5, 51))
    ) %>%
    dplyr::rename(
      s0  = s0_hisp,
      i0w = i0w_hisp,
      i0r = i0r_hisp,
      s1  = s1_hisp,
      i1r = i1r_hisp
    )

  m_black <- df_all %>%
    dplyr::select(c(time, s0_black, i0w_black, i0r_black, s1_black, i1r_black)) %>%
    dplyr::mutate(
      race = "NH Black",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51), rep(4, 51), rep(5, 51))
    ) %>%
    dplyr::rename(
      s0  = s0_black,
      i0w = i0w_black,
      i0r = i0r_black,
      s1  = s1_black,
      i1r = i1r_black
    )

  m_white <- df_all %>%
    dplyr::select(c(time, s0_white, i0w_white, i0r_white, s1_white, i1r_white)) %>%
    dplyr::mutate(
      race = "NH White",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51), rep(4, 51), rep(5, 51))
    ) %>%
    dplyr::rename(
      s0  = s0_white,
      i0w = i0w_white,
      i0r = i0r_white,
      s1  = s1_white,
      i1r = i1r_white
    )

  m_other <- df_all %>%
    dplyr::select(c(time, s0_other, i0w_other, i0r_other, s1_other, i1r_other)) %>%
    dplyr::mutate(
      race = "Other",
      scenario = c(rep(1, 51), rep(2, 51), rep(3, 51), rep(4, 51), rep(5, 51))
    ) %>%
    dplyr::rename(
      s0  = s0_other,
      i0w = i0w_other,
      i0r = i0r_other,
      s1  = s1_other,
      i1r = i1r_other
    )

  df_all_long <- rbind(m_hisp, m_black, m_white, m_other) %>%
    dplyr::mutate(I = i0w + i0r + i1r) %>%
    dplyr::mutate(Ir = i0r + i1r) %>%
    dplyr::mutate(Iw = i0w) %>%
    dplyr::mutate(scen = ifelse(scenario > 1, 1, 0)) %>%
    dplyr::mutate(scen = factor(scen, levels = c(0, 1))) %>%
    dplyr::mutate(scenario = factor(scenario, levels = c(1, 2, 3, 4, 5),
                                    labels = c("No Treatment",
                                               "All Ages: 90% Test Sensitivity",
                                               "All Ages: 60% Test Sensitivity",
                                               "18+: 90% Test Sensitivity",
                                               "18+: 60% Test Sensitivity")))

  return(df_all_long)
}

graph_results_grid <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  # df_notrt <- df_all %>%
  #   filter(scen == 0) %>%
  #   mutate(I_notrt = I) %>%
  #   dplyr::select(time, I_notrt, race) %>%
  #   mutate(scenario_name = "Status Quo")

  df_scen  <- df_all %>%
    filter(scen == 1) %>%
    mutate(scenario_name = "Treatment") %>%
    gather(strain, prev, I:Iw)


  # df_plot  <- merge(df_scen, df_notrt, by = "race")

  plot_1 <- ggplot(data = df_scen, aes(color = scenario)) +
    # geom_line(data = df_notrt,
    #           mapping = aes(x = time, y = I_notrt, size = scenario_name, linetype = scenario_name)) +
    geom_line(mapping = aes(x = time, y = prev, linetype = strain)) +
    facet_grid(race ~ .,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
    # scale_linetype_manual(values = c(2, 1, 3),
    #                       guide = guide_legend(override.aes=list(linetype=c(2,1,3)))) +
    scale_y_continuous(limits = c(0, 0.6) ,
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    xlab("Year") + ylab("Prevalence by Strain") +
    guides(color = guide_legend(title = element_blank(),
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

scenario_results_assort <- get_sis_abr_test_treat(
  m_waifw = beta_assort,
  m_waifw_burn = beta_rp,
  v_race = c("Hispanic", "NH White", "NH Black", "Other"),
  v_race_names = c("hisp", "white", "black", "other"))
plot_results_grid_assort <- graph_results_grid(scenario_results_assort)
plot_results_grid_assort
ggsave(plot = plot_results_grid_assort,
       filename = "results/abr_tnt_strain_assort.pdf",
       width = 10, height = 8)

scenario_results_rp <- get_sis_abr_test_treat(
  m_waifw = beta_rp,
  m_waifw_burn = beta_rp,
  v_race = c("Hispanic", "NH White", "NH Black", "Other"),
  v_race_names = c("hisp", "white", "black", "other"))
plot_results_grid_rp <- graph_results_grid(scenario_results_rp)
plot_results_grid_rp
ggsave(plot = plot_results_grid_rp,
       filename = "results/abr_tnt_strain_rp.pdf",
       width = 10, height = 8)

## Combine Results
scenario_results_assort$waifw <- "Assortative"
scenario_results_rp$waifw <- "Random Partnership"
scenario_results_all <- rbind(scenario_results_assort, scenario_results_rp)
graph_results_grid_all <- function(df_all) {
  require(ggplot2)
  require(ggsci)
  require(ggpubr)

  df_notrt <- df_all %>%
    filter(scen == 0) %>%
    # mutate(I_notrt = I) %>%
    mutate(scenario_name = "Status Quo") %>%
    gather(strain, prev, I:Iw)

  df_scen  <- df_all %>%
    filter(scen == 1) %>%
    mutate(scenario_name = "Treatment") %>%
    gather(strain, prev, I:Iw)


  # df_plot  <- merge(df_scen, df_notrt, by = "race")

  plot_1 <- ggplot(data = df_scen, aes(color = scenario)) +
    geom_line(data = df_notrt,
              mapping = aes(x = time, y = prev, linetype = strain)) +
    geom_line(mapping = aes(x = time, y = prev, linetype = strain)) +
    facet_grid(race ~ waifw,
               labeller = labeller(scenario = label_wrap_gen(21))) +
    scale_color_nejm() +
    theme_bw() +
    scale_linetype_manual(values = c(1, 2, 3),
                          guide = guide_legend(override.aes=list(linetype=c(1,2,3)))) +
    scale_y_continuous(limits = c(0, 0.6) ,
                       breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) +
    xlab("Year") + ylab("Prevalence by Strain") +
    guides(color = guide_legend(title = "Scenario:"),
           linetype = guide_legend(title = "Strain:")) +
    theme(legend.position = "bottom",
          legend.box="vertical",
          strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12))

  return(plot_1)
}

plot_results_all <- graph_results_grid_all(scenario_results_all)
plot_results_all
ggsave(plot = plot_results_all,
       filename = "results/abr_tnt_strain_rp_assort.pdf",
       width = 12, height = 10)

pct_impact_abr <- function(df_all) {
  require(ggplot2)
  require(dplyr)
  require(ggsci)

  df_all <- df_all %>%
    filter(time == 0 | time == 50)

  df_wide <- dcast(setDT(df_all), scenario + waifw + race ~ time,
                   value.var = "I")
  colnames(df_wide) <- c("scenario", "waifw", "race", "t_0", "t_50")

  df_wide <- df_wide %>%
    filter(scenario != "No Treatment") %>%
    dplyr::mutate(pct_diff = abs(t_50 - t_0) / t_0) %>%
    dplyr::mutate(net_diff = t_0 - t_50 )

  plot1 <- ggplot(data = df_wide, aes(x = waifw, y = net_diff, fill = waifw)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste(round(net_diff * 100, 1), "%"), vjust = -0.2)) +
    facet_grid(scenario ~ race,
               labeller = labeller(scenario = label_wrap_gen(width = 22))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Mixing Pattern") + ylab("Net Reduction in Prevalence (%)") +
    scale_y_continuous(limits = c(0, 0.6),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")

  plot2 <- ggplot(data = df_wide, aes(x = waifw, y = pct_diff, fill = waifw)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste(round(pct_diff * 100, 1), "%"), vjust = -0.2)) +
    facet_grid(scenario ~ race,
               labeller = labeller(scenario = label_wrap_gen(width = 22))) +
    scale_fill_nejm() +
    theme_bw() +
    xlab("Mixing Pattern") + ylab("Relative Reduction in Prevalence (%)") +
    scale_y_continuous(limits = c(0, 0.6),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
    theme(strip.text = element_text(size = 13),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          legend.position = "none")

  return(list(plot1, plot2))
}

pct_plots_abr <- pct_impact_abr(scenario_results_all)
ggsave(plot = pct_plots_abr[[1]], filename = "results/abr_assort_rp_net.pdf",
       height = 8, width = 8)
ggsave(plot = pct_plots_abr[[2]], filename = "results/abr_assort_rp_rel.pdf",
       height = 8, width = 8)


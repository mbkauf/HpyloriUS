scale_pct <- function(x) (x*100)

split_data_clean <- function(df_all) {
  df_all <- df_all %>%
    filter(time == 0 | time == 50)

  df_wide_I <- df_all %>%
    select(c(time, scenario, waifw, race, id, I)) %>%
    mutate(I = ifelse(race == "Hispanic", I / (608.095 / 3297.255),
                      ifelse(race == "NH White", I / (1960.185 / 3297.255),
                             ifelse(race == "NH Black", I / (402.225 / 3297.255),
                                    I / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = I) %>%
    mutate(assort = ifelse(waifw == "Assortative", 1, 0)) %>%
    mutate(rp = ifelse(waifw == "Random Partnership", 1, 0)) %>%
    mutate(pm = ifelse(waifw == "Proportional Mixing", 1, 0))

  df_wide_Ir <- df_all %>%
    select(c(time, scenario, waifw, race, id, Ir)) %>%
    mutate(Ir = ifelse(race == "Hispanic", Ir / (608.095 / 3297.255),
                       ifelse(race == "NH White", Ir / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Ir / (402.225 / 3297.255),
                                     Ir / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = Ir) %>%
    mutate(assort = ifelse(waifw == "Assortative", 1, 0)) %>%
    mutate(rp = ifelse(waifw == "Random Partnership", 1, 0)) %>%
    mutate(pm = ifelse(waifw == "Proportional Mixing", 1, 0))

  df_wide_Iw <- df_all %>%
    select(c(time, scenario, waifw, race, id, Iw)) %>%
    mutate(Iw = ifelse(race == "Hispanic", Iw / (608.095 / 3297.255),
                       ifelse(race == "NH White", Iw / (1960.185 / 3297.255),
                              ifelse(race == "NH Black", Iw / (402.225 / 3297.255),
                                     Iw / (326.750 / 3297.255))))) %>%
    pivot_wider(names_from = c(time, scenario),
                values_from = Iw) %>%
    mutate(assort = ifelse(waifw == "Assortative", 1, 0)) %>%
    mutate(rp = ifelse(waifw == "Random Partnership", 1, 0)) %>%
    mutate(pm = ifelse(waifw == "Proportional Mixing", 1, 0))

  df_wide_I$strain  <- "All Strains"
  df_wide_Ir$strain <- "Treatment Resistant"
  df_wide_Iw$strain <- "Treatment Sensitive"

  colnames(df_wide_I) <- c("waifw", "race", "id", "t_0_sq",
                           "t_50_sq", "t_0_p1", "t_50_p1",
                           "assort", "rp", "pm",  "strain")
  colnames(df_wide_Ir) <- c("waifw", "race", "id", "t_0_sq",
                            "t_50_sq", "t_0_p1", "t_50_p1",
                            "assort", "rp", "pm",  "strain")
  colnames(df_wide_Iw) <- c("waifw", "race", "id", "t_0_sq",
                            "t_50_sq", "t_0_p1", "t_50_p1",
                            "assort", "rp", "pm",  "strain")

  df_wide_I <- df_wide_I %>%
    mutate(across(c("t_0_sq", "t_50_sq", "t_0_p1", "t_50_p1"), scale_pct))
  df_wide_Ir <- df_wide_Ir %>%
    mutate(across(c("t_0_sq", "t_50_sq", "t_0_p1", "t_50_p1"), scale_pct))
  df_wide_Iw <- df_wide_Iw %>%
    mutate(across(c("t_0_sq", "t_50_sq", "t_0_p1", "t_50_p1"), scale_pct))


  df_wide_I <- df_wide_I %>%
    mutate(net_diff = (t_50_p1 - t_50_sq))
  df_wide_Ir <- df_wide_Ir %>%
    mutate(net_diff = (t_50_p1 - t_50_sq))
  df_wide_Iw <- df_wide_Iw %>%
    mutate(net_diff = (t_50_p1 - t_50_sq))

  l_wide_I  <- split(df_wide_I, as.factor(df_wide_I$race))
  l_wide_Ir <- split(df_wide_Ir, as.factor(df_wide_Ir$race))
  l_wide_Iw <- split(df_wide_Iw, as.factor(df_wide_Iw$race))

  l_df_all <- list(l_wide_I[[1]], l_wide_I[[2]], l_wide_I[[3]],
                   l_wide_I[[4]], l_wide_Ir[[1]], l_wide_Ir[[2]],
                   l_wide_Ir[[3]], l_wide_Ir[[4]], l_wide_Iw[[1]],
                   l_wide_Iw[[2]], l_wide_Iw[[3]], l_wide_Iw[[4]])

  return(l_df_all)
}



reg_and_tex <- function(l_df) {
  library(texreg)
  l_df_all <- list
  l_tex <- list()
  for (i in 1:12) {
    tmp_reg <- lm(lm(net_diff ~ t_0_sq + rp + pm, data = l_df[[i]]))
    tmp_tex <- createTexreg(coef.names = c("Intercept", "Prev. SQ - T0",
                                           "Random Partnership",
                                           "Proportional Mixing"),
                            coef = tmp_reg$coefficients,
                            se = summary(tmp_reg)$coefficients[, 2],
                            pvalues = summary(tmp_reg)$coefficients[, 4])
    l_tex[[i]] <- tmp_tex
  }
  return(l_tex)
}

l_df_all      <- split_data_clean(scenario_results_all)
l_df_all_race <- split_data_clean(scenario_results_all_race)

l_tex_all  <- reg_and_tex(l_df_all)
l_tex_race <- reg_and_tex(l_df_all_race)

all_latex <- texreg(l_tex_all, digits = 4, stars = c(0.001, 0.01, 0.05),
                    booktabs = TRUE,
                    custom.model.names = c("I_hisp", "I_black", "I_white",
                                           "I_other", "Ir_hisp", "Ir_black",
                                           "Ir_white", "Ir_other", "Iw_hisp",
                                           "Iw_black", "Iw_white", "Iw_other"))
write.table(all_latex, file = "results/reg_latex_all.txt")

all_latex_race <- texreg(l_tex_race, digits = 4, stars = c(0.001, 0.01, 0.05),
                    booktabs = TRUE,
                    custom.model.names = c("I_hisp", "I_black", "I_white",
                                           "I_other", "Ir_hisp", "Ir_black",
                                           "Ir_white", "Ir_other", "Iw_hisp",
                                           "Iw_black", "Iw_white", "Iw_other"))
write.table(all_latex_race, file = "results/reg_latex_all_race.txt")


library(modelsummary)
library(ggplot2)
library(ggsci)

results_all <- lapply(1:12, function(x) lm(net_diff ~ t_0_sq + rp + pm, data = l_df_all[[x]])) |>
  modelplot(draw = FALSE) |>
  transform("Policy" = "Test-and-Treat All Adults")
results_target <- lapply(1:12, function(x) lm(net_diff ~ t_0_sq + rp + pm, data = l_df_all_race[[x]])) |>
  modelplot(draw = FALSE) |>
  transform("Policy" = "Targeted Test-and-Treat")
results <- rbind(results_all, results_target) %>%
  filter(term == "rp" | term == "pm") %>%
  mutate(term = ifelse(term == "rp", "Random Partnership", "Proportional Mixing" ))
results$model <- factor(rep(c("Hispanic - All Strains", "NH Black - All Strains",
                   "NH White - All Strains", "NH Other - All Strains",
                   "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                   "NH White - Resistant Strains", "NH Other - Resistant Strains",
                   "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                   "NH White - Sensitive Strains", "NH Other - Sensitive Strains"), 2),
                   levels = c("Hispanic - All Strains", "NH Black - All Strains",
                                  "NH White - All Strains", "NH Other - All Strains",
                                  "Hispanic - Resistant Strain", "NH Black - Resistant Strains",
                                  "NH White - Resistant Strains", "NH Other - Resistant Strains",
                                  "Hispanic - Sensitive Strain", "NH Black - Sensitive Strains",
                                  "NH White - Sensitive Strains", "NH Other - Sensitive Strains"))

coef_plot <- ggplot(results, aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_pointrange(aes(color = Policy), position = position_dodge(width = .5)) +
  facet_wrap(~model) +
  theme_bw() +
  scale_color_nejm() +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 7)) +
  xlab("Percentage Point Difference in Prevalence Compared to Assortative WAIFW") +
  ylab("Mixing Matrix") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")
ggsave("results/coef_plot.png", coef_plot, height = 7, width = 10)


burn_and_current <- function(v_params) {
  require(rootSolve)
  require(deSolve)
  # Make copy of v_parameter and update for burn-in
  v_parameter_burn         <- v_params
  v_parameter_burn$sigma   <- 0
  v_parameter_burn$v_alpha <- rep(0, length(v_parameter_burn$v_alpha))
  v_parameter_burn$v_psi   <- 0
  v_parameter_burn$v_eta   <- rep(0, length(v_parameter_burn$v_eta))
  v_parameter_burn$burn    <- TRUE

  # Run burn-in period and store starting states
  burn_results <- runsteady(y = v_parameter_burn$v_state,
                            times = c(0, 1E5), func = sis_abr_model,
                            parms = v_parameter_burn,
                            verbose = FALSE)
  v_burn_state <- burn_results$y

  v_parameter_2024            <- v_params
  v_parameter_2024$v_state    <- v_burn_state
  v_parameter_2024$v_eta      <- rep(0, length(v_parameter_2024$v_eta))
  v_parameter_2024$to_present <- TRUE
  v_parameter_2024$v_time     <- seq(0, 33, by = 1)

  tmp_results <- get_sis_abr_model_results_all(
    v_params = v_parameter_2024,
    v_race_names = v_race_names
  )

  v_current_state <- as.numeric(tmp_results[nrow(tmp_results),
                                            2:1601])
  return(list(v_burn_state, v_current_state))
}

burn_time <- function(v_params) {
  require(rootSolve)
  require(deSolve)
  # Make copy of v_parameter and update for burn-in
  v_parameter_burn         <- v_params
  v_parameter_burn$sigma   <- 0
  v_parameter_burn$v_alpha <- rep(0, length(v_parameter_burn$v_alpha))
  v_parameter_burn$v_psi   <- 0
  v_parameter_burn$v_eta   <- rep(0, length(v_parameter_burn$v_eta))
  v_parameter_burn$burn    <- TRUE


  # Run burn-in period and store starting states (runsteady)
  start_time <- Sys.time()
  burn_results <- runsteady(y = v_parameter_burn$v_state,
                            times = c(0, 1E5), func = sis_abr_model,
                            parms = v_parameter_burn,
                            verbose = FALSE)
  end_time <- Sys.time()
  time_runsteady <- end_time - start_time
  # Run burn-in period and store starting states (stodes)
  start_time <- Sys.time()
  burn_results <- stodes(y = v_parameter_burn$v_state,
                         func = sis_abr_model,
                         parms = v_parameter_burn,
                         verbose = FALSE)
  end_time <- Sys.time()
  time_stodes <- end_time - start_time

  return(list(time_runsteady, time_stodes))
}

assort_resamp <- imis_resample_5_assort
rp_resamp <- imis_resample_5_rp
alphas_resamp <- imis_resample_5_alphas

w2 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
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
groups  <- c(1:80)
n_ages  <- length(groups)
v_p <- seq(from = 0.005, to = 0.05, by = 0.005)

v_race       <- c("Hispanic", "NH White", "NH Black", "Other")
v_race_names <- c("hisp", "white", "black", "other")

v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race,
  waifw = NULL,
  trt_year = 50000000,
  end_t = 1000,
  ages = groups,
  prob = FALSE,
  p_trans = v_p[5]
)

v_race_prop <- v_parameter$v_pop / sum(v_parameter$v_pop)
v_race_prop <- rep(v_race_prop, each = length(groups))
v_race_prop <- rep(v_race_prop, 5)

v_parameter$v_state <- v_parameter$v_state * v_race_prop
v_parameter$b <- v_parameter$b * (v_parameter$v_pop / sum(v_parameter$v_pop))

waifw_assort <- foi_transition_matrix_assort(assort_resamp[i, ])
waifw_rp     <- foi_transition_matrix_rp(rp_resamp[i, ])
v_alphas_within     <- as.numeric(alphas_resamp[i, c(1, 6, 11, 16)])
v_alphas_within_off <- as.numeric(1 - alphas_resamp[i, c(1, 6, 11, 16)])
v_alphas_between    <- as.numeric(alphas_resamp[i, c(2:5, 7:10, 12:15)])
m_within      <- alpha_matrix_within(v_alphas_within)
m_within_off  <- alpha_matrix_within(v_alphas_within_off)
m_between     <- alpha_matrix_between(v_alphas_between)
waifw_pm  <- (m_within * waifw_assort) + (m_within_off * waifw_rp) +
  (m_between * waifw_rp)

v_parameter_assort <- v_parameter
v_parameter_assort$m_waifw <- waifw_assort

v_parameter_rp <- v_parameter
v_parameter_rp$m_waifw <- waifw_rp

v_parameter_pm <- v_parameter
v_parameter_pm$m_waifw <- waifw_pm

results_assort <- burn_and_current(v_parameter_assort)
results_rp     <- burn_and_current(v_parameter_rp)
results_pm     <- burn_and_current(v_parameter_pm)

plot(v_parameter$v_state[321:640] - results_assort[[1]][321:640])
plot(v_parameter$v_state[321:640] - results_rp[[1]][321:640])
plot(v_parameter$v_state[321:640] - results_pm[[1]][321:640])

sum(v_parameter$v_state[321:640] - results_assort[[1]][321:640])
sum(v_parameter$v_state[321:640] - results_rp[[1]][321:640])
sum(v_parameter$v_state[321:640] - results_pm[[1]][321:640])

sum(results_assort[[2]][321:960])
sum(results_rp[[2]][321:960])
sum(results_pm[[2]][321:960])

burn_time(v_parameter_assort)
burn_time(v_parameter_rp)
burn_time(v_parameter_pm)

## WAIFW parameters
waifw_breaks      <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
beta0             <- c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64",
                       "65-69", "70-80")

## Age groups
groups <- c(1:80)

## Number of age groups
n_ages <- length(groups)

## Population growth
q <- 0

## WAIFW matrix
beta1 <- get_transmission_matrix(betas = calibrate_white$m_beta_hat[1, ],
                                 breaks = waifw_breaks,
                                 w = w1,
                                 group_names = groups)
beta2 <- get_transmission_matrix(betas = calibrate_white$m_beta_hat[2, ],
                                 breaks = waifw_breaks,
                                 w = w2,
                                 group_names = groups)
beta3 <- get_transmission_matrix(betas = calibrate_white$m_beta_hat[3, ],
                                 breaks = waifw_breaks,
                                 w = w3,
                                 group_names = groups)



## SI model with betas from WAIFW structures 1-3 (NH White)
v_demo <- get_demographic_vars(spec_groups = groups,
                               race = "NH White")
v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                   race = "NH White")

v_parameter <- load_si_model_params(waifw = beta1, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_si_1 <- get_si_model_results(si_model)

v_parameter <- load_si_model_params(waifw = beta2, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_si_2 <- get_si_model_results(si_model)

v_parameter <- load_si_model_params(waifw = beta3, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_si_3 <- get_si_model_results(si_model)

## SIS model with betas from WAIFW structure 3
v_demo <- get_demographic_vars(spec_groups = groups,
                               race = "NH White")
v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                   race = "NH White")
v_parameter <- load_sis_model_params(waifw = beta1,
                                     demography_vars = v_demo,
                                     prevalence_vars = v_prev_vars,
                                     ages = groups,
                                     trt_year = 500,
                                     end_t = 500)

model_results_sis <- get_sis_model_results(v_params = v_parameter)


### Find Steady State
beta_hisp <- get_transmission_matrix(betas = calibrate_hispanic$m_beta_hat[2, ],
                                     breaks = waifw_breaks,
                                     w = w2,
                                     group_names = groups)
beta_white <- get_transmission_matrix(betas = calibrate_white$m_beta_hat[2, ],
                                      breaks = waifw_breaks,
                                      w = w2,
                                      group_names = groups)
beta_black <- get_transmission_matrix(betas = calibrate_black$m_beta_hat[2, ],
                                      breaks = waifw_breaks,
                                      w = w2,
                                      group_names = groups)
beta_other <- get_transmission_matrix(betas = calibrate_other$m_beta_hat[2, ],
                                      breaks = waifw_breaks,
                                      w = w2,
                                      group_names = groups)

v_demo <- get_demographic_vars(spec_groups = groups,
                               race = "Hispanic")
v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                   race = "Hispanic")
v_parameter <- load_sis_model_params(waifw = beta_hisp,
                                     demography_vars = v_demo,
                                     prevalence_vars = v_prev_vars,
                                     ages = groups,
                                     trt_year = 1000,
                                     end_t = 500)

model_results_sis <- get_sis_model_results(v_params = v_parameter)


# Check deaths = births
hisp_mu <- sum(v_parameter$v_mu * v_demo$v_age_prop)

# Run mutli-race load params
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
beta_all <- rbind(cbind(beta_hisp,  beta_mix_0, beta_mix_0, beta_mix_0),
                  cbind(beta_mix_0, beta_black, beta_mix_0, beta_mix_0),
                  cbind(beta_mix_0, beta_mix_0, beta_white, beta_mix_0),
                  cbind(beta_mix_0, beta_mix_0, beta_mix_0, beta_other))

v_parameter <- load_sis_model_params_all(
  v_race = c("Hispanic", "NH Black", "NH White", "Other"),
  waifw = beta_all,
  ages = groups)

model_results_sis_all <- get_sis_model_results_all(
  v_params = v_parameter,
  v_race_names = c("hisp", "black", "white", "other"))


## Test ABR model
# Age groups
groups <- c(1:80)
# Number of age groups
n_ages <- length(groups)
# Population growth
q <- 0
# Vector of Races
v_race <- c("Hispanic", "NH Black", "NH White", "Other")

# Assortative
beta_assort <- as.matrix(read.csv(file = "results/waifw/assort_Beta_hat_w2_1.csv"))[,-1]

v_parameter <- load_sis_abr_model_params_all(
  v_race = v_race,
  waifw = rand_waifw,
  trt_year = 1,
  end_t = 50,
  ages = groups,
  prob = TRUE)

results <- deSolve::ode(y = v_parameter$v_state, times = v_parameter$v_time,
                        func = sis_abr_model, parms = v_parameter,
                        method = "ode45")


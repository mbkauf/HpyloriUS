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
                                 group.names = groups)
beta2 <- get_transmission_matrix(betas = calibrate_white$m_beta_hat[2, ],
                                 breaks = waifw_breaks,
                                 w = w2,
                                 group.names = groups)
beta3 <- get_transmission_matrix(betas = calibrate_white$m_beta_hat[3, ],
                                 breaks = waifw_breaks,
                                 w = w3,
                                 group.names = groups)



## SI model with betas from WAIFW structures 1-3 (NH White)
v_demo <- get_demographic_vars(spec_groups = groups,
                               race = "NH White")
v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                   race = "NH White")

v_parameter <- load_si_model_params(waifw = beta1, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_si_1 <- get_model_results(si_model)

v_parameter <- load_si_model_params(waifw = beta2, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_si_2 <- get_model_results(si_model)

v_parameter <- load_si_model_params(waifw = beta3, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_si_3 <- get_model_results(si_model)

## SIS model with betas from WAIFW structure 3
v_demo <- get_demographic_vars(spec_groups = groups,
                               Race = "NH White")
v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                   Race = "NH White")
v_parameter <- load_sis_model_params(waifw = beta3, demography_vars = v_demo,
                                     prevalence_vars = v_prev_vars,
                                     ages = groups)

model_results_sis <- get_sis_model_results()
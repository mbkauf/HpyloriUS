parameters <- list(beta = beta0, # beta.hat
                   #c(1.631379e-08, 7.865094e-09, 6.735575e-09, 8.119121e-09, 3.556421e-09, 1.207767e-09)*150000,
                   # c(0.0117301564, 0.0103771644, 0.0084906907, 0.0051029252, 0.0007120402, 0.0000000000),
                   W = W,
                   # alpha = alpha,
                   mu = mu, b = b, d = d,
                   waifw.breaks = c(0, 2, 6, 12, 19, 45, 71),
                   groups = groups, age = age, n.ages = n.ages,
                   # age.groups.inf = age.groups.inf,
                   n.strains = n.strains,
                   prev = prev, prev_SD = prev_SD)


## WAIFW parameters
waifw.breaks <- c(1, 5, 15, 25, 45, 55, 65, 70, 81)
beta0 <- c(0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
v_age_group_names <- c("1-4", "5-14", "15-24", "25-44", "45-54", "55-64", "65-69", "70-80")

## Age groups
groups <- c(1:80)

## Number of age groups
n.ages <- length(groups)

## Population growth
q <- 0

## WAIFW matrix
test.w <- matrix(c(1, 9, 9, 9, 9, 9, 9, 9,
                   0, 2, 9, 9, 9, 9, 9, 9,
                   0, 0, 3, 9, 9, 9, 9, 9,
                   0, 0, 0, 4, 9, 9, 9, 9,
                   0, 0, 0, 0, 5, 9, 9, 9,
                   0, 0, 0, 0, 0, 6, 9, 9,
                   0, 0, 0, 0, 0, 0, 7, 9,
                   0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = T)

index.upper.w2 <- matrix(c(1, 1, 3, 4, 5, 6, 7, 8,
                           0, 2, 3, 4, 5, 6, 7, 8,
                           0, 0, 3, 4, 5, 6, 7, 8,
                           0, 0, 0, 4, 5, 6, 7, 8,
                           0, 0, 0, 0, 5, 6, 7, 8,
                           0, 0, 0, 0, 0, 6, 7, 8,
                           0, 0, 0, 0, 0, 0, 7, 8,
                           0, 0, 0, 0, 0, 0, 0, 8), ncol = 8, byrow = T)

Beta1 <- get_transmission_matrix(betas = calibrateWhite$m_beta_hat[1,],
                                 breaks = waifw.breaks,
                                 W = W1,
                                 group.names = groups)
Beta2 <- get_transmission_matrix(betas = calibrateWhite$m_beta_hat[2,],
                                 breaks = waifw.breaks,
                                 W = W2,
                                 group.names = groups)
Beta3 <- get_transmission_matrix(betas = calibrateWhite$m_beta_hat[3,],
                                 breaks = waifw.breaks,
                                 W = W3,
                                 group.names = groups)

# v_state <- v_parameter$v_state

# older_upper_lim <- 101 # max(groups) + 1

v_demo <- get_demographic_vars(spec_groups = groups,
                               Race = "Non-Hispanic White")
v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                   Race = "Non-Hispanic White")


v_parameter <- load_SI_model_params(waifw = Beta1, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_SI_1 <- get_model_results(SI_model)

v_parameter <- load_SI_model_params(waifw = Beta2, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_SI_2 <- get_model_results(SI_model)

v_parameter <- load_SI_model_params(waifw = Beta3, demography_vars = v_demo,
                                    prevalence_vars = v_prev_vars,
                                    ages = groups)
model_results_SI_3 <- get_model_results(SI_model)



v_parameter <- load_SIS_model_params(waifw = Beta3, demography_vars = v_demo,
                                     prevalence_vars = v_prev_vars,
                                     ages = groups)

model_results_SIS <- get_SIS_model_results()
View(model_results_SIS[[1]])

sum(model_results_SIS[[2]]$S0, model_results_SIS[[2]]$I0,
    model_results_SIS[[2]]$S1, model_results_SIS[[2]]$I1)

steady_results <- find_steady()
steady_results

sum(model_results_SI_1[[2]]$S, model_results_SI_1[[2]]$I)


## Run SIS Model for NH White
start_time <- Sys.time()

v_demo <- get_demographic_vars(spec_groups = groups,
                               Race = "Non-Hispanic White")
v_prev_vars <- get_prevalence_vars(spec_groups = groups,
                                   Race = "Non-Hispanic White")
Beta3 <- get_transmission_matrix(betas = calibrateWhite$m_beta_hat[3,],
                                 breaks = waifw.breaks,
                                 W = W3,
                                 group.names = groups)
v_parameter <- load_SIS_model_params(waifw = Beta3, demography_vars = v_demo,
                                     prevalence_vars = v_prev_vars,
                                     ages = groups)


start_time <- Sys.time()
model_results_SIS <- get_SIS_model_results()
end_time <- Sys.time()
model_results_SIS[[1]]


end_time - start_time

## Testing parts of SIS model
m_state <- matrix(
  v_parameter$v_state,
  nrow = v_parameter$n_age_groups,
  ncol = v_parameter$n_id_states,
  dimnames = list(v_parameter$v_names_age_groups, v_parameter$v_names_id_states)
)

load_SIS_model_params <- function(waifw, demography_vars, prevalence_vars,
                                  ages) {

  ## Get initial state
  sus <- (1-prevalence_vars$prevalence) * demography_vars$v_age_prop # susceptibles
  inf <- prevalence_vars$prevalence * demography_vars$v_age_prop  # infected
  a.init <- cbind(sus, inf,
                  rep(0, length(ages)),
                  rep(0, length(ages)))  # combine susceptibles and infected
  dimnames(a.init) <- list(groups, c("S0","I0","I1","S1"))  # add names
  v.init <- as.vector(a.init)  # transform into vector to be used by solver

  v_parameter <- list(
    v_time             = seq (0, 500, by=1),   # time steps to test
    N                  = 1,                    # population size
    b                  = demography_vars$b,    # birth rate
    v_d                = demography_vars$v_d,  # aging rates
    v_mu               = demography_vars$v_mu, # death rates
    m_waifw            = waifw,                # transmission rates
    v_state            = v.init,               # list of initial state values
    v_names_age_groups = as.character(ages),   # ages as character vector
    n_id_states        = length(c("S0","I0","I1","S1")),  # number of states
    v_names_id_states  = c("S0","I0","I1","S1"),  # character vector of states
    n_age_groups       = length(ages),         # number of age groups
    v_test             = 0, # probability of testing infected
    v_sens             = 0, # test sensitivity
    v_trt_success      = 1, # probability of treatment working
    v_f                = 0, # 1 / duration of infection untreated
    v_f_prime          = 0, # 1 / duration of infection treated
    v_gamma            = 1/14, # 1 / treatment length
    v_trtS             = 0, # probability of treatment when susceptible
    v_trtI             = 0  # probability of treatment when infected
  )

  return(v_parameter)
}


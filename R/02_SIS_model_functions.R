SIS_model <- function(v_time,
                      v_state,
                      v_params) {

  # Assign values from parameters
  m_waifw            <- v_params$m_waifw
  n_age_groups       <- v_params$n_age_groups
  n_id_states        <- v_params$n_id_states
  v_names_age_groups <- v_params$v_names_age_groups
  v_names_id_states  <- v_params$v_names_id_states
  b                  <- v_params$b
  v_d                <- v_params$v_d
  v_mu               <- v_params$v_mu
  v_test             <- v_params$v_test  # probability of testing infected
  v_sens             <- v_params$v_sens  # test sensitivity
  v_trt_success      <- v_params$v_trt_success # probability of treatment working
  v_f                <- v_params$v_f
  v_f_prime          <- v_params$v_f_prime
  v_gamma            <- v_params$v_gamma
  v_trtS             <- v_params$v_trtS
  v_trtI             <- v_params$v_trtI

  # One-time intervention
  # v_trtS <- ifelse(
  #   (v_time==0), 0.05, v_trtS
  # )
  # v_trtI <- ifelse(
  #   (v_time==0), 0.05, v_trtI
  # )

  # Generate matrix of proportion in each compartment
  m_state <- matrix(
    v_state,
    nrow = n_age_groups,
    ncol = n_id_states,
    dimnames = list(v_names_age_groups, v_names_id_states)
  )

  # Split m_state matrix into vectors of susceptibles and infected
  v_S0 <- m_state[, "S0"]
  v_I0 <- m_state[, "I0"]
  v_S1 <- m_state[, "S1"]
  v_I1 <- m_state[, "I1"]

  # Force of infection lambda values
  v_lambda       <- m_waifw %*% (v_I0 + v_I1)
  v_lambda_prime <- m_waifw %*% (v_I0 + v_I1)

  # Growth rates (Birth rate and aging rate)
  v_S0g <- c(b, (v_d * v_S0)[-n_age_groups])
  v_I0g <- c(0, (v_d * v_I0)[-n_age_groups])
  v_S1g <- c(0, (v_d * v_S1)[-n_age_groups])
  v_I1g <- c(0, (v_d * v_I1)[-n_age_groups])

  # Transitions
  dS0_dt = v_S0g + (v_f * v_I0) + (v_gamma * v_S1) - (v_mu + v_d + v_lambda + v_trtS) * v_S0
  dI0_dt = v_I0g + (v_lambda * v_S0) - (v_mu + v_d + v_trtI + v_f) * v_I0
  dS1_dt = v_S1g + (v_trtS * v_S0) + (v_f_prime * v_I1) + (v_trtI * v_trt_success) * v_I0 - (v_mu + v_d + v_gamma + v_lambda_prime) * v_S1
  dI1_dt = v_I1g + (v_trtI * (1-v_trt_success))* v_I0 + (v_lambda_prime * v_S1) - (v_mu + v_d + v_f_prime) * v_I1

  # combine results
  return(list(c(dS0_dt, dI0_dt, dS1_dt, dI1_dt)))

}


#' Get and plot model results
#'
#' @return list(plot, desolver))
#' @export
#'
#' @examples
get_SIS_model_results <- function() {
  require(deSolve)
  require(ggplot2)

  desolver <- deSolve::lsoda(y=v_parameter$v_state, times = v_parameter$v_time,
                             func = SIS_model, parms = v_parameter)

  model_results <- as.data.frame(desolver) %>%
    mutate(S0 = rowSums(.[2:82])) %>%
    mutate(I0 = rowSums(.[83:163])) %>%
    mutate(S1 = rowSums(.[164:244])) %>%
    mutate(I1 = rowSums(.[245:325]))

  plot <- ggplot() +
    geom_line(data = model_results, aes(x = time, y = S0), color = "blue") +
    geom_line(data = model_results, aes(x = time, y = I0), color = "red") +
    ylab("Proportion")

  return(list(plot, model_results))
}




### This file contains functions to generate a discrete time version of the model
n_cycles <- 52  # number of cycles per year

convert_params <- function(t, v_params) {

  # v_params$b      <- 1 - exp(-v_params$b*(1/t))
  v_params$v_d    <- 1 - exp(-v_params$v_d*(1/t))
  v_params$v_mu   <- 1 - exp(-v_params$v_mu*(1/t))
  v_params$gamma  <- 1 - exp(-v_params$gamma*(1/t))

  return(v_params)
}

v_parameter_test <- convert_params(t = n_cycles, v_params = v_parameter)


## Initialize model
m_state <- matrix(
  v_parameter_test$v_state,
  nrow = v_parameter_test$n_age_groups,
  ncol = v_parameter_test$n_id_states,
  dimnames = list(v_parameter_test$v_names_age_groups, v_parameter_test$v_names_id_states)
)

v_S0 <- m_state[, "S0"]
v_I0 <- m_state[, "I0"]
v_S1 <- m_state[, "S1"]
v_I1 <- m_state[, "I1"]

v_S0g <- c(b, (v_d * v_S0)[-n_age_groups])

for(i in 1:520) {
  v_S0g <- c(v_parameter_test$b,
             (v_parameter_test$v_d * v_S0)[-v_parameter_test$n_age_groups])

  v_S0_new <- v_S0 + v_S0g + (gamma * v_S1) - v


  v_S0g + (gamma * v_S1) - (v_mu + v_d + v_lambda + p_trtS) * v_S0

}

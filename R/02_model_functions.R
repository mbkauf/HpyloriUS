#' Function for SI model
#'
#' @param v_time  vector of times to use in model
#' @param v_state  vector of propotion in each compartment
#' @param v_params  vector of input parameters
#'
#' @return list(c(ds_dt, di_dt))
#' @export
#'
#' @examples
#' deSolve::lsoda(y=v_parameter$v_state, times = v_parameter$v_time,
#'                func = SI_model, parms = v_parameter)
si_model <- function(v_time,
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

  # Generate matrix of proportion in each compartment
  m_state <- matrix(
    v_state,
    nrow = n_age_groups,
    ncol = n_id_states,
    dimnames = list(v_names_age_groups, v_names_id_states)
  )

  # Split m_state matrix into vectors of susceptibles and infected
  v_s <- m_state[, "S"]
  v_i <- m_state[, "I"]

  # Force of infection lambda values
  v_lambda <- m_waifw %*% v_i

  # Growth rates (Birth rate and aging rate)
  v_sg <- c(b, (v_d * v_s)[-n_age_groups])
  v_ig <- c(0, (v_d * v_i)[-n_age_groups])

  # Transitions
  ds_dt <- v_sg - (v_mu + v_d + v_lambda) * v_s
  di_dt <- v_ig + (v_lambda * v_s) - (v_mu + v_d) * v_i
  # ds_dt = v_sg - v_s
  # di_dt = v_ig + (v_lambda * v_s) - v_i

  # combine results
  return(list(c(ds_dt, di_dt)))

}

#' Function for SIS Model
#'
#' @param v_time  vector of times to use in model
#' @param v_state  vector of propotion in each compartment
#' @param v_params  vector of input parameters
#'
#' @return list(c(ds_dt, di_dt))
#'
#' @examples
sis_model <- function(v_time,
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
  trt_eff            <- v_params$trt_eff    # probability of treatment working
  gamma              <- v_params$gamma
  p_trt_s             <- v_params$p_trt_s
  p_trt_i             <- v_params$p_trt_i

  # One-time intervention
  # p_trtS <- ifelse(
  #   floor(v_time)==10, 0.1, 0
  # )
  # p_trtI <- ifelse(
  #   floor(v_time)==10, 0.1, 0
  # )

  # Generate matrix of proportion in each compartment
  m_state <- matrix(
    v_state,
    nrow = n_age_groups,
    ncol = n_id_states,
    dimnames = list(v_names_age_groups, v_names_id_states)
  )

  # Split m_state matrix into vectors of susceptibles and infected
  v_s0 <- m_state[, "S0"]
  v_i0 <- m_state[, "I0"]
  v_s1 <- m_state[, "S1"]
  v_i1 <- m_state[, "I1"]

  # Force of infection lambda values
  v_lambda       <- m_waifw %*% (v_i0 + v_i1)
  # v_lambda_prime <- m_waifw %*% (v_i0 + v_i1)
  # FIXME - figure out what to do with model not working
  # TODO - test forced non-negative growth rate and lambda

  # Growth rates (Birth rate and aging rate)
  v_s0g <- c(b, (v_d * v_s0)[-n_age_groups])
  v_i0g <- c(0, (v_d * v_i0)[-n_age_groups])
  v_s1g <- c(0, (v_d * v_s1)[-n_age_groups])
  v_i1g <- c(0, (v_d * v_i1)[-n_age_groups])

  # Transitions
  ds0_dt <- v_s0g + (gamma * v_s1) -
    (v_mu + v_d + v_lambda + p_trt_s) * v_s0
  di0_dt <- v_i0g + (v_lambda * v_s0) + (gamma * v_i1) -
    (v_mu + v_d + p_trt_i) * v_i0
  ds1_dt <- v_s1g + (p_trt_s * v_s0) + (p_trt_i * trt_eff) * v_i0 -
    (v_mu + v_d + gamma + v_lambda) * v_s1
  di1_dt <- v_i1g + (p_trt_i * (1 - trt_eff)) * v_i0 + (v_lambda * v_s1) -
    (v_mu + v_d + gamma) * v_i1

  # combine results
  return(list(c(ds0_dt, di0_dt, ds1_dt, di1_dt)))

}


sis_abr_model <- function(v_time,
                          v_state,
                          v_params) {

  ### Assign values from parameters
  m_waifw            <- v_params$m_waifw
  m_waifw_prime      <- v_params$m_waifw_prime
  m_waifw_xi         <- v_params$m_waifw_xi
  n_age_groups       <- v_params$n_age_groups
  n_id_states        <- v_params$n_id_states
  v_names_age_groups <- v_params$v_names_age_groups
  v_names_id_states  <- v_params$v_names_id_states
  b                  <- v_params$b
  v_d                <- v_params$v_d
  v_mu               <- v_params$v_mu
  gamma              <- v_params$gamma
  v_alpha            <- v_params$v_alpha
  v_psi              <- v_params$v_psi

  m_state <- matrix(
    v_state,
    nrow = n_age_groups,
    ncol = n_id_states,
    dimnames = list(v_names_age_groups, v_names_id_states)
  )

  ### Antibiotic policy
  # AB policy 1: No uptake
  # AB policy 2: Treat only in year x
  # if (ab_policy == 1) {
  #   v_psi <- rep(0, n_age_groups)
  # } else if (ab_policy == 2) {
  #   if(floor(v_time) == 0){
  #     v_psi <- rep(0.05, n_age_groups)
  #   } else {
  #     v_psi <- rep(0, n_age_groups)
  #   }
  # }

  ### Split m_state matrix into vectors of susceptibles and infected
  v_s0  <- m_state[, "S0"]
  v_i0w <- m_state[, "I0w"]
  v_i0r <- m_state[, "I0r"]
  v_s1  <- m_state[, "S1"]
  v_i1r <- m_state[, "I1r"]

  ### Force of infection lambda values
  v_lambda       <- m_waifw %*% v_i0w
  v_lambda_prime <- m_waifw_prime %*% (v_i0r + v_i1r)
  v_lambda_xi    <- m_waifw_xi %*% (v_i0r + v_i1r)

  ### Growth rates (Birth rate and aging rate)
  v_s0g  <- c(b, (v_d * v_s0)[-n_age_groups])
  v_i0wg <- c(0, (v_d * v_i0w)[-n_age_groups])
  v_i0rg <- c(0, (v_d * v_i0r)[-n_age_groups])
  v_s1g  <- c(0, (v_d * v_s1)[-n_age_groups])
  v_i1rg <- c(0, (v_d * v_i1r)[-n_age_groups])

  ### Transitions
  ds0_dt  <- v_s0g +
    (gamma * v_s1) -
    ((v_mu + v_d + v_lambda + v_lambda_prime + v_alpha + v_psi) * v_s0)
  di0w_dt <- v_i0wg +
    (v_lambda * v_s0) -
    ((v_mu + v_d + v_alpha + v_psi + v_lambda_xi) * v_i0w)
  di0r_dt <- v_i0rg +
    (v_lambda_prime * v_s0) +
    (v_lambda_xi * v_i0w) +
    (gamma * v_i1r) -
    ((v_alpha + v_psi + v_mu + v_d) * v_i0r)
  ds1_dt  <- v_s1g +
    ((v_alpha + v_psi) * v_s0) +
    (((v_alpha + v_psi) * (1 - sigma)) * v_i0w) -
    ((v_lambda_prime + gamma + v_mu + v_d) * v_s1)
  di1r_dt <- v_i1rg +
    ((v_alpha + v_psi) * v_i0r) +
    ((v_alpha + v_psi) * sigma * v_i0w) +
    (v_lambda_prime * v_s1) -
    ((gamma + v_mu + v_d) * v_i1r)

  ### combine results
  return(list(c(ds0_dt, di0w_dt, di0r_dt, ds1_dt, di1r_dt)))
}


#' Get and plot model results
#'
#' @return list(plot, desolver))
#' @export
#'
#' @examples
get_sis_model_results <- function(n_ages = length(groups)) {
  require(deSolve)
  require(ggplot2)
  require(dtplyr)
  require(dplyr, warn.conflicts = FALSE)

  desolver <- deSolve::ode(y = v_parameter$v_state, times = v_parameter$v_time,
                           func = sis_model, parms = v_parameter,
                           method = "ode45")
  desolver2 <- as.data.frame(desolver)

  model_results <- lazy_dt(desolver2) %>%
    group_by(time) %>%
    mutate(S0 = sum(across(1:as.numeric(n_ages)))) %>%
    mutate(I0 = sum(across(as.numeric(n_ages + 1):
                             as.numeric(2 * n_ages)))) %>%
    mutate(S1 = sum(across(as.numeric(2 * n_ages + 1):
                             as.numeric(3 * n_ages)))) %>%
    mutate(I1 = sum(across(as.numeric(3 * n_ages + 1):
                             as.numeric(4 * n_ages)))) %>%
    as_tibble()

  return(model_results)
}


#' Get and plot model results
#'
#' @return list(plot, desolver))
#' @export
#'
#' @examples
get_model_results <- function(model) {
  require(deSolve)
  require(ggplot2)

  desolver <- deSolve::lsoda(y = v_parameter$v_state,
                             times = v_parameter$v_time,
                             func = model, parms = v_parameter)

  model_results <- as.data.frame(desolver) %>%
    mutate(S = rowSums(.data[2:81])) %>%
    mutate(I = rowSums(.data[82:161]))

  plot <- ggplot() +
    geom_line(data = model_results, aes(x = .data$time, y = .data$S),
              color = "blue") +
    geom_line(data = model_results, aes(x = .data$time, y = .data$I),
              color = "red") +
    ylab("Proportion")

  return(list(plot, model_results))
}


#' Find the steady state of the SI model
#'
#' @param initial_state
#' @param params
#'
#' @return
#' @export
#'
#' @examples
find_steady <- function(initial_state = v_parameter$v_state,
                        params = v_parameter) {
  equil <- rootSolve::runsteady(y = initial_state,
                                times = c(0, 500),
                                func = SI_model,
                                parms = params)


  return(equil$y)  # return rounded steady states
}
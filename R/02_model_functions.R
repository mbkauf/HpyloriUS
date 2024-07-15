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

  # combine results
  return(list(c(ds_dt, di_dt)))

}

#' Function for SIS Model
#'
#' @param v_time  vector of times to use in model
#' @param v_state  vector of propotion in each compartment
#' @param v_params  vector of input parameters
#'
#' @return
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
  p_trt_s            <- v_params$p_trt_s
  p_trt_i            <- v_params$p_trt_i
  trt_year           <- v_params$trt_year

  # One-time intervention
  p_trt_s <- ifelse(
    floor(v_time) == trt_year, p_trt_s, 0
  )
  p_trt_i <- ifelse(
    floor(v_time) == trt_year, p_trt_i, 0
  )

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
  v_psi_0            <- v_params$v_psi_0
  v_psi_s            <- v_params$v_psi_s
  sigma              <- v_params$sigma
  phi                <- v_params$phi
  n_race             <- v_params$n_race
  v_eta              <- v_params$v_eta
  test_sens          <- v_params$test_sens
  test_spec          <- v_params$test_spec
  delta              <- v_params$delta
  age_based_test     <- v_params$age_based_test
  rho                <- v_params$rho
  trt_year           <- v_params$trt_year
  rho                <- v_params$rho
  sus_test_sens      <- v_params$sus_test_sens
  sus_test_spec      <- v_params$sus_test_spec
  sus_test_policy    <- v_params$sus_test_policy
  trt_sus_policy     <- v_params$trt_sus_policy
  one_time_policy    <- v_params$one_time_policy
  to_present         <- v_params$to_present
  burn               <- v_params$burn
  p_trans            <- v_params$p_trans

  m_state <- matrix(
    v_state,
    nrow = n_age_groups,
    ncol = n_id_states,
    dimnames = list(v_names_age_groups, v_names_id_states)
  )

  ### Policy specific variables
  # If policy is time specific (i.e. one-time mass treatment)
  if (one_time_policy == TRUE) {
    v_eta <- rep(0, length(v_eta))
    if (floor(v_time) == trt_year) {
      # If policy is to treat everyone regardless of infection status
      if (trt_sus_policy == TRUE) {
        v_psi_s <- v_psi
      }
      v_psi <- v_psi
    } else {
      v_psi <- v_psi_0
    }
  } else {
    v_psi <- v_psi_0
  }

  ### Background antibiotics usage
  if (to_present == TRUE) {
    if (floor(v_time) < 2) {
      v_alpha <- rep(0, length(v_alpha))
    } else if (floor(v_time) < 8) {
      v_alpha <- rep((0.0456 / 6) * v_time, length(v_alpha))
    } else if (floor(v_time) >= 8) {
      v_alpha <- rep(0.0456 * exp(-0.125 * v_time), length(v_alpha))
    }
    v_eta <- rep(0, length(v_eta))
  } else if (to_present == FALSE && burn == FALSE) {
    # Constant usage after 2023
    v_alpha <- rep(0.0456 * exp(-0.125 * 34), length(v_alpha))
    # Continue along exponential
    # v_alpha <- rep(0.0456 * exp(-0.125 * (v_time + 34)), length(v_alpha))
    if (floor(v_time) <= 1 && floor(v_time) > 0) {
      v_eta <- v_eta
    }
  } else if (burn == TRUE) {
    v_alpha <- rep(0, length(v_alpha))
    v_eta   <- rep(0, length(v_eta))
  }


  ### Split m_state matrix into vectors of susceptibles and infected
  v_s0  <- m_state[, "S0"]
  v_i0w <- m_state[, "I0w"]
  v_i0r <- m_state[, "I0r"]
  v_s1  <- m_state[, "S1"]
  v_i1r <- m_state[, "I1r"]

  ### Force of infection lambda values
  v_lambda       <- (p_trans * m_waifw) %*% v_i0w
  v_lambda_prime <- (p_trans * m_waifw * phi) %*% (v_i0r + v_i1r)

  ### Growth rates (Birth rate and aging rate)
  # Growth rates (Birth rate and aging rate)
  n      <- n_age_groups / n_race
  v_s0d  <- v_d * v_s0
  v_i0wd <- v_d * v_i0w
  v_i0rd <- v_d * v_i0r
  v_s1d  <- v_d * v_s1
  v_i1rd <- v_d * v_i1r

  l_v_s0d  <- split(v_s0d, ceiling(seq_along(v_s0d) / n))
  l_v_i0wd <- split(v_i0wd, ceiling(seq_along(v_i0wd) / n))
  l_v_i0rd <- split(v_i0rd, ceiling(seq_along(v_i0rd) / n))
  l_v_s1d  <- split(v_s1d, ceiling(seq_along(v_s1d) / n))
  l_v_i1rd <- split(v_i1rd, ceiling(seq_along(v_i1rd) / n))

  v_s0g  <- c()
  v_i0wg <- c()
  v_i0rg <- c()
  v_s1g  <- c()
  v_i1rg <- c()
  for (i in 1:n_race) {
    v_s0g  <- c(v_s0g, b[i], l_v_s0d[[i]][-n])
    v_i0wg <- c(v_i0wg, 0, l_v_i0wd[[i]][-n])
    v_i0rg <- c(v_i0rg, 0, l_v_i0rd[[i]][-n])
    v_s1g  <- c(v_s1g, 0, l_v_s1d[[i]][-n])
    v_i1rg <- c(v_i1rg, 0, l_v_i1rd[[i]][-n])
  }

  ### Transitions
  ds0_dt  <- v_s0g +
    (gamma * v_s1) -
    ((v_mu + v_d + v_lambda + v_lambda_prime + v_alpha + v_psi_s +
        (v_eta * (1 - test_spec))) * v_s0)
  di0w_dt <- v_i0wg +
    (v_lambda * v_s0) +
    (rho * v_i0r) -
    ((v_mu + v_d + v_psi + v_alpha + (v_eta * test_sens)) * v_i0w)
  di0r_dt <- v_i0rg +
    (v_lambda_prime * v_s0) +
    (gamma * v_i1r) -
    ((v_psi + v_alpha + v_mu + v_d + (v_eta * test_sens)) * v_i0r) #+ rho
  ds1_dt  <- v_s1g +
    ((v_alpha + v_psi_s + (v_eta * (1 - test_spec))) * v_s0) +
    (((v_psi + v_alpha) * (1 - sigma)) * v_i0w) -
    ((v_lambda_prime + gamma + v_mu + v_d) * v_s1)
  di1r_dt <- v_i1rg +
    ((v_psi + v_alpha + (v_eta * test_sens)) * v_i0r) +
    ((v_psi + v_alpha + (v_eta * test_sens)) * sigma * v_i0w) +
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
get_sis_model_results <- function(n_ages = length(groups),
                                  v_params) {
  require(deSolve)
  require(ggplot2)
  require(dtplyr)
  require(dplyr, warn.conflicts = FALSE)

  desolver <- deSolve::ode(y = v_params$v_state, times = v_params$v_time,
                           func = sis_model, parms = v_params,
                           method = "ode45")
  desolver2 <- as.data.frame(desolver)

  v_col_names <- c("time", paste0("s0_", v_params$v_names_age_groups),
                   paste0("i0_", v_params$v_names_age_groups),
                   paste0("s1_", v_params$v_names_age_groups),
                   paste0("i1_", v_params$v_names_age_groups))
  colnames(desolver2) <- v_col_names

  model_results <- lazy_dt(desolver2) %>%
    dplyr::mutate(
      s0 = dplyr::select(., starts_with("s0")) %>% rowSums(),  #nolint
      i0 = dplyr::select(., starts_with("i0")) %>% rowSums(),  #nolint
      s1 = dplyr::select(., starts_with("s1")) %>% rowSums(),  #nolint
      i1 = dplyr::select(., starts_with("i1")) %>% rowSums()  #nolint
    ) %>%
    dplyr::mutate(all = rowSums(dplyr::select(., .data$s0:.data$i1))) %>%
    as_tibble()

  return(model_results)
  # return(desolver2)
}


#' Get and plot model results
#'
#' @return list(plot, desolver))
#' @export
#'
#' @examples
get_si_model_results <- function(model, v_params) {
  require(deSolve)
  require(ggplot2)
  require(dtplyr)
  require(dplyr, warn.conflicts = FALSE)

  desolver <- deSolve::lsoda(y = v_params$v_state,
                             times = v_params$v_time,
                             func = model, parms = v_params)

  model_results <- as.data.frame(desolver)

  v_col_names <- c("time", paste0("s_", v_params$v_names_age_groups),
                   paste0("i_", v_params$v_names_age_groups))
  colnames(model_results) <- v_col_names

  model_results <- lazy_dt(model_results) %>%
    dplyr::mutate(
      s = dplyr::select(., starts_with("s")) %>% rowSums(),  #nolint
      i = dplyr::select(., starts_with("i")) %>% rowSums()   #nolint
    ) %>%
    dplyr::mutate(all = rowSums(dplyr::select(., .data$s:.data$i))) %>%
    as_tibble()

  return(model_results)
}


#' Function for SIS Model
#'
#' @param v_time  vector of times to use in model
#' @param v_state  vector of propotion in each compartment
#' @param v_params  vector of input parameters
#'
#' @return
#'
#' @examples
sis_model_all <- function(v_time,
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
  p_trt_s            <- v_params$p_trt_s
  p_trt_i            <- v_params$p_trt_i
  trt_year           <- v_params$trt_year
  n_race             <- v_params$n_race

  # One-time intervention
  p_trt_s <- ifelse(
    floor(v_time) == trt_year, p_trt_s, 0
  )
  p_trt_i <- ifelse(
    floor(v_time) == trt_year, p_trt_i, 0
  )

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

  # Growth rates (Birth rate and aging rate)
  n <- n_age_groups / n_race
  v_s0d <- v_d * v_s0
  v_i0d <- v_d * v_i0
  v_s1d <- v_d * v_s1
  v_i1d <- v_d * v_i1

  l_v_s0d <- split(v_s0d, ceiling(seq_along(v_s0d) / n))
  l_v_i0d <- split(v_i0d, ceiling(seq_along(v_i0d) / n))
  l_v_s1d <- split(v_s1d, ceiling(seq_along(v_s1d) / n))
  l_v_i1d <- split(v_i1d, ceiling(seq_along(v_i1d) / n))

  v_s0g <- c()
  v_i0g <- c()
  v_s1g <- c()
  v_i1g <- c()
  for (i in 1:n_race) {
    v_s0g <- c(v_s0g, b[i], l_v_s0d[[i]][-n])
    v_i0g <- c(v_i0g, 0, l_v_i0d[[i]][-n])
    v_s1g <- c(v_s1g, 0, l_v_s1d[[i]][-n])
    v_i1g <- c(v_i1g, 0, l_v_i1d[[i]][-n])
  }

  # v_s0g <- c(b, (v_d * v_s0)[-n_age_groups])
  # v_i0g <- c(0, (v_d * v_i0)[-n_age_groups])
  # v_s1g <- c(0, (v_d * v_s1)[-n_age_groups])
  # v_i1g <- c(0, (v_d * v_i1)[-n_age_groups])

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

get_sis_model_results_all <- function(n_ages = length(groups),
                                      v_params, v_race_names) {
  require(deSolve)
  require(ggplot2)
  require(dtplyr)
  require(dplyr, warn.conflicts = FALSE)

  desolver <- deSolve::ode(y = v_params$v_state, times = v_params$v_time,
                           func = sis_model_all, parms = v_params,
                           method = "ode45")
  desolver2 <- as.data.frame(desolver)

  v_state_names <- c("s0", "i0", "s1", "i1")
  n <- v_params$n_age_groups / v_params$n_race
  v_names <- as.character(1:n)
  v_col_names <- c("time")
  for (state in v_state_names) {
    for (race in v_race_names) {
      v_col_names <- c(v_col_names, paste0(state, "_", v_names, "_", race))
    }
  }

  colnames(desolver2) <- v_col_names

  model_results <- lazy_dt(desolver2) %>%
    dplyr::mutate(
      s0_hisp  = dplyr::select(., intersect(starts_with("s0"), ends_with("hisp")))  %>% rowSums(),  #nolint
      i0_hisp  = dplyr::select(., intersect(starts_with("i0"), ends_with("hisp")))  %>% rowSums(),  #nolint
      s1_hisp  = dplyr::select(., intersect(starts_with("s1"), ends_with("hisp")))  %>% rowSums(),  #nolint
      i1_hisp  = dplyr::select(., intersect(starts_with("i1"), ends_with("hisp")))  %>% rowSums(),  #nolint
      s0_black = dplyr::select(., intersect(starts_with("s0"), ends_with("black"))) %>% rowSums(),  #nolint
      i0_black = dplyr::select(., intersect(starts_with("i0"), ends_with("black"))) %>% rowSums(),  #nolint
      s1_black = dplyr::select(., intersect(starts_with("s1"), ends_with("black"))) %>% rowSums(),  #nolint
      i1_black = dplyr::select(., intersect(starts_with("i1"), ends_with("black"))) %>% rowSums(),  #nolint
      s0_white = dplyr::select(., intersect(starts_with("s0"), ends_with("white"))) %>% rowSums(),  #nolint
      i0_white = dplyr::select(., intersect(starts_with("i0"), ends_with("white"))) %>% rowSums(),  #nolint
      s1_white = dplyr::select(., intersect(starts_with("s1"), ends_with("white"))) %>% rowSums(),  #nolint
      i1_white = dplyr::select(., intersect(starts_with("i1"), ends_with("white"))) %>% rowSums(),  #nolint
      s0_other = dplyr::select(., intersect(starts_with("s0"), ends_with("other"))) %>% rowSums(),  #nolint
      i0_other = dplyr::select(., intersect(starts_with("i0"), ends_with("other"))) %>% rowSums(),  #nolint
      s1_other = dplyr::select(., intersect(starts_with("s1"), ends_with("other"))) %>% rowSums(),  #nolint
      i1_other = dplyr::select(., intersect(starts_with("i1"), ends_with("other"))) %>% rowSums()  #nolint
    ) %>%
    dplyr::mutate(
      all_hisp  = rowSums(dplyr::select(., .data$s0_hisp:.data$i1_hisp)),
      all_black = rowSums(dplyr::select(., .data$s0_black:.data$i1_black)),
      all_white = rowSums(dplyr::select(., .data$s0_white:.data$i1_white)),
      all_other = rowSums(dplyr::select(., .data$s0_other:.data$i1_other)),
    ) %>%
    as_tibble()

  return(model_results)
}

get_sis_abr_model_results_all <- function(n_ages = length(groups),
                                          v_params, v_race_names) {
  require(deSolve)
  require(ggplot2)
  require(dtplyr)
  require(dplyr, warn.conflicts = FALSE)

  desolver <- deSolve::ode(y = v_params$v_state, times = v_params$v_time,
                           func = sis_abr_model, parms = v_params,
                           method = "ode45")
  desolver2 <- as.data.frame(desolver)

  v_state_names <- v_params$v_names_id_states
  n <- v_params$n_age_groups / v_params$n_race
  v_names <- as.character(1:n)
  v_col_names <- c("time")
  for (state in v_state_names) {
    for (race in v_race_names) {
      v_col_names <- c(v_col_names, paste0(state, "_", v_names, "_", race))
    }
  }

  colnames(desolver2) <- v_col_names

  model_results <- lazy_dt(desolver2) %>%
    dplyr::mutate(
      s0_hisp   = dplyr::select(., intersect(starts_with("s0"),
                                             ends_with("hisp")))  %>% rowSums(),  #nolint
      i0w_hisp  = dplyr::select(., intersect(starts_with("i0w"),
                                             ends_with("hisp")))  %>% rowSums(),  #nolint
      i0r_hisp  = dplyr::select(., intersect(starts_with("i0r"),
                                             ends_with("hisp")))  %>% rowSums(),  #nolint
      s1_hisp   = dplyr::select(., intersect(starts_with("s1"),
                                             ends_with("hisp")))  %>% rowSums(),  #nolint
      i1r_hisp  = dplyr::select(., intersect(starts_with("i1r"),
                                             ends_with("hisp")))  %>% rowSums(),  #nolint
      s0_black  = dplyr::select(., intersect(starts_with("s0"),
                                             ends_with("black"))) %>% rowSums(),  #nolint
      i0w_black = dplyr::select(., intersect(starts_with("i0w"),
                                             ends_with("black"))) %>% rowSums(),  #nolint
      i0r_black = dplyr::select(., intersect(starts_with("i0r"),
                                             ends_with("black"))) %>% rowSums(),  #nolint
      s1_black  = dplyr::select(., intersect(starts_with("s1"),
                                             ends_with("black"))) %>% rowSums(),  #nolint
      i1r_black = dplyr::select(., intersect(starts_with("i1r"),
                                             ends_with("black"))) %>% rowSums(),  #nolint
      s0_white  = dplyr::select(., intersect(starts_with("s0"),
                                             ends_with("white"))) %>% rowSums(),  #nolint
      i0w_white = dplyr::select(., intersect(starts_with("i0w"),
                                             ends_with("white"))) %>% rowSums(),  #nolint
      i0r_white = dplyr::select(., intersect(starts_with("i0r"),
                                             ends_with("white"))) %>% rowSums(),  #nolint
      s1_white  = dplyr::select(., intersect(starts_with("s1"),
                                             ends_with("white"))) %>% rowSums(),  #nolint
      i1r_white = dplyr::select(., intersect(starts_with("i1r"),
                                             ends_with("white"))) %>% rowSums(),  #nolint
      s0_other  = dplyr::select(., intersect(starts_with("s0"),
                                             ends_with("other"))) %>% rowSums(),  #nolint
      i0w_other = dplyr::select(., intersect(starts_with("i0w"),
                                             ends_with("other"))) %>% rowSums(),  #nolint
      i0r_other = dplyr::select(., intersect(starts_with("i0r"),
                                             ends_with("other"))) %>% rowSums(),  #nolint
      s1_other  = dplyr::select(., intersect(starts_with("s1"),
                                             ends_with("other"))) %>% rowSums(),  #nolint
      i1r_other = dplyr::select(., intersect(starts_with("i1r"),
                                             ends_with("other"))) %>% rowSums()  #nolint
    ) %>%
    dplyr::mutate(
      all_hisp  = rowSums(dplyr::select(., .data$s0_hisp:.data$i1r_hisp)),
      all_black = rowSums(dplyr::select(., .data$s0_black:.data$i1r_black)),
      all_white = rowSums(dplyr::select(., .data$s0_white:.data$i1r_white)),
      all_other = rowSums(dplyr::select(., .data$s0_other:.data$i1r_other)),
    ) %>%
    as_tibble()

  return(model_results)
}


sis_abr_sus_model <- function(v_time,
                              v_state,
                              v_params) {

  ### Assign values from parameters
  m_waifw            <- v_params$m_waifw
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
  v_psi_0            <- v_params$v_psi_0
  v_psi_s            <- v_params$v_psi_s
  sigma              <- v_params$sigma
  phi                <- v_params$phi
  n_race             <- v_params$n_race
  v_eta              <- v_params$v_eta
  test_sens          <- v_params$test_sens
  test_spec          <- v_params$test_spec
  delta              <- v_params$delta
  age_based_test     <- v_params$age_based_test
  rho                <- v_params$rho
  trt_year           <- v_params$trt_year
  rho                <- v_params$rho
  sus_test_sens      <- v_params$sus_test_sens
  sus_test_spec      <- v_params$sus_test_spec
  sus_test_policy    <- v_params$sus_test_policy
  trt_sus_policy     <- v_params$trt_sus_policy
  one_time_policy    <- v_params$one_time_policy
  to_present         <- v_params$to_present
  burn               <- v_params$burn
  p_trans            <- v_params$p_trans

  m_state <- matrix(
    v_state,
    nrow = n_age_groups,
    ncol = n_id_states,
    dimnames = list(v_names_age_groups, v_names_id_states)
  )

  ### Policy specific variables
  # If policy is time specific (i.e. one-time mass treatment)
  if (one_time_policy == TRUE) {
    v_eta <- rep(0, length(v_eta))
    if (floor(v_time) == trt_year) {
      # If policy is to treat everyone regardless of infection status
      if (trt_sus_policy == TRUE) {
        v_psi_s <- v_psi
      }
      v_psi <- v_psi
    } else {
      v_psi <- v_psi_0
    }
  } else {
    v_psi <- v_psi_0
  }

  ### Background antibiotics usage
  if (to_present == TRUE) {
    if (floor(v_time) < 2) {
      v_alpha <- rep(0, length(v_alpha))
    } else if (floor(v_time) < 8) {
      v_alpha <- rep((0.0456 / 6) * v_time, length(v_alpha))
    } else if (floor(v_time) >= 8) {
      v_alpha <- rep(0.0456 * exp(-0.125 * v_time), length(v_alpha))
    }
    v_eta <- rep(0, length(v_eta))
  } else if (to_present == FALSE && burn == FALSE) {
    # Constant usage after 2023
    v_alpha <- rep(0.0456 * exp(-0.125 * 34), length(v_alpha))
    # Continue along exponential
    # v_alpha <- rep(0.0456 * exp(-0.125 * (v_time + 34)), length(v_alpha))
    if (floor(v_time) <= 1 && floor(v_time) > 0) {
      v_eta <- v_eta
    }
  } else if (burn == TRUE) {
    v_alpha <- rep(0, length(v_alpha))
    v_eta   <- rep(0, length(v_eta))
  }


  ### Split m_state matrix into vectors of susceptibles and infected
  v_s0  <- m_state[, "S0"]
  v_i0w <- m_state[, "I0w"]
  v_i0r <- m_state[, "I0r"]
  v_s1  <- m_state[, "S1"]
  v_i1r <- m_state[, "I1r"]

  ### Force of infection lambda values
  v_lambda       <- (p_trans * m_waifw) %*% v_i0w
  v_lambda_prime <- (p_trans * m_waifw * phi) %*% (v_i0r + v_i1r)

  ### Growth rates (Birth rate and aging rate)
  # Growth rates (Birth rate and aging rate)
  n      <- n_age_groups / n_race
  v_s0d  <- v_d * v_s0
  v_i0wd <- v_d * v_i0w
  v_i0rd <- v_d * v_i0r
  v_s1d  <- v_d * v_s1
  v_i1rd <- v_d * v_i1r

  l_v_s0d  <- split(v_s0d, ceiling(seq_along(v_s0d) / n))
  l_v_i0wd <- split(v_i0wd, ceiling(seq_along(v_i0wd) / n))
  l_v_i0rd <- split(v_i0rd, ceiling(seq_along(v_i0rd) / n))
  l_v_s1d  <- split(v_s1d, ceiling(seq_along(v_s1d) / n))
  l_v_i1rd <- split(v_i1rd, ceiling(seq_along(v_i1rd) / n))

  v_s0g  <- c()
  v_i0wg <- c()
  v_i0rg <- c()
  v_s1g  <- c()
  v_i1rg <- c()
  for (i in 1:n_race) {
    v_s0g  <- c(v_s0g, b[i], l_v_s0d[[i]][-n])
    v_i0wg <- c(v_i0wg, 0, l_v_i0wd[[i]][-n])
    v_i0rg <- c(v_i0rg, 0, l_v_i0rd[[i]][-n])
    v_s1g  <- c(v_s1g, 0, l_v_s1d[[i]][-n])
    v_i1rg <- c(v_i1rg, 0, l_v_i1rd[[i]][-n])
  }

  ### Transitions
  ds0_dt  <- v_s0g +
    (gamma * v_s1) -
    ((v_mu + v_d + v_lambda + v_lambda_prime + v_alpha + v_psi_s +
        (v_eta * (1 - test_spec))) * v_s0)
  di0w_dt <- v_i0wg +
    (v_lambda * v_s0) +
    (rho * v_i0r) -
    ((v_mu + v_d + v_psi + v_alpha + (v_eta * test_sens)) * v_i0w)
  di0r_dt <- v_i0rg +
    (v_lambda_prime * v_s0) +
    (gamma * v_i1r) -
    ((v_psi + v_alpha + v_mu + v_d + (v_eta * test_sens)) * v_i0r) #+ rho
  ds1_dt  <- v_s1g +
    ((v_alpha + v_psi_s + (v_eta * (1 - test_spec))) * v_s0) +
    (((v_psi + v_alpha) * (1 - sigma)) * v_i0w) -
    ((v_lambda_prime + gamma + v_mu + v_d) * v_s1)
  di1r_dt <- v_i1rg +
    ((v_psi + v_alpha + (v_eta * test_sens)) * v_i0r) +
    ((v_psi + v_alpha + (v_eta * test_sens)) * sigma * v_i0w) +
    (v_lambda_prime * v_s1) -
    ((gamma + v_mu + v_d) * v_i1r)

  ### combine results
  return(list(c(ds0_dt, di0w_dt, di0r_dt, ds1_dt, di1r_dt)))
}

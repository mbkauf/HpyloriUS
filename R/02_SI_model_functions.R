#' Function for SI model
#'
#' @param v_time  vector of times to use in model
#' @param v_state  vector of propotion in each compartment
#' @param v_params  vector of input parameters
#'
#' @return list(c(dS_dt, dI_dt))
#' @export
#'
#' @examples
#' deSolve::lsoda(y=v_parameter$v_state, times = v_parameter$v_time,
#'                func = SI_model, parms = v_parameter)
SI_model <- function(v_time,
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
  v_S <- m_state[, "S"]
  v_I <- m_state[, "I"]

  # Force of infection lambda values
  v_lambda <- m_waifw %*% v_I

  # Growth rates (Birth rate and aging rate)
  v_Sg <- c(b, (v_d * v_S)[-n_age_groups])
  v_Ig <- c(0, (v_d * v_I)[-n_age_groups])

  # Transitions
  dS_dt = v_Sg - (v_mu + v_d + v_lambda) * v_S
  dI_dt = v_Ig + (v_lambda * v_S) - (v_mu + v_d) * v_I

  # combine results
  return(list(c(dS_dt, dI_dt)))

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

  desolver <- deSolve::lsoda(y=v_parameter$v_state, times = v_parameter$v_time,
                           func = model, parms = v_parameter)

  model_results <- as.data.frame(desolver) %>%
    mutate(S = rowSums(.[2:82])) %>%
    mutate(I = rowSums(.[83:163]))

  plot <- ggplot() +
    geom_line(data = model_results, aes(x = time, y = S), color = "blue") +
    geom_line(data = model_results, aes(x = time, y = I), color = "red") +
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
  equil <- rootSolve::runsteady(y=initial_state,
                                times=c(0,1E8),
                                func=SI_model,
                                parms=params)


  return(equil$y)  # return rounded steady states
}



#' Title
#'
#' @return
#' @export
#'
#' @examples
get_time_solver <- function() {
  timer <- microbenchmark::microbenchmark(
    desolve <- deSolve::ode(l_state, time, SI_model, l_parameter),
    # odin    <-
  )
}



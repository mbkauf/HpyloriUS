#' Generate WAIFW matrix
#'
#' @param betas       vector of beta values for waifw matrix
#' @param breaks      breaks (should be length of age groups + 1) ?
#' @param w           contact matrix shell for betas (frequently upper
#'                    triangular index matrix).
#' @param upper       TRUE/FALSE about w being an upper triangular index matrix
#' @param group_names vector of numbers for age range (used for naming matrix)
#'
#' @return waifw
#' @export
#'
#' @examples Beta <- get_transmission_matrix(betas = beta0,
#'                                           breaks = v_waifw_breaks,
#'                                           w = test_w,
#'                                           group_names = groups)
get_transmission_matrix <- function(betas, breaks, w, upper = TRUE,
                                    group_names = groups) {

  ### Number of age groups = number of data to estimate
  n_params <- n_a <- (length(breaks) - 1)

  ### Max age
  # age_max <- breaks[length(breaks)] - 1 nolint

  ### Age group labels
  age_names <- paste(breaks[-length(breaks)], breaks[-1] - 1, sep = "-")

  ### Check that we have 8 betas. For complete assortativity, add a zero
  # if (n_params != 8) {
  #   stop("Number of beta parameters is different than 8")
  # }
  # if (max(w) > n_params) {
  #   betas <- c(betas, 0)
  # }
  ## Complete index matrix from upper triangular index matrix
  complete_index <- function(index_upper) {
    index <- index_upper + t(index_upper)
    diag(index) <- diag(index_upper)
    return(index)
  }

  if (upper == TRUE) {
    index <- complete_index(index_upper = w)
    ### Generate 1-year age-specific waifw
    index_full <- index[rep(1:n_a, diff(breaks)),
                        rep(1:n_a, diff(breaks))]
    colnames(index_full) <- rep(age_names, diff(breaks))
    m_waifw <- matrix(betas[index_full], ncol = (length(group_names)))
  } else {
    index <- w
    ### Generate 1-year age-specific waifw
    index_full <- index[rep(1:n_a, diff(breaks)),
                        rep(1:n_a, diff(breaks))]
    colnames(index_full) <- rep(age_names, diff(breaks))
    m_waifw <- matrix(betas[index_full], ncol = (length(group_names)))
  }

  m_waifw[is.na(m_waifw)] <- 0  # replace NAs with 0
  colnames(m_waifw) <- rownames(m_waifw) <- group_names  # add names
  return(m_waifw = m_waifw)  # return WAIFW matrix
}


#' Function to generate demographic variables needed for SI model
#'
#' @param life_table_data  life tables that contain mortality rates
#' @param spec_groups  vector of ages used in model
#' @param r_pop_growth  rate of population growth
#'
#' @return list(b, v_age, v_age_prop, v_mu, v_d)
#' @export
#'
#' @examples v_demo <- get_demographic_vars(spec_groups = groups)
get_demographic_vars <- function(spec_groups = NULL,
                                 r_pop_growth = 0,
                                 race = c("Hispanic",
                                          "NH White",
                                          "NH Black",
                                          "NH AAPI",
                                          "NH AIAN",
                                          "Other")) {
  ### Check that an appropriate vector of ages were entered --------------------
  if (!is.null(spec_groups)) {
    n_groups <- length(spec_groups)
    if (!is.numeric(spec_groups)) {
      spec_groups <- as.numeric(spec_groups)
      warning("spec_groups was not entered as a numeric vector, so it was
              converted to a numeric vector. Results may not be correct")
    }
    older_upper_lim <- max(spec_groups)
    lower_age_lim   <- min(spec_groups)
  }

  ### Calculate death rates ----------------------------------------------------
  # Select appropriate life-table
  if (race == "Hispanic") {
    life_table_data <- "LifeTable_US_Hispanic_All.csv"
  } else if (race == "NH White") {
    life_table_data <- "LifeTable_US_White_All.csv"
  } else if (race == "NH Black") {
    life_table_data <- "LifeTable_US_Black_All.csv"
  } else if (race == "NH AAPI") {
    life_table_data <- "LifeTable_US_Asian_All.csv"
  } else if (race == "NH AIAN") {
    life_table_data <- "LifeTable_US_AIAN_All.csv"
  } else if (race == "Other") {
    life_table_data <- "LifeTable_US_All.csv"
  }

  # Import life-table
  temp <- read.csv(paste0("data/", life_table_data))

  # Rates of death for ages of interest
  v_mu        <- as.numeric(matrix(temp[lower_age_lim:older_upper_lim, 2]))
  names(v_mu) <- as.character(temp[lower_age_lim:older_upper_lim, 1])

  ### Calculate aging rate -----------------------------------------------------
  # Get age character vector for naming purposes
  v_age <- as.character(spec_groups)

  #initialize age_proportions vector
  v_age_prop <- rep(0, n_groups)

  v_d        <- (v_mu + r_pop_growth) / (exp((v_mu + r_pop_growth)) - 1)
  v_d        <- c(v_d[-length(v_d)], 0)
  names(v_d) <- names(v_age_prop) <- names(v_mu) <- v_age

  ### Calculate proportion of population in each age class ---------------------
  # Using Merck pg. 10
  # Explicit translation of first formula on pg. 10
  # Calculating youngest age group
  acc <- 0
  inner_acc <- 1
  # inner product/sum
  for (rr in 2:n_groups) {
    for (jj in 2:rr) {
      inner_acc <- (v_d[v_age[jj - 1]] / (v_d[v_age[jj]] + v_mu[v_age[jj]])) *
        inner_acc
    }
    acc <- acc + inner_acc
    inner_acc <- 1
  }

  #calculate proportion in first age group
  v_age_prop[1] <- 1 / (1 + acc)

  for (AG in v_age[-1]) {
    prev_age <- v_age[which(v_age == AG) - 1]
    v_age_prop[AG] <- v_d[prev_age] * v_age_prop[prev_age] /
      (v_d[AG] + v_mu[AG] + r_pop_growth)
  }


  ### Calculate birth rate -----------------------------------------------------
  b <- (v_d[1] + v_mu[1] + r_pop_growth) * v_age_prop[1]

  return(list(b = b,
              v_age = v_age,
              v_age_prop = v_age_prop,
              v_mu = v_mu,
              v_d = v_d))
}

#' Function to generate prevalence varaiables needed for SI model
#'
#' @param data  prevalence data set
#' @param spec_groups  vector of age groups
#'
#' @return list(prevalence = prevalence, prevalence_se = prevalence_se)
#' @export
#'
#' @examples v_prev_vars <- get_prevalence_vars(spec_groups = groups)
get_prevalence_vars <- function(data_foi = "df_foi_period.csv",
                                data_prev = "df_prev_period.csv",
                                spec_groups = NULL,
                                race = c("Hispanic",
                                         "NH White",
                                         "NH Black",
                                         "NH AAPI",
                                         "NH AIAN",
                                         "Other")) {
  require(dtplyr)
  require(dplyr)
  ### Check that an appropriate vector of ages were entered --------------------
  if (!is.numeric(spec_groups)) {
    spec_groups <- as.numeric(spec_groups)
    warning("spec_groups was not entered as a numeric vector, so it was
            converted to a numeric vector. Results may not be correct")
  }

  ### Load prevalence/foi dataset ----------------------------------------------
  df_foi  <- read.csv(paste0("data/", data_foi))
  df_prev <- read.csv(paste0("data/", data_prev))

  # Perform the following actions to the dataset:
  # 1) filter to ages used in model
  # 2) filter to specified birth cohort
  # 3) filter to specified race group
  # 4) get prevalence value [0, 1]
  df_foi <- lazy_dt(df_foi) %>%
    dplyr::rename(race_eth = race) %>%
    dplyr::filter(.data$age %in% spec_groups) %>%
    dplyr::filter(.data$race_eth == race) %>%
    dplyr::filter(.data$period == 1) %>%
    as_tibble()

  foi    <- as.numeric(df_foi$Mean)
  foi_sd <- as.numeric(df_foi$SD)
  foi_lb <- as.numeric(df_foi$LB)
  foi_ub <- as.numeric(df_foi$UB)

  df_prev <- lazy_dt(df_prev) %>%
    dplyr::rename(race_eth = race) %>%
    dplyr::filter(.data$age %in% spec_groups) %>%
    dplyr::filter(.data$race_eth == race) %>%
    dplyr::filter(.data$period == 1) %>%
    as_tibble()

  prevalence    <- as.numeric(df_prev$Mean)
  prevalence_sd <- as.numeric(df_prev$SD)

  return(list(prevalence = prevalence,
              prevalence_sd = prevalence_sd,
              foi = foi,
              foi_sd = foi_sd,
              foi_lb = foi_lb,
              foi_ub = foi_ub))
}

#' Load the SI model parameters
#'
#' @param waifw  WAIFW matrix
#' @param demography_vars  vector of demography variables
#' @param prevalence_vars  vector of prevalence varaiables
#' @param ages  vector of ages used in model
#'
#' @return list(v_parameter)
#' @export
#'
#' @examples  v_parameter <- load_SI_model_params(waifw = Beta,
#'                                                demography_vars = v_demo,
#'                                                prevalence_vars = v_prev_vars,
#'                                                ages = groups)
load_si_model_params <- function(waifw, demography_vars,
                                 prevalence_vars, ages) {

  ## Get initial states
  # Susceptibles
  sus <- (1 - prevalence_vars$prevalence) * demography_vars$v_age_prop
  # Infected
  inf <- prevalence_vars$prevalence * demography_vars$v_age_prop
  v_init <- c(sus, inf)  # combine susceptibles and infected

  v_parameter <- list(
    v_time             = seq(0, 100, by = 1),  # time steps to test
    N                  = 100000,               # population size
    b                  = demography_vars$b,    # birth rate
    v_d                = demography_vars$v_d,  # aging rates
    v_mu               = demography_vars$v_mu, # death rates
    m_waifw            = waifw,                # transmission rates
    v_state            = v_init,               # list of initial state values
    v_names_age_groups = as.character(ages),   # ages as character vector
    n_id_states        = length(c("S", "I")),  # number of states
    v_names_id_states  = c("S", "I"),          # character vector of states
    n_age_groups       = length(ages)          # number of age groups
  )

  return(v_parameter)
}

#' Load SIS model parameters
#'
#' @param waifw  WAIFW matrix
#' @param demography_vars  vector of demography variables
#' @param prevalence_vars  vector of prevalence varaiables
#' @param ages  vector of ages used in model
#'
#' @return list(v_parameter)
#' @export
#'
#' @examples v_parameter <- load_SIS_model_params(waifw = Beta,
#'                                                demography_vars = v_demo,
#'                                                prevalence_vars = v_prev_vars,
#'                                                ages = groups)
load_sis_model_params <- function(waifw, demography_vars, prevalence_vars,
                                  ages, trt_year = 1, trt_eff = 1,
                                  end_t = 100, burn_state = FALSE,
                                  v_init_state = NULL) {
  ## Get initial state
  if (burn_state == FALSE) {
    # Susceptibles
    sus <- (1 - prevalence_vars$prevalence) * demography_vars$v_age_prop
    # Infected
    inf <- prevalence_vars$prevalence * demography_vars$v_age_prop
    # Combine susceptibles and infected
    v_init <- c(sus, inf, sus * 0, inf * 0)
  } else {
    v_init <- v_init_state
  }

  v_parameter <- list(
    v_time             = seq(0, end_t, by = 1),  # time steps to test
    N                  = 1,                    # population size
    b                  = demography_vars$b,    # birth rate
    v_d                = demography_vars$v_d,  # aging rates
    v_mu               = demography_vars$v_mu, # death rates
    m_waifw            = waifw,                # transmission rates
    v_state            = v_init,               # list of initial state values
    v_names_age_groups = as.character(ages),   # ages as character vector
    n_id_states        = length(c("S0", "I0", "S1", "I1")),  # number of states
    v_names_id_states  = c("S0", "I0", "S1", "I1"), # character vector of states
    n_age_groups       = length(ages),         # number of age groups
    p_test             = 0, # probability of testing infected
    test_sens          = 0, # test sensitivity
    trt_eff            = trt_eff, # probability of treatment working
    gamma              = 1 / (14 / 365), # 1 / treatment length
    p_trt_s            = 1, # probability of treatment when susceptible
    p_trt_i            = 1, # probability of treatment when infected
    trt_year           = trt_year  # year of treatment
  )

  return(v_parameter)
}


load_sis_abr_model_params <- function(waifw, waifw_prime, waifw_xi,
                                      demography_vars, prevalence_vars,
                                      ab_vars, ages) {

  ## Get initial state
  # Susceptibles
  sus <- (1 - prevalence_vars$prevalence) * demography_vars$v_age_prop
  # Infected
  inf <- prevalence_vars$prevalence * demography_vars$v_age_prop
  # Combine susceptibles and infected
  v_init <- cbind(sus * 0.9, inf * 0.5, inf * 0.4, sus * 0.1, inf * 0.1)
  # rep(0, length(ages)), rep(0, length(ages)))

  v_parameter <- list(
    v_time             = seq(0, 100, by = 1),  # time steps to test
    N                  = 1,                    # population size
    v_state            = v_init,               # list of initial state values
    v_names_age_groups = as.character(ages),   # ages as character vector
    # number of states
    n_id_states        = length(c("S0", "I0w", "I0r", "S1", "I1r")),
    # character vector of states
    v_names_id_states  = c("S0", "I0w", "I0r", "S1", "I1r"),
    n_age_groups       = length(ages),         # number of age groups
    b                  = demography_vars$b,    # birth rate
    v_d                = demography_vars$v_d,  # aging rates
    v_mu               = demography_vars$v_mu, # death rates
    # transmission rates (sensitive strains)
    m_waifw            = waifw,
    gamma              = 1 / (14 / 365),       # 1 / treatment length (14 days)
    # background antibiotic use
    v_alpha            = rep(0.05, n_age_groups), # ab_vars$v_alpha
    v_psi              = rep(0.05, n_age_groups), # assume no antibiotic policy
    # probability treatment induces sensitive strains
    sigma              = 0.1,
    phi                = 1  #
  )

  return(v_parameter)
}

get_ab_vars <- function(spec_groups = NULL,
                        race = c("Hispanic",
                                 "NH White",
                                 "NH Black",
                                 "NH AAPI",
                                 "NH AIAN",
                                 "Other")) {

  # TODO - write antibiotic vars function
  v_alpha <- rep()
}


get_demographic_vars_all <- function(ages,
                                     v_race = c("Hispanic",
                                                "NH White",
                                                "NH Black",
                                                "NH AAPI",
                                                "NH AIAN",
                                                "Other")) {
  n_race <- length(v_race)
  v_demo <- list()
  i <- 0
  for (race in v_race) {
    i <- i + 1
    if (race == "NH White") {
      race_text <- "white"
    } else if (race == "NH Black") {
      race_text <- "black"
    } else if (race == "Hispanic") {
      race_text <- "hispanic"
    } else if (race == "Other") {
      race_text <- "other"
    } else if (race == "NH AAPI") {
      race_text <- "aapi"
    } else if (race == "NH AIAN") {
      race_text <- "aian"
    }

    tmp_demo <- get_demographic_vars(spec_groups = ages,
                                     race = race)
    if (i == 1) {
      v_demo <- tmp_demo
    } else {
      v_demo <- lapply(seq_along(tmp_demo), function(x) c(v_demo[[x]],
                                                          tmp_demo[[x]]))
    }
  }
  names(v_demo) <- c("b", "v_age", "v_age_prop", "v_mu", "v_d")
  return(v_demo)
}

get_prevalence_vars_all <- function(ages,
                                    v_race = c("Hispanic",
                                               "NH White",
                                               "NH Black",
                                               "NH AAPI",
                                               "NH AIAN",
                                               "Other")) {
  n_race <- length(v_race)
  v_prev <- list()
  i <- 0
  for (race in v_race) {
    i <- i + 1
    if (race == "NH White") {
      race_text <- "white"
    } else if (race == "NH Black") {
      race_text <- "black"
    } else if (race == "Hispanic") {
      race_text <- "hispanic"
    } else if (race == "Other") {
      race_text <- "other"
    } else if (race == "NH AAPI") {
      race_text <- "aapi"
    } else if (race == "NH AIAN") {
      race_text <- "aian"
    }

    tmp_prev <- get_prevalence_vars(spec_groups = ages,
                                    race = race)
    if (i == 1) {
      v_prev <- tmp_prev
    } else {
      v_prev <- lapply(seq_along(tmp_prev), function(x) c(v_prev[[x]],
                                                          tmp_prev[[x]]))
    }
  }
  names(v_prev) <- c("prevalence", "prevalence_sd", "foi",
                     "foi_sd", "foi_lb", "foi_ub")
  return(v_prev)
}

load_sis_model_params_all <- function(
    v_race = c("Hispanic", "NH White", "NH Black",
               "NH AAPI", "NH AIAN", "Other"),
    waifw, ages, trt_year = 1, trt_eff = 1,
    end_t = 100, burn_state = FALSE,
    v_init_state = NULL) {

  n_race <- length(v_race)
  v_demo <- list()
  v_prev <- list()
  i <- 0
  for (race in v_race) {
    i <- i + 1
    if (race == "NH White") {
      race_text <- "white"
      pop <- 197.6
    } else if (race == "NH Black") {
      race_text <- "black"
      pop <- 40.1
    } else if (race == "Hispanic") {
      race_text <- "hispanic"
      pop <- 62.6
    } else if (race == "Other") {
      race_text <- "other"
    } else if (race == "NH AAPI") {
      race_text <- "aapi"
      pop <- 22.4
    } else if (race == "NH AIAN") {
      race_text <- "aian"
      pop <- 9.7
    }

    tmp_demo <- get_demographic_vars(spec_groups = ages,
                                     race = race)
    tmp_prev <- get_prevalence_vars(spec_groups = ages,
                                    race = race)
    if (i == 1) {
      v_demo <- tmp_demo
      v_prev <- tmp_prev
    } else {
      v_demo <- lapply(seq_along(tmp_demo), function(x) c(v_demo[[x]],
                                                          tmp_demo[[x]]))
      v_prev <- lapply(seq_along(tmp_prev), function(x) c(v_prev[[x]],
                                                          tmp_prev[[x]]))
    }

  }

  names(v_demo) <- c("b", "v_age", "v_age_prop", "v_mu", "v_d")
  names(v_prev) <- c("prevalence", "prevalence_sd", "foi",
                     "foi_sd", "foi_lb", "foi_ub")
  ## Get initial state
  if (burn_state == FALSE) {
    # Susceptibles
    sus <- (1 - v_prev$prevalence) * v_demo$v_age_prop
    # Infected
    inf <- v_prev$prevalence * v_demo$v_age_prop
    # Combine susceptibles and infected
    v_init <- c(sus, inf, sus * 0, inf * 0)
  } else {
    v_init <- v_init_state
  }

  v_parameter <- list(
    v_time             = seq(0, end_t, by = 1),  # time steps to test
    N                  = 1,                    # population size
    b                  = v_demo$b,    # birth rate
    v_d                = v_demo$v_d,  # aging rates
    v_mu               = v_demo$v_mu, # death rates
    m_waifw            = waifw,                # transmission rates
    v_state            = v_init,               # list of initial state values
    v_names_age_groups = as.character(rep(ages, n_race)), # ages as character vector
    n_id_states        = length(c("S0", "I0", "S1", "I1")),  # number of states
    v_names_id_states  = c("S0", "I0", "S1", "I1"), # character vector of states
    n_age_groups       = length(ages) * n_race,         # number of age groups
    p_test             = 0, # probability of testing infected
    test_sens          = 0, # test sensitivity
    trt_eff            = trt_eff, # probability of treatment working
    gamma              = 1 / (14 / 365), # 1 / treatment length
    p_trt_s            = 1, # probability of treatment when susceptible
    p_trt_i            = 1, # probability of treatment when infected
    trt_year           = trt_year,  # year of treatment
    n_race             = n_race  # number of race groups
  )

  return(v_parameter)
}

load_sis_abr_model_params_all <- function(
    v_race = c("Hispanic", "NH White", "NH Black",
               "NH AAPI", "NH AIAN", "Other"),
    waifw, ages, trt_year = 1,
    end_t = 100, burn_state = FALSE,
    v_init_state = NULL,
    prob = FALSE,
    p_trans) {

  n_race <- length(v_race)
  v_demo <- list()
  v_prev <- list()
  v_pop  <- c()
  i <- 0
  for (race in v_race) {
    i <- i + 1
    if (race == "NH White") {
      race_text <- "white"
      pop       <- 1960.185
    } else if (race == "NH Black") {
      race_text <- "black"
      pop       <- 402.225
    } else if (race == "Hispanic") {
      race_text <- "hispanic"
      pop       <- 608.095
    } else if (race == "Other") {
      race_text <- "other"
      pop       <- 326.750
    } else if (race == "NH AAPI") {
      race_text <- "aapi"
      pop       <- 224
    } else if (race == "NH AIAN") {
      race_text <- "aian"
      pop       <- 97
    }

    tmp_demo <- get_demographic_vars(spec_groups = ages,
                                     race = race)
    tmp_prev <- get_prevalence_vars(spec_groups = ages,
                                    race = race)
    tmp_demo$v_age_prop <- tmp_demo$v_age_prop # *
      # (pop / (1960.185 + 402.225 + 608.095 + 326.750))
    v_pop <- c(v_pop, pop)
    if (i == 1) {
      v_demo <- tmp_demo
      v_prev <- tmp_prev
    } else {
      v_demo <- lapply(seq_along(tmp_demo), function(x) c(v_demo[[x]],
                                                          tmp_demo[[x]]))
      v_prev <- lapply(seq_along(tmp_prev), function(x) c(v_prev[[x]],
                                                          tmp_prev[[x]]))
    }

  }

  names(v_demo) <- c("b", "v_age", "v_age_prop", "v_mu", "v_d")
  names(v_prev) <- c("prevalence", "prevalence_sd", "foi",
                     "foi_sd", "foi_lb", "foi_ub")
  ## Get initial state
  if (burn_state == FALSE) {
    # Susceptibles
    sus <- (1 - v_prev$prevalence) * v_demo$v_age_prop
    # Infected
    inf <- v_prev$prevalence * v_demo$v_age_prop
    # Combine susceptibles and infected
    v_init <- c(sus, inf, inf * 0, sus * 0, inf * 0)
  } else {
    v_init <- v_init_state
  }

  v_parameter <- list(
    v_time             = seq(0, end_t, by = 1),  # time steps to test
    N                  = 1,                    # population size
    v_state            = v_init,               # list of initial state values
    v_prev             = v_prev$prevalence,
    v_prev_sd          = v_prev$prevalence_sd,
    v_foi              = v_prev$foi,
    v_foi_sd           = v_prev$foi_sd,
    v_foi_lb           = v_prev$foi_lb,
    v_foi_ub           = v_prev$foi_ub,
    v_age_prop         = v_demo$v_age_prop,
    v_pop              = v_pop,
    n_race             = n_race,
    # ages as character vector
    v_names_age_groups = as.character(rep(ages, n_race)),
    # number of states
    n_id_states        = length(c("S0", "I0w", "I0r", "S1", "I1r")),
    # character vector of states
    v_names_id_states  = c("S0", "I0w", "I0r", "S1", "I1r"),
    n_age_groups       = length(ages) * n_race,  # number of age groups
    b                  = v_demo$b,    # birth rate
    v_d                = v_demo$v_d,  # aging rates
    v_mu               = v_demo$v_mu, # death rates
    # transmission rates (sensitive strains)
    m_waifw            = waifw,
    gamma              = 1 / (14 / 365),       # 1 / treatment length (14 days)
    # background antibiotic use
    v_alpha            = rep(0.00023, length(ages) * n_race), # ab_vars$v_alpha
    v_psi              = rep(c(rep(0, 17), rep(1, max(ages) - 17)), n_race),
    v_psi_0            = rep(0, length(ages) * n_race), # testing this
    v_psi_s            = rep(0, length(ages) * n_race), # testing this
    # probability treatment induces sensitive strains
    sigma              = 0.375,
    phi                = 1,  # multiplier to waifw for increased transmission
    # testing rate
    v_eta              = rep(c(rep(0, 17), rep(1, max(ages) - 17)), n_race),
    test_sens          = 0.96, # test sensitivity
    test_spec          = 0.93, # test specificity
    delta              = 0.924, # probability of treatment working
    age_based_test     = FALSE,
    trt_year           = trt_year,  # year of treatment
    rho                = 0,  # probability of resitant strain becoming sensitive
    sus_test_sens      = 0.846,  # susceptibility test sensitivity
    sus_test_spec      = 0.932,   # susceptibility test specificity
    sus_test_policy    = FALSE,  # TRUE if policy uses susceptibility test
    trt_sus_policy     = FALSE,  # TRUE if policy treats everyone
    one_time_policy    = FALSE,  # TRUE if one-time policy
    to_present         = FALSE,  # TRUE if simulating to 2024
    burn               = FALSE,   # TRUE if running burn in
    p_trans            = p_trans
  )
  if (prob == TRUE) {
    v_sens_spec <- sens_spec_copula(corr = 0.1615, sens_params = c(2432, 100),
                                    spec_params = c(1471, 112), n_draws = 1)
    v_sus_sens_spec <- sens_spec_copula(corr = -0.381, sens_params = c(170, 31),
                                        spec_params = c(426, 31), n_draws = 1)
    v_parameter$sigma     <- rbeta(n = 1, shape1 = 3, shape2 = 8)
    v_parameter$delta     <- rbeta(n = 1, shape1 = 110, shape2 = 9)
    v_parameter$test_sens <- v_sens_spec[1, 1]
    v_parameter$test_spec <- v_sens_spec[1, 2]
    v_parameter$sus_test_sens <- v_sus_sens_spec[1, 1]
    v_parameter$sus_test_spec <- v_sus_sens_spec[1, 2]
  }

  return(v_parameter)
}


#' Get random draws of correlated sensitivity and specificity parameters
#'
#' @param corr correlation coefficient
#' @param sens_params alpha and beta for test sensitivity
#' @param spec_params alpha and beta for test specificity
#' @param n_draws number of draws from copula
#'
#' @return sens_spec
#' @export
#'
#' @examples v_sens_spec <- sens_spec_copula(corr = 0.1615,
#'                                           sens_params = c(2432, 100),
#'                                           spec_params = c(1471, 112),
#'                                           n_draws = 1)
sens_spec_copula <- function(corr, sens_params, spec_params, n_draws) {
  cop <- copula::normalCopula(param = corr, dim = 2, dispstr = "un")
  mvd <- copula::mvdc(copula = cop, margins = c("beta", "beta"),
                      paramMargins = list(list(shape1 = sens_params[1],
                                               shape2 = sens_params[2]),
                                          list(shape1 = spec_params[1],
                                               shape2 = spec_params[2])))
  sens_spec <- copula::rMvdc(n = n_draws, mvdc = mvd)
  return(sens_spec)
}

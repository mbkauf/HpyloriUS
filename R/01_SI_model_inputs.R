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
load_SI_model_params <- function(waifw, demography_vars, prevalence_vars,
                                 ages) {

  ## Get initial state
  sus <- (1-prevalence_vars$prevalence) * demography_vars$v_age_prop # susceptibles
  inf <- prevalence_vars$prevalence * demography_vars$v_age_prop  # infected
  a.init <- cbind(sus, inf)  # combine susceptibles and infected
  dimnames(a.init) <- list(groups, c("S", "I"))  # add names
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
    n_id_states        = length(c("S", "I")),  # number of states
    v_names_id_states  = c("S", "I"),          # character vector of states
    n_age_groups       = length(ages)          # number of age groups

  )

  return(v_parameter)
}



#' Generate WAIFW matrix
#'
#' @param betas       vector of beta values for waifw matrix
#' @param breaks      breaks (should be length of age groups + 1) ?
#' @param W           contact matrix shell for betas (frequently upper triangular index matrix)
#' @param upper       TRUE/FALSE about W being an upper triangular index matrix
#' @param group.names vector of numbers for age range (used for naming matrix)
#'
#' @return waifw
#' @export
#'
#' @examples Beta <- get_transmission_matrix(betas = beta0,
#'                                           breaks = waifw.breaks,
#'                                           W = test.w,
#'                                           group.names = groups)
get_transmission_matrix <- function(betas, breaks, W, upper = TRUE, group.names = groups) {

  ### Number of age groups = number of data to estimate
  n.params <- n.a <- (length(breaks)-1)

  ### Age.max
  age.max <- breaks[length(breaks)]-1

  ### Age group labels
  age.names <- paste(breaks[-9], breaks[-1]-1, sep = "-")

  if (n.params!=8){stop("Number of beta parameters is different than 8")}
  if (max(W) > n.params ) {betas <- c(betas, 0)}
  ## Complete index matrix from upper triangular index matrix
  complete.index <- function(index.upper){
    index <- index.upper + t(index.upper)
    diag(index) <- diag(index.upper)
    return(index)
  }

  if (upper == TRUE) {
    index <- complete.index(index.upper = W)
    ### Generate age-group waifw
    bij <- matrix(betas[index], ncol = n.params, dimnames = list(age.names, age.names))
    ### Generate 1-year age-specific waifw
    index.full <- index[rep(1:n.a, diff(breaks)),
                        rep(1:n.a, diff(breaks))]
    colnames(index.full) <- rep(age.names, diff(breaks))
    m_waifw <- matrix(betas[index.full], ncol = (age.max+1))

  } else {
    index <- W
    ### Generate age-group waifw
    bij <- matrix(betas[index], ncol = n.params, dimnames = list(age.names, age.names))
    ### Generate 1-year age-specific waifw
    index.full <- index[rep(1:n.a, diff(breaks)),
                        rep(1:n.a, diff(breaks))]
    colnames(index.full) <- rep(age.names, diff(breaks))
    m_waifw <- matrix(betas[index.full], ncol = (age.max+1))
  }

  colnames(m_waifw) <- rownames(m_waifw) <- group.names
  return(m_waifw)
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
get_demographic_vars <- function(life_table_data = "LifeTable_US_All.csv",
                                 spec_groups = NULL,
                                 r_pop_growth = 0,
                                 by_race_sex = TRUE,
                                 Race = c("Non-Hispanic White",
                                          "Non-Hispanic Black",
                                          "Other Hispanic",
                                          "Mexican-American",
                                          "Other"),
                                 Sex = c("Male", "Female")) {
  if (!is.null(spec_groups)) {
    ngroups <- length(spec_groups)
    if (!is.numeric(spec_groups)) {
      spec_groups <- as.numeric(spec_groups)
      warning("spec_groups was not entered as a numeric vector, so it was converted to a numeric vector.
              Results may not be correct")
    }
    older_upper_lim <- max(spec_groups) + 1
    lower_age_lim   <- min(spec_groups)
  }

  #========================
  # calculate death rates
  #========================
  # import life-table data here
  if(by_race_sex==TRUE) {
    if(sex=="Male"){
      if(race=="Non-Hispanic White") {
        life_table_data <- "LifeTable_US_White_Male.csv"
      } else if(race=="Non-Hispanic Black") {
        life_table_data <- "LifeTable_US_Black_Male.csv"
      } else if(race=="Other Hispanic") {
        life_table_data <- "LifeTable_US_Hispanic_Male.csv"
      } else if(race=="Mexican-American") {
        life_table_data <- "LifeTable_US_Hispanic_Male.csv"
      } else if(race=="Other") {
        life_table_data <- "LifeTable_US_All.csv"
      }
    } else {
      if(race=="Non-Hispanic White") {
        life_table_data <- "LifeTable_US_White_Female.csv"
      } else if(race=="Non-Hispanic Black") {
        life_table_data <- "LifeTable_US_Black_Female.csv"
      } else if(race=="Other Hispanic") {
        life_table_data <- "LifeTable_US_Hispanic_Female.csv"
      } else if(race=="Mexican-American") {
        life_table_data <- "LifeTable_US_Hispanic_Female.csv"
      } else if(race=="Other") {
        life_table_data <- "LifeTable_US_All.csv"
      }
    }
  } else {
    life_table_data <- "LifeTable_US_All.csv"
  }

  temp <- read.csv(paste0("data/", life_table_data))

  # Rates of death for ages of interest
  v_mu <- as.numeric(matrix(temp[lower_age_lim:older_upper_lim, 2]))
  names(v_mu) <- as.character(temp[lower_age_lim:older_upper_lim, 1])

  #=========
  # calculate proportion of population in each age class
  #=========
  #age are the all-purpose names, which contain the beginning year of each age group
  v_age <- as.character(spec_groups)

  #initialize age_proportions vector
  v_age_prop <- rep(0, ngroups)

  #combine younger and older groups for proportion calculations and for outputting later
  # v_d  <- (v_mu + r_pop_growth)/(exp((v_mu + r_pop_growth)) - 1)
  v_d <- c(rep(1, (ngroups-1)), 0)
  names(v_d) <- names(v_mu) <- names(v_age_prop) <- v_age

  #using Merck pg. 10
  #explicit translation of first formula on pg. 10
  #calculating youngest age group
  acc <- 0
  inner_acc <- 1
  # inner product/sum
  for (rr in 2:ngroups) {
    for (jj in 2:rr) {
      inner_acc <- ( v_d[v_age[jj - 1]] / (v_d[v_age[jj]] + v_mu[v_age[jj]]) ) * inner_acc
    }
    acc <- acc + inner_acc
    inner_acc <- 1
  }

  #calculate proportion in first age group
  v_age_prop[1] <- 1 / (1 + acc)

  for (AG in v_age[-1]) {
    prev_age <- v_age[which(v_age == AG) - 1]
    v_age_prop[AG] <- v_d[prev_age] * v_age_prop[prev_age]  / (v_d[AG] + v_mu[AG] + r_pop_growth)
  }

  #====================
  #calculate entry rate
  #====================
  b <- (v_d[1] + v_mu[1] + r_pop_growth)*v_age_prop[1]

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
get_prevalence_vars <- function(data_foi = "prev_foi.csv",
                                data_prev = "prev_foi.csv",
                                spec_groups = NULL,
                                Race = c("Non-Hispanic White",
                                         "Non-Hispanic Black",
                                         "Other Hispanic",
                                         "Mexican-American",
                                         "Other"),
                                Sex = c("Male", "Female")) {
  require(dtplyr)
  require(dplyr)

  if (!is.null(spec_groups)) {
    n_groups <- length(spec_groups)
    if (!is.numeric(spec_groups)) {
      spec_groups <- as.numeric(spec_groups)
      warning("spec_groups was not entered as a numeric vector, so it was converted to a numeric vector.
              Results may not be correct")
    }
    keep_ages <- c(spec_groups, (max(spec_groups)+1))
  }


  #========================
  # Load NHANES Hp dataset
  #========================
  df_foi <- read.csv(paste0("data/", data_foi))

  df_foi <- df_foi %>%
    filter(age %in% spec_groups) %>%
    filter(NHANES_year == "1991") %>%
    filter(race_f == Race) %>%
    filter(sex == Sex)

  foi    <- as.numeric(df_foi$Mean)
  foi_sd <- as.numeric(df_foi$SD)

  df_prev <- read.csv(paste0("data/", data_prev))

  df_prev <- df_prev %>%
    filter(age %in% spec_groups) %>%
    filter(NHANES_year == "1991") %>%
    filter(race_f == "Non-Hispanic White") %>%
    filter(sex == Sex) %>%
    mutate(se_fit = (prev_fit - LB_scl)/1.96)

  prevalence    <- as.numeric(df_prev$prev_fit)
  prevalence_se <- as.numeric(df_prev$se_fit)

  return(list(prevalence = prevalence, prevalence_se = prevalence_se,
              foi = foi, foi_sd = foi_sd))
}

# Age groups: [0, 4), [5, 14), [15, 24), [25, 44), [45, 54), [55, 64), [65, 69), [70,100)

# index.upper.w2 <- matrix(c(1, 1, 3, 4, 5, 6,
#                            0, 2, 3, 4, 5, 6,
#                            0, 0, 3, 4, 5, 6,
#                            0, 0, 0, 4, 5, 6,
#                            0, 0, 0, 0, 5, 6,
#                            0, 0, 0, 0, 0, 6), ncol = n.params, byrow = T)
#
# index.upper.w3 <- matrix(c(1, 1, 1, 4, 5, 6,
#                            0, 2, 3, 4, 5, 6,
#                            0, 0, 3, 4, 5, 6,
#                            0, 0, 0, 4, 5, 6,
#                            0, 0, 0, 0, 5, 6,
#                            0, 0, 0, 0, 0, 6), ncol = n.params, byrow = T)
#
# index.w4 <- matrix(c(1, 1, 1, 1, 1, 1,
#                      2, 2, 2, 2, 2, 2,
#                      3, 3, 3, 3, 3, 3,
#                      4, 4, 4, 4, 4, 4,
#                      5, 5, 5, 5, 5, 5,
#                      6, 6, 6, 6, 6, 6), ncol = n.params, byrow = T)
#
# betas.w6 <- c(betas, 0) # needed for w6
# index.upper.w6 <- matrix(c(1, 7, 7, 7, 7, 7,
#                            0, 2, 7, 7, 7, 7,
#                            0, 0, 3, 7, 7, 7,
#                            0, 0, 0, 4, 7, 7,
#                            0, 0, 0, 0, 5, 7,
#                            0, 0, 0, 0, 0, 6), ncol = n.params, byrow = T)



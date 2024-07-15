library(data.table)
library(dplyr)
library(dtplyr)
library(survey)
library(ipumsr)

# Load Data
ddi <- read_ipums_ddi("data/usa_00002.xml")
data <- read_ipums_micro(ddi)


ipums <- fread("data/usa_00002.csv")


ipums_rev <- lazy_dt(ipums) %>%
  mutate(
    NH_WHITE = ifelse(RACE == 1 & HISPAN == 0, 1, 0),
    NH_BLACK = ifelse(RACE == 2 & HISPAN == 0, 1, 0),
    HISPANIC = ifelse(HISPAN != 0, 1, 0),
    NH_OTHER = ifelse(RACE > 2 & HISPAN == 0, 1, 0)
    ) %>%
  group_by(SERIAL) %>%
  summarise(any_nhwhite = max(NH_WHITE),
            any_nhblack = max(NH_BLACK),
            any_hispanic = max(HISPANIC),
            any_nhother = max(NH_OTHER),
            hwt = max(HHWT),
            pop = max(row_number())) %>%
  mutate(house_race = ifelse(any_hispanic == 1, 1,
                             ifelse(any_nhblack == 1, 2,
                                    ifelse(any_nhwhite == 1, 3,
                                           ifelse(any_nhother == 1, 4, NaN))))) %>%
  # filter(pop > 1) %>%
  as_tibble()

w_ipums <- svydesign(id=~SERIAL, weights=~hwt, data=ipums_rev)
tb_hisp_black  <- svytable(~any_hispanic+any_nhblack, w_ipums)
tb_hisp_white  <- svytable(~any_hispanic+any_nhwhite, w_ipums)
tb_hisp_other  <- svytable(~any_hispanic+any_nhother, w_ipums)
tb_black_white <- svytable(~any_nhblack+any_nhwhite, w_ipums)
tb_black_other <- svytable(~any_nhblack+any_nhother, w_ipums)
tb_white_other <- svytable(~any_nhwhite+any_nhother, w_ipums)

get_prop <- function(x) {x[2,2] / sum(x[1,2], x[2,])}
house_hisp_black <- get_prop(tb_hisp_black)
house_hisp_white <- get_prop(tb_hisp_white)
house_hisp_other <- get_prop(tb_hisp_other)
house_black_white <- get_prop(tb_black_white)
house_black_other <- get_prop(tb_black_other)
house_white_other <- get_prop(tb_white_other)

# Look at PUMA level
# ipums_pumas <- lazy_dt(ipums) %>%
#   mutate(
#     NH_WHITE = ifelse(RACE == 1 & HISPAN == 0, 1, 0),
#     NH_BLACK = ifelse(RACE == 2 & HISPAN == 0, 1, 0),
#     HISPANIC = ifelse(HISPAN != 0, 1, 0),
#     NH_OTHER = ifelse(RACE > 2 & HISPAN == 0, 1, 0)
#   ) %>%
#   group_by(PUMA) %>%
#     mutate(pop_nhwhite = mean(NH_WHITE),
#            pop_nhblack = mean(NH_BLACK),
#            pop_hispanic = mean(HISPANIC),
#            pop_nhother = mean(NH_OTHER),
#            pwt = mean(PERWT),
#            pop = max(row_number()),
#            fraction_group = estimate / pop,
#            max_fraction_group = max(fraction_group)) %>%
#   filter(fraction_group == max_fraction_group) %>%
#   ungroup() %>%
#   as_tibble()
  # summarise(pop_nhwhite = mean(NH_WHITE),
  #           pop_nhblack = mean(NH_BLACK),
  #           pop_hispanic = mean(HISPANIC),
  #           pop_nhother = mean(NH_OTHER),
  #           pwt = mean(PERWT),
  #           pop = max(row_number())) %>%


# ipums_pumas <- lazy_dt(ipums_pumas) %>%
#   mutate(hisp_nhwhite = pop_hispanic + pop_nhwhite) %>%
#   mutate(hisp_nhblack = pop_hispanic + pop_nhblack) %>%
#   mutate(hisp_nhother = pop_hispanic + pop_nhother) %>%
#   mutate(nhwhite_nhblack = pop_nhwhite + pop_nhblack) %>%
#   mutate(nhwhite_nhother = pop_nhwhite + pop_nhother) %>%
#   mutate(nhblack_nhother = pop_nhblack + pop_nhother) %>%
#   as_tibble()

ipums_pumas <- lazy_dt(ipums) %>%
  mutate(
    NH_WHITE = ifelse(RACE == 1 & HISPAN == 0, 1, 0),
    NH_BLACK = ifelse(RACE == 2 & HISPAN == 0, 1, 0),
    HISPANIC = ifelse(HISPAN != 0, 1, 0),
    NH_OTHER = ifelse(RACE > 2 & HISPAN == 0, 1, 0),
    NHANES_RACE = ifelse(HISPANIC == 1, 1,
                         ifelse(NH_WHITE == 1, 2,
                                ifelse(NH_BLACK == 1, 3, 4)))
  ) %>%
  group_by(PUMA, NHANES_RACE) %>%
  summarize(race_pop = sum(PERWT)) %>%
  mutate(
    fraction_group = race_pop / sum(race_pop),
    max_fraction_group = max(fraction_group)
  ) %>%
  as_tibble()

df_pumas_max <-  ipums_pumas %>%
  filter(fraction_group == max_fraction_group)


### Within race parameters
# Look at households with more than 1 person and all the same race
ipums_house <- lazy_dt(ipums) %>%
  mutate(
    NH_WHITE = ifelse(RACE == 1 & HISPAN == 0, 1, 0),
    NH_BLACK = ifelse(RACE == 2 & HISPAN == 0, 1, 0),
    HISPANIC = ifelse(HISPAN != 0, 1, 0),
    NH_OTHER = ifelse(RACE > 2 & HISPAN == 0, 1, 0)
  ) %>%
  group_by(SERIAL) %>%
  summarise(only_nhwhite  = +(n_distinct(NH_WHITE) == 1),
            only_nhblack  = +(n_distinct(NH_BLACK) == 1),
            only_hispanic = +(n_distinct(HISPANIC) == 1),
            only_nhother  = +(n_distinct(NH_OTHER) == 1),
            hwt = max(HHWT),
            pop = max(row_number())) %>%
  mutate(alone = ifelse(pop == 1, 1, 0)) %>%
  filter(pop > 1) %>%
  as_tibble()

w_ipums_house <- svydesign(id=~SERIAL, weights=~hwt, data=ipums_house)

tb_hisp_only  <- svytable(~only_hispanic + alone, w_ipums_house)
tb_white_only  <- svytable(~only_nhwhite + alone, w_ipums_house)
tb_black_only  <- svytable(~only_nhblack + alone, w_ipums_house)
tb_other_only <- svytable(~only_nhother + alone, w_ipums_house)

# % of households > 1 person and contain only that race/ethnicity
house_hisp  <- tb_hisp_only[2,1]/sum(tb_hisp_only)
house_white <- tb_white_only[2,1]/sum(tb_white_only)
house_black <- tb_black_only[2,1]/sum(tb_black_only)
house_other <- tb_other_only[2,1]/sum(tb_other_only)
1- house_hisp
1- house_white
1- house_black
1- house_other
# % of PUMAs where race makes less than 75% of population
hisp_75  <- sum(df_pumas_max$NHANES_RACE == 1 &
                  df_pumas_max$max_fraction_group >= 0.75) / 982
white_75 <- sum(df_pumas_max$NHANES_RACE == 2 &
                  df_pumas_max$max_fraction_group >= 0.75) / 982
black_75 <- sum(df_pumas_max$NHANES_RACE == 3 &
                  df_pumas_max$max_fraction_group >= 0.75) / 982
other_75 <- sum(df_pumas_max$NHANES_RACE == 4 &
                  df_pumas_max$max_fraction_group >= 0.75) / 982

### alpha_ij
## Function to get bounds
get_ij_bounds <- function(race_i, race_j, x, y, z) {
  df_ij <- ipums_pumas %>%
    filter(NHANES_RACE == race_i | NHANES_RACE == race_j) %>%
    group_by(PUMA) %>%
    mutate(min_fraction_group = min(fraction_group)) %>%
    mutate(max_fraction_group = max(fraction_group)) %>%
    mutate(other_fraction_group = ifelse(min_fraction_group == fraction_group,
                                         max_fraction_group,
                                         min_fraction_group)) %>%
    ungroup() %>%
    filter(NHANES_RACE == race_i)
  lb <- sum(df_ij$fraction_group >= x &
              df_ij$other_fraction_group >= z)/982
  ub <- sum(df_ij$fraction_group >= y &
              df_ij$other_fraction_group >= z)/982
  return(c(lb, ub))
}

# Get bounds
pct_x <- 0.6
pct_y <- 0.2
pct_z <- 0.1

alpha_12_bounds <- round(get_ij_bounds(race_i = 1, race_j = 2,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_13_bounds <- round(get_ij_bounds(race_i = 1, race_j = 3,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_14_bounds <- round(get_ij_bounds(race_i = 1, race_j = 4,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_21_bounds <- round(get_ij_bounds(race_i = 2, race_j = 1,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_23_bounds <- round(get_ij_bounds(race_i = 2, race_j = 3,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_24_bounds <- round(get_ij_bounds(race_i = 2, race_j = 4,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_31_bounds <- round(get_ij_bounds(race_i = 3, race_j = 1,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_32_bounds <- round(get_ij_bounds(race_i = 3, race_j = 2,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_34_bounds <- round(get_ij_bounds(race_i = 3, race_j = 4,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_41_bounds <- round(get_ij_bounds(race_i = 4, race_j = 1,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_42_bounds <- round(get_ij_bounds(race_i = 4, race_j = 2,
                                 x = pct_x, y = pct_y, z = pct_z), 3)
alpha_43_bounds <- round(get_ij_bounds(race_i = 4, race_j = 3,
                                 x = pct_x, y = pct_y, z = pct_z), 3)

df_ij <- ipums_pumas %>%
  filter(NHANES_RACE == 4 | NHANES_RACE == 1) %>%
  group_by(PUMA) %>%
  mutate(min_fraction_group = min(fraction_group)) %>%
  mutate(max_fraction_group = max(fraction_group)) %>%
  mutate(other_fraction_group = ifelse(min_fraction_group == fraction_group,
                                       max_fraction_group,
                                       min_fraction_group)) %>%
  ungroup() %>%
  filter(NHANES_RACE == 4)

## alpha_aa values
alpha_upper_hisp_hisp   <- house_hisp + hisp_80
alpha_upper_white_white <- house_white + white_80
alpha_upper_black_black <- house_black + black_80
alpha_upper_other_other <- house_other + other_80

alpha_lower_hisp_hisp   <- sum(house_hisp_black,
                               house_hisp_white,
                               house_hisp_other)
alpha_lower_white_white <- sum(house_black_white,
                               house_hisp_white,
                               house_white_other)
alpha_lower_black_black <- sum(house_hisp_black,
                               house_black_white,
                               house_black_other)
alpha_lower_other_other <- sum(house_hisp_other,
                               house_black_other,
                               house_white_other)

## alpha_ab values
alpha_upper_hisp_white <- sum(df_puma_hisp_white$min_fraction_group >= 0.2)/982
alpha_upper_hisp_black <- sum(df_puma_hisp_black$min_fraction_group >= 0.2)/982
alpha_upper_hisp_other <- sum(df_puma_hisp_other$min_fraction_group >= 0.2)/982


## Check population sizes
ipums_pop <- lazy_dt(ipums) %>%
  mutate(
    NH_WHITE = ifelse(RACE == 1 & HISPAN == 0, 1, 0),
    NH_BLACK = ifelse(RACE == 2 & HISPAN == 0, 1, 0),
    HISPANIC = ifelse(HISPAN != 0, 1, 0),
    NH_OTHER = ifelse(RACE > 2 & HISPAN == 0, 1, 0),
    NHANES_RACE = ifelse(HISPANIC == 1, 1,
                         ifelse(NH_WHITE == 1, 2,
                                ifelse(NH_BLACK == 1, 3, 4)))
  ) %>%
  as_tibble()

sum(ipums_pop$HISPANIC * ipums_pop$PERWT)/100000
sum(ipums_pop$NH_WHITE * ipums_pop$PERWT)/100000
sum(ipums_pop$NH_BLACK * ipums_pop$PERWT)/100000
sum(ipums_pop$NH_OTHER * ipums_pop$PERWT)/100000


### Create scatter of race/ethnicity combinations ###
df_pumas_long <- ipums_pumas %>%
  select(c("PUMA", "NHANES_RACE", "fraction_group")) %>%
  pivot_longer(!PUMA,
               names_to = "race",
               values_to = "pop_frac")

df_pumas_wide <- ipums_pumas %>%
  select(c("PUMA", "NHANES_RACE", "fraction_group")) %>%
  pivot_wider(names_from = NHANES_RACE,
              values_from = fraction_group)
colnames(df_pumas_wide) <- c("PUMA", "HISPANIC", "NH_WHITE", "NH_BLACK", "NH_OTHER")


plotfun_ub = function(x,y){
  points(x,y,pch=18);abline(h=0.1, v = 0.2, lty=8,col="red")
  rect(xleft = 0.2, xright = 1, ybottom = 0.1, ytop = 1, col = rgb(1,0,0,alpha=0.5))
}
plotfun_lb = function(x,y){
  points(x,y,pch=18);abline(h=0.1, v = 0.6, lty=8,col="red")
  rect(xleft = 0.6, xright = 1, ybottom = 0.1, ytop = 1, col = rgb(1,0,0,alpha=0.5))
}
png("results/PUMAS_ub.png", height = 1000, width = 1000)
pairs(df_pumas_wide[2:5], panel = plotfun_ub)
dev.off()
png("results/PUMAS_lb.png", height = 1000, width = 1000)
pairs(df_pumas_wide[2:5], panel = plotfun_lb)
dev.off()


library(GGally)
ggpairs(df_pumas_wide[2:5])

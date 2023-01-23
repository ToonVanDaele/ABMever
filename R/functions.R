# Functions for (ABM) population models


#-----------------------------------------------
# Set initial ages

set_init_pop <- function(init_agecl, birth_month, Sm, max_age = 180){

  # This function sets distribution of initial ages for the ABM boar model close
  # to the expected 'stable age distribution'. The distribution over the
  # ages could also be derived from an extended matrix population model.
  # For a more accurate initial age distribution the ABM should be run for a
  # few (4/5) years.
  #
  # - Sex ratio is fixed to 50/50 for the initial population.
  # - The initial age for classes 0 (juvenile) and 1 (yearling) are set
  # following the distribution of births per month. This is then corrected
  # for the survival during the months between birth and the start of the model.
  # This function only works for models that start at 1st of January!!
  #
  # @param init_agecl Initial distribution in each age class (vector of length 3)
  # @param birth_month Relative distribution of births during the year (vector of length 12)
  # @param Sm Monthly survival for each age class (vector of length 3)
  # @param max_age Maximum age (default = 180 months)
  #
  # @return Dataframe with 2 columns (age, sex) and one row per individual
  #
  df <- birth_month %>%

    # The age distribution of the juveniles and yearlings is according
    # the birth distribution and the respective survival till the end of the month
    rownames_to_column(var = "m") %>%   # add column with number of the month
    mutate(m = as.numeric(m)) %>%
    mutate(perc_0s = perc * Sm[1] ^ (13 - m), # birth+survival juv. till end of year
           perc_1s = perc * Sm[2] ^ (13 - m)) %>%   # idem for yearling
    map_at(vars(perc_0s, perc_1s), ~./sum(.)) %>%   # calculate proportions
    as_tibble() %>%
    mutate(cl_0 = init_agecl[1] * perc_0s,               # juveniles
           cl_1 = init_agecl[2] * perc_1s) %>%           # yearling
    mutate(cl_0_F = round(cl_0 / 2, 0),      # juveniles female
           cl_1_F = round(cl_1 / 2, 0)) %>%   # yearling female
    mutate(cl_0_M = round(cl_0, 0) - cl_0_F, # juveniles male
           cl_1_M = round(cl_1, 0) - cl_1_F) # yearling male

  # We put everything in separate data frames
  agecl_0_F <- data.frame(age = rep(1:12, df$cl_0_F), sex = "F")
  agecl_1_F <- data.frame(age = rep(1:12, df$cl_1_F) + 12, sex = "F")
  agecl_0_M <- data.frame(age = rep(1:12, df$cl_0_M), sex = "M")
  agecl_1_M <- data.frame(age = rep(1:12, df$cl_1_M) + 12, sex = "M")

  # The initial age for the adults is set with a geometric distribution to
  # simulate the expected adult age distribution (shape and rate to be tuned!!)
  # The age distribution of adults has no impact on the population projection
  agecl_2_nb <- round(init_agecl[3], 0)

  # Age from a geometric distribution - maximum is max_age
  agecl_2_age <- rgamma(n = agecl_2_nb, shape = 2, rate = 0.8) * 12 + 24
  # While values > max age exist -> choose new values from the distribution
  while (length(agecl_2_age[agecl_2_age > max_age]) > 0) {
    new_values <- rgamma(n = length(agecl_2_age[agecl_2_age > max_age]),
                         shape = 2, rate = 0.8) * 12 + 24
    agecl_2_age[agecl_2_age > max_age] <- new_values
  }

  agecl_2_sex <- c(rep("F", round(agecl_2_nb / 2, 0)),
                   rep("M", agecl_2_nb - round(agecl_2_nb / 2, 0)))
  agecl_2 <- data.frame(age = agecl_2_age,
                        sex = agecl_2_sex)

  # rbind all data frames together
  init_pop <- bind_rows(agecl_0_F, agecl_0_M, agecl_1_F, agecl_1_M, agecl_2)

  return(init_pop)
}

#-------------------------------------------------------------------------
# Set initial ages - uniform distribution for class 0 and 1
set_init_pop2 <- function(init_agecl){

  # We set a fixed (50/50) sex ratio for the initial population.
  # This makes comparision with the (female only) matrix population model easier

  # Initial sex can also be set randomly. Although this can initially create
  # extra variability in the population. Stable sex ratio will only be reached
  # after several years.

  nb_f <- round(init_agecl / 2, 0)
  nb_m <- round(init_agecl - nb_f, 0)
  # first two age classes - uniform distributed
  age_0 <- runif(n = round(init_agecl[1], 0),min = 1, max = 12)
  sex_0 <- c(rep("F", nb_f[1]), rep("M", nb_m[1]))

  age_1 <- runif(n = round(init_agecl[2], 0),min = 13, max = 24)
  sex_1 <- c(rep("F", nb_f[2]), rep("M", nb_m[2]))

  # adult geometric distribution
  age_2 <- rgamma(n = round(init_agecl[3], 0), shape = 2, rate = 0.8) * 12 + 24
  sex_2 <- c(rep("F", nb_f[3]), rep("M", nb_m[3]))

  ages <- c(age_0, age_1, age_2)
  sexes <- c(sex_0, sex_1, sex_2)

  init_pop <- data.frame(age = ages, sex = sexes)

  return(init_pop)
}

#----------------------------------------------------------------
# Get summary of boar (population by age class, sex)
#
# @param turtles
#
# @return data frame with number of turtles (n) by sex & age class
#
get_boar <- function(turtles){

  df <- turtles@.Data %>%
    as.data.frame() %>%
    filter(agecl >= 0) %>%
    group_by(sex, agecl) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(sex = turtles@levels$sex[sex],
           agecl = ageclasses[agecl + 1])
  return(df)
}

#----------------------------------------------------------------
# Get summary of harvested boar (population by age class, sex)
#
# @param turtles
#
# @return data frame with number of turtles (n) by sex & age class
#
get_boar_harvest <- function(turtles){

  df <- turtles@.Data %>%
    as.data.frame() %>%
    filter(agecl >= 0) %>%      # change to age >=0
    mutate(dj = ifelse(newb < 5, offspring, 0)) %>%
    group_by(sex, agecl) %>%
    summarise(n = n(),
              dep_juv = sum(dj), .groups = "drop") %>%
    mutate(sex = turtles@levels$sex[sex],
           agecl = ageclasses[agecl + 1])
  return(df)
}

#-------------------------------------------------------------
# Set Fertility
#
# Yearly fertility is distributed over months according
# monthly birth rates.
#
# @param F yearly fertility by age class (vector length 3)
# @param birth_month ratio of births by month (vector length 12)
#
# @return matrix monthly fertility by class (columns) and month(rows)
#
set_F <- function(F, birth_month){

  Fm <- matrix(data = c(rep(birth_month$perc / 100, 3)),
               nrow = 12, ncol = 3,
               dimnames = list(NULL, ageclasses))
  Fm <- t(t(Fm) * F)

  return(Fm)
}

#----------------------------------------------------------------------
# Get geboortepiek
#
# Load geboortepiek as csv
#
# @param csv_filename (string)
#
# @return data frame with percentage of births by month (12 rows)
#
get_birth_month <- function(csv_filename){

#  filename <- "https://raw.githubusercontent.com/inbo/fis-projecten/305-rerun-geboortepiek/Grofwild/Populatiemodel-Everzwijn/Output/Geboortepiek/bp_results_1month.csv?token=GHSAT0AAAAAABYGSVSY725WYIFIZ4K4DKM6YZ23W6Q"

# df <- read.csv(file = filename)
# write.csv(df, file = csv_filename")

  df <- read.csv(file = csv_filename)

  if (round(sum(df$per),2) != 100)
    warning("Sum of monthly birth percentages is not 100%")

  df <- df %>%
    dplyr::select(perc = per)
  return(df)
}

#-----------------------------------------------------------------------

# Get hunting scenario's
#
# @param path filename of excelsheet
#
# @return list with matrices for each hunting scenario
#
get_hunting_scen <- function(path){

  f <- function(path, sheet){
    df <- readxl::read_excel(path = path, sheet = sheet, range = "A1:G13")
    return(as.matrix(df))
  }

  Hscen <- path %>%
    readxl::excel_sheets() %>%
    set_names() %>%
    map(f, path = path)
  return(Hscen)
}


#-----------------------------------------------------------------------
# Get survival table
#
# @param path filename of excelsheet
#
# @return list with matrices for survival
#
# for each age class, month, female and male

get_survival <- function(path){

  f <- function(path, sheet){
    df <- readxl::read_excel(path = path, sheet = sheet, range = "A1:F13")
    return(as.matrix(df))
  }

  surv <- path %>%
    readxl::excel_sheets() %>%
    set_names() %>%
    map(f, path = path)
  return(surv)
}

#-------------------------------------------------------------
# process output - get df_numboar
#
# @param mytb simulation output from scen_sim_boar function (tibble)
#
# @return data frame with number of boars by age class, sex, hunting scenario and simulation
#
get_numboar <- function(mytb){

  df_num <- mytb$result %>%
    map_dfr("df_numboar", .id = "rowname") %>%
    left_join(mytb %>%
              dplyr::select(Hs, sim) %>%
              rownames_to_column(),
            by = "rowname") %>%
  mutate(sex = as.factor(sex),
         agecl = as.factor(agecl),
         Hs = as.factor(Hs))
  return(df_num)
}

#-------------------------------------------------------------
# process output - get df_harvest
#
# @param mytb simulation output from scen_sim_boar function (tibble)
#
# @return data frame with number of harvest boars by age class, sex, hunting scenario and simulation
#

get_harvest <- function(mytb){

  df_num <- mytb$result %>%
    map_dfr("df_harvest", .id = "rowname") %>%
    left_join(mytb %>%
                dplyr::select(Hs, sim) %>%
                rownames_to_column(),
              by = "rowname") %>%
    mutate(sex = as.factor(sex),
           agecl = as.factor(agecl),
           Hs = as.factor(Hs))
  return(df_num)
}

#-------------------------------------------------------------
# process output - get age distribution
#
# @param mytb simulation output from scen_sim_boar function (tibble)
#
# @return data frame with end population of each simulation
#
get_pop <- function(mytb){

  df_pop <- mytb$result %>%
    map_dfr("df_pop", .id = "rowname") %>%
    left_join(mytb %>%
                dplyr::select(Hs, sim) %>%
                rownames_to_column(),
              by = "rowname") %>%
    as_tibble()
  return(df_pop)
}

#---------------------------------------------------------
# Calculate lambda based on the output of the 'get_numboar' function
#
# @param df_num output from the get_numboar function (data frame)
#
# @return lambda (yearly) for each simulation, hunting scenario and time step

calc_lambda <- function(df_num, burnin = 0){

  df <- df_num %>%
    filter(time > burnin) %>%
    group_by(sim, Hs, time) %>%
    summarise(ntot = sum(n), .groups = "drop_last") %>%
    mutate(lambda = ntot / lag(ntot, 1)) %>%
    # Calculate geometric mean
    summarise(gm_lambda_m = exp(mean(log(lambda), na.rm = TRUE)), .groups = "drop") %>%
    mutate(gm_lambda_y = (gm_lambda_m ^ 12 ))
  return(df)
}


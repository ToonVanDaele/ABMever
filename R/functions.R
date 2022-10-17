# Function


#-----------------------------------------------
# Set initial ages

set_init_pop <- function(init_agecl, birth_month, Sm){

  # The initial population should be as close to the stable stage population
  # of the (female only) matrix population model. This enables optimal
  # comparision between BM and matrix models.
  #
  # Sex ratio is fixed to 50/50 for the initial population.
  # Initial sex ratio could also be set randomly. This would introduce extra
  # variability in the population.
  #
  # Initial age for classes 0 (juvenile) and 1 (yearling) are set according
  # the distribution of births per month. This is corrected for the survival
  # during the months between birth and the start of the model.

  df <- birth_month %>%
    rownames_to_column(var = "a") %>%
    mutate(a = as.numeric(a)) %>%
    mutate(perc_0s = perc * Sm[1] ^ (13 - a),
           perc_1s = perc * Sm[2] ^ (13 - a)) %>%
    mutate_at(vars(perc_0s, perc_1s), funs(./ sum(.))) %>%
    mutate(cl_0 = init_agecl[1] * perc_0s,
           cl_1 = init_agecl[2] * perc_1s) %>%
    mutate(cl_0_F = round(cl_0 / 2, 0),
           cl_1_F = round(cl_1 / 2, 0)) %>%
    mutate(cl_0_M = round(cl_0, 0) - cl_0_F,
           cl_1_M = round(cl_1, 0) - cl_1_F)

  agecl_0_F <- data.frame(age = rep(1:12, df[,"cl_0_F"]), sex = "F")
  agecl_1_F <- data.frame(age = rep(1:12, df[,"cl_1_F"]) + 12, sex = "F")
  agecl_0_M <- data.frame(age = rep(1:12, df[,"cl_0_M"]), sex = "M")
  agecl_1_M <- data.frame(age = rep(1:12, df[,"cl_1_M"]) + 12, sex = "M")

  # The inital age for the adults is set with a geometric distribution
  # which simulates the expected age distribution for adults (to be tuned!!)
  # The age distribution of adults has no impact on the population projection
  sex_2_F <- round(init_agecl[3] / 2, 0)
  sex_2_M <- round(init_agecl[3], 0) - sex_2_F
  agecl_2_F <- data.frame(age = rgamma(n = sex_2_F, shape = 2, rate = 0.8) * 12 + 24,
                          sex = "F")
  agecl_2_M <- data.frame(age = rgamma(n = sex_2_M, shape = 2, rate = 0.8) * 12 + 24,
                          sex = "M")

  init_pop <- bind_rows(agecl_0_F, agecl_0_M, agecl_1_F, agecl_1_M,
                        agecl_2_F, agecl_2_M)

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
# Get summary of boar (population by age clas, sex)

get_boar <- function(turtles){

  df <- turtles@.Data %>%
    as.data.frame() %>%
    mutate(agecl = ifelse(agecl == -1, 0, agecl)) %>%   # aanpassen newborns = juvenile !!!!
    group_by(sex, agecl) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(sex = turtles@levels$sex[sex],
           agecl = ageclasses[agecl + 1])
  return(df)
}


#-------------------------------------------------------------
# Set Fertility
#
# Yearly fertility is distributed over months according
# monthly birth rates.
set_F <- function(F = F, birth_month){

  Fm <- matrix(data = c(rep(birth_month$perc / 100, 3)),
               nrow = 12, ncol = 3,
               dimnames = list(NULL, ageclasses))
  Fm <- t(t(Fm) * F)

  return(Fm)
}



#----------------------------------------------------------------------
# Get geboortepiek

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

get_hunting_scen <- function(path){

  f <- function(path, sheet){
    df <- readxl::read_excel(path = path, sheet = sheet)
    return(as.matrix(df))
  }

  Hscen <- path %>%
    readxl::excel_sheets() %>%
    set_names() %>%
    map(f, path = path)
  return(Hscen)
}

#------------------------------------------------------------------

# Run a single simulation to estimate the total simulation time (seconds)

checktime <- function(mylist){

  koffie <- system.time({ sim_boar(init_pop = mylist$init_pop,
                                   max_year = mylist$max_year,
                                   Sm = mylist$Sm,
                                   Fm = mylist$Fm,
                                   world = mylist$world,
                                   Hm = mylist$Hs[[1]]) })
  return(as.double(koffie[3] * nsim * length(Hs)))
}

#-------------------------------------------------------------
# process output - get df_numboar

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
}

#-------------------------------------------------------------
# process output - get df_harvest

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
}


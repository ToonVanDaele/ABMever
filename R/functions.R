### ABM functions

#----------------------------------------------------------------
# initialisation (time base = year)
abm_init_y <- function(init_agecls){
  boar <- createTurtles(n = length(init_agecls), world = dummy,
                        breed = "wildboar",
                        color = rep("red", length(init_agecls)))
  boar <- turtlesOwn(turtles = boar, tVar = "age",
                     tVal = init_agecls)
  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_agecls)
}

#--------------------------------------------------------------
# initialisation (time base = month)

# Initial age = lowest possible age for each class

# Dit geeft in de eerste tijdstappn een scheve leeftijdsdistributie
# Welke alternatieve initiële leeftijdsdistributies zijn mogelijk?

abm_init_m <- function(init_age, world){
  boar <- createTurtles(n = length(init_age), world = dummy,
                        breed = "wildboar",
                        color = rep("red", length(init_age)))
  boar <- turtlesOwn(turtles = boar, tVar = "sex",
                     tVal = sample(c("F", "M"), size = length(init_age), replace = TRUE))
  boar <- turtlesOwn(turtles = boar, tVar = "age",
                     tVal = init_age)

  init_age_cl <- init_age %>%
    as.data.frame() %>%
    mutate(agecl = case_when(. > 24 ~ 2,
                             . > 12 ~ 1,
                             TRUE ~ 0)) %>%
    pull(agecl)

  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_age_cl)
  return(boar)
}

# Numboar init
numboar_init <- function(){
  numboar <- matrix(0, ncol = 6, nrow = max_year * 12 + 1,
                    dimnames = list(NULL,
                                    c("time","year","month", ageclasses)))
  return(numboar)
}

# Numboar init with sex
numboar_init2 <- function(){
  numboar <- matrix(0, ncol = 9, nrow = max_year * 12 + 1,
                    dimnames = list(NULL,
                                    c("time","year","month",
                                      paste0(ageclasses, "F"),
                                      paste0(ageclasses, "M"))))
  return(numboar)
}

# tracker

get_track <- function(turtles){

  df <- turtles@.Data %>%
    as.data.frame() %>%
    mutate(agecl = ifelse(agecl == -1, 0, agecl)) %>%   # aanpassen newborns = juvenile !!!!
    group_by(sex, agecl) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(sex = turtles@levels$sex[sex],
           agecl = ageclasses[agecl + 1]) %>%
    pivot_wider(names_from = c("agecl", "sex"), names_sep = "", values_from = n) %>%
    as.matrix()
  return(df)
}

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



#------------------------------------------------------------
# init H
set_H <- function(){
  Hm <- matrix(data = c(seq(from = 1, to = 12),
                        rep(0.0, 5), rep(0.02, 7),
                        rep(0.05, 3), rep(0.00, 9),
                        rep(0.02, 12),
                        rep(0.01, 5), rep(0.02, 7),
                        rep(0.05, 3), rep(0.00, 9),
                        rep(0.05, 12)),
               nrow = 12, ncol = 7,
               dimnames = list(NULL, c("month", paste0(ageclasses, "F"),
                                       paste0(ageclasses, "M"))))
  return(Hm)
}
#-------------------------------------------------------------
# Set Fertility
#
# Yearly fertility is distributed over months according
# monthly birth rates.
set_F <- function(F = F, csv_filename){

  birth_month <- get_birth_month(csv_filename)

  Fm <- matrix(data = c(rep(birth_month$perc / 100, 3)),
               nrow = 12, ncol = 3,
               dimnames = list(NULL, ageclasses))
  Fm <- t(t(Fm) * F)

  return(Fm)
}


#-------------------------------------------------------------
## natural death (independent of time base)
#
# @param turtles
# @param S vector with survival per age class

# Make sure survival is correct for the given time base

mortality <- function(turtles, S) {
  # Select wildboars (newborns don't die -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # don't include newborns
  age_t2 <- of(agents = t2, var = "agecl")
  tdie <- rbinom(n = NLcount(t2), size = 1, prob = 1 - S[age_t2 + 1])
  who_dies <- who_t2[tdie == 1]    # ID's of dead turtles
  turtles <- die(turtles = turtles, who = who_dies) # remove from list
  return(turtles)
}

#-------------------------------------------------------------------
# hunting (independent of time base)
#
# @param turtles
# @param H vector with hunting probabilities for each age class
#
# Make sure the hunting probability is correct for the given time base
#
hunt <- function(turtles, H, time) {
  # Select wildboars only (newborns are not hunted -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # newborns not included
  age_t2 <- of(agents = t2, var = "agecl") # get ages
  sex_t2 <- of(agents = t2, var = "sex") # get sex
  s <- ifelse(sex_t2 == "F", 1, 4)  # select female or male columns
  tdie <- rbinom(n = NLcount(t2), size = 1, prob = H[age_t2 + s])
  who_dies <- who_t2[tdie == 1]    # ID's of hunted turtles

  # track hunted individuals
  harvest <- NLwith(agents = turtles, var = "who", val = who_dies)
  trackhunt[[time]] <<- get_boar(harvest)

  turtles <- die(turtles = turtles, who = who_dies) # remove from list
  return(turtles)
}

#--------------------------------------------------------------------
## reproduction (independent of time step)
#
# @param turtles turtles object
# @param F vector with reproduction for each class

# Only turtles of age > 10 reproduce
# Newborns have 50/50 sex ratio
# Make sure reproduction values are correct for the given time base

reproduce <- function(turtles = boar, F) {
  # Select female turtles of age > 10 months
  # reproduction (n) according to F value of respective age class
  rturtle <- turtles@.Data %>%
    as.data.frame() %>%
    filter(sex == 1 & age > 10) %>%
    mutate(n = rpois(n = nrow(.), lambda = F[.$agecl + 1])) %>%
    filter(n > 0) %>%
    dplyr::select(who, n)

  # Hatch (add offspring to the population)
  turtles <- hatch(turtles = turtles, who = rturtle$who, n = rturtle$n,
                   breed = "newborn")

  # Set some variable values for the newborns
  newborn <- NLwith(agents = turtles, var = "breed", val = "newborn")
  turtles <- NLset(turtles = turtles, agents = newborn,
                   var = c("age", "agecl", "sex"),
                   val = data.frame(age = -1,
                                    agecl = -1,
                                    sex = sample(c("F", "M"),
                                                 NLcount(newborn), replace = TRUE)))
  return(turtles)
}

#---------------------------------------------------------
## aging (time step = one year)
aging_y <- function(turtles){
  # Newbborns become wildboars
  newborn <- NLwith(agents = turtles, var = "breed", val = "newborn")
  turtles <- NLset(turtles = turtles, agents = newborn,
                   var = c("breed", "agecl"),
                   val = data.frame(breed = "wildboar",
                                    agecl = 0))
  # All age + 1
  turtles@.Data[, "age"] <- turtles@.Data[, "age"] + 1

  # Set age class
  turtles@.Data[turtles@.Data[, "age"] > 12, "agecl"] <- 1
  turtles@.Data[turtles@.Data[, "age"] > 24, "agecl"] <- 2
  return(turtles)
}

#--------------------------------------------------------------------
## aging (time step = one month)
aging_m <- function(turtles){
  # Newbborns become wildboars
  newborn <- NLwith(agents = turtles, var = "breed", val = "newborn")
  turtles <- NLset(turtles = turtles, agents = newborn,
                   var = "breed", val = "wildboar")
  # All age + 1 (month)
  turtles@.Data[, "age"] <- turtles@.Data[, "age"] + 1
  # Set agecl
  turtles@.Data[turtles@.Data[, "age"] <= 12, "agecl"] <- 0
  turtles@.Data[turtles@.Data[, "age"] > 12, "agecl"] <- 1
  turtles@.Data[turtles@.Data[, "age"] > 24, "agecl"] <- 2
  return(turtles)
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

# Run a single simulation to estimate the total simulation time

checktime <- function(mylist){

  koffie <- system.time({ sim_boar(init_age = mylist$init_age,
                                   max_year = mylist$max_year,
                                   S = mylist$S,
                                   Fm = mylist$Fm,
                                   world = mylist$world,
                                   Hm = mylist$Hs[[1]]) })
  return(as.double(koffie[3] * nsim * length(Hs)))
}

#-------------------------------------------------------------
# process output - get df_numboar

get_numboar <- function(mytb){

  df_num <- mytb$result %>%
  map_dfr("df_numboar", .id = "run") %>%
  left_join(mytb %>%
              dplyr::select(Hs, sim, run),
            by = "run") %>%
  mutate(sex = as.factor(sex),
         agecl = as.factor(agecl),
         Hs = as.factor(Hs))
}


# process output - get df_harvest

get_harvest <- function(mytb){

  df_num <- mytb$result %>%
    map_dfr("df_harvest", .id = "run") %>%
    left_join(mytb %>%
                dplyr::select(Hs, sim, run),
              by = "run") %>%
    mutate(sex = as.factor(sex),
           agecl = as.factor(agecl),
           Hs = as.factor(Hs))
}

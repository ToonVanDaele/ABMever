# ABM functions
#
# require NetlogoR
#
#----------------------------------------------------------------------
# Run a single ABM simulation
#
# - time step is months
# - three age classes (juvenile, yearling, adult)
#
# @param Initial population (vector 2 columns (age (months) & sex))
# @param max_month number of months to simulate
# @param start_month month to start the simulation (1..12)
# @param Sm monthly survival for each age class, sex and month (matrix)
# @param Fm monthly fertility for each age class and month (matrix)
# @param Hm monthly hunting ratio or absolute hunting for each age class & sex (matrix)
#
# @return list with 3 data frames:
#              df_numboar = number of individuals (start of time step)
#              df_harvest = number of individuals harvested (end of time step)
#              df_pop = individuals in population at the end of the simulation
#
sim_boar <- function(init_pop, max_month, start_month, Sm, Fm, Hm){

  require(NetLogoR)

  # initialisation
  boar <- abm_init_m(init_pop = init_pop)

  hunt_type <- substr(names(Hm), 1, 1)  # Hunting type can be 'N', 'P', 'R' or 'A'

  tracknum <- trackhunt <- NULL
  time <- 1
  month <- start_month

  while (NLany(boar) & NLcount(boar) < 5000 & time <= max_month) {

    tracknum[[time]] <- get_boar(boar) # track individuals in each age class

    # natural mortality
    boar <- mortality(turtles = boar, S = Sm[month,])

    # reproduction
    boar <- reproduce(turtles = boar, F = Fm[month,])

    # hunting
    who_dies <- hunt(turtles = boar, H = Hm[[1]][month,], hunt_type = hunt_type)
    harvest <- NLwith(agents = boar, var = "who", val = who_dies)
    trackhunt[[time]] <- get_boar_harvest(turtles = harvest) # track hunted
    boar <- die(turtles = boar, who = who_dies) # remove hunted from population

    # aging
    boar <- aging_m(boar)

    # time
    time <- time + 1
    month <- month + 1
    if (month > 12) month <- 1
  }

  if (!NLany(boar)) warning("aborted early: reached 0 boar")
  if (!(NLcount(boar) < 5000)) warning("aborted early: reached max boar (5000)")

  # Process tracking data
  # Set start date
  g <- lubridate::make_date(year = 0, month = start_month, day = 1)

  # Number of individuals
  df_numboar <- tracknum %>%
    map_dfr(rbind, .id = "time") %>%
    mutate(time = as.integer(time)) %>%
    mutate(date = lubridate::add_with_rollback(g, months(time - 1)))

    # mutate(month = as.integer(lubridate::month(date)),
    #        year = as.integer(lubridate::year(date)))

  # harvested individuals
  df_harvest <- trackhunt %>%
    map_dfr(rbind, .id = "time") %>%
    mutate(time = as.integer(time))

  if (!nrow(df_harvest) == 0) {
      df_harvest <- df_harvest %>%
      mutate(date = lubridate::add_with_rollback(g, months(time - 1)))
      # mutate(month = as.integer(lubridate::month(date)),
      #        year = as.integer(lubridate::year(date)))
  }

  # store the whole population after the final simulation time step
  df_pop <- boar@.Data %>%
    as.data.frame() %>%
    dplyr::select(who, breed, sex, age, newb, offspring) %>%
    mutate(sex = boar@levels$sex[sex],
           breed = boar@levels$breed[breed]) %>%
    mutate(sex = as.factor(sex),
           breed = as.factor(breed))

  return(list(df_numboar = df_numboar, df_harvest = df_harvest,
              df_pop = df_pop))
}


#---------------------------------------------------------
# run simulation for 1 or more scenarios
#
# @param max_month Number of months to simulate
# @param init_pop Vector with initial population 2 columns: age (months) & sex
# @param Sm Monthly survival: age class (vector) or ageclass, sex and month (matrix)
# @param Fm Monthly fertility for each age class (vector length 3)
# @param Hs Monthly hunting scenarios (list of matrices)

sim_scen_boar <- function(init_pop, max_month, start_month = 1,
                          Sm, Fm, Hs, nsim){

  # Check if Sm is a vector (length 3) -> make it a matrix (12 x 3)
  if (!is.matrix(Sm)) {
    Smm <- matrix(data = Sm, nrow = 12, ncol = 6, byrow = TRUE)
    colnames(Smm) <- c(paste0(ageclasses, "F"), paste0(ageclasses, "M"))
    Sm <- Smm
  }

  # Create dataframe with simulations to run
  df <- expand.grid(Hs = names(Hs),
                    sim = seq(from = 1, to = nsim),
                    result = list(NULL)) %>%
    mutate(run = row.names(.)) %>%
    as_tibble()

  pb <- txtProgressBar(max = nrow(df), style = 3, width = 50)

    for (i in 1:nrow(df)) {

    outsim <- sim_boar(init_pop = init_pop,
                       max_month = max_month,
                       start_month = start_month,
                       Sm = Sm,
                       Fm = Fm,
                       Hm = Hs[df$Hs[i]])
    df$result[i] <- list(outsim)

    setTxtProgressBar(pb, value = i)
  }
  close(pb)

  return(df)
}

#--------------------------------------------------------------
# Initialisation (time base = Month)
#
# Initialisation of a new NetLogo population
#
# @param init_pop initial population. matrix with 2 columns (age & sex)
#
# @return netlogo population with age, sex, newb and offspring
#
abm_init_m <- function(init_pop){

  # Create world (required, but not used)
  world <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)

  nb <- nrow(init_pop)  # number of individuals
  boar <- createTurtles(n = nb,
                        world = world,
                        breed = "wildboar",
                        color = rep("red", nb))
  boar <- turtlesOwn(turtles = boar,
                     tVar = "sex",
                     tVal = init_pop$sex)
  boar <- turtlesOwn(turtles = boar,
                     tVar = "age",
                     tVal = init_pop$age)
  boar <- turtlesOwn(turtles = boar,
                     tVar = "newb",
                     tVal = 99)
  boar <- turtlesOwn(turtles = boar,
                     tVar = "offspring",
                     tVal = 0)

  init_age_cl <- init_pop$age %>%
    as.data.frame() %>%
    mutate(agecl = case_when(. > 24 ~ 2,
                             . > 12 ~ 1,
                             TRUE ~ 0)) %>%
    pull(agecl)

  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_age_cl)
  return(boar)
}

#-------------------------------------------------------------
## natural death (independent of time base)
#
# @param turtles
# @param S vector with survival per age class

# Make sure survival is correct for the given time step!

mortality <- function(turtles, S) {
  # Select wildboars (newborns don't die -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # don't include newborns
  age_t2 <- of(agents = t2, var = "agecl") # get age class
  sex_t2 <- of(agents = t2, var = "sex") # get sex
  s <- ifelse(sex_t2 == "F", 1, 4)  # select female or male columns

  tdie <- rbinom(n = NLcount(t2), size = 1, prob = 1 - S[age_t2 + s])
  who_dies <- who_t2[tdie == 1]    # ID's of dead turtles
  turtles <- die(turtles = turtles, who = who_dies) # remove from list
  return(turtles)
}

#-------------------------------------------------------------------
# hunting (monthly time step)
#
# @param turtles netlogo population
# @param H vector with hunting for each age class and sex (vector length 8)
# @param time timestep of simulation. Needed to track the harvest
# @param hunt_abs hunting in absolute numbers (TRUE), or probabilities (FALSE) (default)
#
# @return turtles population + extra element in global variabel trackhunt.
#
hunt <- function(turtles, H, hunt_type = "P") {
  # Select wildboars only (newborns are not hunted -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # newborns of current month not included
  age_t2 <- of(agents = t2, var = "agecl") # get age class
  agb_t2 <- of(agents = t2, var = "age") # get age
  sex_t2 <- of(agents = t2, var = "sex") # get sex
  s <- ifelse(sex_t2 == "F", 2, 5)  # select female or male columns

  if (hunt_type == "P"){  # hunting proportions
    tdie <- rbinom(n = NLcount(t2), size = 1, prob = H[age_t2 + s])
    who_dies <- who_t2[tdie == 1]    # ID's of hunted turtles
  }

  if (hunt_type == "A") {  # hunting absolute numbers
    # Get the population with specified age class and sex
    jF <- who_t2[age_t2 == 0 & sex_t2 == "F"]
    jM <- who_t2[age_t2 == 0 & sex_t2 == "M"]
    yF <- who_t2[age_t2 == 1 & sex_t2 == "F"]
    yM <- who_t2[age_t2 == 1 & sex_t2 == "M"]
    aF <- who_t2[age_t2 == 2 & sex_t2 == "F"]
    aM <- who_t2[age_t2 == 2 & sex_t2 == "M"]

    # sample the required number or remove all (if H > population)
    who_dies <- c(sample(jF, size = min(H["juvenileF"], length(jF))),
                  sample(jM, size = min(H["juvenileM"], length(jM))),
                  sample(yF, size = min(H["yearlingF"], length(yF))),
                  sample(yM, size = min(H["yearlingM"], length(yM))),
                  sample(aF, size = min(H["adultF"], length(aF))),
                  sample(aM, size = min(H["adultM"], length(aM))))
  }

  if (hunt_type == "R") {   # Double relative hunting strategy

    # Here comes the double relative hunting strategy!!

  }


  if (hunt_type == "N") {    # No hunting

    who_dies <- numeric()

  }

  return(who_dies)
}

#--------------------------------------------------------------------
## reproduction
#
# Only turtles of age > 10 reproduce
# Newborns have 50/50 sex ratio
# Make sure reproduction values are correct for the given time base
#
# @param turtles turtles object
# @param F vector with reproduction for each class
#
# @return turtles

reproduce <- function(turtles = boar, F) {
  # Select female turtles of age > 10 months
  # reproduction (n) according to F value of respective age class
  who <- turtles[turtles$sex == "F" & turtles$age > 10,
                 c("who", "agecl"), drop = FALSE]

  # If any females of age > 10 in population
  if (nrow(who) > 0) {

    n <- rpois(n = nrow(who), lambda = F[who[, "agecl"] + 1])
    who_repr <- who[n > 0, "who"]

    # Months since newborns = 0 -> to track dependent turtles

    turtles[turtles$who %in% who_repr, "newb"] <- 0
    turtles[turtles$who %in% who_repr, "offspring"] <- n[n > 0]

    # Hatch (add offspring to the population)
    turtles <- hatch(turtles = turtles, who = who_repr, n = n[n > 0],
                   breed = "newborn")
  }

  # Set some variable values for the newborns
  newborn <- NLwith(agents = turtles, var = "breed", val = "newborn")
  # sex 50/50 ratio
  nb_sex <- sample(c("F", "M"), NLcount(newborn), replace = TRUE)

  turtles <- NLset(turtles = turtles, agents = newborn,
                   var = c("age", "agecl", "sex", "offspring", "newb"),
                   val = data.frame(age = -1,
                                    agecl = -1,
                                    sex = nb_sex,
                                    offspring = 0,
                                    newb = 99))
  return(turtles)
}

#--------------------------------------------------------------------
# aging
#
# Increase age, change ageclass and breed type
#
# @param turtles
#
# @return turtles
#
aging_m <- function(turtles){
  # Newbborns become wildboars
  newborn <- NLwith(agents = turtles, var = "breed", val = "newborn")
  turtles <- NLset(turtles = turtles, agents = newborn,
                   var = "breed", val = "wildboar")
  # All age + 1 (month)
  turtles@.Data[, "age"] <- turtles@.Data[, "age"] + 1
  # Set agecl
  turtles@.Data[turtles$age <= 12, "agecl"] <- 0
  turtles@.Data[turtles$age > 12, "agecl"] <- 1
  turtles@.Data[turtles$age > 24, "agecl"] <- 2

  # All months since new newborns + 1
  turtles@.Data[, "newb"] <- turtles@.Data[, "newb"] + 1
  return(turtles)
}

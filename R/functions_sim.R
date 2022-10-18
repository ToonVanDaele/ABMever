# ABM functions
#
# uses Netlogo (NetlogoR)
#
#----------------------------------------------------------------------
# A ABM run
#
#
sim_boar <- function(max_year = max_year, init_pop = init_pop,
                     Hm = Hm, Sm = Sm, Fm = Fm, world = world){

  # initialisation
  boar <- abm_init_m(init_pop = init_pop, world = world)
  tracknum <- NULL
  trackhunt <<- NULL

  time <- 1
  year <- 1
  month <- 1   # In welke maand starten?

  while (NLany(boar) & NLcount(boar) < 5000 & year <= max_year) {
    #  print(paste(year, month))

    # track number of individuals in each age class
    d <- get_boar(boar)
    tracknum[[time]] <- d

    # events
    boar <- hunt(turtles = boar, H = Hm[month,], time)
    boar <- reproduce(turtles = boar, F = Fm[month,])
    boar <- mortality(boar, Sm)
    boar <- aging_m(boar)

    # time
    time <- time + 1
    month <- month + 1
    if (month > 12) {
      month <- 1
      year <- year + 1
    }
  }

  # Process tracking data
  # Number of individuals
  df_numboar <- tracknum %>%
    map_dfr(rbind, .id = "time") %>%
    mutate(time = as.integer(time))
  # harvested individuals
  df_harvest <- trackhunt %>%
    map_dfr(rbind, .id = "time") %>%
    mutate(time = as.integer(time))
  trackhunt <<- NULL

  # store age distribution at the end of the simulation
  age_distr <- boar@.Data %>%
    as.data.frame() %>%
    dplyr::select(age, agecl)

  return(list(df_numboar = df_numboar, df_harvest = df_harvest,
              age_distr = age_distr))
}


#---------------------------------------------------------
# run simulation for 1 or more scenarios

sim_scen_boar <- function(scenlist){

  df <- expand.grid(Hs = names(scenlist$Hs),
                    sim = seq(from = 1, to = nsim),
                    result = list(NULL))
  df <- df %>%
    mutate(run = row.names(.)) %>%
    as_tibble()

  for (i in 1:nrow(df)) {

  outsim <- sim_boar(init_pop = scenlist$init_pop,
                     max_year = scenlist$max_year,
                     Sm = scenlist$Sm,
                     Fm = scenlist$Fm,
                     world = scenlist$world,
                     Hm = scenlist$Hs[[df$Hs[i]]])
  df$result[i] <- list(outsim)
  }
  return(df)
}

#--------------------------------------------------------------
# Initialisation (time base = Month)

abm_init_m <- function(init_pop, world){

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
  who <- turtles@.Data[turtles@.Data[,"sex"] == 1 &
                  turtles@.Data[,"age"] > 10, c("who", "agecl"), drop = FALSE]
  n <- rpois(n = nrow(who), lambda = F[who[, "agecl"] + 1])

  # Hatch (add offspring to the population)
  turtles <- hatch(turtles = turtles, who = who[,"who"], n = n,
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



#----------------------------------------------------------------
# initialisation (time base = year)
abm_init_y <- function(init_pop, world){
  boar <- createTurtles(n = length(init_agecls), world = world,
                        breed = "wildboar",
                        color = rep("red", length(init_agecls)))
  boar <- turtlesOwn(turtles = boar, tVar = "age",
                     tVal = init_agecls)
  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_agecls)
}



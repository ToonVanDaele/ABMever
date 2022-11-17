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
# @param max_month number of months to simulate
# @param init_pop vector with initial population 2 columns: age (months) & sex
# @param Sm monthly survival for each age class (vector)
# @param Fm monthly fertility for each age class (vector)
# @param Hm monthly hunting ratio or absolute hunting for each age class & sex (matrix)
# @param hunt_abs hunting Hm in absolute numbers (default FALSE)
#
sim_boar <- function(max_month = max_month, init_pop = init_pop,
                     Sm = Sm, Fm = Fm, Hm = Hm, hunt_abs = FALSE){

  require(NetLogoR)

  # initialisation
  boar <- abm_init_m(init_pop = init_pop)
  tracknum <- NULL
  trackhunt <<- NULL

  time <- 1
  month <- 1   # Always start 1ste of January?

  while (NLany(boar) & NLcount(boar) < 5000 & time <= max_month) {

    d <- get_boar(boar)   # get number of individuals in each age class
    tracknum[[time]] <- d # store in a list

    # events
    boar <- hunt(turtles = boar, H = Hm[month,], time, hunt_abs = hunt_abs)
    boar <- reproduce(turtles = boar, F = Fm[month,])
    boar <- mortality(boar, Sm)
    boar <- aging_m(boar)

    # time
    time <- time + 1
    month <- month + 1
    if (month > 12) month <- 1
  }

  if (!NLany(boar)) warning("aborted early: reached 0 boar")
  if (!(NLcount(boar) < 5000)) warning("aborted early: reached max boar (5000)")

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

  # store the whole population after the final simulation time step
  df_pop <- boar@.Data %>%
    as.data.frame() %>%
    mutate(sex = boar@levels$sex[sex],
           breed = boar@levels$breed[breed]) %>%
    mutate(sex = as.factor(sex),
           breed = as.factor(breed))

  return(list(df_numboar = df_numboar, df_harvest = df_harvest,
              df_pop = df_pop))
}


#---------------------------------------------------------
# run simulation for 1 or more scenarios

sim_scen_boar <- function(init_pop = init_pop,
                          max_year = max_year,
                          nsim = nsim,
                          Sm = Sm,
                          Fm = Fm,
                          Hs = Hs,
                          hunt_abs = FALSE,
                          dochecktime = FALSE){

  # Estimate runtime before full simulation
  if (dochecktime == TRUE & hunt_abs == FALSE) {
    cat("estimating runtime... \n")
    est_time <- checktime(init_pop = init_pop, max_year = max_year, Sm = Sm,
              Fm = Fm, Hs = Hs, nsim = nsim)
    cat("Estimated runtime (seconds): ", est_time)
  }

  # Create dataframe with simulations to run
  df <- expand.grid(Hs = names(Hs),
                    sim = seq(from = 1, to = nsim),
                    result = list(NULL))
  df <- df %>%
    mutate(run = row.names(.)) %>%
    as_tibble()

  for (i in 1:nrow(df)) {

  outsim <- sim_boar(init_pop = init_pop,
                     max_month = max_year * 12 + 1,
                     Sm = Sm,
                     Fm = Fm,
                     Hm = Hs[[df$Hs[i]]],
                     hunt_abs = hunt_abs)
  df$result[i] <- list(outsim)
  }
  return(df)
}

#--------------------------------------------------------------
# Initialisation (time base = Month)

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
# turtles
# H vector with hunting probabilities of absolute number for each age class
# hunt_abs = TRUE for absolute numbers, FALSE for probabilities
#
hunt <- function(turtles, H, time, hunt_abs = FALSE) {
  # Select wildboars only (newborns are not hunted -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # newborns of current month not included
  age_t2 <- of(agents = t2, var = "agecl") # get age class
  agb_t2 <- of(agents = t2, var = "age") # get age
  sex_t2 <- of(agents = t2, var = "sex") # get sex
  s <- ifelse(sex_t2 == "F", 1, 5)  # select female or male columns
  a <- ifelse(agb_t2 > 60, 1, 0)    # select adult 5 column (instead of adult 3)

  if (hunt_abs == FALSE){  # hunting proportions
    tdie <- rbinom(n = NLcount(t2), size = 1, prob = H[age_t2 + s + a])
    who_dies <- who_t2[tdie == 1]    # ID's of hunted turtles
  }else{    # hunting absolute numbers
    # Get the population with specified age class and sex
    jF <- who_t2[age_t2 == 0 & sex_t2 == "F"]
    jM <- who_t2[age_t2 == 0 & sex_t2 == "M"]
    yF <- who_t2[age_t2 == 1 & sex_t2 == "F"]
    yM <- who_t2[age_t2 == 1 & sex_t2 == "M"]
    a3F <-who_t2[age_t2 == 2 & sex_t2 == "F"]
    a3M <-who_t2[age_t2 == 2 & sex_t2 == "M"]
    a5F <-who_t2[age_t2 == 2 & sex_t2 == "F"]
    a5M <-who_t2[age_t2 == 2 & sex_t2 == "M"]

    # sample the required number or remove all (if H > population)
    who_dies <- c(sample(jF, size = min(H["juvenileF"], length(jF))),
                  sample(jM, size = min(H["juvenileM"], length(jM))),
                  sample(yF, size = min(H["yearlingF"], length(yF))),
                  sample(yM, size = min(H["yearlingM"], length(yM))),
                  sample(a3F, size = min(H["adult3F"], length(a3F))),
                  sample(a3M, size = min(H["adult3M"], length(a3M))),
                  sample(a5F, size = min(H["adult5F"], length(a5F))),
                  sample(a5M, size = min(H["adult5M"], length(a5M))))
  }

  # track hunted individuals
  harvest <- NLwith(agents = turtles, var = "who", val = who_dies)
  trackhunt[[time]] <<- get_boar_harvest(turtles = harvest)

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

  if (nrow(who) > 0) {
    n <- rpois(n = nrow(who), lambda = F[who[, "agecl"] + 1])
    # Months since new juveniles = 0 -> used to track turtles with dependent juveniles
    # whohasnewborn <- NLwith(agents = turtles, var = "who", val = who[n > 0, "who"])
    # NLset(turtles = turtles, agents = whohasnewborn, var = "newb", val = 0)
    # NLset(turtles = turtles, agents = whohasnewborn, var = "offspring", val = n[n>0])
    turtles@.Data[turtles@.Data[,"who"] %in% who[n > 0, "who"],"newb"] <- 0
    turtles@.Data[turtles@.Data[,"who"] %in% who[n > 0, "who"],"offspring"] <- n[n > 0]

    # Hatch (add offspring to the population)
    turtles <- hatch(turtles = turtles, who = who[,"who"], n = n,
                   breed = "newborn")
  }

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

  # All months since new newborns + 1
  turtles@.Data[, "newb"] <- turtles@.Data[, "newb"] + 1
  return(turtles)
}

#----------------------------------------------------------------
# initialisation (time base = year)

abm_init_y <- function(init_pop){

  # Create world (required, but not used)
  world <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)

  boar <- createTurtles(n = length(init_agecls), world = world,
                        breed = "wildboar",
                        color = rep("red", length(init_agecls)))
  boar <- turtlesOwn(turtles = boar, tVar = "age",
                     tVal = init_agecls)
  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_agecls)
}



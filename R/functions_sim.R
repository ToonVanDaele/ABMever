sim_boar <- function(mylist){

  max_year <- mylist[["max_year"]]
  init_age <- mylist[["init_age"]]
  Hm <- mylist[["Hm"]]
  S <- mylist[["S"]]
  Fm <- mylist[["Fm"]]
  world <- mylist[["world"]]

  # initialisation
  boar <- abm_init_m(init_age = init_age, world = world)
  tracknum <- NULL
  trackhunt <<- NULL

  time <- 1
  year <- 1
  month <- 1   # In welke maand starten?

  while (NLany(boar) & NLcount(boar) < 5000 & year <= max_year) {
    #  print(paste(year, month))

    boar <- hunt(turtles = boar, H = Hm[month,], time)
    boar <- reproduce(turtles = boar, F = Fm[month,])
    boar <- mortality(boar, S^(1/12))
    boar <- aging_m(boar)

    # track number of individuals in each age class
    d <- get_boar(boar)
    tracknum[[time]] <- d

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
    as_tibble() %>%
    dplyr::select(age, agecl)

  return(list(df_numboar = df_numboar, df_harvest = df_harvest,
              age_distr = age_distr))
}


#---------------------------------------------------------
# run simulation

sim_scen_boar <- function(scenlist){

  df <- expand.grid(Hs = names(scenlist$Hs),
                    sim = seq(from = 1, to = nsim),
                    result = list(NULL))
  df <- df %>%
    mutate(run = row.names(.)) %>%
    as_tibble()

  for (i in 1:nrow(df)){

  simlist <- list(init_age = scenlist$init_age,
                  max_year = scenlist$max_year,
                  S = scenlist$S,
                  Fm = scenlist$Fm,
                  world = scenlist$world,
                  Hm = scenlist$Hs[[df$Hs[i]]])

  outsim <- sim_boar(simlist)
  df$result[i] <- list(outsim)
  }
  return(df)
}



# sim_scen_boar <- function(scenlist){
#
#   mylist <- NULL
#   for (i in 1:length(scenlist$Hs)){
#
#   simlist <- list(init_age = scenlist$init_age,
#                   max_year = scenlist$max_year,
#                   S = scenlist$S,
#                   Fm = scenlist$Fm,
#                   world = scenlist$world,
#                   Hm = scenlist$Hs[[1]])
#
#   outsim <- map(1:nsim, ~ sim_boar(simlist))
#   mylist[[names(scenlist$Hs[i])]] <- outsim
#   }
#   return(mylist)
# }

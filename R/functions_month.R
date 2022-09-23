### ABM functions

# initialisation
abm_init <- function(init_agecls){
  boar <- createTurtles(n = length(init_agecls), world = dummy,
                        breed = "wildboar",
                        color = rep("red", length(init_agecls)))
  boar <- turtlesOwn(turtles = boar, tVar = "age",
                     tVal = init_agecls * 12)
  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_agecls)
}


## natural death
death <- function(turtles, S) {
  # Select wildboars (newborns don't die -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # don't include newborns
  age_t2 <- of(agents = t2, var = "agecl")
  tdie <- rbinom(n = NLcount(t2), size = 1, prob = 1 - S[age_t2 + 1])
  who_dies <- who_t2[tdie == 1]    # ID's of dead turtles
  turtles <- die(turtles = turtles, who = who_dies) # remove from list
  return(turtles)
}

# hunting
hunt <- function(turtles, H) {
  # Select wildboars (newborns don't die -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # newborns not included
  age_t2 <- of(agents = t2, var = "agecl")
  tdie <- rbinom(n = NLcount(t2), size = 1, prob = H[age_t2 + 1])
  who_dies <- who_t2[tdie == 1]    # ID's of hunted turtles
  turtles <- die(turtles = turtles, who = who_dies) # remove from list
  return(turtles)
}

## reproduction
reproduce <- function(turtles, F) {
  whoTurtles <- of(agents = turtles, var = "who") # get all turtles
  ageTurtle <- of(agents = turtles, var = "agecl")
  # Some reproduce (poisson distribution) -> newborns
  repro <- rpois(n = NLcount(turtles), lambda = F[ageTurtle + 1])

  who_repro <- whoTurtles[repro > 0] # ID's of reproducing turtles
  turtles <- hatch(turtles = turtles, who = who_repro, n = repro[repro > 0],
                   breed = "newborn") # add offspring
  newborn <- NLwith(agents = turtles, var = "breed", val = "newborn")
  turtles <- NLset(turtles = turtles, agents = newborn,
                   var = c("age","agecl"),
                   val = data.frame(age = -1,
                                    agecl = -1))
  return(turtles)
}

## aging
aging <- function(turtles){
  # Newbborns become wildboars
  newborn <- NLwith(agents = turtles, var = "breed", val = "newborn")
  turtles <- NLset(turtles = turtles, agents = newborn,
                   var = "breed", val = "wildboar")
  # All age + 1 (month)
  turtles@.Data[, "age"] <- turtles@.Data[, "age"] + 1
  # Set agecl
  turtles@.Data[turtles@.Data[, "age"] > 12, "agecl"] <- 1
  turtles@.Data[turtles@.Data[, "age"] > 24, "agecl"] <- 2
  return(turtles)
}


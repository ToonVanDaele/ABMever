# Tryout ever met NetlogoR
# Vergelijking van de resultaten met een matrix model en ABM

library(tidyverse)
library(NetLogoR)
library(popbio)

# Population
ageclasses <- c("Juvenile", "Yearling", "Adult")

S <- c(0.6, 0.7, 0.8) # Survival probability
F <- c(0.0, 0.1, 0.5) # Fertility
H <- c(0.0, 0.0, 0.2) # Hunting mortality

nboar <- 1000  # total population at time = 0
tmax <- 20     # time steps (years)

#---------------------------------
# Matrix model
mat <- matrix(c(F[1], F[2],	F[3],
                S[1], 0,  	0,
                0,    S[2],	S[3]), ncol = 3, byrow = TRUE)

mat_h <- t(t(mat) * (1 - H))

init_agecl <- stable.stage(mat_h) * nboar   # initial age class = stable stage
mm <- pop.projection(A = mat_h, n = init_agecl, iterations = tmax)
rownames(mm$stage.vectors) <- ageclasses

#---------------------------------
# ABM model

init_agecls <- c(rep(0, init_agecl[1]),
                 rep(1, init_agecl[2]),
                 rep(2, init_agecl[3]))  # initial ages for all individuals

# Create world (required, but not used)
dummy <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)

# ABM functions

## natural death
death <- function(turtles) {
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
hunt <- function(turtles) {
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
reproduce <- function(turtles) {
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
  # All age + 1
  turtles@.Data[, "age"] <- turtles@.Data[, "age"] + 1
  # All agecl + 1 (max = 2)
  turtles@.Data[, "agecl"] <- turtles@.Data[, "agecl"] + 1
  turtles@.Data[turtles@.Data[,"agecl"] > 2, "agecl"] <- 2
  return(turtles)
}

## The ABM packed in a function
sim <- function(){

  # initialisation
  boar <- createTurtles(n = length(init_agecls), world = dummy,
                        breed = "wildboar",
                        color = rep("red", length(init_agecls)))
  boar <- turtlesOwn(turtles = boar, tVar = "age",
                     tVal = init_agecls)
  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_agecls)

  numboar <- matrix(0, ncol = 4, nrow = tmax,
                    dimnames = list(NULL, c("time", ageclasses)))
  numboar[1,2:4] <- round(init_agecl)

  time <- 1

  while (NLany(boar) & NLcount(boar < 5000) & time < tmax) {

    boar <- hunt(boar)
    boar <- reproduce(boar)
    boar <- death(boar)
    boar <- aging(boar)

    # track number of individuals in each age class
    numboar[time + 1,] <- c(time,
                            sum(boar@.Data[,"agecl"] == 0),
                            sum(boar@.Data[,"agecl"] == 1),
                            sum(boar@.Data[,"agecl"] == 2))

    time <- time + 1
  }

  age_distr <- boar@.Data %>%
    as.data.frame() %>%
    dplyr::select(age, agecl)

  return(list(numboar = as.data.frame(numboar),
              age_distr = age_distr))
}

sim()

# Run simulation 20x
outsim <- rerun(.n = 20, sim())

# Process output ABM
df <- outsim %>%
  map_dfr("numboar", .id = "sim") %>%
  pivot_longer(cols = all_of(ageclasses), names_to = "agecl", values_to = "n")

# Process output matrix model
mms <- mm$stage.vectors %>%
  t() %>%
  as.data.frame() %>%
  mutate(time = row_number() - 1,
         sim = "mm") %>%
  pivot_longer(cols = ageclasses, names_to = "ageclass", values_to = "n" )

# Plot all simulations (thick line = matrix model)
ggplot(df, aes(x = time, y = n, color = agecl, group = paste(sim, agecl))) +
  geom_line() +
  geom_line(data = mms, aes(x = time, y = n,
                            group = ageclass, color = ageclass), size = 1)

# Plot mean of simulations (thick line = matrix model)
df %>%
  group_by(time, agecl) %>%
  summarise(meann = mean(n)) %>%
  ggplot(aes(x = time, y = meann, color = agecl)) + geom_line() +
  geom_line(data = mms, aes(x = time, y = n, color = ageclass), size = 1)

# Plot age distribution ABM
df_age_distr <- outsim %>%
  map_dfr("age_distr", .id = "sim")
ggplot(df_age_distr, aes(x = age, group = sim)) + geom_density()

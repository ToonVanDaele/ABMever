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

abm_init_m <- function(init_agecls){
  boar <- createTurtles(n = length(init_agecls), world = dummy,
                        breed = "wildboar",
                        color = rep("red", length(init_agecls)))
  boar <- turtlesOwn(turtles = boar, tVar = "sex",
                     tVal = sample(c("F", "M"), size = length(init_agecls), replace = TRUE))
  boar <- turtlesOwn(turtles = boar, tVar = "age",
                     tVal = init_agecls * 12)
  boar <- turtlesOwn(turtles = boar, tVar = "agecl",
                     tVal = init_agecls)
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

# init F
set_F <- function(){
  Fm <- matrix(data = c(seq(from = 1, to = 12),
                        rep(0, 12),
                        rep(0.05, 12),
                        rep(0.1, 12)),
               nrow = 12, ncol = 4,
               dimnames = list(NULL, c("month", ageclasses)))
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
hunt <- function(turtles, H) {
  # Select wildboars (newborns don't die -> analoog aan matrix model)
  t2 <- NLwith(agents = turtles, var = "breed", val = "wildboar")
  who_t2 <- of(agents = t2, var = "who") # newborns not included
  age_t2 <- of(agents = t2, var = "agecl") # get ages
  sex_t2 <- of(agents = t2, var = "sex") # get sex
  s <- ifelse(sex_t2 == "F", 2, 5)  # select female of male column
  tdie <- rbinom(n = NLcount(t2), size = 1, prob = H[age_t2 + s])
  who_dies <- who_t2[tdie == 1]    # ID's of hunted turtles
  turtles <- die(turtles = turtles, who = who_dies) # remove from list
  return(turtles)
}

#--------------------------------------------------------------------
## reproduction (independent of time step)
#
# @param turtles turtles object
# @param F vector with reproduction for each class

# Make sure reproduction values are correct for the given time base

reproduce <- function(turtles, F) {
  FTurtles <- NLwith(agents = turtles, var = "sex", val = "F")
  whoTurtles <- of(agents = FTurtles, var = "who") # get all female turtles ID
  ageTurtle <- of(agents = FTurtles, var = "agecl") # get age class
  # Some reproduce (poisson distribution) -> newborns
  repro <- rpois(n = NLcount(FTurtles), lambda = F[ageTurtle + 2])

  who_repro <- whoTurtles[repro > 0] # ID's of reproducing turtles
  turtles <- hatch(turtles = turtles, who = who_repro, n = repro[repro > 0],
                   breed = "newborn") # add offspring
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
                   var = "breed", val = "wildboar")
  # All age + 1
  turtles@.Data[, "age"] <- turtles@.Data[, "age"] + 1
  # All agecl + 1 (max = 2)
  turtles@.Data[, "agecl"] <- turtles@.Data[, "agecl"] + 1
  turtles@.Data[turtles@.Data[,"agecl"] > 2, "agecl"] <- 2
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


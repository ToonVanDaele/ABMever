# Tryout ever met NetlogoR
# Vergelijking van de resultaten met een matrix model en ABM

library(tidyverse)
library(NetLogoR)
library(popbio)
source("R/functions.R")

# Population
ageclasses <- c("Juvenile", "Yearling", "Adult")

S <- c(0.6, 0.7, 0.8) # yearly survival probability
F <- c(0.0, 0.2, 0.5) # yearly fertility
H <- c(0.0, 0.0, 0.0) # yearly hunting mortality

nboar <- 1000    # total population at time = 0
max_year <- 15   # number of years to simulate
nsim <- 20       # number of simulations

#---------------------------------
# Matrix model (yearly)
mat <- matrix(c(F[1], F[2],	F[3],
                S[1], 0,  	0,
                0,    S[2],	S[3]), ncol = 3, byrow = TRUE)

mat_h <- t(t(mat) * (1 - H))

init_agecl <- stable.stage(mat_h) * nboar   # initial age class = stable stage
mm <- pop.projection(A = mat_h, n = init_agecl, iterations = max_year)
rownames(mm$stage.vectors) <- ageclasses

# Process output
mms <- mm$stage.vectors %>%
  t() %>%
  as.data.frame() %>%
  mutate(time = row_number() - 1,
         sim = "mm") %>%
  pivot_longer(cols = all_of(ageclasses),
               names_to = "ageclass", values_to = "n" )

#---------------------------------
# Matrix model (monthly)
Sm <- S^(1/12)
Fm <- F / 12
p <- 11/12
matm <- matrix(c(Fm[1] + Sm[1] * p, Fm[2],	         Fm[3],
                 Sm[1] * (1 - p),   Sm[2] * p,       0,
                 0,                 Sm[2] * (1 - p), Sm[3]), ncol = 3, byrow = TRUE)

matm_h <- t(t(matm) * (1 - H^(1/12)))

#init_agecl <- stable.stage(matm_h) * nboar   # initial age class = stable stage
mmm <- pop.projection(A = matm_h, n = init_agecl, iterations = max_year * 11)
rownames(mmm$stage.vectors) <- ageclasses

# Process output
mmms <- mmm$stage.vectors %>%
  t() %>%
  as.data.frame() %>%
  mutate(time = row_number() - 1,
         sim = "mm") %>%
  pivot_longer(cols = all_of(ageclasses),
               names_to = "ageclass", values_to = "n" )

# Plot matrix yearly and by month
ggplot(mms, aes(x = time, y = n, group = ageclass, color = ageclass)) +
  geom_line(size = 1) +
  geom_line(data = mmms, aes(x = time / 12), size = 0.3)


#---------------------------------
# ABM model - monthly

init_agecls <- c(rep(0, init_agecl[1]),
                 rep(1, init_agecl[2]),
                 rep(2, init_agecl[3]))  # initial ages for all individuals

# Create world (required, but not used)
dummy <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)


## The ABM packed in a function
sim <- function(){

  # initialisation
  boar <- abm_init_m(init_agecls = init_agecls)

  numboar <- matrix(0, ncol = 6, nrow = max_year * 12 + 1,
                    dimnames = list(NULL,
                      c("time","year","month", ageclasses)))
  numboar[1,4:6] <- round(init_agecl)

  time <- 1
  year <- 1
  month <- 1

  while (NLany(boar) & NLcount(boar < 5000) & year <= max_year) {

    #boar <- hunt(boar, H)
    boar <- reproduce(boar, F/12)
    boar <- mortality(boar, S^(1/12))
    boar <- aging_m(boar)

    # track number of individuals in each age class
    numboar[time + 1,] <- c(time, year, month,
                            sum(boar@.Data[,"agecl"] == 0),
                            sum(boar@.Data[,"agecl"] == 1),
                            sum(boar@.Data[,"agecl"] == 2))

    time <- time + 1
    month <- month + 1
    if (month > 12) {
      month <- 1
      year <- year + 1
    }
  }

  # store age distribution at the end of the simulation
  age_distr <- boar@.Data %>%
    as.data.frame() %>%
    dplyr::select(age, agecl)

  return(list(numboar = as.data.frame(numboar),
              age_distr = age_distr))
  cat("-")  # dit werkt niet
}

temp <- sim()

#--------------------------------------------------------------
# Run simulation nsim times
outsim <- rerun(.n = nsim, sim())

# Process output ABM
df <- outsim %>%
  map_dfr("numboar", .id = "sim") %>%
  pivot_longer(cols = all_of(ageclasses), names_to = "agecl", values_to = "n")

# Plot all simulations (thick line = matrix model)
# ggplot(df, aes(x = time, y = n, color = agecl, group = paste(sim, agecl))) +
#   geom_line() +
#   geom_line(data = mms, aes(x = time, y = n,
#                             group = ageclass, color = ageclass), size = 1)

# Plot mean of simulations (thick line = matrix model)
df %>%
  group_by(time, agecl) %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = time, color = agecl)) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  geom_line(data = mms, aes(x = time * 12, y = n,
                            color = ageclass), size = 1)

# Plot age distribution ABM
df_age_distr <- outsim %>%
  map_dfr("age_distr", .id = "sim")
ggplot(df_age_distr, aes(x = age, group = sim)) + geom_density()


#-------------------------------------------------------------
# Hunting and Fertility differentiated in time

Hm <- set_H()
Hm

Fm <- set_F()
Fm

sim_h <- function(){

  # initialisation
  boar <- abm_init_m(init_agecls = init_agecls)
  numboar <- numboar_init2()
  numboar[1, 4:9] <- get_track(boar)

  time <- 1
  year <- 1
  month <- 1   # In welke maand starten?

  while (NLany(boar) & NLcount(boar) < 5000 & year <= max_year) {

    boar <- hunt(turtles = boar, H = Hm[month,])
    boar <- reproduce(boar, Fm[month,])
    boar <- mortality(boar, S^(1/12))
    boar <- aging_m(boar)

    # track number of individuals in each age class
    d <- get_track(boar)
    numboar[time + 1,1:3] <- c(time, year, month)
    for (i in 1:length(d)) numboar[time + 1, colnames(d)[i]] <- d[i]

    time <- time + 1
    month <- month + 1
    if (month > 12) {
      month <- 1
      year <- year + 1
    }
    print(year)
  }

  # store age distribution at the end of the simulation
  age_distr <- boar@.Data %>%
    as.data.frame() %>%
    dplyr::select(age, agecl)

  return(list(numboar = as.data.frame(numboar),
              age_distr = age_distr))
}

out <- sim_h()

df_out <- out$numboar %>%
  pivot_longer(cols = -c("time", "year", "month"), names_to = "agecl", values_to = "n")

ggplot(df_out, aes(x = time, y = n, colour = agecl)) + geom_line()


###---------------------------------

# Vragen

#- In welke maand start het model? Beter jaarkalender of modelkalender?
#- Initiële verdeling van de leeftijden en geslacht
#- Jacht per maand -> hoe intensiteit uitdrukken? % afschot per maand? verschil male female
#- Jacht veranderlijk per gewicht -> groeifunctie?
#- periode tussen twee drachten? een jaar? gemiddeld aantal newborns?
#- realistische periode (maanden) voor F>0
#-

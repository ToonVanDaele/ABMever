# Base model for wild boar with NetlogoR
#
# time base: monthly

library(tidyverse)
library(NetLogoR)
library(popbio)
source("R/functions.R")

# Global population properties
ageclasses <- c("Juvenile", "Yearling", "Adult")

nboar <- 1000    # total population at time = 0
max_year <- 10   # number of years to simulate
nsim <- 10       # number of simulations

# Initial population matrix (yearly)
S <- c(0.6, 0.8, 0.9) # yearly survival probability
F <- c(0.0, 0.2, 0.5) # yearly fertility
H <- c(0.0, 0.0, 0.0) # yearly hunting mortality

# Matrix model
# Used to asses the initial age distribution (based on stable stage)
mat <- matrix(c(F[1], F[2],	F[3],
                S[1], 0,  	0,
                0,    S[2],	S[3]), ncol = 3, byrow = TRUE)
mat_h <- t(t(mat) * (1 - H))

init_agecl <- stable.stage(mat_h) * nboar   # initial number by age class


#----------------------------------------------------------------------
# Initial ABM population

# init_agecls <- c(rep(0, init_agecl[1]),
#                  rep(1, init_agecl[2]),
#                  rep(2, init_agecl[3]))  # initial ages for all individuals

# nog te verbeteren met een fit
init_age <- rgamma(n = nboar, shape = 2, rate = 0.8) * 12

# Create world (required, but not used)
dummy <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)

# Hunting differentiated in time
Hm <- set_H()
Hm

# Jachtdruk bepaald
# Bijkomend differentiëren met leefijd (vb juveniel <3 maand en > 3maand)
# Alternatieve jachtdruk in aantallen per jaar -> iteratief -> quota
# Vooral belangrijk om verschillende scenario's te kunnen vergelijken
# Wel cijfers over verhoudingen bij afschot / niet de percentage tov van wat
# er rondloopt.

# Fertility differentiated in time
Fm <- set_F(F = F)
Fm

# Reproductie op basis van effectieve leeftijd ipv leeftijdsklasse

# vb. reproductie minimum 10 maanden (minimum gewicht) + geboortepiek in maart



sim_h <- function(){

  # initialisation
  boar <- abm_init_m(init_age = init_age)
  tracknum <- NULL
  trackhunt <<- NULL

  time <- 1
  year <- 1
  month <- 1   # In welke maand starten?

  while (NLany(boar) & NLcount(boar) < 5000 & year <= max_year) {
  #  print(paste(year, month))

    boar <- hunt(turtles = boar, H = Hm[month,], time)
    boar <- reproduce(boar, Fm[month,])
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
    as.data.frame() %>%
    dplyr::select(age, agecl)

  return(list(df_numboar = df_numboar, df_harvest = df_harvest,
              age_distr = age_distr))
}

starttime <- Sys.time()
set.seed(5)
out <- sim_h()
endtime <- Sys.time()
endtime - starttime

# Number of individuals
df_num <- out$df_numboar
ggplot(df_num, aes(x = as.integer(time), y = n, color = paste(agecl, sex))) + geom_line()

# Harvested individuals
df_hunt <- out$df_harvest
ggplot(df_hunt, aes(x = as.integer(time), y = n, color = paste(agecl, sex))) + geom_line()

# Total harvested
df_hunt %>%
  group_by(sex, agecl) %>%
  summarise(n = sum(n))


#--------------------------------------------------------------
# Run simulation nsim times
outsim <- rerun(.n = nsim, sim_h())

# Process output ABM
df <- outsim %>%
  map_dfr("df_numboar", .id = "sim")
  #pivot_longer(cols = all_of(ageclasses), names_to = "agecl", values_to = "n")

# Plot all simulations (thick line = matrix model)
# ggplot(df, aes(x = time, y = n, color = agecl, group = paste(sim, agecl))) +
#   geom_line() +
#   geom_line(data = mms, aes(x = time, y = n,
#                             group = ageclass, color = ageclass), size = 1)

# Plot mean of simulations (thick line = matrix model)
# df %>%
#   group_by(time, agecl) %>%
#   summarise(mean = mean(n),
#             p90 = quantile(n, prob = 0.9),
#             p10 = quantile(n, prob = 0.1)) %>%
#   ggplot(aes(x = time, color = agecl)) +
#   geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
#   geom_line(data = mms, aes(x = time * 12, y = n,
#                             color = ageclass), size = 1)

#Plot mean of simulations
df %>%
  group_by(time, agecl, sex) %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = time, color = paste(agecl, sex))) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Plot age distribution ABM
df_age_distr <- outsim %>%
  map_dfr("age_distr", .id = "sim")
ggplot(df_age_distr, aes(x = age, group = sim)) + geom_density()


#-------------------------------------------------------

# Set (hunting) scenario's

# H1 - all equal - all year
Hscen[[1]] <-
# H2 - all equal - hunting season (sep-march)
# H3 - adult and yearling only - all year
# H4 - adult and yearling only - hunting season
# H5 - male adults all year - yearling during hunting season (sep-march)






# run scenario's

# process results


###---------------------------------

# Vragen

#- In welke maand start het model? Beter jaarkalender of modelkalender? Variabel?
#- Initiële verdeling van de leeftijden en geslacht
#- Jacht per maand -> hoe intensiteit uitdrukken? % afschot per maand? verschil male female
#- Jacht veranderlijk per gewicht -> groeifunctie?
#- periode tussen twee drachten? een jaar? gemiddeld aantal newborns?
#- realistische periode (maanden) voor F>0
#- scenarios in excel of google sheets


# later (eventueel) nog toe te voegen:
# F proportie reproductie * worpgrootte *
#? Kans op overleven van moederdier is dit met post-natal survival
# ? cijfers over dracht bij jacht op basis van embryo's

# Scenario's example - no hunting

library(tidyverse)
library(NetLogoR)
source("R/functions.R")
source("R/functions_sim.R")

#---------------------------------
# Set population parameters
ageclasses <- c("Juvenile", "Yearling", "Adult")

# Survival
S <- c(0.6, 0.8, 0.9) # yearly survival probability

# Fertility
F <- c(0, 0.1, 0.5) # yearly fertility
Fm <- set_F(F = F, csv_filename = "./data/input/birth_month.csv") # by month

# Hunting (by month)
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")

# We only use H0 - no hunting
Hs <- Hscen[c("H0")]

# ----------------------------------------------------
# ABM related parameters

nboar0 <- 1000    # initial population size
max_year <- 5    # number of years to simulate
nsim <- 4        # number of simulations per scenario

# Set initial age distribution  (nog aan te passen!!)
init_age <- rgamma(n = nboar0, shape = 2, rate = 0.8) * 12

# Create world (required, but not used)
dummy <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)


#----------------------------------------------------
# Put everything together in a list for multiple scenarios
mypop <- list(init_age = init_age,
              max_year = max_year,
              nsim = nsim,
              S = S,
              Fm = Fm,
              Hs = Hs,
              world = dummy)

# --------------------------------------------------
# run a single simulation to estimate required time (seconds)
checktime(mypop)

#-----------------------------------------------------
# run full simulation
scen_H0 <- sim_scen_boar(mypop)
# store output
saveRDS(scen_H0, file = "./data/interim/scen_H0.RDS")

#-----------------------------------------------------
# process results
#
# Alle output (list) wordt bewaard in een tibble
# met rbind_rows kunnen resultaten van vroegere simulaties worden samengevoegd.
# Zo moeten simulaties niet steeds opnieuw uitgevoerd worden.

df_num <- get_numboar(scen_H0)
df_har <- get_harvest(scen_H0)

#----------------------------------------------------
# plot

# Time series number of individuals by age class and hunting scenario
df_num %>%
  filter(sex == "F") %>%
  group_by(time, agecl, Hs) %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = time, color = paste(agecl, Hs))) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity")


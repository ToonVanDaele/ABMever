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
Sm <- S^(1/12)        # monthly survival probability

# Fertility
F <- c(0, 0.1, 0.5) # yearly fertility
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting (by month)
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")

# We only use scenario "N" - no hunting
Hs <- Hscen[c("N")]

# ----------------------------------------------------
# ABM related parameters

nboar0 <- 1000    # initial population size
max_year <- 5    # number of years to simulate
nsim <- 4        # number of simulations per scenario

# Set initial age distribution  (nog aan te passen!!)
init_pop <- set_init_pop(init_agecl = c(200, 100, 500),
                         birth_month = birth_month, Sm = Sm)

#-----------------------------------------------------
# run full simulation
scen_N <- sim_scen_boar(init_pop = init_pop,
                        max_year = max_year,
                        nsim = nsim,
                        Sm = Sm,
                        Fm = Fm,
                        Hs = Hs)
# store output
saveRDS(scen_N, file = "./data/interim/scen_N.RDS")

#-----------------------------------------------------
# process results
#
# Alle output (list) wordt bewaard in een tibble
# met rbind_rows kunnen resultaten van vroegere simulaties worden samengevoegd.
# Zo moeten simulaties niet steeds opnieuw uitgevoerd worden.

df_num <- get_numboar(scen_N)

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


# Scenario's example

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

# H0 - no hunting
# H1 - all equal - all year
# H2 - all equal - hunting season (sep-march)
# H3 - adult and yearling only - all year
# H4 - adult and yearling only - hunting season
# H5 - male adults all year - yearling during hunting season (sep-march)
# ...

# Hunting scenarios are stored in a list of dataframes
# Here we only use H1 & H3
Hs <- Hscen[c("H1", "H3")]

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
scen1 <- sim_scen_boar(mypop)
# store output
saveRDS(scen1, file = "./data/interim/scen1.RDS")

#-----------------------------------------------------
# process results
#
# Alle output (list) wordt bewaard in een tibble
# met rbind_rows kunnen resultaten van vorige simulaties toevoegen.
scen0 <- readRDS(file = "./data/interim/scen_H0.RDS") #scenario no hunting

scen <- scen1 %>%
  bind_rows(scen0)

df_num <- get_numboar(scen)
df_har <- get_harvest(scen)

#----------------------------------------------------
# plot

# Time series: number of female individuals by age class and hunting scenario
df_num %>%
  filter(sex == "F") %>%
  group_by(time, agecl, Hs) %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, group = paste(agecl, Hs), color = agecl, linetype = Hs)) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Total harvest by age class and scenario
df_har %>%
  group_by(agecl, Hs, sim) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = paste(agecl, Hs), y = mean)) + geom_point() +
      geom_errorbar(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Cumulative harvest by age class and scenario
df_har %>%
  group_by(agecl, Hs, time) %>%
  summarise(mean = mean(n), .groups = "drop_last") %>%
  ungroup() %>%
  complete(time, Hs, agecl, fill = list(mean = 0)) %>%
  group_by(agecl, Hs) %>%
  mutate(n = cumsum(mean)) %>%
  ggplot(aes(x = time, y = n, color = agecl, linetype = Hs)) +
  geom_line()


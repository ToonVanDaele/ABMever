# Optimalisation examples

# This script searches for a parameter combinations (for hunting) with
# an expectance value for lambda of 1.

# Scenario's example

library(tidyverse)
library(NetLogoR)
source("R/functions.R")
source("R/functions_sim.R")

#---------------------------------
# Set population parameters
ageclasses <- c("Juvenile", "Yearling", "Adult")

# Survival
S <- c(0.81, 0.876, 0.876) # yearly survival probability (Toigo 2008)
Sm <- S^(1/12)        # monthly survival probability

# Fertility
F <- c(0, 0.1, 0.5) # yearly fertility
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting (by month)
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")

# H0 - no hunting
Hm <- Hscen[[c("H0")]]

# Projection matrix (for stable stage)
mat <- set_projmat_post(S = S, F = F, H = H)

# ----------------------------------------------------
# ABM related parameters

nboar0 <- 100    # initial population size
max_year <- 5    # number of years to simulate
nsim <- 4        # number of simulations per scenario

# Set initial age distribution
init_agecl <- stable.stage(mat$mat) * nboar0
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# Create world (required, but not used)
world <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)

#----------------------------------------------------
# Put everything together in a list for multiple scenarios
mypop <- list(init_pop = init_pop,
              max_year = max_year,
              nsim = nsim,
              Sm = Sm,
              Fm = Fm,
              Hs = Hs,
              world = world)

# --------------------------------------------------
# run a single simulation to estimate required time (seconds)
checktime(mypop)

#-----------------------------------------------------
# run full simulation and store results
scen1 <- sim_scen_boar(mypop)
saveRDS(scen1, file = "./data/interim/scen1.RDS")

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

df_num %>%
  group_by(sim, Hs, time) %>%
  summarise(ntot = sum(n), .groups = "drop_last") %>%
  view()
  ggplot(aes(x = time, y = ntot, colour = as.factor(sim))) + geom_line()

calc_lambda(df_num)

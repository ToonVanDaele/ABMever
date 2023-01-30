# Scenario's example - no hunting
# Example

library(tidyverse)
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

Hscen
# We only use scenario "N" - no hunting
Hs <- Hscen[c("N")]

# ----------------------------------------------------
# ABM related parameters

nboar0 <- 500    # initial population size

# Set initial age distribution  (nog aan te passen!!)
init_pop <- set_init_pop(init_agecl = c(0.3, 0.3, 0.4) * nboar0,
                         birth_month = birth_month, Sm = Sm)

#-----------------------------------------------------
# run full simulation
scen_N <- sim_scen_boar(init_pop = init_pop,
                        max_month = 5 * 12 + 1,
                        nsim = 4,
                        Sm = Sm,
                        Fm = Fm,
                        Hs = Hs)
# store output
saveRDS(scen_N, file = "./data/interim/scen_N.RDS")

scen_N$result[[1]]$df_numboar
scen_N$result[[1]]$df_harvest
scen_N$result[[1]]$df_pop

#-----------------------------------------------------
# process results
#
# Alle output (list) wordt bewaard in een tibble
# met rbind_rows kunnen resultaten van vroegere simulaties worden samengevoegd.
# Zo moeten simulaties niet steeds opnieuw uitgevoerd worden.

df_num <- get_numboar(scen_N, df = "df_numboar")

#----------------------------------------------------
# plot

# Time series number of individuals by age class and hunting scenario
df_num %>%
  filter(sex == "F") %>%
  group_by(time, agecl, Hs) %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, color = paste(agecl, Hs))) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity")


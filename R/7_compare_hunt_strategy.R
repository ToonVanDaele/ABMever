# Met welk scenario is het totaal afschot het kleinst en wat is het verschil
# tussen jacht het hele jaar door en een beperkt jachtseizoen.

library(tidyverse)
library(NetLogoR)
library(popbio)
source("R/functions.R")
source("R/functions_sim.R")
source("R/functions_matrix.R")

#---------------------------------
# Overall demographic parameters
ageclasses <- c("juvenile", "yearling", "adult")
nboar0 <- 1000    # initial population size

# Survival
S <- c(0.81, 0.876, 0.876) # yearly survival probability (Toigo 2008)
Sm <- S^(1/12)             # monthly survival probability

# Fertility
e <- c(4.21, 5.44, 6.12)    # average number of embryos  (tabel 3)
r <- c(0.5, 0.9, 0.95)      # proportion females reproducing    (tabel 3)
F <- r * e                  # Fertility (yearly)
# load percentage births by month
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
# Caclulate monthly fertility
Fm <- set_F(F = F, birth_month = birth_month)

# Hunting
# Load hunting scenario's from excel sheet
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")
names(Hscen)

# Select "H0"
H0 <- Hscen[c("H0")]  # Select hunting scenario H0

# Run an ABM simulation for 4 years, required to find stable stage distribution

# Create projection matrix without hunting (for stable stage as initial guess)
mat <- set_projmat_post(S = S, F = F, H = c(0, 0, 0))
init_agecl <- stable.stage(mat$mat) * nboar0     # Stable stage distribution
# Initial population (age, sex)
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# run simulation and store results
scen_7a <- sim_scen_boar(init_pop = init_pop,
                         max_year = 4,
                         nsim = 5,
                         Sm = Sm,
                         Fm = Fm,
                         Hs = H0)

saveRDS(scen_7a, file = "./data/interim/scen_7a.RDS")
#scen_7a <- readRDS(file = "./data/interim/scen_7a.RDS")

# Get the end population of the ABM simulation
# and sample 1000 individuals as initial population
# Sample is taken from all 5 repeated simulations

df_pop <- get_pop(scen_7a)

df_init_pop <- df_pop %>%
  sample_n(size = 1000, replace = TRUE) %>%
  dplyr::select(age, sex)

# We now run simulations with new hunting strategies
Hnew <- Hscen[c("H0", "H1", "H2", "H3")]

scen_7b <- sim_scen_boar(init_pop = df_init_pop,
                         max_year = 10,
                         nsim = 5,
                         Sm = Sm,
                         Fm = Fm,
                         Hs = Hnew,
                         dochecktime = TRUE)

saveRDS(scen_7b, file = "./data/interim/scen_7b.RDS")
#scen_7b <- readRDS(file = "./data/interim/scen_7b.RDS")
df_num <- get_numboar(scen_7b)
df_har <- get_harvest(scen_7b)

view(df_har)

df_har <- scen_7b$result[[1]]$df_harvest

# Population
df_num %>%
  group_by(time, Hs, sim) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  summarise(mean_n = mean(tot),
            p90 = quantile(tot, prob = 0.9),
            p10 = quantile(tot, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_n, colour = Hs)) +
  geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Age class distribution by year (1st of January)
df_num %>%
  filter(time %in% seq(from = 1, to = 121, by = 12)) %>%
  group_by(time, Hs, sim, agecl) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  mutate(rel_n = tot / sum(tot)) %>%
  group_by(time, Hs, agecl) %>%
  summarise(mean_rel_n = mean(rel_n),
            p90 = quantile(rel_n, prob = 0.9),
            p10 = quantile(rel_n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_rel_n, colour = agecl, shape = Hs)) + geom_line() +
  geom_point()

# Impact on harvest
df_har %>%
  group_by(agecl, Hs, sim) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = paste(agecl, Hs), y = mean)) + geom_point() +
  geom_errorbar(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Harvest during the first 5 years
df_har %>%
  filter(time < 61) %>%
  mutate(year = floor((time - 1) / 12) + 1) %>%
  group_by(agecl, Hs, year, sim) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = year, y = mean, colour = agecl)) + geom_point() + geom_line() +
  facet_wrap(~Hs)

# Cumulative harvest by age class and scenario
# ! H0 max 5000 is reached at +/- 60 months
df_har %>%
  group_by(Hs, time, sim) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  summarise(mean_n = mean(n), .groups = "drop_last") %>%
  mutate(n = cumsum(mean_n)) %>%
  ggplot(aes(x = time, y = n, color = Hs)) + geom_line()


## # Beginnen we met de stable stage zonder jacht of stable stage met jacht?
#
# Een afschot van een bepaald jaar -> hoe hetzelfde afschot verdelen in het jaar
# ttz een jachtseizoen introduceren.
#
# Heeft het zin om alle afschot voor de geboortepiek te concentreren?
#
# zomerafschot is om schade te vermijden, minder om de populatie te reguleren
# ? wat is de impact van het zomerafschot op

# Waar je controle over hebt (qua jacht) zijn de absolute aantallen,
# de relatieve verdeling over het jaar en tussen de leeftijdsklassen.
# Er is geen controle over de ratios afschot t.a.v. de werkelijke populatie.



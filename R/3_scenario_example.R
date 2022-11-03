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
Sm <- S^(1/12)        # monthly survival probability

# Fertility
F <- c(0, 0.1, 0.5) # yearly fertility
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting (by month)
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")

# N  - no hunting
# T1 - all equal - all year
# T2 - all equal - hunting season (sep-march)
# T3 - adult and yearling only - all year

# Hunting scenarios are stored in a list of dataframes
# Here we only use H1 & H3
Hs <- Hscen[c("N", "T1", "T3")]

# ----------------------------------------------------
# ABM related parameters

nboar0 <- 500    # initial population size
max_year <- 5    # number of years to simulate
nsim <- 4        # number of simulations per scenario

# Set initial population with age distribution
init_pop <- set_init_pop(init_agecl = c(150, 150, 200),
                         birth_month = birth_month, Sm = Sm)

#-----------------------------------------------------
# run full simulation and store results
scen1 <- sim_scen_boar(init_pop = init_pop,
                       max_year = max_year,
                       nsim = nsim,
                       Sm = Sm,
                       Fm = Fm,
                       Hs = Hs,
                       dochecktime = TRUE)
saveRDS(scen1, file = "./data/interim/scen1.RDS")
#scen1 <- readRDS(file = "./data/interim/scen1.RDS")

scen1
head(scen1$result[[1]][[1]])

#-----------------------------------------------------
# process results
#
# Alle output (list) wordt bewaard in een tibble.
# Met rbind_rows voegen we resultaten van vorige simulaties toe.
scenN <- readRDS(file = "./data/interim/scen_N.RDS") #scenario no hunting

scen <- scen1 %>%
  bind_rows(scenN)
scen

df_num <- get_numboar(scen)
df_har <- get_harvest(scen1)


df_har %>%
  filter(time > 23 & time < 36 & Hs == "T1" & sim == 1 & agecl == "Adult" & sex == "F")


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

# Number of dependent juveniles by scenario
df_har %>%
  mutate(year = floor(time / 12) + 1) %>%
  group_by(Hs, year, sim) %>%
  summarise(dep_juv = sum(dep_juv)) %>%
  summarise(m_dep_juv = mean(dep_juv)) %>%
  ggplot(aes(x = year, y = m_dep_juv, colour = Hs)) + geom_line()


#----------------------------------------------------
# Scenario example absolute numbers

# Load hunting scenarios with absolute numbers
Hscen_abs <- get_hunting_scen(path = "./data/input/hunting_scenarios_abs.xlsx")

names(Hscen_abs)

#-----------------------------------------------------
# run full simulation and store results
scen_abs <- sim_scen_boar(init_pop = init_pop,
                       max_year = max_year,
                       nsim = nsim,
                       Sm = Sm,
                       Fm = Fm,
                       Hs = Hscen_abs,
                       hunt_abs = TRUE,     # Hunting is in absolute numbers!!
                       dochecktime = TRUE)
saveRDS(scen_abs, file = "./data/interim/scen_abs.RDS")
#scen_abs <- readRDS(file = "./data/interim/scen_abs.RDS")

#-----------------------------------------------------
# process results
#
df_num <- get_numboar(scen_abs)

# Total for a specific time, sim and scenario
df_num %>%
  filter(time == 24 & sim == 1 & Hs == "abs1") %>%
  summarise(n = sum(n))


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


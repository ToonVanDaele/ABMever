# Scenario's example

library(tidyverse)
source("R/functions.R")
source("R/functions_sim.R")

#---------------------------------
# Set population parameters

ageclasses <- c("Juvenile", "Yearling", "Adult")

# Survival
S <- c(0.6, 0.8, 0.9) # yearly survival probability
Sm <- S^(1/12)        # monthly survival probability

# Fertility (yearly and distribution by month)
F <- c(0, 0.1, 0.5) # yearly fertility by age class
# Distribution of fertility by month
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
# Distribution by month and age class
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting (by month)
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")

# N  - no hunting
# P_ = proportional by month, age class & sex
# A_ = absolute number by month, age class & sex
# R_ = 2 step proportional first by population, then age class & sex

# P_T1 - all equal - all year
# P_T2 - all equal - hunting season (sep-march)
# P_T3 - adult and yearling only - all year

# Hunting scenarios are stored in a list of dataframes
# Here we select some
Hs <- Hscen[c("N", "P_T1", "P_T2", "R_1", "A_1")]

# ----------------------------------------------------
# ABM related parameters
nboar0 <- 500    # initial population size

# Set initial population with age distribution
init_pop <- set_init_pop(init_agecl = c(0.3, 0.3, 0.4) * nboar0,
                         birth_month = birth_month, Sm = Sm)

#-----------------------------------------------------
# run full simulation and store results
scen1 <- sim_scen_boar(init_pop = init_pop,
                       max_month = 5 * 12 + 1,  # 5 years
                       start_month = 1,          # start in January
                       nsim = 4,
                       Sm = Sm,
                       Fm = Fm,
                       Hs = Hs)
saveRDS(scen1, file = "./data/interim/scen1.RDS")
#scen1 <- readRDS(file = "./data/interim/scen1.RDS")

scen1
head(scen1$result[[1]]$df_numboar)
head(scen1$result[[2]]$df_harvest)

#-----------------------------------------------------
# process results
#
# Alle output (list) wordt bewaard in een tibble.

df_num <- get_numboar(scen1, df = "df_numboar")
df_har <- get_numboar(scen1, df = "df_harvest")

#----------------------------------------------------
# Number of individuals by hunting scenario
df_num %>%
  group_by(date, Hs, sim) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = date, group = Hs, linetype = Hs)) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Total harvest by age class and scenario
df_har %>%
  group_by(agecl, Hs, sim) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = paste(agecl, Hs), y = mean)) + geom_point() +
      geom_errorbar(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Harvest by age class and time for R_1
df_har %>%
  filter(Hs == "R_1") %>%
  group_by(date, agecl, sim) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = date, y = mean, group = agecl, color = agecl)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  scale_x_date(date_labels =  "%Y")

df_har %>%
  filter(Hs == "P_T1") %>%
  group_by(date, agecl, sim) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = date, y = mean, group = agecl, color = agecl)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  scale_x_date(date_labels =  "%Y")


# Cumulative harvest by scenario
df_har %>%
  group_by(Hs, time, sim) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  summarise(mean = mean(n), .groups = "drop_last") %>%
  ungroup() %>%
  complete(time, Hs, fill = list(mean = 0)) %>%
  group_by(Hs) %>%
  mutate(n = cumsum(mean)) %>%
  ggplot(aes(x = time, y = n, linetype = Hs)) +
  geom_line()

# Mean number of dependent juveniles by scenario
df_har %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(Hs, year, sim) %>%
  summarise(dep_juv = sum(offspring),.groups = "drop_last") %>%
  summarise(m_dep_juv = mean(dep_juv), .groups = "drop") %>%
  ggplot(aes(x = year, y = m_dep_juv, colour = Hs)) + geom_line()


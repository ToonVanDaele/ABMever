# Assessement of time till stable stage for a hunting scenario

library(tidyverse)
source("R/functions.R")
source("R/functions_sim.R")
source("R/functions_matrix.R")

#---------------------------------
# Overall demographic parameters
ageclasses <- c("juvenile", "yearling", "adult")
nboar0 <- 100    # initial population size
max_year <- 15    # number of years to simulate

# Survival
S <- c(0.81, 0.876, 0.876) # yearly survival probability (Toigo 2008)
Sm <- S^(1/12)             # monthly survival probability

# Fertility
e <- c(4.21, 5.44, 6.12)    # average number of embryos  (tabel 3)
r <- c(0.5, 0.9, 0.95)      # proportion females reproducing    (tabel 3)
F <- r * e
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # Fertility by Month

# Hunting
# Load hunting scenario's from excel sheet
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")
names(Hscen)

Hs <- Hscen[c("N", "P_1", "P_2")]  # Select 3 different hunting scenarios

# We run a ABM simulation to see how long it takes till a stable stage
# distribution is reached for this hunting scenario.

# Create projection matrix (for stable stage as initial population)
mat <- set_projmat_post(S = S, F = F, H = c(0, 0, 0))

# Set initial age distribution
init_pop <- set_init_pop2(n = nboar0)

# run simulation and store results
scen_h1 <- sim_scen_boar(init_pop = init_pop,
                          n_month = max_year * 12 + 1,
                          nsim = 5,
                          Sm = Sm,
                          Fm = Fm,
                          Hs = Hs)

saveRDS(scen_h1, file = "./data/interim/scen_h1.RDS")
#scen_h1 <- readRDS(file = "./data/interim/scen_h1.RDS")

df_num <- get_numboar(scen_h1, df = "df_numboar")

# Population
df_num %>%
  group_by(time, Hs, sim) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  summarise(mean_n = mean(tot),
            p90 = quantile(tot, prob = 0.9),
            p10 = quantile(tot, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_n, colour = Hs)) +
    geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Age class distribution
df_num %>%
  group_by(time, sim, Hs, agecl) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  mutate(rel_n = tot / sum(tot)) %>%
  group_by(time, Hs, agecl) %>%
  summarise(mean_rel_n = mean(rel_n),
            p90 = quantile(rel_n, prob = 0.9),
            p10 = quantile(rel_n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_rel_n, colour = agecl, shape = Hs)) + geom_line() +
  geom_point() + geom_errorbar(aes(ymax = p90, ymin = p10))

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
  geom_point() + geom_errorbar(aes(ymax = p90, ymin = p10), width = 0.8)

# After a few years the age distribution stabilizes. We run the model again
# for 4 years (till January 1st) for hunting scenario 'N' and use a sample of
# the end population as an initial population for the further analysis.

N <- Hscen["N"]

# run simulation and store results
scen_h2 <- sim_scen_boar(init_pop = init_pop,
                         n_month = 4 * 12 + 1,
                         nsim = 5,
                         Sm = Sm,
                         Fm = Fm,
                         Hs = N)

saveRDS(scen_h2, file = "./data/interim/scen_h2.RDS")
#scen_h2 <- readRDS(file = "./data/interim/scen_h2.RDS")

df_num <- get_numboar(scen_h2, df = "df_numboar")

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
  geom_point() + geom_errorbar(aes(ymax = p90, ymin = p10))

# Get the population and sample 1000 individuals

df_pop <- get_pop(scen_h2)

df_init_pop <- df_pop %>%
  sample_n(size = 1000, replace = TRUE) %>%
  dplyr::select(age, sex)

# We run again with the new initial population. The age class distribution
# should be stable from the beginning

scen_h3 <- sim_scen_boar(init_pop = df_init_pop,
                         n_month = 10 * 12 + 1,
                         nsim = 5,
                         Sm = Sm,
                         Fm = Fm,
                         Hs = N)

saveRDS(scen_h3, file = "./data/interim/scen_h3.RDS")
#scen_h3 <- readRDS(file = "./data/interim/scen_h3.RDS")

df_num <- get_numboar(scen_h3, df = "df_numboar")

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
  geom_point() + geom_errorbar(aes(ymax = p90, ymin = p10))

# Comparision of initial and end population age distribution

df_plot_init <- init_pop %>%
  group_by(age) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(p = n / sum(n))

ggplot(data = df_plot_init, aes(x = age, y = p)) + geom_point()

df_pop <- get_pop(mytb = scen_h3)

df_pop %>%
  group_by(sim, age) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(p = n / sum(n)) %>%
  group_by(age) %>%
  summarise(p_mean = mean(p),
            p90 = quantile(p, prob = 0.9),
            p10 = quantile(p, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = age, y = p_mean)) + geom_point() +
    geom_errorbar(aes(ymax = p90, ymin = p10)) +
    geom_point(data = df_plot_init, aes(x = age, y = p), color = "blue")


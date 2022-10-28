# Introducing variation of hunting effort during the year.
# Assessement of impact on total hunted individuals

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
max_year <- 5    # number of years to simulate

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

H1 <- Hscen[c("H1")]  # Select current hunting scenario ("H1")

# We run a ABM simulation to define the stable stage distribution with this
# hunting scenario.

nsim <- 5        # number of simulations per scenario

# Create projection matrix (for stable stage as initial population)
mat <- set_projmat_post(S = S, F = F, H = c(0, 0, 0))

# Set initial age distribution
init_agecl <- stable.stage(mat$mat) * nboar0
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# run simulation and store results
scen_h1 <- sim_scen_boar(init_pop = init_pop,
                          max_year = max_year,
                          nsim = nsim,
                          Sm = Sm,
                          Fm = Fm,
                          Hs = H1,
                          dochecktime = TRUE)

saveRDS(scen_h1, file = "./data/interim/scen_h1.RDS")
#scen_h1 <- readRDS(file = "./data/interim/scen_h1.RDS")

df_num <- get_numboar(scen_h1)

# Population
df_num %>%
  group_by(time, sim) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  summarise(mean_n = mean(tot),
            p90 = quantile(tot, prob = 0.9),
            p10 = quantile(tot, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_n)) +
    geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Age class distribution
df_num %>%
  group_by(sim, time, agecl) %>%
  summarise(tot = sum(n),
            rel_n = tot / sum(tot), .groups = "drop_last")

  summarise(rel_n = tot / sum(tot), .groups = "drop_last") %>%
  summarise(mean_n = mean(tot),
            p90 = quantile(tot, prob = 0.9),
            p10 = quantile(tot, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_n, colour = Hs)) +
  geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

#----------------------------------------------------
# 6.3.2.  Exclusieve jacht

# We define scenarios with increasing hunting pressure 0 -> 0.9 (yearly basis)
# on each age class  (3 * 10 scenarios)

df <- expand.grid(agecl = ageclasses,
                  Hy = c(seq(from = 0, to = 0.9, by = 0.1), 0.95))

# On a monthly basis this means...
df$Hm <- 1 - ((1 - df$Hy)^(1/12))

# Helpfunction to set values in a list of matrices
fdd2 <- function(value, agecl){
  mm <- N                                              # "N" scenario (as template)
  mycol <- colnames(N)[str_detect(colnames(N), agecl)] # select ageclass columns
  mm[,mycol] <- value                                  # assign value
  return(mm)
}

Hsel <- map2(.x = df$Hm, .y = as.character(df$agecl), .f = fdd2) # apply function
names(Hsel) <- paste(df$agecl, df$Hy, sep = "_")  # give list elements a name

# run simulation and store results
scen_sel <- sim_scen_boar(init_pop = init_pop,
                          max_year = max_year,
                          nsim = nsim,
                          Sm = Sm,
                          Fm = Fm,
                          Hs = Hsel,
                          dochecktime = TRUE)
saveRDS(scen_sel, file = "./data/interim/scen_sel.RDS")
#scen_sel <- readRDS(file = "./data/interim/scen_sel.RDS")

df_num <- get_numboar(scen_sel)

#----------------------------------------------------
# Figuur 15 Exclusieve jacht
df_num %>%
  calc_lambda(burnin = 24) %>%
  mutate(temp = str_split_fixed(string = Hs, pattern = "_", n = 2)) %>%
  mutate(agecl = temp[,1],
         Hscen = temp[,2]) %>%
  dplyr::select(-"temp") %>%
  group_by(Hscen, agecl) %>%
  summarise(lambda = mean(gm_lambda_y),
            p90 = quantile(gm_lambda_y, prob = 0.9),
            p10 = quantile(gm_lambda_y, prob = 0.1), .groups = "drop") %>%
  mutate(Hscen = as.numeric(Hscen)) %>%
  ggplot(aes(x = Hscen, y = lambda, color = agecl)) +
  geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  scale_y_continuous(limits = c(0, 2.5)) +
  geom_hline(yintercept = 1)

#-------------------------------------------------
# Hoofdstuk 7 - Afschotregimes en absolute aantallen

# Select 'afschotscenarios' (tabel 15, p. 46)
Hs <- Hscen[c("N", "H0", "H1", "H2", "H3", "H4", "H5", "H6")]

# ABM related parameters
nboar0 <- 1000    # initial population size
max_year <- 10    # number of years to simulate
nsim <- 5         # number of simulations per scenario

# run simulation and store results
scen_ch7 <- sim_scen_boar(init_pop = init_pop,
                          max_year = max_year,
                          nsim = nsim,
                          Sm = Sm,
                          Fm = Fm,
                          Hs = Hs,
                          dochecktime = TRUE)

saveRDS(scen_ch7, file = "./data/interim/scen_ch7.RDS")
#scen_ch7 <- readRDS(file = "./data/interim/scen_ch7.RDS")

df_num <- get_numboar(scen_ch7)

# Population in time for each scenario
df_num %>%
  group_by(time, Hs, sim) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  summarise(mean_n = mean(tot), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_n, colour = Hs)) + geom_line()

# Show lambda for each hunting scenario
df_num %>%
  calc_lambda(burnin = 24) %>%
  group_by(Hs) %>%
  summarise(lambda = mean(gm_lambda_y, na.rm = TRUE),
            p90 = quantile(gm_lambda_y, prob = 0.9, na.rm = TRUE),
            p10 = quantile(gm_lambda_y, prob = 0.1, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Hs, y = lambda)) + geom_point() +
  geom_errorbar(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  scale_y_continuous(limits = c(0, 2.5)) +
  geom_hline(yintercept = 1)

# Number of individuals after 2, 5 & 10 years
df_num %>%
  filter(time %in% c(24, 60, 120)) %>%
  group_by(time, Hs, sim) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  summarise(n_mean = median(n), .groups = "drop") %>%
  mutate(year = time / 12) %>%
  dplyr::select(year, Hs, n_mean) %>%
  pivot_wider(names_from = Hs, values_from = n_mean)

df_num %>%
  filter(time %in% c(24, 60, 120)) %>%
  group_by(time, Hs, sim) %>%
  summarise(n = sum(n), .groups = "drop_last") %>%
  summarise(n_mean = median(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1), .groups = "drop") %>%
  mutate(year = as.factor(time / 12)) %>%
  ggplot(aes(x = Hs, y = n_mean, fill = year)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymax = p90, ymin = p10), width = 0.5, position = position_dodge(.9))

# Maximum age after 10 years
df_age <- get_agedistr(scen_ch7)

df_age %>%
  mutate(age_y = round(age / 12, 0)) %>%
  group_by(Hs, age_y, sim) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  summarise(n = mean(n), .groups = "drop_last") %>%
  mutate(rel_n = n / sum(n)) %>%
  ggplot(aes(x = age_y, y = rel_n, group = Hs, colour = Hs)) + geom_line() +
  scale_x_continuous(limits = c(0, 15))

# max age
df_age %>%
  group_by(Hs, sim) %>%
  summarise(max_age = max(age), .groups = "drop_last") %>%
  summarise(max_age = mean(max_age)) %>%
  mutate(age_y = round(max_age / 12, 0))


# Beginnen we met de stable stage zonder jacht of stable stage met jacht?
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



# Reconstruction of Pallemaerts 2022 with ABM

# Comparision of the ABM models with results of population matrix

library(tidyverse)
library(popbio)
source("R/functions.R")
source("R/functions_sim.R")
source("R/functions_matrix.R")

#---------------------------------
# Set overall demographic parameters
ageclasses <- c("juvenile", "yearling", "adult")
nboar0 <- 100    # initial population size
max_year <- 10    # number of years to simulate

# Demographic parameters (table 10, p.30)

S <- c(0.697, 0.353, 0.367) # yearly survival (including hunting)

e <- c(4.21, 5.44, 6.12)    # average number of embryos
r <- c(0.5, 0.9, 0.95)      # proportion of females reproducing
F <- r * e                  # yearly fertility

H <- c(0,0,0)               # no hunting. (hunting is included in S)
#-----------------------------------------------------------
# Set projection matrix (yearly)
mat <- set_projmat_post(S = S, F = F, H = H)

# Set initial ageclasses as stable stage
init_agecl <- stable.stage(mat$mat) * nboar0

# run population model
# init_agecl / 2 as we use a female only model
out_m <- matrix_proj(A = mat$mat_h, n = init_agecl / 2, iterations = max_year + 1)

# Plot
ggplot(data = out_m, aes(x = time, y = n, color = agecl, group = agecl)) +
  geom_line()

lambda(mat$mat)    # ok p.30
stable.stage(mat$mat)  # ok figuur 8

#--------------------------------------------
# Parameters for ABM model

# Surival monthly
Sm <- S^(1/12)
# Fertility monthly according the proportions of births by month
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting monthly
# Load hunting scenario's from excel sheet
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")
# Select scenario "N" (no hunting) (Hunting is already included in S)
Hs <- Hscen[c("N")]

# Create initial population
# According stable stage and births by month
init_pop <- set_init_pop(init_agecl = init_agecl, birth_month = birth_month, Sm = Sm)

init_pop %>%
  as_tibble() %>%
  ggplot(aes(x = age / 12)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(breaks = seq(0, 15,1))

#-----------------------------------------------------
# run ABM model
scen_h4_3 <- sim_scen_boar(init_pop = init_pop,
                           max_month = max_year * 12 + 1,
                           nsim = 5,
                           Sm = Sm,
                           Fm = Fm,
                           Hs = Hs)
# store output
saveRDS(scen_h4_3, file = "./data/interim/scen_h4_3.RDS")
#scen_h4_3 <- readRDS(file = "./data/interim/scen_h4_3.RDS")

#-----------------------------------------------------
# process results

df_num <- get_numboar(scen_h4_3)

#----------------------------------------------------
# plot matrix + abm

df_num %>%
  filter(sex == "F") %>%
  group_by(time, agecl, Hs) %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, color = paste(agecl, Hs))) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  geom_line(data = out_m,
            aes(x = (time - 1) * 12 + 1, y = n, color = agecl, group = agecl)) +
  geom_point(data = out_m,
             aes(x = (time - 1) * 12 + 1, y = n, color = agecl, group = agecl))

# Total population
df_num %>%
  group_by(time, sim) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  summarise(pop = mean(tot), .groups = "drop") %>%
  ggplot(aes(x = time, y = pop)) + geom_line()

#------------------------------------------------------------
# hoofdstuk 6
#
# 6.3.1 Niet selectieve jacht  (figuur 14)

# Survival
S <- c(0.81, 0.876, 0.876) # yearly survival probability (Toigo 2008)
Sm <- S^(1/12)             # monthly survival probability

# Fertility
e <- c(4.21, 5.44, 6.12)    # average number of embryos  (tabel 3)
r <- c(0.5, 0.9, 0.95)      # proportion females reproducing    (tabel 3)
F <- r * e
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting
# Get 0 hunting as a template
N <- Hscen[[c("N")]]

# We define 10 scenarios with increasing hunting pressure 0 -> 0.9 (yearly basis)
Hy <- seq(from = 0, to = 0.9, by = 0.1)

# On a monthly basis this means...
Hm <- 1 - ((1 - Hy)^(1/12))

# Helpfunction to set the values in a list of matrices.
fdd <- function(value){
  mm <- N                # use "N" scenario as template
  mm[,] <- value         # assign value
  return(mm)
}

Hs <- map(.x = Hm, .f = fdd)   # apply help function to each element of Hm
names(Hs) <- Hy  # give the elements in the list a name

# ABM related parameters

nboar0 <- 100    # initial population size
max_year <- 5    # number of years to simulate

# Create projection matrix (for stable stage as initial population)
mat <- set_projmat_post(S = S, F = F, H = c(0, 0, 0))

# Set initial age distribution
init_agecl <- stable.stage(mat$mat) * nboar0
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# run simulation and store results
scen_int <- sim_scen_boar(init_pop = init_pop,
                          max_month = max_year * 12 + 1,
                          nsim = 4,
                          Sm = Sm,
                          Fm = Fm,
                          Hs = Hs,
                          dochecktime = TRUE)

saveRDS(scen_int, file = "./data/interim/scen_int.RDS")
#scen_int <- readRDS(file = "./data/interim/scen_int.RDS")

df_num <- get_numboar(scen_int)

#----------------------------------------------------
# Figuur 6.3.1. Niet selectieve jacht
df_num %>%
  calc_lambda() %>%
  group_by(Hs) %>%
  summarise(mean = mean(gm_lambda_y),
            p90 = quantile(gm_lambda_y, prob = 0.9),
            p10 = quantile(gm_lambda_y, prob = 0.1), .groups = "drop") %>%
  mutate(Hs = as.numeric(levels(Hs))[Hs]) %>%
  ggplot(aes(x = Hs, y = mean)) +
  geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  geom_hline(yintercept = 1)

#----------------------------------------------------
# 6.3.2.  Exclusieve jacht

# We define scenarios with increasing hunting pressure 0 -> 0.9 (yearly basis)
# in each age class  (3 * 10 scenarios)

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
                          max_month = max_year * 12 + 1,
                          nsim = 4,
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

# Set initial age distribution
init_agecl <- stable.stage(mat$mat) * nboar0
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# run simulation and store results
scen_ch7 <- sim_scen_boar(init_pop = init_pop,
                          max_month = max_year * 12 + 1,
                          nsim = 5,
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
df_pop <- get_pop(scen_ch7)

df_pop %>%
  mutate(age_y = round(age / 12, 0)) %>%
  group_by(Hs, age_y, sim) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  summarise(n = mean(n), .groups = "drop_last") %>%
  mutate(rel_n = n / sum(n)) %>%
  ggplot(aes(x = age_y, y = rel_n, group = Hs, colour = Hs)) + geom_line() +
  scale_x_continuous(limits = c(0, 15))

# max age
df_pop %>%
  group_by(Hs, sim) %>%
  summarise(max_age = max(age), .groups = "drop_last") %>%
  summarise(max_age = mean(max_age)) %>%
  mutate(age_y = round(max_age / 12, 0))



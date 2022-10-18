# Reconstruction of Pallemaerts 2022

# Comparision of models with results in matrix and ABM

library(tidyverse)
library(NetLogoR)
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

S <- c(0.697, 0.353, 0.367) # jaarlijkse overleving (hunting included)

e <- c(4.21, 5.44, 6.12)    # gemiddeld aantal embryo's
r <- c(0.5, 0.9, 0.95)      # proportie reproducerend
F <- r * e                  # jaarlijkse fertiliteit

H <- c(0,0,0)               # geen hunting. (Zit mee in S)
#-----------------------------------------------------------
# Set projection matrix (yearly)
mat <- set_projmat_post(S = S, F = F, H = H)

# Set initial ageclasses as stable stage
init_agecl <- stable.stage(mat$mat) * nboar0

# run population model
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
# Fertility monthly --  nog te bekijken!!
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting monthly
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")
# We only use H0 - no hunting (Hunting is here included in the survival parameter)
Hs <- Hscen[c("H0")]

nsim <- 5   # number of simulations

# Initial population
init_agecl <- stable.stage(mat$mat) * nboar0
# Set initial age distribution  (nog aan te passen!!)
init_pop <- set_init_pop(init_agecl = init_agecl, birth_month = birth_month, Sm = Sm)

hist(init_pop$age / 12)

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
# run full simulation
scen_h4_3 <- sim_scen_boar(mypop)
# store output
saveRDS(scen_h4_3, file = "./data/interim/scen_h4_3.RDS")

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
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = time, color = paste(agecl, Hs))) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  geom_line(data = out_m, aes(x = (time - 1) * 12 + 1, y = n, color = agecl, group = agecl)) +
  geom_point(data = out_m, aes(x = (time - 1) * 12 + 1, y = n, color = agecl, group = agecl))


# Total population
df_num %>%
  group_by(time, sim) %>%
  summarise(tot = sum(n)) %>%
  ggplot(aes(x = time, y = tot, colour = as.factor(sim))) + geom_line()


#------------------------------------------------------------
# hoofdstuk 6
#
# 6.3.1 Niet selectieve jacht  (figuur 14)

# Survival
S <- c(0.81, 0.876, 0.876) # yearly survival probability (Toigo 2008)
Sm <- S^(1/12)        # monthly survival probability

# Fertility
e <- c(4.21, 5.44, 6.12)    # gemiddeld aantal embryo's  (tabel 3)
r <- c(0.5, 0.9, 0.95)      # proportie reproducerend    (tabel 3)
F <- r * e
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month


# Hunting
# Get 0 hunting as a template
H0 <- Hscen[[c("H0")]]

# We define 10 scenarios with increasing hunting pressure 0 -> 0.9 (yearly basis)
Hy <- seq(from = 0, to = 0.9, by = 0.1)

# On a monthly basis this means...
Hm <- 1 - ((1 - Hy)^(1/12))

# Helpfunctie om de waarden in een lijst van matrices te zetten
fdd <- function(value){
  mm <- H0               # use H0 scenario as template
  mm[,] <- value         # assign value
  return(mm)
}

Hs <- map(.x = Hm, .f = fdd)
names(Hs) <- Hy  # give the elements in the list a name

# ----------------------------------------------------
# ABM related parameters

nboar0 <- 100    # initial population size
max_year <- 5    # number of years to simulate
nsim <- 4        # number of simulations per scenario

# Create projection matrix (for stable stage as initial population)
mat <- set_projmat_post(S = S, F = F, H = c(0, 0, 0))

# Set initial age distribution
init_agecl <- stable.stage(mat$mat) * nboar0
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# Put everything together in a list for multiple scenarios
mypop <- list(init_pop = init_pop,
              max_year = max_year,
              nsim = nsim,
              Sm = Sm,
              Fm = Fm,
              Hs = Hs,
              world = world)

checktime(mypop)

# run full simulation and store results
scen_int <- sim_scen_boar(mypop)
saveRDS(scen_int, file = "./data/interim/scen_int.RDS")

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
#
# 6.3.2.  Exclusieve jacht

# We define 3 * 10 scenarios with increasing hunting pressure 0 -> 0.9 (yearly basis)
# On each age class

df <- expand.grid(agecl = ageclasses,
                  Hy = seq(from = 0, to = 0.9, by = 0.1))

# On a monthly basis this means...
df$Hm <- 1 - ((1 - df$Hy)^(1/12))

# Helpfunctie om de waarden in een lijst van matrices te zetten
fdd2 <- function(value, agecl){
  mm <- H0                                 # use H0 scenario as template
  mycol <- paste0(agecl, c("F", "M"))  # in the ageclass column for M & F
  mm[,mycol] <- value               # assign value
  return(mm)
}

Hsel <- map2(.x = df$Hm, .y = df$agecl, .f = fdd2)
names(Hsel) <- paste(df$agecl, df$Hy, sep = "_")  # give the elements in the list a name

# ----------------------------------------------------
# ABM related parameters

nboar0 <- 100    # initial population size
max_year <- 5    # number of years to simulate
nsim <- 4        # number of simulations per scenario

# Create projection matrix (for stable stage as initial population)
mat <- set_projmat_post(S = S, F = F, H = c(0, 0, 0))

# Set initial age distribution
init_agecl <- stable.stage(mat$mat) * nboar0
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# Put everything together in a list for multiple scenarios
mypop <- list(init_pop = init_pop,
              max_year = max_year,
              nsim = nsim,
              Sm = Sm,
              Fm = Fm,
              Hs = Hsel,
              world = world)

checktime(mypop)

# run full simulation and store results
scen_sel <- sim_scen_boar(mypop)
saveRDS(scen_sel, file = "./data/interim/scen_sel.RDS")

df_num <- get_numboar(scen_sel)

#----------------------------------------------------

# Figuur 15 Exclusieve jacht

df_num %>%
  calc_lambda() %>%
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


#

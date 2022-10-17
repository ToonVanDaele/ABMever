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
ageclasses <- c("Juvenile", "Yearling", "Adult")
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
  ggplot(aes(x = time, y = tot, group = sim)) + geom_line()



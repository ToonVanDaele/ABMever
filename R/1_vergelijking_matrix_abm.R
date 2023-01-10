# Compare population matrix and individual based model

# This script compares the results of both approaches with
# similar population parameters

library(tidyverse)
library(popbio)
source("R/functions.R")
source("R/functions_sim.R")
source("R/functions_matrix.R")

#---------------------------------
# Set overall demographic parameters
ageclasses <- c("Juvenile", "Yearling", "Adult")
nboar0 <- 1000    # initial population size
max_year <- 15    # number of years to simulate

# Survival rate by age class
S <- c(0.6, 0.8, 0.9) # yearly survival probability
# Fertility rate by age class
F <- c(0.1, 0.2, 0.5) # yearly fertility
# Harvest rate by age class
H <- c(0, 0, 0)

#-----------------------------------------------------------
# Set projection matrix (yearly)
mat <- set_projmat_post(S = S, F = F, H = H)

# Set initial ageclasses as stable stage without hunting
# in order to be able to compare different huntig scenarios later on.
init_agecl <- stable.stage(mat$mat) * nboar0

# run population model
out_m <- matrix_proj(A = mat$mat_h, n = init_agecl / 2, iterations = max_year + 1)

# Plot
ggplot(data = out_m, aes(x = time, y = n, color = agecl, group = agecl)) +
  geom_line()

#--------------------------------------------------
# Parameters for ABM model

# Surival monthly
Sm <- S^(1/12)
# Fertility monthly --  nog te bekijken!!
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Hunting monthly
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")
# We only use scenario "N" - no hunting
Hs <- Hscen[c("N")]

# Set initial age distribution  (nog aan te passen!!)
init_pop <- set_init_pop(init_agecl = init_agecl, birth_month = birth_month, Sm = Sm)

hist(init_pop$age / 12)

#-----------------------------------------------------
# run full simulation
scen_comp <- sim_scen_boar(init_pop = init_pop,
                           max_year = max_year,
                           nsim = 5,
                           Sm = Sm,
                           Fm = Fm,
                           Hs = Hs,
                           dochecktime = TRUE)

# store output
saveRDS(scen_comp, file = "./data/interim/scen_comp.RDS")
#scen_comp <- readRDS("./data/interim/scen_comp.RDS")

#-----------------------------------------------------
# process results
#
df_num <- get_numboar(scen_comp)

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
  geom_line(data = out_m, aes(x = (time - 1) * 12 + 1, y = n, color = agecl, group = agecl)) +
  geom_point(data = out_m, aes(x = (time - 1) * 12 + 1, y = n, color = agecl, group = agecl))


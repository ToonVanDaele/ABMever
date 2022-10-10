# Compare population matrix and individual based model

# This script compares the results of both approaches with
# similar population parameters

library(tidyverse)
library(NetLogoR)
library(popbio)
source("R/functions.R")
source("R/functions_sim.R")
source("R/functions_matrix.R")

#---------------------------------
# Set overall demographic parameters
ageclasses <- c("Juvenile", "Yearling", "Adult")
nboar0 <- 1000    # initial population size
max_year <- 5    # number of years to simulate

# Survival rate
S <- c(0.6, 0.8, 0.9) # yearly survival probability
# Fertility rate
F <- c(0, 0, 0.5) # yearly fertility
# Harvest rate
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
Fm <- set_F(F = F, csv_filename = "./data/input/birth_month.csv") # by month
# reproduction only in month 12
# Fm <- matrix(data = 0, nrow = 12, ncol = 3, dimnames = list(NULL, ageclasses))
# Fm[12,] <- F
Fm

# geometrische gemiddelde
exp(mean(log(Fm[,"Adult"]))) * 12
rpois(n = 12, lambda = Fm[,"Adult"])

# Hunting monthly
# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")
# We only use H0 - no hunting
Hs <- Hscen[c("H0")]

nsim <- 5   # number of simulations

# Set initial age distribution  (nog aan te passen!!)
#init_pop <- set_init_pop(init_agecl = init_agecl)
init_pop <- set_init_pop2(init_agecl = init_agecl)

hist(init_pop$age / 12)

# Create world (required, but not used)
world <- createWorld(minPxcor = -5, maxPxcor = 5, minPycor = -5, maxPycor = 5)

#----------------------------------------------------
# Put everything together in a list for multiple scenarios
mypop <- list(init_pop = init_pop,
              max_year = max_year,
              nsim = nsim,
              S = S,
              Fm = Fm,
              Hs = Hs,
              world = world)

# --------------------------------------------------
# run a single simulation to estimate required time (seconds)
checktime(mypop)

#-----------------------------------------------------
# run full simulation
scen_comp <- sim_scen_boar(mypop)
# store output
saveRDS(scen_comp, file = "./data/interim/scen_comp.RDS")

#-----------------------------------------------------
# process results
#
df_num <- get_numboar(scen_comp)
#df_har <- get_harvest(scen_comp)

#----------------------------------------------------
# plot

# Time series number of individuals by age class and hunting scenario
df_num %>%
  filter(sex == "F") %>%
  group_by(time, agecl, Hs) %>%
  summarise(mean = mean(n),
            p90 = quantile(n, prob = 0.9),
            p10 = quantile(n, prob = 0.1)) %>%
  ggplot(aes(x = time, color = paste(agecl, Hs))) +
  geom_smooth(aes(y = mean, ymax = p90, ymin = p10), size = 0.5, stat = "identity")

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



df_num %>%
  filter(sim == 1) %>%
  arrange(agecl, sex) %>%
  view()


100 * 0.25

mpois <- function(a){
  sum(rpois(n = 100, lambda = 0.25))
}

df <- data.frame(n = 1:1000)
df <- df %>%
  rowwise() %>%
  mutate(pp = mpois())
hist(df$pp)
mean(df$pp)
median(df$pp)



# Optimalisation examples

# This script searches for a parameter combinations (for hunting) with
# an expectance value for lambda of 1.

# Scenario's example

library(tidyverse)
library(NetLogoR)
source("R/functions.R")
source("R/functions_sim.R")
source("R/functions_matrix.R")

#---------------------------------
# Set population parameters
ageclasses <- c("Juvenile", "Yearling", "Adult")

# Survival
S <- c(0.81, 0.876, 0.876) # yearly survival probability (Toigo 2008)
Sm <- S^(1/12)        # monthly survival probability

# Fertility
e <- c(4.21, 5.44, 6.12)    # gemiddeld aantal embryo's
r <- c(0.5, 0.9, 0.95)      # proportie reproducerend
F <- r * e
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
Fm <- set_F(F = F, birth_month = birth_month) # by month

# Load all hunting scenario's from excel sheets
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")

# N - no hunting
N <- Hscen[[c("N")]]

# We define 10 scenarios with increasing hunting pressure 0 -> 0.9
# 0 till 0.9 on a yearly basis
Hy <- seq(from = 0, to = 0.9, by = 0.1)
# This means on a monthly basis...
Hm <- 1 - ((1 - Hy)^(1/12))

# Helpfunctie om de waarden in een lijst van matrices te zetten
fdd <- function(value){

  mm <- N                  # use N (no hunting) scenario as template
  mm[,] <- round(value, 4) # assign value
  return(mm)
}

Hs <- map(.x = Hm, .f = fdd)
names(Hs) <- Hy  # give the list a name

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

#-----------------------------------------------------
# run full simulation and store results
scen_int <- sim_scen_boar(init_pop = init_pop,
                          max_year = max_year,
                          nsim = nsim,
                          Sm = Sm,
                          Fm = Fm,
                          Hs = Hs,
                          dochecktime = TRUE)
saveRDS(scen_int, file = "./data/interim/scen_int.RDS")

df_num <- get_numboar(scen_int)
df_har <- get_harvest(scen_int)

#----------------------------------------------------

# Time series: number individuals by hunting scenario

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


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
nboar0 <- 1000    # initial population size
max_year <- 5    # number of years to simulate

# Demographic parameters (table 10, p.30)

S <- c(0.697, 0.353, 0.367) # jaarlijkse overleving (hunting included)

e <- c(4.21, 5.44, 6.12)    # gemiddeld aantal embryo's
r <- c(0.5, 0.9, 0.95)      # proportie reproducerend
F <- r * e                  # jaarlijkse fertiliteit

H <- c(0,0,0)               # geen hunting. Zit al in S
#-----------------------------------------------------------
# Set projection matrix (yearly)
mat <- set_projmat_post(S = S, F = F, H = H)

lambda(mat$mat)    # ok p.30
stable.stage(mat$mat)  # ok figuur 8

# Set initial ageclasses as stable stage without hunting
# in order to be able to compare different huntig scenarios later on.





#--------------------------------------------------------------
# Vetter et al


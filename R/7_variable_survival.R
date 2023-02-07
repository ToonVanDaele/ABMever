# Met welk scenario is het totaal afschot het kleinst en wat is het verschil
# tussen jacht het hele jaar door en een beperkt jachtseizoen.
#
# Startdatum modellering is 1 april

library(tidyverse)
source("R/functions.R")
source("R/functions_sim.R")
source("R/functions_matrix.R")

#---------------------------------
# Overall demographic parameters
ageclasses <- c("juvenile", "yearling", "adult")
nboar0 <- 100    # initial population size

# Survival
S <- c(0.81, 0.876, 0.876) # yearly survival probability (Toigo 2008)
Sm <- S^(1/12)             # monthly survival probability

S2 <- get_survival("./data/input/survival.xlsx")
S_toigo <- S2[["Toigo"]]

S_toigo

# Fertility
e <- c(4.21, 5.44, 6.12)    # average number of embryos  (tabel 3)
r <- c(0.5, 0.9, 0.95)      # proportion females reproducing    (tabel 3)
F <- r * e                  # Fertility (yearly)
# load percentage births by month
birth_month <- get_birth_month(csv_filename = "./data/input/birth_month.csv")
# Caclulate monthly fertility
Fm <- set_F(F = F, birth_month = birth_month)

# Hunting
# Load hunting scenario's from excel sheet
Hscen <- get_hunting_scen(path = "./data/input/hunting_scenarios.xlsx")
N <- Hscen["N"]  # Select hunting scenario N

# Create projection matrix without hunting (for stable stage as initial guess)
mat <- set_projmat_post(S = S, F = F, H = c(0, 0, 0))
init_agecl <- stable.stage(mat$mat) * nboar0     # Stable stage distribution

# Initial population (age, sex)
init_pop <- set_init_pop(init_agecl = init_agecl,
                         birth_month = birth_month, Sm = Sm)

# run simulation and store results
scen_8 <- sim_scen_boar(init_pop = init_pop,
                        n_month = 4 * 12 + 1,
                        start_month = 4,
                        nsim = 5,
                        Sm = S_toigo,    # variable survival -> matrix!!!
                        Fm = Fm,
                        Hs = N)

saveRDS(scen_8, file = "./data/interim/scen_8.RDS")
#scen_8 <- readRDS(file = "./data/interim/scen_8.RDS")

df_num <- get_numboar(scen_8, df = "df_numboar")

# Population (by model time)
df_num %>%
  group_by(time, Hs, sim) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  summarise(mean_n = mean(tot),
            p90 = quantile(tot, prob = 0.9),
            p10 = quantile(tot, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = time, y = mean_n, colour = Hs)) +
  geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity")

# Population (by calender date)
df_num %>%
  filter(time < 17) %>%
  group_by(date, Hs, sim) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%
  summarise(mean_n = mean(tot),
            p90 = quantile(tot, prob = 0.9),
            p10 = quantile(tot, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = date, y = mean_n, colour = Hs)) +
  geom_smooth(aes(ymax = p90, ymin = p10), size = 0.5, stat = "identity") +
  scale_x_date(date_labels = "%m/%y", date_breaks = "1 month")

# Age class distribution by year (at 1st of January)
df_num %>%
  filter(lubridate::month(date) == 1) %>%
  group_by(date, Hs, sim, agecl) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%  # sum female & male
  mutate(rel_n = tot / sum(tot)) %>%                  # proportions agecl
  group_by(date, Hs, agecl) %>%
  summarise(mean_rel_n = mean(rel_n),
            p90 = quantile(rel_n, prob = 0.9),
            p10 = quantile(rel_n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = date, y = mean_rel_n, colour = agecl, shape = Hs)) +
  geom_line() + geom_point() + xlab("year")

# Age class distribution by year (at 1st of April)
df_num %>%
  filter(lubridate::month(date) == 4) %>%
  group_by(date, Hs, sim, agecl) %>%
  summarise(tot = sum(n), .groups = "drop_last") %>%  # sum female & male
  mutate(rel_n = tot / sum(tot)) %>%                  # proportions agecl
  group_by(date, Hs, agecl) %>%
  summarise(mean_rel_n = mean(rel_n),
            p90 = quantile(rel_n, prob = 0.9),
            p10 = quantile(rel_n, prob = 0.1), .groups = "drop") %>%
  ggplot(aes(x = date, y = mean_rel_n, colour = agecl, shape = Hs)) +
  geom_line() + geom_point() +
  scale_x_date(date_labels = "%m/%y", date_breaks = "1 year")


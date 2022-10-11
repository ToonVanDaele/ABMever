# map purrr simulations

df1 <- data.frame(rb = rbinom(n = 5, size = 1, prob = 0.5),
                  ne = runif(n = 5, min = 5, max = 10))

scen <- c("H1", "H2")
sim <- seq(from = 1, to = 2)

df_sims <- expand.grid(scen = scen, sim = sim)
df_sims <- df_sims %>%
  mutate(name = paste(scen, sim))

out <- df_sims %>%
  rowwise() %>%
  mutate(z = list(list(df1 = df1, df2 = df1)))


out$z[[1]][["df2"]]

map_df(.x = out$z, .f = "df2", .id = "sim")

map_df(.x = out, .f = ~ bind_rows(.x, .id = 'scen'), .id = 'sim')


unnest(out, z)


do.call("rbind", lapply(out$z, "[[", 1) )

lapply(out, "[[", 1) %>% bind_rows


for(i in 1:nrow(out)){

}


mutate(df = map(.x = , .f = , ))


df_sims %>%
  map_dfr("z", .id = "name")

map_df(.x = df_sims, .f = )




df <- data.frame(X1=1:3, X2=4:6, X3=7:9)
FINAL <- list(Result=list(Smalldf1=df, Mat1=as.matrix(df)),
           Result=list(Smalldf1=df+1, Mat1=as.matrix(df+1)))

do.call("rbind", lapply(FINAL, "[[", 1) )


#### # Or doing it with dplyr:

library(dplyr)
lapply(FINAL, "[[", 1) %>% bind_rows

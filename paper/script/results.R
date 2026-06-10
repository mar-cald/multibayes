rm(list=ls())

# load pkgs
pkg <- c("dplyr","tidyr","purrr")
invisible(sapply(pkg, require, character.only = T))

# custom functions
source("paper/script/utl.R")

# Results simulation 

load("paper/script/output/sim_pi0.rda")
sim_tb <- tibble(sim)

df_pi0 <- sim_tb |>
  mutate(metrics = purrr::map2(res, pi0, ~ compute_metrics(.x, pi0 = .y),
                               .progress = "Calculating metrics...")) |>
  dplyr::select(-res) |>
  unnest(metrics)

save(df_pi0, file = "paper/script/output/df_pi0.rda")




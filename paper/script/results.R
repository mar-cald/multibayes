rm(list=ls())

# load pkgs
pkg <- c("dplyr","tidyr","purrr")

invisible(sapply(pkg, require, character.only = T))

# custom functions
source("paper/script/utl.R")

# Results simulation 1------

# Results simulation 3------
load("paper/script/output/sim_pi0_opt.rda")
sim_tb <- tibble(sim)

df_pi0_opt <- sim_tb |>
  mutate(metrics = purrr::map(res, compute_metrics, 
                              .progress = "Calculating metrics...")) |>
  dplyr::select(-res) |> 
  unnest(metrics)

save(df_pi0_opt, file = "paper/script/output/df_pi0_opt.rda")


# Results simulation 2-----

load("paper/script/output/sim_pi0.rda")
sim_tb <- tibble(sim)

df_pi0 <- sim_tb |>
  mutate(metrics = purrr::map(res, compute_metrics, 
                              .progress = "Calculating metrics...")) |>
  dplyr::select(-res) |> 
  unnest(metrics)

save(df_pi0, file = "paper/script/output/df_pi0.rda")




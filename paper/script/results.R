rm(list=ls())

# load pkgs
pkg <- c("dplyr","tidyr","purrr")
invisible(sapply(pkg, require, character.only = T))

# custom functions
source("paper/script/utl.R")

# Results simulation

load("paper/script/output/sim.rda")
sim_tb <- tibble(sim)

# prior-odds rule: every specified pi0 on the same replications
df_pi0 <- sim_tb |>
  tidyr::expand_grid(pi0 = seq(0.5, 0.9, by = 0.1)) |>
  mutate(cstar = sapply(pi0, cstar),
         metrics = purrr::map2(res, pi0, ~ compute_metrics(.x, "prior-odds", pi0 = .y),
                               .progress = "Prior-odds rule...")) |>
  dplyr::select(-res) |>
  unnest(metrics)

save(df_pi0, file = "paper/script/output/df_pi0.rda")

# prior-odds rule with pi0 set through a fixed Pr(H0) = .5, so that pi0 = .5^(1/m) grows with m
df_h0 <- sim_tb |>
  mutate(pi0 = 0.5^(1 / m),
         cstar = sapply(pi0, cstar),
         metrics = purrr::map2(res, pi0, ~ compute_metrics(.x, "prior-odds", pi0 = .y),
                               .progress = "Fixed Pr(H0)...")) |>
  dplyr::select(-res) |>
  unnest(metrics)

save(df_h0, file = "paper/script/output/df_h0.rda")

# comparison: the joint statement and the simultaneous band on the same
# replications as the prior-odds rule, for every n, m and pi0_true
df_compare <- sim_tb |>
  tidyr::expand_grid(method = c("joint", "simultaneous")) |>
  mutate(metrics = purrr::map2(res, method, ~ compute_metrics(.x, .y),
                               .progress = "Joint statement and band...")) |>
  dplyr::select(-res) |>
  unnest(metrics) |>
  bind_rows(df_pi0 |> mutate(method = "prior-odds"))

save(df_compare, file = "paper/script/output/df_compare.rda")

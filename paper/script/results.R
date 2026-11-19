# =====================================================================
# Error rates and power from the simulation
# Reads paper/script/output/sim.rda and writes df_pi0.rda, df_h0.rda and
# df_compare.rda to paper/script/output/; every data frame has the columns
# structure and rho (correlation of the observations)
# =====================================================================

rm(list = ls())

# 1. Packages ----------------------------------------------------------
pkg <- c("dplyr", "tidyr", "purrr")
invisible(lapply(pkg, library, character.only = TRUE))

source("paper/script/utl.R")

# 2. Load --------------------------------------------------------------
load("paper/script/output/sim.rda")
sim_tb <- tibble(sim)

# check
table(sim_tb$structure, sim_tb$n)

# 3. Prior-odds rule ---------------------------------------------------
# every specified pi0 on the same replications
df_pi0 <- sim_tb |>
  tidyr::expand_grid(pi0 = seq(0.5, 0.9, by = 0.1)) |>
  mutate(cstar = sapply(pi0, cstar),
         metrics = purrr::map2(res, pi0, ~ compute_metrics(.x, "prior-odds", pi0 = .y),
                               .progress = "Prior-odds rule...")) |>
  dplyr::select(-res) |>
  unnest(metrics)

save(df_pi0, file = "paper/script/output/df_pi0.rda")

# pi0 set through a fixed Pr(H0) = .5, so that pi0 = .5^(1/m) grows with m
df_h0 <- sim_tb |>
  mutate(pi0 = 0.5^(1 / m),
         cstar = sapply(pi0, cstar),
         metrics = purrr::map2(res, pi0, ~ compute_metrics(.x, "prior-odds", pi0 = .y),
                               .progress = "Fixed Pr(H0)...")) |>
  dplyr::select(-res) |>
  unnest(metrics)

save(df_h0, file = "paper/script/output/df_h0.rda")

# 4. Comparison --------------------------------------------------------
# the joint statement and the simultaneous band on the same replications as
# the prior-odds rule, only where they were computed (n = 50)
df_compare <- sim_tb |>
  filter(compare) |>
  tidyr::expand_grid(method = c("joint", "simultaneous")) |>
  mutate(metrics = purrr::map2(res, method, ~ compute_metrics(.x, .y),
                               .progress = "Joint statement and band...")) |>
  dplyr::select(-res) |>
  unnest(metrics) |>
  bind_rows(df_pi0 |> filter(compare) |> mutate(method = "prior-odds"))

# check
table(df_compare$method, df_compare$structure)

save(df_compare, file = "paper/script/output/df_compare.rda")

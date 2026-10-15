set.seed(150595)

rm(list=ls())

# load pkgs
pkg <- c("dplyr","MASS","tidyr","furrr","tibble","purrr")

invisible(sapply(pkg, require, character.only = T))

# custom functions
source("paper/script/utl.R")

# Simulation -----
m <- c(4, 8, 10, 20)
n <- c(30, 50, 100)
nsim <- 1e4
r <- 0
s <- 2
eff <- 0.3
pi0_true <- seq(0.5, 1, by = 0.1)

sim <- expand.grid(
  n = n,
  m = m,
  s = s,
  eff = eff,
  r = r,
  pi0_true = pi0_true
)


plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(
  sim[, c("n", "m", "s", "eff", "r", "pi0_true")],
  ~ sim_ind(n = ..1, m = ..2, s = ..3, eff = ..4,
            r = ..5, pi0_true = ..6, nsim = nsim,
            compare = TRUE),   # also store the claims of the joint statement and of the band
  .options  = furrr_options(seed = TRUE, packages = "multibayes"),
  .progress = TRUE
)
plan(sequential)

save(sim, file = "paper/script/output/sim.rda")

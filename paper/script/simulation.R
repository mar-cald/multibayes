# =====================================================================
# Simulation of the prior-odds rule and of the two comparison procedures
# Writes paper/script/output/sim.rda (raw pd, favored signs, true effects and,
# for n = 50, the claims of the joint statement and of the simultaneous band)
# Normal model with known covariance, diffuse prior N(0, 2^2), effects 0 or
# 0.3; two correlation structures of the observations: independent and .8
# between every pair
# =====================================================================

rm(list = ls())
set.seed(150595)

# 1. Packages ----------------------------------------------------------
pkg <- c("dplyr", "MASS", "tidyr", "furrr", "tibble", "purrr")
invisible(lapply(pkg, library, character.only = TRUE))

source("paper/script/utl.R")

# 2. Design ------------------------------------------------------------
m <- c(4, 8, 10, 20)
n <- c(30, 50, 100)
nsim <- 1e4
s <- 2                      # prior sd
eff <- 0.3
pi0_true <- seq(0.5, 1, by = 0.1)

# correlation structures of the observations
corr <- data.frame(structure = c("independent", "equi"),
                   rho       = c(0, 0.8))

sim <- expand.grid(n = n, m = m, s = s, eff = eff, pi0_true = pi0_true,
                   structure = corr$structure, stringsAsFactors = FALSE) |>
  left_join(corr, by = "structure") |>
  mutate(compare = n == 50)   # the comparison procedures only at n = 50 (slow)

# check
table(sim$structure, sim$compare)

# 3. Run ---------------------------------------------------------------
plan(multisession(workers = parallel::detectCores() - 2))
sim$res <- future_pmap(
  sim[, c("n", "m", "s", "eff", "structure", "rho", "pi0_true", "compare")],
  \(n, m, s, eff, structure, rho, pi0_true, compare)
    sim_ind(n = n, m = m, s = s, eff = eff, structure = structure, rho = rho,
            pi0_true = pi0_true, nsim = nsim, compare = compare),
  .options  = furrr_options(seed = TRUE, packages = "multibayes"),
  .progress = TRUE
)
plan(sequential)

# 4. Save --------------------------------------------------------------
save(sim, file = "paper/script/output/sim.rda")

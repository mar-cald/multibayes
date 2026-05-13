set.seed(150595)

rm(list=ls())

# load pkgs
pkg <- c("dplyr","MASS","tidyr","shape",
         "furrr","tibble","purrr", "multibayes")

invisible(sapply(pkg, require, character.only = T))

# custom functions
source("paper/script/utl.R")

# Simulation 1-----
m <- c(4, 6, 8, 10, 20)
n <- c(30, 50, 100, 500)
nsim <- 1e4
r <- 0
s <- 2
eff <- 0.3
p <- seq(0.05, 0.5, by = 0.1)

sim_ind <- function(n, m, r, s, eff, nsim, p) {
  replicate(nsim, {
    # cov
    S  = r + diag(1 - r, m)
    # --- null/non-null assignment ---
    n0 <- rbinom(1, size = m, prob = p^(1/m))
    is_null <- logical(m)
    is_null[sample(m, n0)] <- TRUE
    effsim <- rep(eff, m)
    effsim[is_null] <- 0
    
    # simulate data
    X <- MASS::mvrnorm(n, effsim, Sigma = S)
    
    # bayes
    pd = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    
    pd_adj = pd.adjust(pd = pd, p = p)$pd.adj
    
    list(pd = pd, pd_adj = pd_adj, effsim = effsim)
  }, simplify = FALSE)
} 

sim = expand.grid(
  n = n,
  m = m,
  s = s,
  eff = eff,
  r = r,
  p = p
)

plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(
  sim[, c("n", "m", "s", "eff", "r", "p")],
  ~ sim_ind(n = ..1, m = ..2, s = ..3, eff = ..4,
            r = ..5, p = ..6, nsim = nsim),
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
plan(sequential)

save(sim, file = "paper/script/output/sim_pi0_opt.rda")



# Simulation 2-----
m <- c(4, 8, 10, 20)
n <- c(30, 50, 100)
nsim <- 1e4
r <- 0
s <- 2
eff <- 0.3
pi0 <- seq(0.5, 0.9, by = 0.1) # implied per-test H0 probability
pi0_true <- seq(0.5, 1, by = 0.1)

sim_ind <- function(n, m, r, s, eff, nsim, p, pi0_true) {
  replicate(nsim, {
    # cov
    S  = r + diag(1 - r, m)
    # --- null/non-null assignment ---
    n0 <- rbinom(1, size = m, prob = pi0_true)
    is_null <- logical(m)
    is_null[sample(m, n0)] <- TRUE
    effsim <- rep(eff, m)
    effsim[is_null] <- 0
    
    # simulate data
    X <- MASS::mvrnorm(n, effsim, Sigma = S)
    
    # bayes
    pd = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    
    pd_adj = pd.adjust(pd = pd, p = p)$pd.adj
    
    list(pd = pd, pd_adj = pd_adj, effsim = effsim)
  }, simplify = FALSE)
} 


sim <- expand.grid(
  n        = n,
  m        = m,
  s        = s,
  eff      = eff,
  r        = r,
  pi0      = pi0,
  pi0_true = pi0_true
) |>
  dplyr::mutate(p = pi0^m)# obtain p


plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(
  sim[, c("n", "m", "s", "eff", "r", "p", "pi0_true")],
  ~ sim_ind(n = ..1, m = ..2, s = ..3, eff = ..4,
            r = ..5, p = ..6, pi0_true = ..7, nsim = nsim),
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
plan(sequential)

save(sim, file = "paper/script/output/sim_pi0.rda")



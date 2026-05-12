set.seed(150595)

rm(list=ls())

# load pkgs
pkg <- c("dplyr","MASS","tidyr","shape",
         "furrr","tibble","purrr", "multibayes")

invisible(sapply(pkg, require, character.only = T))

# custom functions
source("paper/script/utl.R")

# Simulation 1 ---------------------------------
m <- c(1,2,3,4,5,10,20) # number of tests
n <- c(30,50,100) # number of subjects
nsim <- 1e4 # number of simulations
r <- 0 # correlation
s <- c(.2, .5, 1, 1.5, 2) # prior sd
eff <- c(0,0.3)
q <- c(0.5, 0.4, 0.3,0.2, 0.1)  # prior all null

sim_ind <- function(n, m, r, s, eff, nsim, q){
  replicate(nsim, {
    
    r = r
    S = r + diag(1 - r, m) 
    R = S
    
    effsim = rep(0, m)
    
    if(eff != 0) {
      effsim[m] = eff
    }
    
    X = MASS::mvrnorm(n, effsim, Sigma = R) 
    
    # frequentist
    pval = apply(X, 2, function(x) z_test(x))
    #  adjustment
    pval_bonf = stats::p.adjust(pval, method = "bonferroni")
    pval_fdr = stats::p.adjust(pval, method = "fdr")
    
    # bayes
    pd = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    
    # adjustment
    pd_adj = pd.adjust(pd = pd, q = q)$pd.adj
    
    # save data
    data.frame(pval, pval_bonf,pval_holm,pval_fdr, pd, pd_adj)
  }, simplify = FALSE)
}


sim = expand.grid(
  n = n,
  m = m,
  s = s,
  eff = eff,
  r = r,
  q = q
)

plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(sim, 
                      ~sim_ind(n = ..1, m = ..2, s = ..3,eff = ..4,
                               r = ..5,  q = ..6,nsim = nsim), 
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)
plan(sequential)
save(sim, file = "paper/script/output/sim1.rda")


# Simulation 2-----
m <- c(4, 8, 10, 20)
n <- c(30, 50, 100)
nsim <- 1e4
r <- 0
s <- 2
eff <- 0.3
pi0 <- seq(0.5, 0.9, by = 0.1) # implied per-test H0 probability
pi0_true <- seq(0.5, 1, by = 0.1)

sim_ind <- function(n, m, r, s, eff, nsim, q, pi0_true) {
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
    
    pd_adj = pd.adjust(pd = pd, q = q)$pd.adj
    
    list(pd = pd, pd_adj = pd_adj, effsim = effsim)
  }, simplify = FALSE)
} 


sim <- expand.grid(
  n        = n,
  m        = m,
  s        = s,
  eff      = eff,
  r        = r,
  pi0       = pi0,
  pi0_true = pi0_true
) |>
  dplyr::mutate(q = pi0^m)   # q derivato, non direttamente variato


plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(
  sim[, c("n", "m", "s", "eff", "r", "q", "pi0_true")],
  ~ sim_ind(n = ..1, m = ..2, s = ..3, eff = ..4,
            r = ..5, q = ..6, pi0_true = ..7, nsim = nsim),
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
plan(sequential)

save(sim, file = "paper/script/output/sim_pi0q.rda")


# Relationship between Pr(H0) and true nulls----------
m <- c(4, 6, 8, 10, 20)
n <- c(30, 50, 100, 500)
nsim <- 1e4
r <- 0
s <- 2
eff <- 0.3
q <- seq(0.05, 0.5, by = 0.1)

sim_ind <- function(n, m, r, s, eff, nsim, q) {
  replicate(nsim, {
    # cov
    S  = r + diag(1 - r, m)
    # --- null/non-null assignment ---
    n0 <- rbinom(1, size = m, prob = q^(1/m))
    is_null <- logical(m)
    is_null[sample(m, n0)] <- TRUE
    effsim <- rep(eff, m)
    effsim[is_null] <- 0
    
    # simulate data
    X <- MASS::mvrnorm(n, effsim, Sigma = S)
    
    # bayes
    pd = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    
    pd_adj = pd.adjust(pd = pd, q = q)$pd.adj
    
    list(pd = pd, pd_adj = pd_adj, effsim = effsim)
  }, simplify = FALSE)
} 

sim = expand.grid(
  n = n,
  m = m,
  s = s,
  eff = eff,
  r = r,
  q = q
)

plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(
  sim[, c("n", "m", "s", "eff", "r", "q")],
  ~ sim_ind(n = ..1, m = ..2, s = ..3, eff = ..4,
            r = ..5, q = ..6, nsim = nsim),
  .options  = furrr_options(seed = TRUE),
  .progress = TRUE
)
plan(sequential)

save(sim, file = "paper/script/output/sim_pi0q_opt.rda")


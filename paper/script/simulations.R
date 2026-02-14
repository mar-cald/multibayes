# migliorabile

set.seed(150595)

rm(list=ls())

# load pkgs
pkg = c("dplyr","MASS","tidyr","shape","jointest",
        "furrr","tibble","purrr", "multibayes")

invisible(sapply(pkg, require, character.only = T))

# custom functions
source("paper/script/utl.R")


# Simulation 1 ---------------------------
m = c(1,2,5,10,20) # number of tests
n = c(30,50,100) # number of subjects
nsim = 1e4 # number of simulations
r = 0 # correlation
s = c(.2, .5, 1, 1.5, 2) # prior sd
eff = 0 # simulation under the null 
q = 0.5 # prior all null

sim_ind_1 = function(n, m, r, s, eff,nsim, q){
  replicate(nsim, {
    r =  r
    S =  r + diag(1 - r, m) 
    R =  S
    X =  MASS::mvrnorm(n, rep(eff, m), Sigma = R) 
    
    # frequentist
    pval = apply(X, 2, function(x) z_test(x))
    # adjustment
    pval_bonf= pval*m #bonf
    
    # bayes
    pd = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    # adjustment
    pd_adj = pd.adjust(pd = pd, q = q)
    
    # save data 
    data.frame(pval,pval_bonf,pd, pd_adj)}, 
    simplify = FALSE)
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
                      ~sim_ind_1(n = ..1, m = ..2, s = ..3,eff = ..4,
                                 r = ..5,  q = ..6, nsim = nsim), 
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)
plan(sequential)
save(sim, file = "paper/script/output/sim1.rda")


# Simuation 2 ---------------------------------
m = c(1,2,3,4,5,10,20,100) # number of tests
n = 50 # number of subjects
nsim = 1e4 # number of simulations
r = 0 # correlation
s = 2 # prior sd
eff = c(0,0.3)
q = c(0.5, 0.4, 0.3,0.2, 0.1)  # prior all null

sim_ind_2 = function(n, m, r, s, eff, nsim, q){
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
    # bonf adjustment
    pval_bonf = pval*m
    
    # bayes
    pd = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    # adjustment
    pd_adj = pd.adjust(pd = pd, q = q)
    
    # save data
    data.frame(pval, pval_bonf, pd, pd_adj)
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
                      ~sim_ind_2(n = ..1, m = ..2, s = ..3,eff = ..4,
                                 r = ..5,  q = ..6,nsim = nsim), 
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)
plan(sequential)
save(sim, file = "paper/script/output/sim2.rda")


# Simulation 3 ---------------------------------
m = c(2,4,6,10,20) # number of tests
r = c(0.2,0.4,0.6,0.8)

sim_corr_1 = function(n = 50, m, eff = 0.3,r = r, s = 2, nsim = 1e4, q = 0.4){
  replicate(nsim, {
    
    r = r
    S = r + diag(1 - r, m)
    R = S
    
    effsim = rep(0, m)
    effsim[m/2:m] = eff
    
    X = MASS::mvrnorm(n, effsim, Sigma = R)
    
    # bayes
    out = bayes_posterior_multivariate(X, tau0 = s)
    pd = out[[1]]
    post_cor = out[[2]]
    # adjustment corr
    pd_adj = pd.adjust(pd = pd, q = q, post_corr = post_cor)
    pd_meff = pd_adj > 0.975
    # adjustment no corr
    pd_adj = pd.adjust(pd = pd, q = q, post_corr = NULL)
    pd_m = pd_adj > 0.975
    
    # save data
    data.frame(pd_meff, pd_m,effsim)
  }, simplify = FALSE)
}

sim = expand.grid(
  m = m,
  r = r
)

plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(sim, 
                      ~sim_corr_1(m = ..1, r = ..2), 
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)
plan(sequential)
save(sim, file = "paper/script/output/sim3.rda")

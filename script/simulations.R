
rm(list=ls())

# load pkgs
pkg = c("dplyr","MASS","tidyr","shape",
        "dplyr","furrr","tibble","purrr")

invisible(sapply(pkg, require, character.only = T))
set.seed(150595)

# Functions------------------------

## Normal normal conjugate model
bayes_posterior_analytical = function(x, sigma = 1, tau0 = 1, mu0=0) {
  n = length(x)
  bar_y  = mean(x)
  post_mean  = (bar_y * (n / sigma^2) + mu0 * (1 / tau0^2)) / (n / sigma^2 + 1 / tau0^2)
  post_sd  = sqrt(1 / (n / sigma^2 + 1 / tau0^2))
  posterior_samples   =  rnorm(10000, mean = post_mean, sd = post_sd)
  return(posterior_samples)
}


## Prior adjustment
prior_adj = function(pd, m, q) {
  H0_i = q^(1 / m)
  H1_i = 1 - H0_i
  
  if (H0_i >= 0.5) {
    pd_adj = (pd * H1_i) / (pd * H1_i + (1 - pd) * H0_i)
  } else {
    warning("Pr(H0) < 0.5, pd will not be corrected")
    pd_adj = pd
  }
  
  return(pd_adj)
}
  

# Simulation 1 ---------------------------
m = c(1,2,5,10,15) # number of tests
n = c(30,50,100) # number of subjects
nsim = 10000 # number of simulations
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
    pval = apply(X, 2, function(x) t.test(x)$p.value)
    # sidak adjustment
    pval_sidak = 1-(1-pval)^m #sidak
    
    # bayes
    post = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    pd = apply(post, 2, function(x) max(mean(x > 0), mean(x < 0))) #two one sided tests
    # Jeffreys adjustment
    pd_adj = prior_adj(pd = pd, m = m, q = q)
    
    # save data 
    data.frame(pval,pval_sidak,pd, pd_adj)}, 
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

save(sim, file = "script/output/sim1.rda")


# Simuation 2 ---------------------------------
m = c(1,2,5,10,15) # number of tests
n = 100 # number of subjects
nsim = 1e4 # number of simulations
r = 0 # correlation
s = 2 # prior sd
eff = 0.25 # theta under H1 
q = c(0.3,0.1)  # prior all null
q_sim = 0.3 # true all null

sim_ind_2 = function(n, m, r, s, eff, nsim, q, q_sim){
  replicate(nsim, {
    
    # chance all tests in family are null
    all_null = runif(1) < q_sim
    
    if(all_null){
      # All null
      effect = rep(0, m)
    } else {
      # Otherwise: one effect, rest null
      which_effect = sample(1:m, 1)  # Randomly pick 1 test
      effect = rep(0, m)             # Start with all null
      effect[which_effect] = eff     # Set exactly one to eff
    }
    
    r = r
    S = r + diag(1 - r, m) 
    R = S
    X = MASS::mvrnorm(n, effect, Sigma = R) 
    
    # frequentist
    pval = apply(X, 2, function(x) t.test(x)$p.value)
    # sidak adjustment
    pval_sidak = 1-(1-pval)^m #sidak
    
    # bayes
    post = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    pd = apply(post, 2, function(x) max(mean(x > 0), mean(x < 0)))
    # Jeffreys adjustment
    pd_adj = prior_adj(pd = pd, m = m, q = q)
    
    # save data 
    data.frame(pval, pval_sidak, pd, pd_adj, is_null = all_null)
  }, simplify = FALSE)
}

sim = expand.grid(
  n = n,
  m = m,
  s = s,
  eff = eff,
  r = r,
  q = q,
  q_sim = q_sim
)

plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(sim, 
                      ~sim_ind_2(n = ..1, m = ..2, s = ..3,eff = ..4,
                               r = ..5,  q = ..6,q_sim = ..7, nsim = nsim), 
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)
plan(sequential)
save(sim, file = "script/output/sim2.rda")



rm(list=ls())

# load pkgs
pkg = c("dplyr","MASS","tidyr","shape",
        "furrr","tibble","purrr")

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

# z.test
z_test=function(x, mu0 = 0, sigma = 1){
  z = (mean(x) - mu0) / (sigma / sqrt(length(x)))
  pval = 2 * (1 - pnorm(abs(z)))
  return(pval)
} 


## Prior adjustment
#' Prior-odds adjustment for Probability of Direction (pd)
#'
#' Adjusts a vector of Probability of Direction (*pd*) values using a global
#' prior probability that all tested hypotheses are null, \eqn{q}. The adjustment
#' converts \eqn{q} into a per-hypothesis prior probability \eqn{H0 = q^{1/m}},
#' where \eqn{m} is the family size, and then reweights each *pd* accordingly.
#'
#' @param pd Numeric vector of *pd* values (typically \eqn{\in [0.5, 1]}).
#'   Each entry should represent the posterior mass on the favored direction,
#'   e.g., \eqn{\max\{\Pr(\theta>0\mid y), \Pr(\theta<0\mid y)\}}.
#' @param q Numeric scalar in (0, 1). Interpreted as the prior probability that
#'   **all** hypotheses in the family are null.
#'
#' @details
#' Let \eqn{m = \mathrm{length}(pd)}. The function computes:
#' \deqn{H0 = q^{1/m}, \quad H1 = 1 - H0.}
#' If \eqn{H0 \ge 0.5}, each element of `pd` is adjusted as:
#' \deqn{pd_{\mathrm{adj}} = \frac{pd \cdot H1}{pd \cdot H1 + (1-pd)\cdot H0}.}
#' If \eqn{H0 < 0.5}, the function returns the unadjusted `pd` values and emits a
#' warning (because the “adjustment” would no longer be conservative under that
#' setting).
#'
#' @return A numeric vector of the same length as `pd`, containing adjusted *pd*
#'   values (or the original `pd` if \eqn{H0 < 0.5}).
#'
#' @examples
#' pd = c(0.55, 0.80, 0.97)
#' prior_adj(pd, q = 0.5)
#'
#' # Larger family -> stronger per-hypothesis skepticism (larger H0):
#' prior_adj(rep(0.97, 20), q = 0.5)
#'

prior_adj = function(pd, q = 0.5, m = length(pd)) {

  stopifnot(
    "`pd` must be numeric"        = is.numeric(pd),
    "`pd` must be in [0.5, 1]"    = all(is.finite(pd)) && all(pd >= 0.5 & pd <= 1),
    "`q` must be a single number" = length(q) == 1L && is.finite(q),
    "`q` must be in (0, 1)"       = (q > 0) && (q < 1),
    "`m` must be >= 1"            = length(m) == 1L && is.finite(m) && (m >= 1)
  )
  
  H0 = q^(1 / m)
  H1 = 1 - H0
  
  if (H0 >= 0.5) {
    pd_adj = (pd * H1) / (pd * H1 + (1 - pd) * H0)
    
  } else {
    warning("Pr(H0) < 0.5; returning unadjusted pd")
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
    pval = apply(X, 2, function(x) z_test(x))
    # sidak adjustment
    pval_sidak = 1-(1-pval)^m #sidak
    
    # bayes
    post = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    pd = apply(post, 2, function(x) max(mean(x > 0), mean(x < 0))) #two one sided tests
    # Jeffreys adjustment
    pd_adj = prior_adj(pd = pd, q = q)
    
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
plan(sequential)
save(sim, file = "script/output/sim1.rda")


# Simuation 2 ---------------------------------
m = c(1,2,3,4,5,10,15) # number of tests
n = 100 # number of subjects
nsim = 1e4 # number of simulations
r = 0 # correlation
s = 2 # prior sd
eff = 0.25 # theta under H1 
q = c(0.5, 0.4, 0.3,0.2, 0.1)  # prior all null
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
    pval = apply(X, 2, function(x) z_test(x))
    # sidak adjustment
    pval_sidak = 1-(1-pval)^m #sidak
    
    # bayes
    post = apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
    pd = apply(post, 2, function(x) max(mean(x > 0), mean(x < 0)))
    # Jeffreys adjustment
    pd_adj = prior_adj(pd = pd, q = q)
    
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


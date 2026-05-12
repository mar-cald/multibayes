# Simulation 2 ---------------------------------
m = c(2,4,6,10,20) # number of tests
r = c(0.2,0.4,0.6,0.8)
n = c(30,50,100)

sim_corr = function(n = n, m, eff = 0.3, r = r, s = 1.5, nsim = 1e4, q = 0.4){
  replicate(nsim, {
    
    r = r
    S = r + diag(1 - r, m)
    R = S
    
    effsim = rep(0, m)
    effsim[m/2:m] = eff
    
    X = MASS::mvrnorm(n, effsim, Sigma = R)
    
    # frequentist
    pval = apply(X, 2, function(x) z_test(x))
    #  adjustment
    pval_holm = stats::p.adjust(pval, method = "holm")
    pval_holm = pval_holm < 0.05
    
    # bayes
    out = bayes_posterior_multivariate(X, tau0 = s)
    pd = out[[1]]
    post_cor = out[[2]]
    
    # adjustment corr
    pd_adj = pd.adjust(pd = pd, q = q, R = post_cor)
    pd_meff = pd_adj[,2] > 0.975
    # save data
    data.frame(pval_holm, pd_meff = pd_meff, 
               effsim)
  }, simplify = FALSE)
}

sim = expand.grid(
  m = m,
  r = r,
  n = n
)

plan(multisession(workers = parallel::detectCores() - 2))
sim$res = future_pmap(sim, 
                      ~sim_corr(m = ..1, r = ..2 , n = ..3), 
                      .options = furrr_options(seed = TRUE),
                      .progress = TRUE)
plan(sequential)
save(sim, file = "paper/script/output/sim2.rda")


# Results simulation 2--------------------------
load(file = "paper/script/output/sim2.rda")


sim2  <-  tibble(sim)

# Compute FWER
sign_err  <-  function(data){
  data |> 
    group_by(m,r,n,effsim,nsim) |> 
    reframe(sign = any(value))  
}
# Compute power (avg)
sign_pw  <-  function(data){
  data |> 
    filter(effsim != 0) |> 
    group_by(m, r,n, nsim) |> 
    reframe(sign = mean(value)) 
}


sim_res2  <-  sim2 |>
  unnest(res) |>
  group_by(m,r,n) |>
  mutate(nsim = 1:n()) |>
  ungroup() |>
  unnest(res) |>
  pivot_longer(c("pd_meff","pval_holm")) |>
  group_by(name) |>
  nest()

sim_res2$fwer = vector(mode = "list", length = nrow(sim_res2))
sim_res2$pw = vector(mode = "list", length = nrow(sim_res2))

for(i in 1:nrow(sim_res2)){
  sim_res2$fwer[[i]] = sign_err(sim_res2$data[[i]])
  sim_res2$pw[[i]] = sign_pw(sim_res2$data[[i]])
}

save(sim_res2, file = "paper/script/output/sim_res2.rda")




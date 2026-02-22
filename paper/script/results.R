library(dplyr)
library(tidyr)

# Results simulation 1

# 1. Load Data
load(file = "paper/script/output/sim1.rda")

sim1 <- tibble(sim)
n_grids <- nrow(sim1)

# Pre-allocate empty lists to hold our two new columns
res_value_list <- vector("list", n_grids)
res_sign_list <- vector("list", n_grids)

# 2. Loop through and compress the lists
for(i in seq_len(n_grids)) {
  
  if (i %% 50 == 0) cat("Processed", i, "of", n_grids, "grids...\r")
  
  # Bind the 10,000 individual data frames into ONE clean data frame
  bound_sims <- bind_rows(sim1$res[[i]], .id = "nsim")
  
  # A. Save the raw continuous values
  res_value_list[[i]] <- bound_sims
  
  # B. Compute and save ONLY the TRUE/FALSE significance
  res_sign_list[[i]] <- bound_sims |> 
    transmute(
      nsim      = nsim,
      pval      = pval < 0.05,
      pval_bonf = pval_bonf < 0.05,
      pval_holm = pval_holm < 0.05,
      pval_fdr  = pval_fdr < 0.05,
      pd        = pd > 0.975,
      pd_adj    = pd_adj > 0.975
    )
  
  # Empty the trash to guarantee RAM stays low
  if(i %% 50 == 0) gc() 
}

# 3. Attach our two new lists to the main dataframe
sim_res1 <- sim1 |> 
  select(eff, m, s, q, n) |> 
  mutate(
    res_value = res_value_list,
    res_sign  = res_sign_list
  )

# Clean up 
rm(sim1, res_value_list, res_sign_list, bound_sims)

# Save
save(sim_res1, file = "paper/script/output/sim_res1.rda")


# Compute quantites of interest----------------
## FWER
df_fwer <- sim_res1 |> 
  # 1. Filter for the null effect
  filter(eff == 0) |> 
  # 2. Keep the parameters we need and unnest the TRUE/FALSE significance data
  dplyr::select(n, s, m, q, res_sign) |> 
  unnest(res_sign) |> 
  
  # 3. First, collapse the 'm' tests into 1 logical result per simulation
  group_by(n, s, m, q, nsim) |> 
  summarise(
    pval_sim = any(pval),
    pval_holm_sim = any(pval_holm),
    pval_fdr_sim = any(pval_fdr),
    pval_bonf_sim= any(pval_bonf),
    pd_adj_sim    = any(pd_adj),
    pd_sim    = any(pd),
    .groups       = "drop"
  ) |> 
  # 4. Second, calculate the FWER (mean) across all 10,000 simulations
  group_by(n,m, s, q) |> 
  summarise(
    pval_bonf = mean(pval_bonf_sim),
    pval_holm    = mean(pval_holm_sim),
    pval_fdr    = mean(pval_fdr_sim),
    pval = mean(pval_sim),
    pd    = mean(pd_sim),
    pd_adj    = mean(pd_adj_sim),
    .groups   = "drop"
  ) |> 
  
  # 5. Pivot into the long format 
  pivot_longer(
    cols      = c(pd_adj,pd, pval_bonf,pval_holm,pval_fdr,pval ),
    names_to  = "name",
    values_to = "FWER"
  )

# save
save(df_fwer, file = "paper/script/output/df_fwer_sim1.rda")


## Power
df_pw <- sim_res1 |> 
  # 1. Filter for the null effect
  filter(eff == 0.3) |> 
  # 2. Keep the parameters we need and unnest the TRUE/FALSE significance data
  dplyr::select(n, s, m, q, res_sign) |> 
  unnest(res_sign) |> 
  
  # 3. First, collapse the 'm' tests into 1 logical result per simulation
  group_by(n, s, m, q, nsim) |> 
  summarise(
    pval_sim = any(pval),
    pval_holm_sim = any(pval_holm),
    pval_fdr_sim = any(pval_fdr),
    pval_bonf_sim= any(pval_bonf),
    pd_adj_sim    = any(pd_adj),
    pd_sim    = any(pd),
    .groups       = "drop"
  ) |> 
  # 4. Second, calculate the power (mean) across all 10,000 simulations
  group_by(n,m, s, q) |> 
  summarise(
    pval_bonf = mean(pval_bonf_sim),
    pval_holm    = mean(pval_holm_sim),
    pval_fdr    = mean(pval_fdr_sim),
    pval = mean(pval_sim),
    pd    = mean(pd_sim),
    pd_adj    = mean(pd_adj_sim),
    .groups   = "drop"
  ) |> 
  
  # 5. Pivot into the long format
  pivot_longer(
    cols      = c(pd_adj,pd, pval_bonf,pval_holm,pval_fdr,pval ),
    names_to  = "name",
    values_to = "pw"
  )

# save
save(df_pw, file = "paper/script/output/df_pw_sim1.rda")


# Results simulation 2--------------------------
load(file = "script/output/sim2.rda")

sim2  <-  tibble(sim)

# Compute FWER
sign_err  <-  function(data){
  data |> 
    group_by(m,r,n,effsim,nsim) |> 
    reframe(sign = any(value))  
}
# Compute power
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
  pivot_longer(c("pd_meff","pd_m","pval_holm","pval_fdr")) |>
  group_by(name) |>
  nest()

sim_res2$fwer = vector(mode = "list", length = nrow(sim_res2))
sim_res2$pw = vector(mode = "list", length = nrow(sim_res2))

for(i in 1:nrow(sim_res2)){
  sim_res2$fwer[[i]] = sign_err(sim_res2$data[[i]])
  sim_res2$pw[[i]] = sign_pw(sim_res2$data[[i]])
}

save(sim_res2, file = "paper/script/output/sim_res2.rda")




# Results simulation 3--------------------------
load(file = "script/output/sim3.rda")


sim3  <-  tibble(sim)

# Compute FWER
sign_err  <-  function(data){
  data |> 
    group_by(m,r,n,effsim,nsim) |> 
    reframe(sign = any(value))  
}
# Compute power
sign_pw  <-  function(data){
  data |> 
    filter(effsim != 0) |> 
    group_by(m, r,n, nsim) |> 
    reframe(sign = mean(value)) 
}


sim_res3  <-  sim3 |>
  unnest(res) |>
  group_by(m,r,n) |>
  mutate(nsim = 1:n()) |>
  ungroup() |>
  unnest(res) |>
  pivot_longer(c("pd_meff","pd_sim","pval_holm","pval_fdr")) |>
  group_by(name) |>
  nest()

sim_res3$fwer = vector(mode = "list", length = nrow(sim_res3))
sim_res3$pw = vector(mode = "list", length = nrow(sim_res3))

for(i in 1:nrow(sim_res3)){
  sim_res3$fwer[[i]] = sign_err(sim_res3$data[[i]])
  sim_res3$pw[[i]] = sign_pw(sim_res3$data[[i]])
}

save(sim_res3, file = "paper/script/output/sim_res3.rda")


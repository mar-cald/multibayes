library(dplyr)
library(tidyr)
library(purrr)


# Results simulation 1------

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
      nsim = nsim,
      pval = pval < 0.05,
      pval_bonf = pval_bonf < 0.05,
      pval_holm = pval_holm < 0.05,
      pval_fdr  = pval_fdr < 0.05,
      pd = pd > 0.975,
      pd_adj = pd_adj > 0.975
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
    pd_adj_sim  = any(pd_adj),
    pd_sim  = any(pd),
    .groups = "drop"
  ) |> 
  # 4. Second, calculate the FWER (mean) across all 10,000 simulations
  group_by(n,m, s, q) |> 
  summarise(
    pval_bonf = mean(pval_bonf_sim),
    pval_holm  = mean(pval_holm_sim),
    pval_fdr= mean(pval_fdr_sim),
    pval = mean(pval_sim),
    pd = mean(pd_sim),
    pd_adj = mean(pd_adj_sim),
    .groups = "drop"
  ) |> 
  
  # 5. Pivot into the long format 
  pivot_longer(
    cols      = c(pd_adj,pd, pval_bonf,pval_holm,pval_fdr,pval ),
    names_to  = "name",
    values_to = "FWER"
  )

# save
save(df_fwer, file = "paper/script/output/df_fwer_sim1.rda")

# Results simulation 2-----
rm(list=ls())

load("paper/script/output/sim_pi0q.rda")
sim_tb <- tibble(sim)

compute_metrics <- function(res_list) {
  # 1. TRASFORMAZIONE IN MATRICI
  # do.call(rbind, ...) estrae i vettori da tutte le iterazioni e 
  # crea matrici con 10.000 righe (iterazioni) e 'm' colonne (test)
  pd_adj_mat <- do.call(rbind, lapply(res_list, `[[`, "pd_adj"))
  pd_mat <- do.call(rbind, lapply(res_list, `[[`, "pd"))
  effsim_mat <- do.call(rbind, lapply(res_list, `[[`, "effsim"))
  
  # Matrice logica: TRUE se in quella posizione c'era un'ipotesi nulla
  is_null_mat <- (effsim_mat == 0)
  
  # ══════════════════════════════════════════════════════════════
  # 2. METRICHE GLOBALI (calcolate su tutte le celle in un colpo solo)
  # ══════════════════════════════════════════════════════════════
  
  # Il check della percentuale reale di H0 simulate
  pi0_check <- mean(is_null_mat) 
  
  gamma <- 0.60
  prob_H0_in_zone <- (gamma - 0.5) / 0.5
  
  # Stima di pi0 (calcoliamo la prop. di pd < gamma e applichiamo i limiti 0-1)
  pi0_est_adj <- min(max(mean(pd_adj_mat < gamma, na.rm = TRUE) / prob_H0_in_zone, 0), 1)
  pi0_est_pd  <- min(max(mean(pd_mat < gamma, na.rm = TRUE) / prob_H0_in_zone, 0), 1)
  
  prop_not_sig_adj <- mean(pd_adj_mat < 0.975, na.rm = TRUE)
  prop_not_sig_pd  <- mean(pd_mat < 0.975,     na.rm = TRUE)
  
  # ══════════════════════════════════════════════════════════════
  # 3. FWER (Almeno un falso positivo nell'iterazione)
  # ══════════════════════════════════════════════════════════════
  
  # Falso positivo: p-value significativo (> 0.975) DOVE is_null_mat è TRUE
  fp_adj_mat <- (pd_adj_mat > 0.975) & is_null_mat
  fp_pd_mat  <- (pd_mat > 0.975)     & is_null_mat
  
  # rowSums somma gli errori per riga (iterazione). Se > 0, c'è un errore FWER.
  fwer_adj <- mean(rowSums(fp_adj_mat, na.rm = TRUE) > 0)
  fwer_pd  <- mean(rowSums(fp_pd_mat,  na.rm = TRUE) > 0)
  
  # ══════════════════════════════════════════════════════════════
  # 4. POWER (Veri Positivi / Totale H1 vere)
  # ══════════════════════════════════════════════════════════════
  
  # Vero positivo: p-value significativo DOVE is_null_mat è FALSE
  tp_adj_mat <- (pd_adj_mat > 0.975) & !is_null_mat
  tp_pd_mat  <- (pd_mat > 0.975)     & !is_null_mat
  
  # Calcolo aggregato della potenza
  pw_adj <- sum(tp_adj_mat, na.rm = TRUE) / sum(!is_null_mat)
  pw_pd  <- sum(tp_pd_mat,  na.rm = TRUE) / sum(!is_null_mat)
  
  # Gestione di sicurezza: se per puro caso in 10.000 sim non c'era Manco una H1
  if (is.nan(pw_adj)) pw_adj <- NA
  if (is.nan(pw_pd))  pw_pd  <- NA
  
  
  # ══════════════════════════════════════════════════════════════
  # FDR (False Discovery Rate)
  # ══════════════════════════════════════════════════════════════
  
  # R: Totale dei rigetti per ogni iterazione (somma per riga)
  R_adj_mat <- rowSums(pd_adj_mat > 0.975, na.rm = TRUE)
  R_pd_mat  <- rowSums(pd_mat > 0.975, na.rm = TRUE)
  
  # V: Totale dei Falsi Positivi per ogni iterazione (li avevamo già calcolati per il FWER)
  V_adj_mat <- rowSums(fp_adj_mat, na.rm = TRUE)
  V_pd_mat  <- rowSums(fp_pd_mat,  na.rm = TRUE)
  
  # FDP (False Discovery Proportion): V / R. Se R = 0, FDP = 0.
  fdp_adj <- ifelse(R_adj_mat == 0, 0, V_adj_mat / R_adj_mat)
  fdp_pd  <- ifelse(R_pd_mat == 0,  0, V_pd_mat / R_pd_mat)
  
  # FDR è il valore atteso (la media) della FDP su tutte le iterazioni
  fdr_adj <- mean(fdp_adj, na.rm = TRUE)
  fdr_pd  <- mean(fdp_pd,  na.rm = TRUE)
  
  
  # ══════════════════════════════════════════════════════════════
  # 5. OUTPUT
  # ══════════════════════════════════════════════════════════════
  
  data.frame(
    method = c("pd_adj", "pd"),
    FWER = c(fwer_adj, fwer_pd),
    pw = c(pw_adj, pw_pd),
    FDR = c(fdr_adj, fdr_pd),   
    pi0_est = c(pi0_est_adj, pi0_est_pd),
    prop_not_sig = c(prop_not_sig_adj, prop_not_sig_pd),
    pi0_check = c(pi0_check, pi0_check)
  )
}


df_pi0q <- sim_tb |>
  mutate(metrics = purrr::map(res, compute_metrics, 
                              .progress = "Calculating metrics...")) |>
  dplyr::select(-res) |> 
  unnest(metrics)

save(df_pi0q, file = "paper/script/output/df_pi0q.rda")


# Results simulation 3------
load("paper/script/output/sim_pi0q_opt.rda")
sim_tb <- tibble(sim)

df_pi0q_opt <- sim_tb |>
  mutate(metrics = purrr::map(res, compute_metrics, 
                              .progress = "Calculating metrics...")) |>
  dplyr::select(-res) |> 
  unnest(metrics)

save(df_pi0q_opt, file = "paper/script/output/df_pi0q_opt.rda")


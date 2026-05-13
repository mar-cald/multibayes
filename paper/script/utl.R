# Functions

## Normal model
bayes_posterior_analytical <- function(x, sigma = 1, tau0 = 1, mu0=0) {
  n = length(x)
  bar_y  = mean(x)
  post_mean  = (bar_y * (n / sigma^2) + mu0 * (1 / tau0^2)) / (n / sigma^2 + 1 / tau0^2)
  post_sd  = sqrt(1 / (n / sigma^2 + 1 / tau0^2))
  pd = pnorm(abs(post_mean - 0) / post_sd)
  return(pd)
}

# Results simulation
compute_metrics <- function(res_list) {
  # 1. TRASFORMAZIONE IN MATRICI
  # do.call(rbind, ...) estrae i vettori da tutte le iterazioni e 
  # crea matrici con 10.000 righe (iterazioni) e 'm' colonne (test)
  pd_adj_mat <- do.call(rbind, lapply(res_list, `[[`, "pd_adj"))
  pd_mat <- do.call(rbind, lapply(res_list, `[[`, "pd"))
  effsim_mat <- do.call(rbind, lapply(res_list, `[[`, "effsim"))
  
  # Matrice logica: TRUE se in quella posizione c'era un'ipotesi nulla
  is_null_mat <- (effsim_mat == 0)
  
  # Il check della percentuale reale di H0 simulate
  pi0_check <- mean(is_null_mat) 
  
  gamma <- 0.60
  prob_H0_in_zone <- (gamma - 0.5) / 0.5
  
  # Stima di pi0 (calcoliamo la prop. di pd < gamma e applichiamo i limiti 0-1)
  pi0_est_adj <- min(max(mean(pd_adj_mat < gamma, na.rm = TRUE) / prob_H0_in_zone, 0), 1)
  pi0_est_pd  <- min(max(mean(pd_mat < gamma, na.rm = TRUE) / prob_H0_in_zone, 0), 1)
  
  prop_not_sig_adj <- mean(pd_adj_mat < 0.975, na.rm = TRUE)
  prop_not_sig_pd  <- mean(pd_mat < 0.975,     na.rm = TRUE)
  
  # Falso positivo: > 0.975 DOVE is_null_mat è TRUE
  fp_adj_mat <- (pd_adj_mat > 0.975) & is_null_mat
  fp_pd_mat  <- (pd_mat > 0.975)     & is_null_mat
  
  # rowSums somma gli errori per riga (iterazione). Se > 0, c'è un errore FWER.
  fwer_adj <- mean(rowSums(fp_adj_mat, na.rm = TRUE) > 0)
  fwer_pd  <- mean(rowSums(fp_pd_mat,  na.rm = TRUE) > 0)
  
  # Vero positivo
  tp_adj_mat <- (pd_adj_mat > 0.975) & !is_null_mat
  tp_pd_mat  <- (pd_mat > 0.975)     & !is_null_mat
  
  # Calcolo aggregato della potenza
  pw_adj <- sum(tp_adj_mat, na.rm = TRUE) / sum(!is_null_mat)
  pw_pd  <- sum(tp_pd_mat,  na.rm = TRUE) / sum(!is_null_mat)
  
  # Gestione di sicurezza
  if (is.nan(pw_adj)) pw_adj <- NA
  if (is.nan(pw_pd))  pw_pd  <- NA
  
  
  # R: Totale dei rigetti per ogni iterazione (somma per riga)
  R_adj_mat <- rowSums(pd_adj_mat > 0.975, na.rm = TRUE)
  R_pd_mat  <- rowSums(pd_mat > 0.975, na.rm = TRUE)
  
  # V: Totale dei Falsi Positivi per ogni iterazione 
  V_adj_mat <- rowSums(fp_adj_mat, na.rm = TRUE)
  V_pd_mat  <- rowSums(fp_pd_mat,  na.rm = TRUE)
  
  # FDP (False Discovery Proportion): V / R. Se R = 0, FDP = 0.
  fdp_adj <- ifelse(R_adj_mat == 0, 0, V_adj_mat / R_adj_mat)
  fdp_pd  <- ifelse(R_pd_mat == 0,  0, V_pd_mat / R_pd_mat)
  
  # FDR è il valore atteso (la media) della FDP su tutte le iterazioni
  fdr_adj <- mean(fdp_adj, na.rm = TRUE)
  fdr_pd  <- mean(fdp_pd,  na.rm = TRUE)
  
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

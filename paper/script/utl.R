# Functions

## Normal model
bayes_posterior_analytical <- function(x, sigma = 1, tau0 = 1, mu0=0) {
  n <- length(x)
  bar_y <- mean(x)
  post_mean <- (bar_y * (n / sigma^2) + mu0 * (1 / tau0^2)) / (n / sigma^2 + 1 / tau0^2)
  post_sd <- sqrt(1 / (n / sigma^2 + 1 / tau0^2))
  pd <- pnorm(abs(post_mean - 0) / post_sd) 
  return(pd)
}

# Results simulation
compute_metrics <- function(res_list) {
  # Matrix transformation
  # from vectors to matrix with nrow = iter  and ncol = m (n test)
  pd_adj_mat <- do.call(rbind, lapply(res_list, `[[`, "pd_adj"))
  pd_mat <- do.call(rbind, lapply(res_list, `[[`, "pd"))
  effsim_mat <- do.call(rbind, lapply(res_list, `[[`, "effsim"))
  
  # Logi matrix
  is_null_mat <- (effsim_mat == 0)
  
  # Check simulation set up
  pi0_check <- mean(is_null_mat) 
  
  # Curiosity-----
  gamma <- 0.60
  prob_H0_in_zone <- (gamma - 0.5) / 0.5
  
  # Estimated null prob given pd in 0.6-0.5
  pi0_est_adj <- min(max(mean(pd_adj_mat < gamma, na.rm = TRUE) / prob_H0_in_zone, 0), 1)
  pi0_est_pd  <- min(max(mean(pd_mat < gamma, na.rm = TRUE) / prob_H0_in_zone, 0), 1)
  # check prop sign adj
  prop_not_sig_adj <- mean(pd_adj_mat < .975, na.rm = TRUE)
  prop_not_sig_pd  <- mean(pd_mat < .975,     na.rm = TRUE)
  
  # False Positive----
  fp_adj_mat <- (pd_adj_mat > .975) & is_null_mat
  fp_pd_mat <- (pd_mat > .975) & is_null_mat
  
  # FWER: sum errors
  fwer_adj <- mean(rowSums(fp_adj_mat, na.rm = TRUE) > 0)
  fwer_pd <- mean(rowSums(fp_pd_mat,  na.rm = TRUE) > 0)
  
  # True positive
  tp_adj_mat <- (pd_adj_mat > .975) & !is_null_mat
  tp_pd_mat <- (pd_mat > .975)     & !is_null_mat
  
  # Avg Power
  pw_adj <- sum(tp_adj_mat, na.rm = TRUE) / sum(!is_null_mat)
  pw_pd <- sum(tp_pd_mat,  na.rm = TRUE) / sum(!is_null_mat)
  
  # Control in case of null
  if (is.nan(pw_adj)) pw_adj <- NA
  if (is.nan(pw_pd))  pw_pd  <- NA
  
  # Sum n rejections for each iter
  R_adj_mat <- rowSums(pd_adj_mat > .975, na.rm = TRUE)
  R_pd_mat <- rowSums(pd_mat > .975, na.rm = TRUE)
  
  # V: total FA
  V_adj_mat <- rowSums(fp_adj_mat, na.rm = TRUE)
  V_pd_mat  <- rowSums(fp_pd_mat,  na.rm = TRUE)
  
  # FDP (False Discovery Proportion): V / R. For R = 0, FDP = 0.
  fdp_adj <- ifelse(R_adj_mat == 0, 0, V_adj_mat / R_adj_mat)
  fdp_pd<- ifelse(R_pd_mat == 0,  0, V_pd_mat / R_pd_mat)
  
  # FDR expected FDP
  fdr_adj <- mean(fdp_adj, na.rm = TRUE)
  fdr_pd<- mean(fdp_pd,  na.rm = TRUE)
  
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

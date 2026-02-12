# Parameters
m_values <- 1:60
alpha_base <- 0.05
target_pd <- 0.95 # We want to maintain this posterior probability

# Initialize data frame to store results
results <- data.frame(
  m = m_values,
  alpha_adj = NA,
  pd_adj1 = NA,
  calculated_q = NA,
  verified_pd = NA
)

for (i in 1:length(m_values)) {
  m <- m_values[i]
  
  # 1. Bonferroni Adjustment
  # The frequentist confidence increases because alpha shrinks
  alpha_adj <- alpha_base / m
  pd_val <- 1 - alpha_adj 
  
  # 2. Calculate q to EQUATE the decrease
  # We solve for q such that the Bayesian Posterior (pd_adj2) stays at 0.95
  # Formula derived from inverting the Bayesian update equation
  numerator <- pd_val * (1 - target_pd)
  denominator <- pd_val * (1 - target_pd) + target_pd * (1 - pd_val)
  H0_needed <- numerator / denominator
  q_calculated <- H0_needed^m
  
  # 3. Verify the result (Reverse check)
  H0_verify <- q_calculated^(1/m)
  H1_verify <- 1 - H0_verify
  # Bayesian Posterior Formula
  pd_adj2_verify <- (pd_val * H1_verify) / 
    (pd_val * H1_verify + (1 - pd_val) * H0_verify)
  
  # Store results
  results$alpha_adj[i] <- alpha_adj
  results$pd_adj1[i] <- pd_val
  results$calculated_q[i] <- q_calculated
  results$verified_pd[i] <- pd_adj2_verify
}

# Print the table
print(round(results, 5))

# Check the asymptotic limit (theoretical q)
cat("\nAsymptotic limit of q (exp(-0.95)):", exp(-0.95))

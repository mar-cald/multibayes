# Functions

## Normal model
bayes_posterior_analytical = function(x, sigma = 1, tau0 = 1, mu0=0) {
  n = length(x)
  bar_y  = mean(x)
  post_mean  = (bar_y * (n / sigma^2) + mu0 * (1 / tau0^2)) / (n / sigma^2 + 1 / tau0^2)
  post_sd  = sqrt(1 / (n / sigma^2 + 1 / tau0^2))
  pd = pnorm(abs(post_mean - 0) / post_sd)
  return(pd)
}

## Multivariate normal  model
bayes_posterior_multivariate = function(X, tau0 = 2) {
  n = nrow(X)
  m = ncol(X)
  
  # Estimate covariance matrix of data
  Sigma = cov(X)
  
  # Precision matrix, matrix inverse
  Sigma_inv = MASS::ginv(Sigma)
  
  # Prior precision matrix, corr =  0
  prior_cov_inv = diag(1 / tau0^2, m)
  
  # Posterior covariance
  post_cov = solve(n * Sigma_inv + prior_cov_inv)
  
  # Posterior mean vector
  post_mean = post_cov %*% (n * Sigma_inv %*% colMeans(X))
  
  # Post SD
  post_sd = sqrt(diag(post_cov))
  
  # extract pd and correlations
  pd = pnorm(abs(post_mean - 0) / post_sd)
  post_cor =  cov2cor(post_cov)
  
  return(list(pd,post_cor))
}

## Z test
z_test=function(x, mu0 = 0, sigma = 1){
  z = (mean(x) - mu0) / (sigma / sqrt(length(x)))
  pval = 2 * (1 - pnorm(abs(z)))
  return(pval)
} 


## MaxT correction Z test
maxZ = function(X){
  n = nrow(X)
  m = ncol(X)
  B = 10000
  
  # observed z
  obs_stats = apply(X, 2, function(x) (mean(x) - 0) / 
                      (1 / sqrt(length(x)))
  )
  
  # init max
  max_null = numeric(B)
  # z null
  for(b in 1:B) {
    signs = sample(c(-1, 1), n, replace = TRUE)
    X_perm = X * signs
    stats_b = apply(X_perm, 2, function(x) (mean(x) - 0) / 
                      (1 / sqrt(length(x))))
    max_null[b] = max(abs(stats_b))  # Store the maximum
  }
  
  # Critical value for two-sided test at alpha = 0.05
  crit_value = quantile(max_null, 0.95)
  
  # Reject if |obs_stat| > crit_value
  rejections = abs(obs_stats) > crit_value
  
  return(rejections)
}


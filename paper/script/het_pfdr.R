# With heterogeneous per-test error rates, pi0 stays a single family-level
# proportion (exchangeable labels), and the family-averaged mixture formula
# gives the ratio of expected counts E[V]/E[R],
#
#   E[V]/E[R] = pi0 * abar / ( pi0 * abar + (1 - pi0) * (1 - bbar) )
#   abar = mean(alpha_i),  bbar = mean(beta_i)
#
# which approximates the pFDR = E[V/R | R>0] (the two coincide for identical
# tests). The script reports pFDR (E[V/R|R>0]), E[V]/E[R], and the formula.

set.seed(2026)

# Closed-form per-test Type I (alpha) and Type II (beta) error rates for the
# normal-normal pd > c rule:  y_i | theta ~ N(theta, 1),  theta ~ N(0, tau^2).
ab_rates <- function(n, tau = 2, eff = 0.3, c = 0.975) {
  kappa <- n / (n + 1 / tau^2)
  t <- qnorm(c) / sqrt(kappa * n)
  alpha <- 2 * (1 - pnorm(t * sqrt(n)))
  power <- (1 - pnorm((t - eff) * sqrt(n))) + pnorm((-t - eff) * sqrt(n))
  c(alpha = alpha, beta = 1 - power)
}

# EXACT pFDR = E[V/R | R>0] for INDEPENDENT heterogeneous tests, via the joint
# law of (V = false positives, S = true positives). Each test is, independently,
# a false positive w.p. pi0*alpha_i, a true positive w.p. (1-pi0)*(1-beta_i), or
# not rejected otherwise; the joint pmf is built by convolution. This is the
# proper pFDR (not the E[V]/E[R] ratio), with no closed form but exact.
pfdr_exact <- function(alpha, beta, pi0) {
  m <- length(alpha)
  a <- pi0 * alpha             # P(false positive)
  b <- (1 - pi0) * (1 - beta)  # P(true positive)
  P <- matrix(0, m + 1, m + 1); P[1, 1] <- 1
  for (i in seq_len(m)) {
    ci <- 1 - a[i] - b[i]; Pn <- matrix(0, m + 1, m + 1)
    for (v in 0:(i - 1)) for (s in 0:(i - 1)) {
      p <- P[v + 1, s + 1]; if (p == 0) next
      Pn[v + 1, s + 1] <- Pn[v + 1, s + 1] + p * ci      # no rejection
      Pn[v + 2, s + 1] <- Pn[v + 2, s + 1] + p * a[i]    # false positive
      Pn[v + 1, s + 2] <- Pn[v + 1, s + 2] + p * b[i]    # true positive
    }
    P <- Pn
  }
  num <- 0; prpos <- 0
  for (v in 0:m) for (s in 0:m) {
    if (v + s == 0) next
    p <- P[v + 1, s + 1]; num <- num + (v / (v + s)) * p; prpos <- prpos + p
  }
  c(pFDR = num / prpos, FDR = num, pr_pos = prpos)       # FDR = pFDR * Pr(R>0) = num
}

# A heterogeneous family: tests with very different sample sizes -> different power
m <- 12
pi0 <- 0.6
tau <- 2; eff <- 0.3; cc <- 0.975
n_i <- rep(c(20, 50, 200, 1000), length.out = m)

rates   <- t(vapply(n_i, ab_rates, numeric(2), tau = tau, eff = eff, c = cc))
alpha_i <- rates[, "alpha"]
beta_i  <- rates[, "beta"]
cat("per-test alpha  :", paste(sprintf("%.3f", alpha_i), collapse = " "), "\n")
cat("per-test 1-beta :", paste(sprintf("%.3f", 1 - beta_i), collapse = " "), "\n\n")

# Simulate the family: each test null with prob pi0 (exchangeable), then the
# normal-normal pd > c decision on simulated data.
B <- 100000
kappa_i <- n_i / (n_i + 1 / tau^2)
post_sd <- sqrt(1 / (n_i + 1 / tau^2))

sim <- replicate(B, {
  is_null <- runif(m) < pi0
  theta <- ifelse(is_null, 0, eff)
  ybar  <- rnorm(m, theta, sqrt(1 / n_i))
  pd <- pnorm(abs(kappa_i * ybar) / post_sd)
  rej <- pd > cc
  c(R = sum(rej), V = sum(rej & is_null), n0 = sum(is_null))
})
R <- sim["R", ]; V <- sim["V", ]

FDR_sim  <- mean(ifelse(R == 0, 0, V / R))
pFDR_sim <- mean((V / R)[R > 0])
EV_ER <- mean(V) / mean(R)            # E[V]/E[R], what the formula equals
pFDR_exact <- pfdr_exact(alpha_i, beta_i, pi0)["pFDR"]   # exact E[V/R|R>0], no simulation

# Predictions
abar <- mean(alpha_i); bbar <- mean(beta_i)
pFDR_avg  <- pi0 * abar / (pi0 * abar + (1 - pi0) * (1 - bbar))               # correct
j         <- which.max(n_i)
pFDR_one  <- pi0 * alpha_i[j] / (pi0 * alpha_i[j] + (1 - pi0) * (1 - beta_i[j]))  # single test
r_nbar    <- ab_rates(mean(n_i), tau, eff, cc)                                # rates at mean n
pFDR_nbar <- pi0 * r_nbar["alpha"] / (pi0 * r_nbar["alpha"] + (1 - pi0) * (1 - r_nbar["beta"]))

cat(sprintf("pi0 set                            = %.3f\n", pi0))
cat(sprintf("realized null proportion           = %.3f   (holds regardless of power)\n", mean(sim["n0", ]) / m))
cat("\n")
cat(sprintf("pFDR  simulated   E[V/R | R>0]     = %.4f\n", pFDR_sim))
cat(sprintf("pFDR  exact       E[V/R | R>0]     = %.4f   <- proper pFDR (joint law of V,S)\n", pFDR_exact))
cat(sprintf("ratio simulated   E[V]/E[R]        = %.4f\n", EV_ER))
cat(sprintf("ratio formula     E[V]/E[R]        = %.4f   <- family-averaged approximation\n", pFDR_avg))
cat("\nWrong single-(alpha,beta) plug-ins:\n")
cat(sprintf("  single test (n = %4d)           = %.4f\n", n_i[j], pFDR_one))
cat(sprintf("  rates evaluated at mean n        = %.4f\n", as.numeric(pFDR_nbar)))


# Per-test prior-odds procedure vs the joint statement vs simultaneous
# credible intervals, all computed with the package functions on posterior draws:
#   * pd.adjust(draws, pi0) -> per-test rule  pd.adj > c  (c = 1 - alpha/2)
#   * joint(draws)-> cumulative joint, accept top-k with joint_cum > 1 - alpha
#   * joint(draws, interval = TRUE)-> 95% simultaneous band, declare when it excludes 0
# alpha = 0.05 throughout. Output saved to paper/script/output/df_compare.rda.

set.seed(2026)

library(multibayes)

n <- 50   # observations per test
s <- 2   # prior sd, theta ~ N(0, s^2)
pi0 <- 0.7   # per-test null prior (and true null proportion: correct spec)
alpha <- 0.05  # nominal level used by all three rules
c <- 1 - alpha / 2   # pd cutoff (per-test), 0.975
jlev <- 1 - alpha # cumulative-joint / simultaneous level, 0.95
eff <- 0.3         # |effect| of the non-null tests
B <- 10000  # replications per cell
S <- 4000   # posterior draws per replication (typical MCMC sample size)
ms <- c(4, 8, 10, 20)

kap <- n / (n + 1 / s^2)
psd <- sqrt(1 / (n + 1 / s^2))

one_rep <- function(m) {
  n0 <- rbinom(1, size = m, prob = pi0)
  is_null <- logical(m)
  is_null[sample(m, n0)] <- TRUE
  theta <- ifelse(is_null, 0, eff)

  ybar  <- rnorm(m, theta, sqrt(1 / n))

  # posterior draws (independent normal-normal posteriors), one column per test
  draws <- matrix(rnorm(S * m, mean = rep(kap * ybar, each = S), sd = psd), S, m)
  colnames(draws) <- paste0("H", seq_len(m))
  names(is_null)  <- colnames(draws)

  # (1) per-test prior-odds procedure 
  adj <- suppressWarnings(pd.adjust(draws = draws, pi0 = pi0, null.value = 0))
  rej_p <- setNames(adj$pd.adj > c, rownames(adj))

  # (2) joint statement: top-k rows with joint_cum >= jlev
  jc <- joint(draws, null.value = 0)
  k <- sum(jc$joint_cum >= jlev)
  rej_j <- setNames(logical(m), colnames(draws))
  if (k > 0) rej_j[rownames(jc)[seq_len(k)]] <- TRUE

  # (3) simultaneous credible intervals
  bd <- joint(draws, interval = TRUE, prob = jlev)
  rej_s <- setNames(bd$lower > 0 | bd$upper < 0, rownames(bd))

  err <- function(r) is_null[names(r)] & r  # Type I: declared true null
  fdp <- function(r) if (sum(r) == 0) 0 else sum(err(r)) / sum(r)
  nn <- sum(!is_null)
  pw <- function(r) if (nn == 0) NA_real_ else sum(r & !is_null[names(r)]) / nn

  c(procedure_FDR = fdp(rej_p), procedure_FWER = as.numeric(any(err(rej_p))), procedure_power = pw(rej_p),
    joint_FDR = fdp(rej_j), joint_FWER = as.numeric(any(err(rej_j))), joint_power = pw(rej_j),
    simultaneous_FDR = fdp(rej_s), simultaneous_FWER = as.numeric(any(err(rej_s))), simultaneous_power = pw(rej_s))
}

df_compare <- do.call(rbind, lapply(ms, function(m) {
  data.frame(m = m, as.list(rowMeans(replicate(B, one_rep(m)), na.rm = TRUE)))
}))

save(df_compare, file = "paper/script/output/df_compare.rda")


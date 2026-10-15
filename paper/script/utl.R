# Functions used by simulation.R and results.R

# Simulation -------------------------------------------------------------------

## Normal model: pd from the conjugate normal-normal posterior
bayes_posterior_analytical <- function(x, sigma = 1, tau0 = 1, mu0=0) {
  n <- length(x)
  bar_y <- mean(x)
  post_mean <- (bar_y * (n / sigma^2) + mu0 * (1 / tau0^2)) / (n / sigma^2 + 1 / tau0^2)
  post_sd <- sqrt(1 / (n / sigma^2 + 1 / tau0^2))
  pd <- pnorm(abs(post_mean - 0) / post_sd)
  return(pd)
}

## nsim replications of one cell; each replication stores the raw pd values,
## the direction favored by each posterior (sgn) and the true parameters. With
## compare = TRUE it also stores the claims of the joint statement and of the
## 95% simultaneous band (+1 / -1 = claimed direction, 0 = no claim),
## computed with multibayes::joint() from 4,000 draws of each posterior.
sim_ind <- function(n, m, r, s, eff, nsim, pi0_true, compare = FALSE) {
  replicate(nsim, {
    # cov
    S  = r + diag(1 - r, m)
    # --- null/non-null assignment ---
    n0 <- rbinom(1, size = m, prob = pi0_true)
    is_null <- logical(m)
    is_null[sample(m, n0)] <- TRUE
    effsim <- rep(eff, m)
    effsim[is_null] <- 0

    # simulate data
    X <- MASS::mvrnorm(n, effsim, Sigma = S)

    # bayes
    pd <- apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))

    # direction favored by the posterior 
    sgn <- sign(colMeans(X))

    out <- list(pd = pd, sgn = sgn, effsim = effsim)

    if (compare) {
      # posterior draws (conjugate normal-normal posterior), one column per test
      post_mean <- colMeans(X) * n / (n + 1 / s^2)
      post_sd   <- sqrt(1 / (n + 1 / s^2))
      draws <- matrix(rnorm(4000 * m, rep(post_mean, each = 4000), post_sd), 4000, m)
      colnames(draws) <- paste0("H", seq_len(m))

      # joint statement: the top-k claims with joint probability >= c = .975,
      # so that, as the other two rules, it reduces to pd >= .975 for one parameter
      jc <- multibayes::joint(draws, null.value = 0)
      k  <- sum(jc$joint_cum >= 0.975)
      jdir <- sign(jc$median.est) * (seq_len(nrow(jc)) <= k)
      out$joint <- jdir[match(colnames(draws), rownames(jc))]

      # 95% simultaneous band: a direction is claimed when the band excludes 0
      bd <- multibayes::joint(draws, interval = TRUE, prob = 0.95)
      bdir <- ifelse(bd$lower > 0, 1, ifelse(bd$upper < 0, -1, 0))
      out$simultaneous <- bdir[match(colnames(draws), rownames(bd))]
    }

    out
  }, simplify = FALSE)
}


# Results ----------------------------------------------------------------------

## threshold c* of the prior-odds rule for a specified pi0, from
## multibayes::pd.adjust() (c = .975; pi0 = .5 gives the unadjusted rule)
cstar <- function(pi0) attr(multibayes::pd.adjust(pd = 0.5, pi0 = pi0), "threshold.adj")

## FWER, average power and FDR of one procedure in one cell: the prior-odds
## rule (pd >= c*) for a given pi0, or the stored claims of the joint statement
## or of the simultaneous band. A claim is false when the effect is null or
## the claimed sign is wrong.
compute_metrics <- function(res_list, method = "prior-odds", pi0 = NA) {
  # matrices with one row per replication and one column per test
  effsim_mat  <- do.call(rbind, lapply(res_list, `[[`, "effsim"))
  is_null_mat <- (effsim_mat == 0)

  # claimed direction: +1 / -1, or 0 when no claim is made
  if (method == "prior-odds") {
    pd_mat  <- do.call(rbind, lapply(res_list, `[[`, "pd"))
    sgn_mat <- do.call(rbind, lapply(res_list, `[[`, "sgn"))
    dir <- (pd_mat >= cstar(pi0)) * sgn_mat
  } else {
    dir <- do.call(rbind, lapply(res_list, `[[`, method))
  }
  claim <- dir != 0

  tp <- claim & !is_null_mat & (dir == sign(effsim_mat))  # correct-sign claims about non-null effects
  fp <- claim & !tp                                       # false claims: null effects or wrong sign

  # FWER: probability of at least one false discovery
  fwer <- mean(rowSums(fp) > 0)

  # average power: proportion of non-null effects detected with the correct sign
  pw <- if (any(!is_null_mat)) sum(tp) / sum(!is_null_mat) else NA_real_

  # FDR: expected V / R, with V / R = 0 when R = 0
  R <- rowSums(claim)
  V <- rowSums(fp)
  fdr <- mean(ifelse(R == 0, 0, V / R))

  data.frame(
    FWER = fwer,
    pw = pw,
    FDR = fdr,
    pi0_check = mean(is_null_mat)   # realized proportion of null effects
  )
}

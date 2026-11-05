# =====================================================================
# Prior-odds rule with priors commonly used in regression models
# Logistic regression with one binary predictor per test, pd computed by
# grid integration under four priors for the slope (normal, Student-t,
# Cauchy, flat). Checks (1) that pd_adj from the formula equals the
# posterior probability under the most skeptical reweighted prior, and
# (2) FDR, FWER, and power of pd >= c*. Writes paper/script/output/sim_glm.rda
# =====================================================================

rm(list = ls())
set.seed(150595)

# 1. Packages ----------------------------------------------------------
pkg <- c("dplyr", "tidyr")
invisible(lapply(pkg, library, character.only = TRUE))

# 2. Settings ----------------------------------------------------------
n_grp <- 50                 # observations per group (two groups per test)
alpha <- qlogis(0.3)        # mean log-odds, x coded -0.5 / 0.5
eff <- 1                    # log-odds ratio of non-null effects
c_thr <- 0.975
pi_tilde <- seq(0.5, 0.9, by = 0.1)
pi0 <- seq(0.5, 1, by = 0.1)  # actual proportion of null effects
m <- c(10, 20)
nsim <- 1e4

# priors on the slope (log density up to a constant); intercept always student_t(3, 0, 2.5)
priors <- list(
  "normal(0, 2.5)"     = function(b) dnorm(b, 0, 2.5, log = TRUE),     # rstanarm default
  "student_t(3, 0, 2.5)" = function(b) dt(b / 2.5, df = 3, log = TRUE),
  "cauchy(0, 2.5)"     = function(b) dcauchy(b, 0, 2.5, log = TRUE),   # Gelman et al. 2008
  "flat"               = function(b) 0 * b                              # brms default
)

# 3. Posterior on a grid -----------------------------------------------
# the likelihood depends on the data only through the successes in each
# group (k0, k1), so pd is computed once for every possible (k0, k1)
a_grid <- seq(-6, 4, by = 0.025)
b_grid <- seq(-6, 6, by = 0.02) + 0.01      # midpoints, none exactly at 0
A <- matrix(a_grid, length(a_grid), length(b_grid))
B <- matrix(b_grid, length(a_grid), length(b_grid), byrow = TRUE)
L1 <- plogis(A + B / 2, log.p = TRUE); M1 <- plogis(A + B / 2, lower.tail = FALSE, log.p = TRUE)
L0 <- plogis(A - B / 2, log.p = TRUE); M0 <- plogis(A - B / 2, lower.tail = FALSE, log.p = TRUE)
lp_a <- dt(A / 2.5, df = 3, log = TRUE)
pos <- B > 0

# Pr(beta > 0 | y) under each prior, and the check of pd_adj at pi_tilde = .8
pplus <- array(NA, c(n_grp + 1, n_grp + 1, length(priors)), dimnames = list(NULL, NULL, names(priors)))
max_diff <- setNames(rep(NA, length(priors)), names(priors))
for (j in seq_along(priors)) {
  lp <- lp_a + priors[[j]](B)
  diffs <- c()
  for (k0 in 0:n_grp) for (k1 in 0:n_grp) {
    ll <- k1 * L1 + (n_grp - k1) * M1 + k0 * L0 + (n_grp - k0) * M0 + lp
    w <- exp(ll - max(ll))
    pp <- sum(w[pos]) / sum(w)
    pplus[k0 + 1, k1 + 1, j] <- pp
    # direct check on outcomes that can occur: reweight the halves of the prior so the
    # favored direction gets .2, recompute its posterior probability, compare with the formula
    if (abs(k1 - k0) <= 30 && min(k0, k1) >= 3) {
      fav_pos <- pp > 0.5
      wr <- w * ifelse(pos == fav_pos, 0.2, 0.8)
      direct <- sum(wr[pos == fav_pos]) / sum(wr)
      pd <- max(pp, 1 - pp)
      formula <- pd * 0.2 / (pd * 0.2 + (1 - pd) * 0.8)
      diffs <- c(diffs, abs(direct - formula))
    }
  }
  max_diff[j] <- max(diffs)
}

# check: largest difference between direct computation and formula
max_diff

# 4. Simulation --------------------------------------------------------
cells <- expand.grid(prior = names(priors), m = m, pi0 = pi0, stringsAsFactors = FALSE)
res <- list()
for (r in seq_len(nrow(cells))) {
  mm <- cells$m[r]
  null <- matrix(runif(nsim * mm) < cells$pi0[r], nsim, mm)   # each effect null with prob pi0
  beta <- ifelse(null, 0, eff)
  k0 <- matrix(rbinom(nsim * mm, n_grp, plogis(alpha - beta / 2)), nsim, mm)
  k1 <- matrix(rbinom(nsim * mm, n_grp, plogis(alpha + beta / 2)), nsim, mm)
  pp <- matrix(pplus[cbind(as.vector(k0) + 1, as.vector(k1) + 1, match(cells$prior[r], names(priors)))], nsim, mm)
  pd <- pmax(pp, 1 - pp)
  wrong <- null | pp < 0.5                  # false if null or claimed with the wrong sign
  for (pt in pi_tilde) {
    c_star <- c_thr * pt / ((1 - c_thr) * (1 - pt) + c_thr * pt)
    claim <- pd >= c_star
    R <- rowSums(claim); V <- rowSums(claim & wrong)
    res[[length(res) + 1]] <- data.frame(
      cells[r, ], pi_tilde = pt,
      FDR = mean(V / pmax(R, 1)),
      FWER = mean(V > 0),
      pw = sum(claim & !null & pp > 0.5) / sum(!null))
  }
}
df_glm <- bind_rows(res)

# check
df_glm |> filter(m == 20, abs(pi_tilde - pmin(pi0, 0.9)) < 1e-9 | pi_tilde == 0.5) |>
  arrange(pi0, pi_tilde, prior)

# 5. Save --------------------------------------------------------------
save(df_glm, max_diff, file = "paper/script/output/sim_glm.rda")

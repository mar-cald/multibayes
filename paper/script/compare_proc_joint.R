# Per-test prior-odds procedure vs the cumulative joint statement vs simultaneous
# credible intervals, scored on the Type I error (declaring a true point null).

set.seed(2026)

source("paper/script/utl.R")
library(multibayes)

n <- 50   # observations per test
s  <- 2   # prior sd
pi0 <- 0.7 # per-test null prior (and true null proportion: correct spec)
c <- 0.975  # pd cutoff (per-test)
jlev <- 0.95  # cumulative-joint / simultaneous level
eff <- 0.3 # |effect| of the non-null tests
B <- 10000  # replications per cell
ms  <- c(4, 8, 10, 20)

one_rep <- function(m) {
  n0 <- rbinom(1, size = m, prob = pi0)
  is_null <- logical(m)
  is_null[sample(m, n0)] <- TRUE
  effsim <- rep(eff, m)
  effsim[is_null] <- 0
  
  # simulate data
  X <- MASS::mvrnorm(n, effsim, Sigma = diag(1, m))
  
  # bayes
  pd <- apply(X, 2, function(x) bayes_posterior_analytical(x, tau0 = s))
  pd_adj <- pd.adjust(pd = pd, pi0 = pi0)$pd.adj
  rej_p <- pd_adj > .975   # per-test procedure

  ord <- order(pd, decreasing = TRUE)  # cumulative joint: top-k with joint >= jlev
  jc <- cumprod(pd[ord])
  kstar <- if (jc[1] > jlev) max(which(jc > jlev)) else 0
  rej_j <- logical(m); if (kstar > 0) rej_j[ord[seq_len(kstar)]] <- TRUE

  rej_s <- pd > 1 - (1 - jlev^(1 / m)) / 2    #simultaneous intervals (Sidak band excludes 0)

  fdp <- function(r) if (sum(r) == 0) 0 else sum(r & is_null) / sum(r)
  nn  <- sum(!is_null)  # average power: correct detections
  pw  <- function(r) if (nn == 0) NA_real_ else sum(r & !is_null) / nn
  c(procedure_FDR  = fdp(rej_p), 
    procedure_FWER = as.numeric(any(rej_p & is_null)), 
    procedure_power = pw(rej_p),
    joint_FDR= fdp(rej_j), joint_FWER = as.numeric(any(rej_j & is_null)), 
    joint_power  = pw(rej_j),
    simultaneous_FDR = fdp(rej_s), simultaneous_FWER = as.numeric(any(rej_s & is_null)), 
    simultaneous_power = pw(rej_s))
}

df_compare <- do.call(rbind, lapply(ms, function(m) {
  data.frame(m = m, as.list(rowMeans(replicate(B, one_rep(m)), na.rm = TRUE)))
}))

save(df_compare, file = "paper/script/output/df_compare.rda")
print(df_compare, digits = 3)

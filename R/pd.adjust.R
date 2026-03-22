#' Prior-odds adjustment for Probability of Direction (pd)
#'
#' The function accepts either a vector of pre-computed \emph{pd} values or
#' a matrix of posterior draws, from which \emph{pd} values are computed
#' internally. The global prior probability that all tested hypotheses are
#' null, \eqn{q}, is decomposed into a per-hypothesis prior
#' \eqn{P(H_0) = q^{1/m}}, where \eqn{m} is the number of hypotheses or,
#' when the correlation structure among parameters is taken into account,
#' the effective number of tests \eqn{m_{\text{eff}}}.
#'
#' @details
#' The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
#' \eqn{P(H_0) = q^{1/m}} and its complement \eqn{P(H_1) = 1 - P(H_0)},
#' the adjusted \emph{pd} is:
#' \deqn{
#'    pd_{adj} = \frac{pd \cdot P(H_1)}{pd \cdot P(H_1) + (1 - pd) \cdot P(H_0)}
#' }
#'
#' Because the prior is conservative (\eqn{P(H_0) > P(H_1)}), the adjustment
#' always moves \emph{pd} toward \eqn{0}. A floor at \eqn{0.5} is applied to
#' \eqn{pd_{adj}}, so the effective result is a shrinkage toward \eqn{0.5}.
#'
#' For direction-agnostic tests (\code{direction = 0}), \emph{pd} is the
#' maximum of the two tail probabilities and is bounded in \eqn{[0.5, 1]}.
#' For directional tests (\code{direction = 1} or \code{-1}), \emph{pd} is
#' the probability mass on the predicted side, \eqn{Pr(\hat\theta > \theta_{null})} or
#' \eqn{Pr(\hat\theta < \theta_{null})}, and is floored at \eqn{0.5} before the adjustment
#' is applied. A value of \eqn{pd = 0.5} should therefore be interpreted as
#' absence of support for the predicted direction; the floor prevents the
#' correction from amplifying contradictory evidence. The raw directional
#' probability before flooring is returned in \code{pd_raw} for directional
#' tests, allowing the researcher to assess whether the floor was triggered
#' and how strongly the data contradicted the predicted direction.
#'
#' When \code{R} is supplied, the effective number of tests \eqn{m_{\text{eff}}}
#' is estimated from the eigenvalues \eqn{\lambda} of the correlation matrix
#' (Cheverud, 2001):
#' \deqn{
#'    m_{\text{eff}} = K \left( 1 - \frac{(K-1)\,\text{Var}(\lambda)}{K^2} \right)
#' }
#' where \eqn{K} is the number of hypotheses.
#'
#' @param pd Numeric vector of \emph{pd} values in \eqn{[0.5, 1]}. For
#'   direction-agnostic tests, this is the maximum of the two tail
#'   probabilities. For directional tests, this is the probability mass
#'   on the predicted side, floored at \eqn{0.5}.
#' @param draws Optional matrix or data frame of posterior draws (columns = parameters).
#'   If provided, \code{pd} is calculated automatically.
#' @param mu0 Numeric scalar or vector. The null (reference) value against which
#'   the posterior is evaluated. A scalar applies the same null to all parameters;
#'   a vector of length equal to \code{ncol(draws)} allows a different null per
#'   parameter. Ignored if \code{pd} is supplied directly. Defaults to \code{0}.
#' @param direction Integer vector of \code{1}, \code{-1}, or \code{0} (or \code{NULL}).
#'   Specifies the predicted direction of each effect: \code{1} for positive
#'   (tests \eqn{Pr(\theta > \theta_{null})}), \code{-1} for negative (tests
#'   \eqn{Pr(\theta < \theta_{null}}), and \code{0} for direction-agnostic (takes the
#'   maximum over both sides). Should be a vector of length \code{ncol(draws)};
#'   a scalar is recycled. Defaults to \code{NULL} (direction-agnostic for all
#'   parameters).
#' @param q Numeric scalar in \eqn{(0, 1)}. The prior probability that
#'   \strong{all} hypotheses are null. Defaults to \code{0.4}.
#' @param R Optional correlation information for computing \eqn{m_{\text{eff}}}.
#'   Accepts \code{TRUE} (correlation estimated from \code{draws}), a numeric
#'   scalar (assumed uniform correlation applied to all parameter pairs), or a
#'   full correlation matrix. When provided, \eqn{m_{\text{eff}}} replaces the
#'   nominal \eqn{m}.
#'
#' @return A \code{data.frame} with one row per hypothesis, containing:
#'   \code{pd} (values used in the adjustment, floored at \eqn{0.5}),
#'   \code{pd_adj} (adjusted values, also floored at \eqn{0.5}), \code{q}
#'   (prior probability of the global null), and \code{m} (number of hypotheses
#'   or effective number of tests). When \code{draws} are supplied,
#'   \code{mean_est} (posterior mean per parameter), \code{mu0} (null reference
#'   values), and \code{direction} are also returned. For directional tests,
#'   \code{pd_raw} is additionally returned, reporting the raw directional
#'   probability before flooring; values below \eqn{0.5} indicate that the
#'   posterior was concentrated opposite to the predicted direction.
#'
#' @references
#' Jeffreys, H. (1938). Significance tests when several degrees of freedom arise
#' simultaneously. \emph{Proceedings of the Royal Society of London. Series A.
#' Mathematical and Physical Sciences, 165}(921), 161--198.
#' <https://doi.org/10.1098/rspa.1938.0052>
#'
#' Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian perspective
#' on the Bonferroni adjustment. \emph{Biometrika, 84}(2), 419--427.
#' <https://doi.org/10.2307/2337467>
#'
#' Cheverud, J. (2001). A simple correction for multiple comparisons in interval
#' mapping genome scans. \emph{Heredity, 87}, 52--58.
#' <https://doi.org/10.1046/j.1365-2540.2001.00901.x>
#'
#' @importFrom stats var cor
#'
#' @export
pd.adjust <- function(pd = NULL, draws = NULL, q = 0.4, mu0 = 0,
                      direction = NULL, R = NULL) {
  
  from_draws <- !is.null(draws)
  
  if (from_draws) {
    draws <- as.matrix(draws)
    p <- ncol(draws)
    
    if (length(mu0) == 1L)       mu0 <- rep(mu0, p)
    if (is.null(direction))      direction <- rep(0L, p)
    if (length(direction) == 1L) direction <- rep(direction, p)
    
    if (length(mu0) != p)       stop("`mu0` must be a scalar or a vector of length `ncol(draws)`.")
    if (length(direction) != p) stop("`direction` must be a scalar or a vector of length `ncol(draws)`.")
    if (!all(direction %in% c(-1L, 0L, 1L))) stop("`direction` must contain only -1, 0, or 1.")
    
    centered <- sweep(draws, 2, mu0, "-")
    
    # Raw directional probability (before flooring)
    pd_raw <- mapply(function(j, d) {
      if      (d ==  1L) mean(centered[, j] > 0)
      else if (d == -1L) mean(centered[, j] < 0)
      else    NA_real_
    }, seq_len(p), direction)
    
    # pd used in adjustment: directional floored at 0.5, agnostic as max
    pd <- mapply(function(j, d) {
      if      (d ==  1L) max(mean(centered[, j] > 0), 0.5)
      else if (d == -1L) max(mean(centered[, j] < 0), 0.5)
      else    max(mean(centered[, j] > 0), mean(centered[, j] < 0))
    }, seq_len(p), direction)
  }
  
  if (is.null(pd)) stop("Either `pd` or `draws` must be provided.")
  
  m <- length(pd)
  
  stopifnot(
    "`pd`: must be numeric in [0.5, 1]" = is.numeric(pd) && all(pd >= 0.5 & pd <= 1, na.rm = TRUE),
    "`q`: must be a single number in (0, 1)" = length(q) == 1L && q > 0 && q < 1
  )
  
  # Effective number of tests (Cheverud, 2001)
  if (!is.null(R)) {
    
    if (isTRUE(R) && is.null(draws))
      stop("`R = TRUE` requires `draws` to be provided.")
    
    if (isTRUE(R) && !is.null(draws))
      R <- cor(draws)
    
    if (is.numeric(R) && length(R) == 1L) {
      R_val <- R
      R <- matrix(R_val, nrow = length(pd), ncol = length(pd))
      diag(R) <- 1
    }
    
    ev   <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    K    <- length(ev)
    v_ev <- if (K == 1L) 0 else var(ev)
    m    <- K * (1 - ((K - 1) * v_ev) / (K^2))
  }
  
  # Prior-odds adjustment
  prior_H0 <- q^(1 / m)
  
  if (prior_H0 > 0.5) {
    prior_H1 <- 1 - prior_H0
    pd_adj <- (pd * prior_H1) / (pd * prior_H1 + (1 - pd) * prior_H0)
  } else {
    warning("Pr(H0_i) <= 0.5 (non-conservative prior); returning unadjusted pd.")
    pd_adj <- pd
  }
  
  # Floor at 0.50 for all tests
  pd_adj <- pmax(pd_adj, 0.5)
  
  if (from_draws) {
    out <- data.frame(
      mean_est  = round(colMeans(draws), 4),
      mu0       = mu0,
      direction = direction,
      pd        = round(pd, 4),
      pd_adj    = round(pd_adj, 4),
      q         = rep(q, length(pd)),
      m         = round(rep(m, length(pd)), 4)
    )
    # Add pd_raw only when at least one directional test is present
    if (any(direction != 0L))
      out$pd_raw <- round(pd_raw, 4)
    out
  } else {
    data.frame(
      pd     = round(pd, 4),
      pd_adj = round(pd_adj, 4),
      q      = rep(q, length(pd)),
      m      = round(rep(m, length(pd)), 4)
    )
  }
}

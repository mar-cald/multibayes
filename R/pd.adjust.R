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
#' When \code{R} is supplied, the effective number of tests \eqn{m_{\text{eff}}}
#' is estimated from the eigenvalues \eqn{\lambda} of the correlation matrix
#' (Cheverud, 2001):
#' \deqn{
#'    m_{\text{eff}} = K \left( 1 - \frac{(K-1)\,\text{Var}(\lambda)}{K^2} \right)
#' }
#' where \eqn{K} is the number of hypotheses.
#'
#' @param pd Numeric vector of \emph{pd} values in \eqn{[0.5, 1]}.
#' @param draws Optional matrix or data frame of posterior draws (columns = parameters).
#'   If provided, \code{pd} is calculated automatically.
#' @param mu0 Numeric scalar or vector. The null (reference) value against which
#'   the posterior is evaluated. A scalar applies the same null to all parameters;
#'   a vector of length equal to \code{ncol(draws)} allows a different null per
#'   parameter. Ignored if \code{pd} is supplied directly. Defaults to \code{0}.
#' @param direction Integer vector of \code{1}, \code{-1}, or \code{0} (or \code{NULL}).
#'   Specifies the expected direction of each effect: \code{1} for positive
#'   (tests \code{draws > mu0}), \code{-1} for negative (tests \code{draws < mu0}),
#'   and \code{0} for direction-agnostic (takes the maximum over both sides).
#'   Should be a vector of length \code{ncol(draws)} to specify a different
#'   direction per parameter; a scalar is recycled across all parameters.
#'   Should be specified when \code{mu0 != 0}. Defaults to \code{NULL} (direction-agnostic
#'   for all parameters).
#' @param q Numeric scalar in \eqn{(0, 1)}. The prior probability that
#'   \strong{all} hypotheses are null. Defaults to \code{0.4}.
#' @param m Positive integer. The number of tested hypotheses.
#'   Defaults to \code{length(pd)}. Overridden if \code{R} is provided.
#' @param R Optional correlation matrix of the posterior draws. Can be provided
#'   directly as a matrix or a single scalar, or computed automatically by the
#'   function when posterior draws are supplied (set \code{R = TRUE}). When
#'   provided, \eqn{m_{\text{eff}}} is calculated from the correlation structure
#'   and used in place of \code{m}.
#'
#' @return A \code{data.frame} with one row per hypothesis, containing the
#'   following columns: \code{pd} (original values), \code{pd_adj} (adjusted
#'   values), \code{q} (prior probability of the global null), and \code{m}
#'   (number of hypotheses or effective number of tests). When \code{draws} are
#'   supplied, \code{mean_est} (posterior mean per parameter), \code{mu0} (null
#'   reference values), and \code{direction} (declared direction per hypothesis;
#'   \code{0} indicates direction-agnostic) are also returned.
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
                      direction = NULL, m = NULL, R = NULL) {
  
  from_draws <- !is.null(draws)
  
  if (from_draws) {
    draws <- as.matrix(draws)
    p <- ncol(draws)
    
    if (length(mu0) == 1L)       mu0       <- rep(mu0, p)
    if (is.null(direction))      direction <- rep(0L, p)
    if (length(direction) == 1L) direction <- rep(direction, p)
    
    if (length(mu0) != p)       stop("`mu0` must be a scalar or a vector of length `ncol(draws)`.")
    if (length(direction) != p) stop("`direction` must be a scalar or a vector of length `ncol(draws)`.")
    if (!all(direction %in% c(-1L, 0L, 1L))) stop("`direction` must contain only -1, 0, or 1.")
    
    centered <- sweep(draws, 2, mu0, "-")
    
    pd <- mapply(function(j, d) {
      if      (d ==  1L) mean(centered[, j] > 0)
      else if (d == -1L) mean(centered[, j] < 0)
      else               max(mean(centered[, j] > 0), mean(centered[, j] < 0))
    }, seq_len(p), direction)
    
    if (is.null(m)) m <- p
  }
  
  if (is.null(pd)) stop("Either `pd` or `draws` must be provided.")
  if (is.null(m))  m  <- length(pd)
  
  stopifnot(
    "`pd`: must be numeric in [0.5, 1]" = is.numeric(pd) && all(pd >= 0.5 & pd <= 1, na.rm = TRUE),
    "`q`: must be a single number (0, 1)" = length(q) == 1L && q > 0 && q < 1,
    "`m`: must be >= 1" = length(m) == 1L && m >= 1
  )
  
  # Effective number of tests (Cheverud, 2001)
  if (!is.null(R)) {
    
    if (isTRUE(R) && is.null(draws)) {
      stop("`R = TRUE` requires `draws` to be provided.")
    }
    
    if (isTRUE(R) && !is.null(draws)) {
      R <- cor(draws)
    }
    
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
    pd_adj   <- (pd * prior_H1) / (pd * prior_H1 + (1 - pd) * prior_H0)
  } else {
    warning("Pr(H0_i) <= 0.5 (Non-conservative prior); returning unadjusted pd.")
    pd_adj <- pd
  }
  
  pd_adj <- ifelse(pd_adj < 0.50, 0.50, pd_adj)
  
  if (from_draws) {
    data.frame(mean_est = round(colMeans(draws),4),
               mu0 = mu0, direction = direction, 
               pd = round(pd,4), 
               pd_adj = round(pd_adj, 4),
               q = rep(q, length(pd)), m = round(rep(m, length(pd)),4))
  } else {
    data.frame(pd = round(pd,4), pd_adj =  round(pd_adj, 4),
               q = rep(q, length(pd)), m = round(rep(m, length(pd)),4))
  }
}

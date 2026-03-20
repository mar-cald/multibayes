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
#' the adjusted *pd* is:
#' \deqn{
#'    pd_{adj} = \frac{pd \cdot P(H_1)}{pd \cdot P(H_1) + (1 - pd) \cdot P(H_0)}
#' }
#'
#' When `R` is supplied, the effective number of tests \eqn{m_{eff}}
#' is estimated from the eigenvalues \eqn{\lambda} of the correlation matrix
#' (Cheverud, 2001):
#' \deqn{
#'    m_{eff} = K \left( 1 - \frac{(K-1) \text{Var}(\lambda)}{K^2} \right)
#' }
#' where \eqn{K} is the number of hypotheses.
#'
#' @param pd Numeric vector of *pd* values in \eqn{[0.5, 1]}.
#' @param draws Optional matrix or data frame of posterior draws (columns = parameters).
#'   If provided, `pd` is calculated automatically.
#' @param q Numeric scalar in \eqn{(0, 1)}: the prior probability that
#'   **all** hypotheses are null. Defaults to `0.4`.
#' @param m Positive integer. The number of tested hypotheses.
#'   Defaults to `length(pd)`. Overridden if `R` is provided.
#' @param R Optional correlation matrix of the posterior draws. Can be provided
#' directly as a matrix or a single scalar, or computed automatically by the function
#' when posterior draws are supplied (set R = TRUE).
#' When provided, \eqn{m_{\text{eff}}} is calculated from the correlation structure and used in place of m.
#'
#' @return A `data.frame` containing original `pd`, `pd_adj`, and the
#'   parameters `q` and `m` used for the correction.
#'
#' @references
#' Jeffreys, H. (1938). Significance tests when several degrees of freedom arise simultaneously.Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences,165(921), 161–198.
#' <https://doi.org/10.1098/rspa.1938.0052>
#' 
#' Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian Perspective on the Bonferroni Adjustment.
#' Biometrika, 84(2), 419–427. <http://www.jstor.org/stable/2337467>
#'
#' Cheverud, J. A simple correction for multiple comparisons in interval mapping genome scans. Heredity 87, 52–58 (2001). <https://doi.org/10.1046/j.1365-2540.2001.00901.x>
#'
#' @importFrom matrixStats colMeans2
#' @importFrom stats var cor
#'
#' @export
pd.adjust <- function(pd = NULL, draws = NULL, q = 0.4, m = NULL, R = NULL) {
  
  # Input handling and validation 
  if (!is.null(draws)) {
    draws <- as.matrix(draws)
    pd <- pmax(
      matrixStats::colMeans2(draws > 0),
      matrixStats::colMeans2(draws < 0)
    )
    if (is.null(m)) m <- ncol(draws)
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
    
    ev  <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    K   <- length(ev)
    v_ev <- var(ev)
    m   <- K * (1 - ((K - 1) * v_ev) / (K^2))
  }
  
  # Prior-odds adjustment 
  prior_H0 <- q^(1 / m)
  
  if (prior_H0 > 0.5) {
    prior_H1 <- 1 - prior_H0
    pd_adj   <- (pd * prior_H1) / (pd * prior_H1 + (1 - pd) * prior_H0)
  } else {
    warning("Pr(H0) <= 0.5 (Non-conservative prior); returning unadjusted pd.")
    pd_adj <- pd
  }
  
  pd_adj <- ifelse(pd_adj < 0.50, 0.50, pd_adj)
  
  data.frame(pd = pd, pd_adj = pd_adj, q = rep(q, length(pd)), 
             m = rep(m, length(pd)))
}

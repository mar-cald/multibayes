#' Prior-odds adjustment for Probability of Direction (pd)
#'
#' Adjusts a vector of Probability of Direction (*pd*) values using a global
#' prior probability that **all** tested hypotheses are null, \eqn{q}. The adjustment
#' converts \eqn{q} into a per-hypothesis prior probability \eqn{H_{0_m} = q^{1/m}},
#' where \eqn{m} is the family size, and then reweights each *pd* accordingly.
#'
#' @param pd Numeric vector of *pd* values (typically \eqn{pd \in [0.5, 1]}).
#' @param q Numeric scalar in \eqn{(0, 1)} giving the prior probability that
#'   **all** hypotheses in the family are null. Defaults to `0.4`.
#' @param m Positive integer giving the number of tested hypotheses. Defaults to
#'   `length(pd)`.
#' @param post_corr Optional correlation matrix of posterior draws (`NULL` by default).
#'   If provided, must be a square matrix with dimensions matching `length(pd)`.
#'   Adjusts the effective number of tests using eigenvalue-based correction 
#'   to account for dependence between hypotheses.
#'   
#' @return A numeric vector of the same length as `pd`, containing adjusted *pd*
#'   values (or the original `pd` if \eqn{H_{0_m} < 0.5}).
#'
#' @export
#' @examples
#' # Without correlation adjustment
#' pd <- c(0.99, 0.98, 0.978)
#' prior_adj(pd, q = 0.4)
#' #> [1] 0.9725000 0.9459554 0.9407567
#' 
#' # With correlation adjustment
#' corr_mat <- matrix(c(1.0, 0.5, 0.3,
#'                      0.5, 1.0, 0.4,
#'                      0.3, 0.4, 1.0), nrow = 3)
#' prior_adj(pd, q = 0.4, post_corr = corr_mat)
#' #>[1] 0.9759573 0.9525872 0.9479914



prior_adj = function(pd, q = 0.4, m = length(pd), post_corr = NULL) {
  
  stopifnot(
    "`pd`: must be numeric"        = is.numeric(pd),
    "`pd`: must be in [0.5, 1]"    = all(is.finite(pd)) && all(pd >= 0.5 & pd <= 1),
    "`q`: must be a single number" = length(q) == 1L && is.finite(q),
    "`q`: must be in (0, 1)"       = (q > 0) && (q < 1),
    "`m`: must be >= 1"            = length(m) == 1L && is.finite(m) && (m >= 1)
  )
  
  # Validate post_corr if provided
  if (!is.null(post_corr)) {
    stopifnot(
      "`post_corr`: must be a matrix" = is.matrix(post_corr),
      "`post_corr`: must be square" = ncol(post_corr) == nrow(post_corr),
      "`post_corr`: dimensions must match length(pd)" = ncol(post_corr) == length(pd)
    )
    
    # Effective number of tests correction
    eigen_vals = eigen(post_corr)$values
    v = sum((eigen_vals - 1)^2 / (m - 1))
    m = m * (1 - (m - 1) * v / (m^2))
  }
  
  H0 = q^(1 / m)
  H1 = 1 - H0
  
  if (H0 >= 0.5) {
    pd_adj = (pd * H1) / (pd * H1 + (1 - pd) * H0)
  } else {
    warning("Pr(H0) < 0.5; returning unadjusted pd")
    pd_adj = pd
  }
  
  return(as.vector(pd_adj))
}

#' Prior-odds adjustment for Probability of Direction (pd)
#'
#' Adjusts a vector of Probability of Direction (*pd*) values using a global
#' prior probability that all tested hypotheses are null, \eqn{q}. The adjustment
#' converts \eqn{q} into a per-hypothesis prior probability \eqn{H0 = q^{1/m}},
#' where \eqn{m} is the family size, and then reweights each *pd* accordingly.
#'
#' @param pd Numeric vector of *pd* values (typically \eqn{\in [0.5, 1]}).
#' @param q Numeric scalar in (0, 1). Interpreted as the prior probability that
#'   **all** hypotheses in the family are null.
#'@param m Numeric scalar (>0). Number of tested hypotheses.
#'   
#' @return A numeric vector of the same length as `pd`, containing adjusted *pd*
#'   values (or the original `pd` if \eqn{H0 < 0.5}).
#'
#' @export
#' @examples
#' pd = c(0.55, 0.80, 0.97)
#' prior_adj(pd, q = 0.5)
#'


prior_adj = function(pd, q = 0.5, m = length(pd)) {
  
  stopifnot(
    "`pd` must be numeric"        = is.numeric(pd),
    "`pd` must be in [0.5, 1]"    = all(is.finite(pd)) && all(pd >= 0.5 & pd <= 1),
    "`q` must be a single number" = length(q) == 1L && is.finite(q),
    "`q` must be in (0, 1)"       = (q > 0) && (q < 1),
    "`m` must be >= 1"            = length(m) == 1L && is.finite(m) && (m >= 1)
  )
  
  H0 = q^(1 / m)
  H1 = 1 - H0
  
  if (H0 >= 0.5) {
    pd_adj = (pd * H1) / (pd * H1 + (1 - pd) * H0)
    
  } else {
    warning("Pr(H0) < 0.5; returning unadjusted pd")
    pd_adj = pd
  }
  
  return(pd_adj)
}



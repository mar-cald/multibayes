#' Prior-odds adjustment for Probability of Direction \emph{pd}
#'
#' The function accepts either a vector of pre-computed \emph{pd} values or
#' a matrix of posterior draws, from which \emph{pd} values are computed
#' internally. Both direction-agnostic and directional tests are supported:
#' the \code{direction} argument controls which formulation is applied per
#' hypothesis. The global prior probability that all tested hypotheses are
#' null, \eqn{\Pi_0}, is decomposed into a per-hypothesis prior
#' \eqn{\pi_0 = \Pi_0^{1/m}}, where \eqn{m} is the number of hypotheses.
#'
#' @details
#' The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
#' \eqn{\pi_0 = \Pi_0^{1/m}} and its complement \eqn{\pi_1 = 1 - \pi_0},
#' the adjusted \emph{pd} is:
#' \deqn{
#'    pd_{adj} = \frac{pd \pi_1}{pd \pi_1 + (1-pd)\pi_0}
#' }
#'
#' Because the prior is conservative (\eqn{\pi_0 > \pi_1}), the adjustment
#' always shrinks \emph{pd} toward its lower bound.
#'
#' \strong{Direction-agnostic tests} (\code{"two.sided"}): \emph{pd} is
#' defined as \eqn{\max\!\big(\Pr(\hat\theta > \theta_\text{null}),\,
#' \Pr(\hat\theta < \theta_\text{null})\big)} and is bounded in
#' \eqn{[0.5, 1]} by construction. \eqn{pd_{adj}} is also floored at
#' \eqn{0.5}, so the adjustment produces shrinkage toward \eqn{0.5}.
#'
#' \strong{Directional tests} (\code{"greater"} or \code{"less"}): \emph{pd}
#' is the raw posterior probability mass on the predicted side,
#' \eqn{\Pr(\hat\theta > \theta_\text{null})} or
#' \eqn{\Pr(\hat\theta < \theta_\text{null})}, and is defined on \eqn{[0, 1]}.
#' Values of \emph{pd} below \eqn{0.5} indicate that the posterior is
#' concentrated opposite to the predicted direction; the adjustment will
#' further shrink such values toward \eqn{0}, reflecting the combined weight
#' of the data and the conservative prior against the hypothesis.
#'
#' Mixed use of directional and direction-agnostic tests within the same call
#' is supported: each element of \code{direction} is handled independently,
#' and the same prior-odds adjustment is applied uniformly across all
#' hypotheses regardless of their directionality.
#'
#' @param pd Numeric vector of \emph{pd} values. For direction-agnostic tests,
#'   values must be in \eqn{[0.5, 1]}. For directional tests, values are raw
#'   one-sided probabilities in \eqn{[0, 1]}. Ignored if \code{draws} is
#'   supplied.
#' @param draws Optional matrix or data frame of posterior draws (columns = parameters).
#'   If provided, \emph{pd} values are computed automatically from the draws
#'   according to \code{direction} and \code{null.value}.
#' @param p Numeric scalar in \eqn{(0, 1)}. The prior probability that
#'   \strong{all} hypotheses are null simultaneously. 
#' @param null.value Numeric scalar or vector. The null (reference) value against 
#'  which the posterior is evaluated, specified on the scale of the posterior. 
#'  A single scalar applies the same null to all parameters; a vector of length \code{ncol(draws)} 
#'  assigns a distinct null to each parameter. Ignored when \code{pd} is supplied directly. 
#'  Defaults to \code{0}.
#' @param direction Character vector of \code{"greater"}, \code{"less"}, or
#'   \code{"two.sided"} (or \code{NULL}).
#'   Specifies the testing mode for each hypothesis: \code{"greater"} for a
#'   positive directional test (\eqn{\Pr(\theta > \theta_\text{null})}),
#'   \code{"less"} for a negative directional test
#'   (\eqn{\Pr(\theta < \theta_\text{null})}), and \code{"two.sided"} for
#'   direction-agnostic testing (maximum over both sides). A scalar is recycled
#'   to match the number of parameters; a mixed vector allows different modes
#'   across hypotheses. Defaults to \code{NULL} (direction-agnostic for all
#'   parameters).
#'
#' @return A \code{data.frame} with one row per hypothesis, containing:
#'   \code{pd} (values used in the adjustment), \code{pd.adj} (adjusted
#'   values), \code{p} (prior probability of the global null), and \code{m}
#'   (number of tests), \code{pm} (prior probability for the null for each test). 
#'   For direction-agnostic tests, both \code{pd} and \code{pd.adj} are bounded in \eqn{[0.5, 1]}; for
#'   directional tests, both are on \eqn{[0, 1]}, with values below \eqn{0.5}
#'   indicating that the data (and the adjustment) favoured the opposite
#'   direction. When \code{draws} are supplied, the output additionally
#'   includes \code{mean.est} (posterior mean per parameter), \code{null.value}
#'   (null reference values), and \code{direction}.
#'
#' @examples
#' # From a vector of pd values 
#' pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763, 
#' H5 = 0.891, H6 = 0.987)
#' pd.adjust(pd = pd_values, p = 0.4)
#'
#' # Simulate draws
#' Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
#' mu    <- c(1, -0.1, 0.8, 0, 2, 3)
#' draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
#' colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")
#'
#' # From posterior draws
#' pd.adjust(draws = draws, p = 0.4, null.value = 0)
#'
#' # Mix of directional and agnostic tests with parameter-specific nulls
#' pd.adjust(draws = draws, p = 0.2, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
#'           direction = c("greater", "two.sided", "greater", 
#'           "two.sided", "greater", "greater"))
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
#' @importFrom stats var cor
#'
#' @export
pd.adjust <- function(pd = NULL, draws = NULL, p = NULL, null.value = 0,
                      direction = NULL) {
  
  stopifnot(
    "`p`: must be a single number in (0, 1)" = length(p) == 1L && p > 0 && p < 1
  )
  
  if(!is.null(pd)){
    if(!is.numeric(pd) || any(pd < 0.5 | pd > 1, na.rm = TRUE)){
      stop("`pd` must be numeric in [0.5, 1].")
    } 
    if(!is.null(direction)){
      warning("`direction` cannot be specified, fixed to `two.sided`")
      direction <- rep("two.sided", length(pd))
    } 
  }
  
  from_draws <- !is.null(draws)
  
  if (from_draws) {
    draws <- as.matrix(draws)
    nc <- ncol(draws)
    
    if (length(null.value) == 1L) null.value <- rep(null.value, nc)
    if (is.null(direction))      direction <- rep("two.sided", nc)
    if (length(direction) == 1L) direction <- rep(direction, nc)
    
    if (length(null.value) != nc) stop("`null.value` must be a scalar or a vector of length `ncol(draws)`.")
    if (length(direction) != nc)  stop("`direction` must be a scalar or a vector of length `ncol(draws)`.")
    if (!all(direction %in% c("two.sided", "less", "greater"))) stop("`direction` must contain only two.sided, less, greater.")
    
    centered <- sweep(draws, 2, null.value, "-")
    
    # pd: raw one-sided probability for directional tests, max for agnostic
    pd <- mapply(function(j, d) {
      if      (d ==  "greater") mean(centered[, j] > 0)
      else if (d == "less") mean(centered[, j] < 0)
      else    max(mean(centered[, j] > 0), mean(centered[, j] < 0))
    }, seq_len(nc), direction)
  }
  
  if (is.null(pd)) stop("Either `pd` or `draws` must be provided.")
  
  m <- length(pd)
  
  # Prior-odds adjustment
  prior_H0 <- p^(1 / m)
  
  if (prior_H0 > 0.5) {
    prior_H1 <- 1 - prior_H0
    pd.adj <- (pd * prior_H1) / (pd * prior_H1 + (1 - pd) * prior_H0)
  } else {
    warning("Pr(H0_i) <= 0.5 (non-conservative prior); returning unadjusted pd.")
    pd.adj <- pd
  }
  
  # Floor pd.adj at 0.5 for agnostic (two-sided) tests
  if (from_draws) {
    if (any(direction == "two.sided" & pd.adj < 0.5)) warning("some pd.adj have been floored to 0.5.")
    pd.adj[direction == "two.sided"] <- pmax(pd.adj[direction == "two.sided"], 0.5)
  } else {
    if (any(pd.adj < 0.5)) warning("some pd.adj have been floored to 0.5.")
    pd.adj <- pmax(pd.adj, 0.5)
  }
  
  if (from_draws) {
    data.frame(
      mean.est  = round(colMeans(draws), 4),
      null.value= null.value,
      pd        = round(pd, 4),
      pd.adj    = round(pd.adj, 4),
      p         = round(rep(p, length(pd)),4),
      m         = round(rep(m, length(pd)), 4),
      pm = round(rep(p^(1/m), length(pd)),4),
      direction = direction
    )
  } else {
    data.frame(
      mean.est  = NA,
      null.value= NA,
      pd        = round(pd, 4),
      pd.adj    = round(pd.adj, 4),
      p         = round(rep(p, length(pd)),4),
      m         = round(rep(m, length(pd)), 4),
      pm = round(rep(p^(1/m), length(pd)),4),
      direction = rep("two.sided",length(pd))
    )
  }
}


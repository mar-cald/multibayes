#' Simultaneous credible intervals from joint posterior draws
#'
#' Computes simultaneous equitailed credible intervals for a set of parameters,
#' guaranteeing that **all** parameters stay within their bounds simultaneously
#' with probability \eqn{1 - \alpha}.
#'
#' @details
#' Simultaneous coverage is calibrated by examining how extreme each draw is
#' across all parameters jointly. For each draw, the minimum tail probability
#' across all \eqn{K} parameters is computed. The simultaneous threshold
#' is the \eqn{(1 - \text{prob})}-quantile of these minima: only
#' \eqn{\alpha} of draws are more extreme than this value in at least one
#' parameter simultaneously. See Goeman et al. (2026) for details.
#'
#' A closely related implementation is \code{\link[credsubs]{sim.cred.band}}
#' in the \pkg{credsubs} package (Schnell et al., 2020), which applies the same
#' procedure to function-valued parameters over a covariate space.
#'
#' @param draws Numeric matrix or data frame of posterior draws
#'   (rows = draws, columns = parameters).
#' @param prob Numeric scalar in \eqn{(0, 1)}. Joint coverage probability
#'   (default \code{0.95}).
#' @param est.FUN Function for point estimates (default \code{median}).
#'
#' @return A named list with:
#' \describe{
#'   \item{lower}{Lower bounds, one per parameter.}
#'   \item{upper}{Upper bounds, one per parameter.}
#'   \item{prob}{Requested joint coverage probability.}
#'   \item{cq}{Simultaneous threshold: the
#'     \eqn{(1 - \text{prob})}-quantile of the worst-case tail probability
#'     distribution. Always smaller than \eqn{(1 - \text{prob}) / 2}, yielding
#'     wider intervals than marginal credible intervals.}
#'   \item{est}{Point estimates via \code{est.FUN}.}
#'   \item{est.FUN}{The estimation function used.}
#' }
#'
#' @note
#' If you use this function in published research, please cite:
#' \itemize{
#'   \item Goeman, J., Calderan, M., & Solari, A. (2026). Bonferroni for
#'     Bayesian: Multiplicity correction based on joint credibility. \emph{(in press)}.
#'   \item The package: \url{https://github.com/yourusername/multibayes}.
#' }
#'
#' @references
#' Goeman, J., Calderan, M., & Solari, A. (2026). Bonferroni for Bayesian:
#' Multiplicity correction based on joint credibility. \emph{(in press)}.
#'
#' Schnell, P. M., Fiecas, M., & Carlin, B. P. (2020). credsubs:
#' Multiplicity-adjusted subset identification.
#' \emph{Journal of Statistical Software}, 94(7), 1--22.
#' \doi{10.18637/jss.v094.i07}
#'
#' @importFrom stats quantile median setNames
#' @importFrom matrixStats colRanks rowMins colQuantiles
#' @export
#' @examples
#' mu    <- c(4, 0, -2)
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
#' draws <- MASS::mvrnorm(2000, mu, Sigma)
#' colnames(draws) <- c("theta1", "theta2", "theta3")
#'
#' joint(draws, prob = 0.95)
joint <- function(draws, prob = 0.95, est.FUN = median) {
  
  if (is.data.frame(draws)) draws <- as.matrix(draws)
  
  stopifnot(
    "`draws`: must be a numeric matrix-like object" = is.numeric(draws),
    "`draws`: must have at least 2 rows"            = nrow(draws) >= 2L,
    "`draws`: must have at least 1 column"          = ncol(draws) >= 1L,
    "`prob`: must be a single number in (0, 1)"     = length(prob) == 1L && prob > 0 && prob < 1,
    "`est.FUN`: must be a function"                 = is.function(est.FUN)
  )
  
  S <- nrow(draws)
  
  # Empirical CDF per parameter:
  # "max", for lower tail; "min", for upper tail
  up <- matrixStats::colRanks(draws, ties.method = "max") / S
  lw <- matrixStats::colRanks(draws, ties.method = "min") / S
  
  # Worst-case tail probability for each draw across all parameters
  # pmin selects the nearer tail per draw per parameter, 
  # rowMins takes the worst parameter
  tp <- matrixStats::rowMins(pmin(up, 1 - lw))
  
  # (1 - prob)-quantile of tp; type = 1 for exact empirical CDF inversion
  cq <- quantile(tp, probs = 1 - prob, type = 1)
  
  # Simultaneous equitailed bounds: cut cp from each tail 
  out <- matrixStats::colQuantiles(draws, probs = c(cq, 1 - cq), type = 1)
  
  # output
  nms <- colnames(draws)
  list(
    lower   = setNames(out[, 1], nms),
    upper   = setNames(out[, 2], nms),
    prob    = prob,
    cq = cq,
    est     = setNames(apply(draws, 2, est.FUN), nms),
    est.FUN = est.FUN
  )
}

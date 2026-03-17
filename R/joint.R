#' Simultaneous credible interval from joint posterior draws
#'
#' Computes simultaneous credible intervals for a set of parameters.
#' Unlike marginal intervals, which guarantee \eqn{1 - \alpha} coverage for
#' each parameter individually, these bands guarantee that **all** parameters
#' stay within their respective bounds simultaneously with probability \eqn{prob}.
#'
#' @details
#' The procedure uses a rank-based approach to find a common quantile level
#' \eqn{W^*} more conservative than the marginal \eqn{(1 - prob)/2}. For each
#' draw \eqn{s}, the "worst-case" tail probability across all \eqn{K} parameters is:
#' \deqn{
#'   WC_s = \min_{k=1,\dots,K} \left(
#'     \min\!\left[\frac{\mathrm{rank}(x_{sk})}{S},\;
#'                 1 - \frac{\mathrm{rank}(x_{sk}) - 1}{S}\right]
#'   \right)
#' }
#' where ranks are computed within each parameter column across all \eqn{S} draws.
#' The simultaneous bands are the \eqn{W^*} and \eqn{1 - W^*} marginal quantiles,
#' where \eqn{W^*} is the \eqn{(1 - prob)}-quantile of the vector \eqn{WC}.
#'
#' @param draws Numeric matrix or data frame of posterior draws
#'   (rows = draws, columns = parameters).
#' @param prob Numeric scalar in \eqn{(0, 1)}. The joint coverage probability
#'   (default `0.95` for 95\% simultaneous intervals).
#' @param est.FUN Function used to compute point estimates for each parameter
#'   (default `median`).
#'
#' @return A named list with components:
#' \describe{
#'   \item{lower}{Numeric vector of lower bounds, one per parameter.}
#'   \item{upper}{Numeric vector of upper bounds, one per parameter.}
#'   \item{prob}{The requested coverage probability.}
#'   \item{WC}{The adjusted tail quantile \eqn{W^*}.}
#'   \item{est}{Point estimates computed via \code{est.FUN}.}
#'   \item{est.FUN}{The estimation function used.}
#' }
#'
#' @importFrom stats quantile
#' @importFrom matrixStats colRanks rowMins rowMaxs colQuantiles
#'
#' @export
#' @examples
#' mu    <- c(4, 0, -2)
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
#' draws <- MASS::mvrnorm(2000, mu, Sigma)
#' colnames(draws) <- c("theta1", "theta2", "theta3")
#'
#' joint(draws, prob = 0.95)
joint <- function(draws, prob = 0.95, est.FUN = median) {
  
  # --- Input validation ----
  if (is.data.frame(draws)) draws <- as.matrix(draws)
  
  stopifnot(
    "`draws`: must be a numeric matrix-like object" = is.numeric(draws),
    "`draws`: must have at least 2 rows (draws)" = nrow(draws) >= 2L,
    "`draws`: must have at least 1 column"  = ncol(draws) >= 1L,
    "`prob`: must be a single number in (0, 1)" = length(prob) == 1L && prob > 0 && prob < 1,
    "`est.FUN`: must be a function" = is.function(est.FUN)
  )
  
  S <- nrow(draws)
  K <- ncol(draws)
  
  # --- Joint interval ----
  
  # 1. Column-wise ranks (rank within each parameter across draws)
  rF <- matrixStats::colRanks(draws, ties.method = "max")
  
  # 2. Worst-case tail probability per draw
  min_rank_prop <- matrixStats::rowMins(rF) / S
  max_rank_prop <- 1 - (matrixStats::rowMaxs(rF) - 1L) / S
  W <- pmin(min_rank_prop, max_rank_prop)
  
  # 3. Empirical threshold W
  WC <- unname(quantile(W, probs = 1 - prob, type = 7))
  
  # 4. Marginal quantiles at adjusted level
  out <- matrixStats::colQuantiles(draws, probs = c(WC, 1 - WC))
  
  # ---- Output ----
  nms <- colnames(draws)
  lower <- setNames(out[, 1], nms)
  upper <- setNames(out[, 2], nms)
  est   <- setNames(apply(draws, 2, est.FUN), nms)
  
  list(
    lower = lower,
    upper = upper,
    prob  = prob,
    est = est,
    est.FUN = est.FUN
  )
}


#' Simultaneous credible bands from joint posterior draws
#'
#' Computes simultaneous (joint) credible bands for a set of parameters. 
#' Unlike marginal intervals, which guarantee \eqn{1-\alpha} coverage for 
#' each parameter individually, these bands guarantee that **all** parameters 
#' stay within their respective bounds simultaneously with \eqn{1-\alpha} probability.
#'
#' @details
#' The procedure utilizes a rank-based approach to find a common quantile 
#' level \eqn{W^*} that is more conservative than the marginal \eqn{\alpha/2}. 
#' For each draw \eqn{s}, we calculate the "extremeness" across all \eqn{K} 
#' parameters:
#' \deqn{
#'    WC_s = \min_{k=1 \dots K} \left( \min \left[ \frac{rank(x_{sk})}{S}, 1 - \frac{rank(x_{sk})-1}{S} \right] \right)
#' }
#' The simultaneous bands are then defined as the \eqn{W^*} and \eqn{1-W^*} 
#' marginal quantiles, where \eqn{W^*} is the \eqn{\alpha}-quantile of the 
#' vector \eqn{WC}.
#'
#' @param draws Numeric matrix or data frame of posterior draws (rows = draws, columns = parameters).
#' @param alpha Numeric scalar in \eqn{(0, 1)}. The joint error rate (default is `0.05` for 95% bands).
#'
#' @return A matrix with two columns, `lower` and `upper`, and one row per parameter.
#'
#' @importFrom stats quantile
#'
#' @export
#' @examples
#' # 3 correlated parameters
#' mu <- c(1, 2, 3)
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
#' draws <- MASS::mvrnorm(2000, mu, Sigma)
#' colnames(draws) <- c("theta1", "theta2", "theta3")
#'
#' joint(draws, alpha = 0.05)
joint <- function(draws, alpha = 0.05) {
  
  # --- Input validation -------------------------------------------------------
  if (is.data.frame(draws)) draws <- as.matrix(draws)
  
  stopifnot(
    "`draws`: must be a numeric matrix-like object" = is.numeric(draws),
    "`draws`: must have at least 2 rows (draws)"    = nrow(draws) >= 2L,
    "`alpha`: must be a single number in (0, 1)"    = length(alpha) == 1L && alpha > 0 && alpha < 1
  )
  
  # --- Efficiency Note: matrixStats is highly recommended for speed ---
  has_matrixStats <- requireNamespace("matrixStats", quietly = TRUE)
  
  # --- Rank-based simultaneous bands ------------------------------------------
  S <- nrow(draws)
  
  # 1. Compute column-wise ranks
  if (has_matrixStats) {
    rF <- matrixStats::colRanks(draws, ties.method = "max")
  } else {
    rF <- apply(draws, 2, rank, ties.method = "max")
  }
  
  # 2. Compute "Worst-Case" tail probability per draw (row)
  # Logic: Find how close each draw is to the edge in ANY dimension
  if (has_matrixStats) {
    min_rank_prop <- matrixStats::rowMins(rF) / S
    max_rank_prop <- 1 - (matrixStats::rowMaxs(rF) - 1L) / S
  } else {
    min_rank_prop <- apply(rF, 1, min) / S
    max_rank_prop <- 1 - (apply(rF, 1, max) - 1L) / S
  }
  
  WC <- pmin(min_rank_prop, max_rank_prop)
  
  # 3. Find the empirical threshold W*
  Wstar <- unname(quantile(WC, probs = alpha, type = 7))
  
  # If alpha is smaller than the probability of a single draw, 
  # we don't have enough data to support that confidence level.
  if (alpha < (1 / S)) {
    warning(sprintf(
      "Alpha (%.4f) is smaller than the data resolution (1/S = %.4f). Results are unstable and default to the most extreme observed draws.",
      alpha, 1/S
    ))
  }

  # 4. Extract marginal quantiles at the adjusted level
  if (has_matrixStats) {
    res <- t(matrixStats::colQuantiles(draws, probs = c(Wstar, 1 - Wstar)))
  } else {
    res <- apply(draws, 2, quantile, probs = c(Wstar, 1 - Wstar))
  }
  
  # Ensure clean output format
  out <- t(res)
  colnames(out) <- c("lower", "upper")
  return(out)
}

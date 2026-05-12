#' Simultaneous credible intervals from joint posterior draws
#' Computes simultaneous equitailed credible intervals for a set of parameters,
#' guaranteeing that **all** parameters stay within their bounds simultaneously
#' with a user specified probability.
#'
#' @details
#' Simultaneous coverage is calibrated by examining how extreme each draw is
#' across all parameters jointly. For each draw, the minimum tail probability
#' across all parameters is computed. The simultaneous threshold
#' is the \eqn{\alpha}-quantile of these minima: only
#' \eqn{\alpha} of draws are more extreme than this value in at least one
#' parameter simultaneously.
#'
#' A closely related implementation is \code{\link[credsubs]{sim.cred.band}}
#' in the \pkg{credsubs} package (Schnell et al., 2020).
#'
#' @param draws Numeric matrix or data frame of posterior draws
#'   (rows = draws, columns = parameters).
#' @param prob Numeric scalar in \eqn{(0, 1)}. Joint coverage probability
#'   (default \code{0.95}).
#' @param est.FUN Function for point estimates (default \code{median}).
#' @param null.value Optional numeric vector of null values to test against
#'   the simultaneous credible intervals, one per parameter. If \code{NULL}
#'   (default), no test is performed. When provided, a logical
#'   column \code{test} is appended: \code{TRUE} means the null value is
#'   **excluded** from the interval.
#'
#' @return A dataframe with:
#' \describe{
#'   \item{lower}{Lower bounds, one per parameter.}
#'   \item{est}{Point estimates via \code{est.FUN}.}
#'   \item{upper}{Upper bounds, one per parameter.}
#'   \item{prob}{Requested joint coverage probability.}
#'   \item{cq}{Critical value.}
#'   \item{sig}{Logical; \code{TRUE} if \code{null.value} is excluded from
#'     the interval. Only present when \code{null.value} is supplied.}
#' }
#'
#' @note
#' If you use this function in published research, please cite:
#' \itemize{
#'   \item The package: \url{https://github.com/mar-cald/multibayes}.
#' }
#'
#' @importFrom stats quantile median setNames
#' @importFrom matrixStats colRanks rowMins colQuantiles
#' @keywords internal
#' @examples
#' mu <- c(4, 0, -2)
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
#' draws <- MASS::mvrnorm(2000, mu, Sigma)
#' colnames(draws) <- c("theta1", "theta2", "theta3")
#'
#' joint(draws, prob = 0.95)
#' joint(draws, prob = 0.95, null.value = c(0, 0, 0))
#' @export 
joint <- function(draws, prob = 0.95, est.FUN = median, null.value = NULL) {
  
  if (is.data.frame(draws)) draws <- as.matrix(draws)
  
  stopifnot(
    "`draws`: must be a numeric matrix-like object" = is.numeric(draws),
    "`draws`: must have at least 1000 rows"         = nrow(draws) >= 1000L,
    "`draws`: must have at least 2 column"          = ncol(draws) >= 2L,
    "`prob`: must be a single number in (0, 1)"     = length(prob) == 1L && prob > 0 && prob < 1,
    "`est.FUN`: must be a function"                 = is.function(est.FUN)
  )
  
  if (!is.null(null.value)) {
    stopifnot(
      "`null.value`: must be a numeric vector"             = is.numeric(null.value),
      "`null.value`: length must match number of parameters" =
        length(null.value) == ncol(draws)
    )
  }
  
  S <- nrow(draws)
  
  # Conservative bounds: upper max and lower min.
  up <- t(matrixStats::colRanks(draws, ties.method = "max") / S)
  lw <- t(matrixStats::colRanks(draws, ties.method = "min") / S)
  
  # Worst-case tail probability for each draw across all parameters
  tp <- matrixStats::rowMins(pmin(up, 1 - lw))
  
  # Critical quantile
  cq <- quantile(tp, probs = 1 - prob, type = 1)
  
  # Simultaneous equitailed bounds
  out <- matrixStats::colQuantiles(draws, probs = c(cq, 1 - cq), type = 1)
  
  nms <- colnames(draws)
  result <- data.frame(
    lower = setNames(out[, 1], nms),
    est   = setNames(apply(draws, 2, est.FUN), nms),
    upper = setNames(out[, 2], nms),
    prob  = prob,
    cq    = cq
  )
  
  if (!is.null(null.value)) {
    result$test <- null.value < result$lower | null.value > result$upper
  }
  
  result
}


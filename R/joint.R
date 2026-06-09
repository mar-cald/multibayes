#' Family-wise summaries from joint posterior draws
#'
#' Two complementary family-wise statements across a set of parameters: the
#' \strong{cumulative joint statement} (default), the posterior probability that
#' the strongest directional claims hold simultaneously, and \strong{simultaneous
#' credible intervals} (\code{interval = TRUE}). Whereas \code{\link{pd.adjust}}
#' targets a per-test (directional false-discovery) criterion, \code{joint}
#' provides the family-wise counterparts.
#'
#' @details
#' \strong{Cumulative joint statement} (\code{interval = FALSE}, the default).
#' Each parameter's claimed direction is its dominant posterior side relative to
#' \code{null.value} (or the side set by \code{direction}); parameters are
#' ordered by decreasing \emph{pd}, and \code{joint_cum} is the cumulative joint
#' posterior probability that the strongest \eqn{k} claims hold simultaneously.
#' It is computed directly from the draws by intersecting the per-draw
#' directional events, so it reflects the posterior dependence among parameters;
#' \eqn{1 - } \code{joint_cum} is the posterior probability that at least one of
#' those \eqn{k} claims is mis-signed (a family-wise, Type S, statement).
#'
#' \strong{Simultaneous credible intervals} (\code{interval = TRUE}). Equitailed
#' intervals calibrated so that \strong{all} parameters lie within their bounds
#' simultaneously with probability \code{prob}, via the quantile-based
#' simultaneous credible band (Besag et al., 1995; Held, 2004). Coverage is set
#' by examining how extreme each draw is across all parameters jointly: for every
#' draw the minimum (nearer) tail probability across parameters is taken, and the
#' threshold is the \eqn{(1 - \texttt{prob})}-quantile of these minima. An effect
#' can be declared when its band excludes the null value. A closely related
#' implementation is \code{sim.cred.band} in the \pkg{credsubs} package
#' (Schnell et al., 2020).
#'
#' @param draws Numeric matrix or data frame of posterior draws
#'   (rows = draws, columns = parameters). At least 1000 rows and 2 columns.
#' @param interval Logical (default \code{FALSE}). If \code{FALSE}, return the
#'   cumulative joint statement; if \code{TRUE}, return simultaneous credible
#'   intervals.
#' @param prob Numeric scalar in \eqn{(0, 1)}. Joint coverage probability for the
#'   simultaneous intervals (default \code{0.95}). Used only when
#'   \code{interval = TRUE}.
#' @param null.value Numeric scalar or vector giving the null reference value(s)
#'   used to define the claimed direction (a scalar is recycled). Used only when
#'   \code{interval = FALSE}. Defaults to \code{0}.
#' @param direction Character vector of \code{"greater"}, \code{"less"}, or
#'   \code{"two.sided"} (or \code{NULL}), the testing mode per parameter. A scalar
#'   is recycled; \code{NULL} (default) uses \code{"two.sided"} (dominant side)
#'   for all. Used only when \code{interval = FALSE}.
#' @param est.FUN Function returning a point estimate from a numeric vector
#'   (default \code{\link[stats]{median}}). Used only when \code{interval = TRUE}.
#'
#' @return A \code{data.frame}. For \code{interval = FALSE}, one row per
#'   parameter ordered by decreasing \emph{pd}, with columns \code{median.est},
#'   \code{null.value}, \code{pd}, \code{direction}, and \code{joint_cum}
#'   (cumulative joint probability that the first \eqn{k} claims hold). For
#'   \code{interval = TRUE}, one row per parameter with columns \code{lower},
#'   \code{est}, \code{upper}, and \code{prob}.
#'
#' @references
#' Box, G. E. P., & Tiao, G. C. (1992). Bayesian inference in statistical
#' analysis (1st ed.). \emph{Wiley}.
#'
#' Besag, J., Green, P., Higdon, D., & Mengersen, K. (1995). Bayesian
#' computation and stochastic systems. \emph{Statistical Science}, 10(1), 3--41.
#'
#' Held, L. (2004). Simultaneous posterior probability statements from Monte
#' Carlo output. \emph{Journal of Computational and Graphical Statistics},
#' 13(1), 20--35.
#'
#' Schnell, P. et al. (2020). \pkg{credsubs}: Credible Subsets. R package.
#' \url{https://CRAN.R-project.org/package=credsubs}.
#'
#' @seealso \code{\link{pd.adjust}} for the per-test directional adjustment.
#'
#' @importFrom stats quantile median setNames
#' @importFrom matrixStats colRanks rowMins colQuantiles
#' @export
#' @examples
#' mu <- c(4, 0, -2)
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
#' draws <- MASS::mvrnorm(2000, mu, Sigma)
#' colnames(draws) <- c("theta1", "theta2", "theta3")
#'
#' joint(draws)                    # cumulative joint statement (default)
#' joint(draws, interval = TRUE)   # simultaneous credible intervals
joint <- function(draws, interval = FALSE, prob = 0.95,
                  null.value = 0, direction = NULL, est.FUN = median) {

  if (is.data.frame(draws)) draws <- as.matrix(draws)
  stopifnot(
    "`draws`: must be a numeric matrix-like object" = is.numeric(draws),
    "`draws`: must have at least 1000 rows" = nrow(draws) >= 1000L,
    "`draws`: must have at least 2 columns" = ncol(draws) >= 2L,
    "`interval`: must be TRUE or FALSE" = length(interval) == 1L && is.logical(interval)
  )
  nms <- colnames(draws)

  ## ---- Simultaneous credible intervals -------------------------------------
  if (isTRUE(interval)) {
    stopifnot(
      "`prob`: must be a single number in (0, 1)" = length(prob) == 1L && prob > 0 && prob < 1,
      "`est.FUN`: must be a function" = is.function(est.FUN)
    )
    S  <- nrow(draws)
    # Conservative empirical CDF per draw per parameter (max/min ties coincide
    # for continuous draws); nearer tail per draw, most extreme across params.
    up <- t(matrixStats::colRanks(draws, ties.method = "max") / S)
    lw <- t(matrixStats::colRanks(draws, ties.method = "min") / S)
    tp <- matrixStats::rowMins(pmin(up, 1 - lw))
    cq <- quantile(tp, probs = 1 - prob, type = 1)
    out <- matrixStats::colQuantiles(draws, probs = c(cq, 1 - cq), type = 1)
    return(data.frame(
      lower = setNames(out[, 1], nms),
      est   = setNames(apply(draws, 2, est.FUN), nms),
      upper = setNames(out[, 2], nms),
      prob  = prob
    ))
  }

  ## ---- Cumulative joint statement ------------------------------------------
  nc <- ncol(draws)
  if (length(null.value) == 1L) null.value <- rep(null.value, nc)
  if (is.null(direction))       direction  <- rep("two.sided", nc)
  if (length(direction) == 1L)  direction  <- rep(direction, nc)
  if (length(null.value) != nc) stop("`null.value` must be a scalar or a vector of length `ncol(draws)`.")
  if (length(direction)  != nc) stop("`direction` must be a scalar or a vector of length `ncol(draws)`.")
  if (!all(direction %in% c("two.sided", "less", "greater")))
    stop("`direction` must contain only two.sided, less, greater.")

  centered <- sweep(draws, 2, null.value, "-")

  # Per-draw indicator that each parameter lies on its claimed side. For
  # two.sided tests the claimed side is the dominant marginal direction.
  side_ind <- mapply(function(j, d) {
    if (d == "greater")   centered[, j] > 0
    else if (d == "less") centered[, j] < 0
    else if (mean(centered[, j] > 0) >= mean(centered[, j] < 0)) centered[, j] > 0
    else                  centered[, j] < 0
  }, seq_len(nc), direction)
  side_ind <- matrix(side_ind, nrow = nrow(draws), ncol = nc)

  pd  <- colMeans(side_ind)
  ord <- order(pd, decreasing = TRUE)         # cumulative only meaningful strongest-first

  joint_cum <- numeric(nc)
  running   <- rep(TRUE, nrow(side_ind))
  for (k in seq_len(nc)) {
    running      <- running & side_ind[, ord[k]]
    joint_cum[k] <- mean(running)
  }

  out <- data.frame(
    median.est = round(apply(draws, 2, median), 4),
    null.value = null.value,
    pd         = round(pd, 4),
    direction  = direction
  )
  if (!is.null(nms) && length(nms) == nc) rownames(out) <- nms
  out <- out[ord, , drop = FALSE]
  out$joint_cum <- round(joint_cum, 4)        # in ord order, aligned with reordered rows
  out
}

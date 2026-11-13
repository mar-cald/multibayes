#' Plot method for \code{pd_adjust} objects
#'
#' Shows, for each test, how the posterior probability of the favored (or
#' predicted) direction changes with the prior. Each bar runs from
#' \code{pd.adj}, the value when that direction has prior probability
#' \eqn{1 - \bar{\pi}}, to the value when it has prior probability
#' \eqn{\bar{\pi}}, which is
#' \eqn{pd\,\bar{\pi} / (pd\,\bar{\pi} + (1 - pd)(1 - \bar{\pi}))}. The point
#' marks \code{pd}, and the two ends of each bar are printed on the right.
#'
#' A direction is claimed when the whole bar lies at or above the threshold
#' \eqn{c} (dashed line), that is, when \code{pd.adj} is at least \eqn{c}. A
#' bar that crosses 0.5 means that the more probable direction depends on the
#' prior.
#'
#' @param x A \code{pd_adjust} object returned by \code{\link{pd.adjust}}.
#' @param labels Labels of the tests, in the order of the rows of \code{x}.
#'   Defaults to the row names.
#' @param xlim Range of the x-axis. By default it covers all bars and \eqn{c}.
#' @param digits Number of digits of the values printed on the right.
#'   Defaults to 3.
#' @param ... Further arguments passed to \code{plot.default}.
#'
#' @return Invisibly, a data frame with one row per test and the columns
#'   \code{lower} (equal to \code{pd.adj}), \code{pd}, \code{upper}, and
#'   \code{claim}.
#'
#' @examples
#' set.seed(1)
#' draws <- MASS::mvrnorm(4000, mu = c(0.4, 0.15, 0.05),
#'                        Sigma = diag(c(0.12, 0.12, 0.12)^2))
#' colnames(draws) <- c("b1", "b2", "b3")
#' res <- pd.adjust(draws = draws, pi.max = 0.8)
#' plot(res)
#'
#' @export
plot.pd_adjust <- function(x, labels = rownames(x), xlim = NULL, digits = 3, ...) {
  df <- as.data.frame(x)
  pd <- df$pd
  p  <- df$pi.max
  cc <- attr(x, "threshold")
  k  <- length(pd)
  if (is.null(labels)) labels <- paste0("H", seq_len(k))

  # posterior probability of the favored direction under the two extreme priors
  lower <- pd * (1 - p) / (pd * (1 - p) + (1 - pd) * p)
  upper <- pd * p / (pd * p + (1 - pd) * (1 - p))
  claim <- lower >= cc

  if (is.null(xlim)) xlim <- c(min(lower, cc), 1)

  # printed ranges, without leading zeros
  # values close to 1 get the digits they need, so .9996 is not shown as 1.000
  fmt <- function(v) {
    d <- pmax(digits, ceiling(-log10(pmax(1 - v, 1e-6))))
    sub("^0\\.", ".", formatC(v, format = "f", digits = d))
  }
  rng <- paste0("[", vapply(lower, fmt, ""), ", ", vapply(upper, fmt, ""), "]")

  # margins in lines: labels on the left, ranges on the right, axis title below
  op <- graphics::par(mar = c(5.5, 1, 2, 1))
  on.exit(graphics::par(op))
  csi <- graphics::par("csi")
  left  <- max(graphics::strwidth(labels, units = "inches", cex = 0.9)) / csi + 1.5
  right <- max(graphics::strwidth(rng, units = "inches", cex = 0.9)) / csi + 1.5
  graphics::par(mar = c(5.5, left, 2, right))

  y <- rev(seq_len(k))  # first hypothesis at the top, as in the printed table
  graphics::plot(NA, xlim = xlim, ylim = c(0.5, k + 0.5), yaxt = "n",
                 xlab = "", ylab = "", ...)
  graphics::title(xlab = "Posterior probability of the favored direction across the class",
                  line = 2.8)
  graphics::axis(2, at = y, labels = labels, las = 1, cex.axis = 0.9, tick = FALSE)
  graphics::axis(4, at = y, labels = rng, las = 1, cex.axis = 0.9, tick = FALSE)
  graphics::abline(v = cc, lty = 2)
  graphics::mtext(bquote(italic(c) == .(fmt(cc))), side = 3, at = cc, line = 0.1, cex = 0.8)
  graphics::segments(lower, y, upper, y, lwd = 6, lend = 1,
                     col = ifelse(claim, "black", "grey70"))
  graphics::points(pd, y, pch = 21, bg = "white", cex = 0.9, lwd = 1.2)

  invisible(data.frame(lower = lower, pd = pd, upper = upper, claim = claim,
                       row.names = labels))
}

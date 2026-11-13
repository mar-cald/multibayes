#' Simultaneous credible intervals and joint directional statements
#'
#' Two ways of making directional claims about several parameters from their
#' joint posterior. By default, \code{joint()} returns simultaneous credible
#' intervals. With \code{interval = FALSE} it returns, for each \eqn{k}, the
#' joint posterior probability that the \eqn{k} strongest directional claims
#' all hold.
#'
#' @details
#' \strong{Simultaneous credible intervals} (\code{interval = TRUE}).
#' Equal-tailed intervals such that all parameters lie inside their interval
#' at the same time with posterior probability \code{prob}. For every draw, the
#' function takes the smallest tail probability across parameters, and the
#' intervals are set at the \eqn{(1 - \texttt{prob})}-quantile of these values.
#' A direction can be claimed for a parameter when its interval excludes the
#' null value.
#'
#' The ends of the intervals lie further in the tails as the number of
#' parameters grows, so more draws are needed. Each end rests on about
#' \eqn{S (1 - \texttt{prob}^{1/m}) / 2} draws, where \eqn{S} is the number of
#' draws and \eqn{m} the number of parameters. With fewer than 10, the
#' intervals are too narrow on average and \code{joint()} gives a warning. With
#' \code{prob = 0.95}, this means at least about 1,600 draws for 4 parameters,
#' 7,800 for 20, and 15,600 for 40. \code{sim.cred.band} in the \pkg{credsubs} package (Schnell et
#' al., 2020) is a closely related implementation.
#'
#' \strong{Joint directional statement} (\code{interval = FALSE}). The
#' direction of each parameter is the side of \code{null.value} favored by its
#' posterior, or the side given in \code{direction}. Parameters are sorted by
#' decreasing \emph{pd}, and \code{joint_cum} gives, for each \eqn{k}, the
#' proportion of draws in which the first \eqn{k} directions all hold. It
#' therefore accounts for the posterior dependence among parameters. This is a
#' probability about the signs of the parameters under the prior of the fitted
#' model, with no adjustment of the prior odds. In families with many null
#' effects it leads to many false claims.
#'
#' @param draws Matrix or data frame of posterior draws, one column per
#'   parameter, with at least 1000 rows and 2 columns.
#' @param interval Logical. If \code{TRUE} (the default), simultaneous
#'   credible intervals are returned. If \code{FALSE}, the joint directional
#'   statement.
#' @param prob Number in \eqn{(0, 1)}. Joint posterior probability of the
#'   intervals. Defaults to \code{0.95}. Used only when \code{interval = TRUE}.
#' @param null.value Number or vector, the null value of each parameter. A
#'   single number applies to all parameters. Defaults to \code{0}. Used only
#'   when \code{interval = FALSE}.
#' @param direction Character vector with elements \code{"two.sided"},
#'   \code{"greater"}, or \code{"less"}, one per parameter, or a single value
#'   for all. \code{NULL} (the default) means \code{"two.sided"} for all, that
#'   is, the side favored by the posterior. Used only when
#'   \code{interval = FALSE}.
#' @param est.FUN Function that computes a point estimate from the draws of
#'   one parameter. Defaults to \code{\link[stats]{median}}. Used only when
#'   \code{interval = TRUE}.
#'
#' @return A data frame with one row per parameter and a print method.
#'   With \code{interval = TRUE}, an object of class \code{joint_interval}
#'   with the columns \code{lower}, \code{est}, \code{upper}, and \code{prob}.
#'   With \code{interval = FALSE}, an object of class \code{joint_directional},
#'   sorted by decreasing \emph{pd}, with the columns \code{median.est},
#'   \code{null.value}, \code{pd}, \code{direction}, and \code{joint_cum}.
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
#' Gelman, A., Tuerlinckx, F. (2000). Type S error rates for classical 
#' and Bayesian single and multiple comparison procedures. 
#' \emph{Computational Statistics} 15, 373–390 (2000). https://doi.org/10.1007/s001800000040
#'
#' Schnell, P. et al. (2020). \pkg{credsubs}: Credible Subsets. R package.
#' \url{https://CRAN.R-project.org/package=credsubs}.
#'
#' @seealso \code{\link{pd.adjust}} for the per-test prior-odds threshold.
#'
#' @importFrom stats quantile median setNames
#' @importFrom matrixStats colRanks rowMins colQuantiles
#' @export
#' @examples
#' mu <- c(1.5, 0, -2)
#' Sigma <- matrix(c(1, 0.8, 0.5, 0.8, 1, 0.3, 0.5, 0.3, 1), 3, 3)
#' draws <- MASS::mvrnorm(2000, mu, Sigma)
#' colnames(draws) <- c("theta1", "theta2", "theta3")
#'
#' joint(draws)                     # simultaneous credible intervals (default)
#' joint(draws, interval = FALSE)   # joint directional statement
joint <- function(draws, interval = TRUE, prob = 0.95,
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
    # Each interval end sits in a tail of probability about (1 - prob^(1/m)) / 2,
    # which shrinks as the number of parameters m grows. With fewer than about
    # 10 draws in that tail the ends are noisy and the joint probability of the
    # intervals falls below prob.
    m_par <- ncol(draws)
    tail_p <- (1 - prob^(1 / m_par)) / 2
    if (S * tail_p < 10) {
      warning(sprintf(paste0("With %d parameters, each interval end rests on about %.0f draws. ",
                             "Use at least %d draws for intervals with joint probability %s."),
                      m_par, S * tail_p, ceiling(10 / tail_p), format(prob)), call. = FALSE)
    }
    # Empirical CDF per draw per parameter; 
    # nearer tail per draw, most extreme across params.
    up <- t(matrixStats::colRanks(draws, ties.method = "max") / S)
    lw <- t(matrixStats::colRanks(draws, ties.method = "min") / S)
    tp <- matrixStats::rowMins(pmin(up, 1 - lw))
    cq <- quantile(tp, probs = 1 - prob, type = 1)
    out <- matrixStats::colQuantiles(draws, probs = c(cq, 1 - cq), type = 1)
    res <- data.frame(
      lower = setNames(out[, 1], nms),
      est   = setNames(apply(draws, 2, est.FUN), nms),
      upper = setNames(out[, 2], nms),
      prob  = prob
    )
    class(res) <- c("joint_interval", "data.frame")
    return(res)
  }

  ## ---- joint directional statement ------------------------------------------
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
  ord <- order(pd, decreasing = TRUE) # cumulative only meaningful strongest-first

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
  out$joint_cum <- round(joint_cum, 4)  # in ord order, aligned with reordered rows
  class(out) <- c("joint_directional", "data.frame")
  out
}

#' Print method for joint directional statements
#'
#' Prints the table returned by \code{\link{joint}} with \code{interval = FALSE}.
#' Columns with the same value in every row (\code{null.value},
#' \code{direction}) are shown once, above the table.
#'
#' @param x A \code{joint_directional} object returned by \code{\link{joint}}.
#' @param digits Number of digits used when printing. Defaults to 4.
#' @param ... Further arguments passed to \code{print.data.frame}.
#'
#' @return \code{x}, invisibly.
#'
#' @keywords internal
#' @export
print.joint_directional <- function(x, digits = 4, ...) {
  df <- as.data.frame(x)
  cat("Joint directional statement\n")
  ch <- const_cols(df, c("null.value", "direction"))
  if (!is.null(ch$text)) cat("  ", ch$text, "\n", sep = "")
  cat("  joint_cum: Pr(the top k claims all hold | y)\n\n")
  print(df[, setdiff(names(df), ch$drop), drop = FALSE], digits = digits, ...)
  invisible(x)
}


#' Print method for simultaneous credible intervals
#'
#' Prints the table returned by \code{\link{joint}} with \code{interval = TRUE}.
#' The joint probability is shown once, above the table.
#'
#' @param x A \code{joint_interval} object returned by \code{\link{joint}}.
#' @param digits Number of decimal places used when printing. Defaults to 4.
#' @param ... Further arguments passed to \code{print.data.frame}.
#'
#' @return \code{x}, invisibly.
#'
#' @keywords internal
#' @export
print.joint_interval <- function(x, digits = 4, ...) {
  df <- as.data.frame(x)
  cat(sprintf("Simultaneous credible intervals (joint probability %.4g)\n\n", df$prob[1]))
  show <- df[, setdiff(names(df), "prob"), drop = FALSE]
  show[] <- lapply(show, round, digits = digits)
  print(show, ...)
  invisible(x)
}

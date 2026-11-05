#' Prior-odds threshold for the Probability of Direction \emph{pd}
#'
#' Converts an upper bound, set by the researcher, on the prior probability
#' of either sign of a parameter into the threshold that \emph{pd}
#' must reach before a direction is claimed. The function
#' accepts either a vector of pre-computed \emph{pd} values or a matrix of
#' posterior draws, from which \emph{pd} values are computed internally. Both
#' direction-agnostic and directional tests are supported: the
#' \code{direction} argument controls how \emph{pd} is computed per
#' hypothesis. The threshold is governed by \eqn{\bar{\pi}} (argument
#' \code{pi.max}), an upper bound, set by the researcher before the data, on
#' the prior probability of each direction of a parameter (either side of
#' \code{null.value}). Equivalently, one may supply \eqn{\Pr(H_0)} (argument
#' \code{p0}), the largest prior probability that no claim in the family is
#' correct, from which \eqn{\bar{\pi} = \Pr(H_0)^{1/m}} is recovered if the
#' parameters are independent a priori. For each hypothesis the function also
#' reports \code{pd.adj}, the posterior probability of the favored (or
#' predicted) direction under the prior in the class most skeptical of it;
#' under the same prior the other side has probability \code{1 - pd.adj}, so
#' the two sum to one. A direction is claimed when \code{pd.adj} is at least
#' \eqn{c}, which is the same as \eqn{pd \ge c^*}.
#'
#' @details
#' Following Jeffreys (1938) and Westfall et al. (1997), the prior odds of each
#' directional claim in the family are penalized. With \code{pi.max} the penalty
#' is the same for any number of tests; with \code{p0},
#' \eqn{\bar{\pi} = \Pr(H_0)^{1/m}} and the penalty grows with the number of
#' tests. Because the
#' direction of a claim is chosen by the data, the penalty is applied over a
#' class of priors: the estimation prior (the prior of the fitted model) is
#' reweighted on the two sides of
#' \code{null.value}, so that the prior probability of each direction lies
#' between \eqn{1-\bar{\pi}} and \eqn{\bar{\pi}}, and a claim is made only if its
#' posterior probability is at least \eqn{c} under every prior in the class that
#' reweights the prior of that parameter. This is equivalent to
#' \deqn{pd \ge c^*, \qquad c^* = \frac{c\,\bar{\pi}}{(1-c)(1-\bar{\pi}) + c\,\bar{\pi}},}
#' and the function reports \eqn{c^*} (\code{threshold.adj}). The smallest
#' posterior probability is reached under the prior that gives the favored
#' direction probability \eqn{1 - \bar{\pi}}; under that prior the posterior
#' probability of the favored direction is
#' \deqn{pd_{adj} = \frac{pd\,(1-\bar{\pi})}{pd\,(1-\bar{\pi}) + (1-pd)\,\bar{\pi}},}
#' reported as \code{pd.adj}, and \eqn{pd_{adj} \ge c} if and only if
#' \eqn{pd \ge c^*}. The posterior used for estimation is not modified. The
#' class is defined for an estimation prior that is symmetric about \code{null.value}, as are priors
#' centred on it and, as a limit, flat priors.
#'
#' \strong{Direction-agnostic tests} (\code{"two.sided"}): \emph{pd} is
#' \eqn{\max\!\big(\Pr(\theta > \theta_\text{null} \mid y),\,
#' \Pr(\theta < \theta_\text{null} \mid y)\big)} on \eqn{[0.5, 1]}, and a
#' claim is made in the favoured direction when \eqn{pd \ge c^*}.
#'
#' \strong{Directional tests} (\code{"greater"} or \code{"less"}): \emph{pd}
#' is the posterior mass on the predicted side, on \eqn{[0, 1]}, and the
#' predicted direction is claimed when \eqn{pd \ge c^*}.
#'
#' Mixed use of directional and direction-agnostic tests within the same call
#' is supported: each element of \code{direction} is handled independently,
#' and the same threshold \eqn{c^*} applies to all hypotheses regardless of
#' their directionality.
#'
#' @param pd Numeric vector of direction-agnostic \emph{pd} values in
#'   \eqn{[0.5, 1]}. With \code{pd} input every test is treated as
#'   direction-agnostic and \code{direction} is ignored; for directional tests,
#'   supply \code{draws}. Ignored if \code{draws} is supplied.
#' @param draws Optional matrix or data frame of posterior draws (columns = parameters).
#'   If provided, \emph{pd} values are computed automatically from the draws
#'   according to \code{direction} and \code{null.value}.
#' @param p0 Numeric scalar in \eqn{(0, 1)}. The largest prior probability
#'   that \strong{no} directional claim in the family is correct
#'   (\eqn{\Pr(H_0)}). Supply \strong{either} \code{pi.max} \strong{or} \code{p0},
#'   not both. When \code{p0} is given, \eqn{\bar{\pi}} is derived, for parameters
#'   independent a priori, as \eqn{\bar{\pi} = \Pr(H_0)^{1/m}}.
#' @param pi.max Numeric scalar in \eqn{(0, 1)}. An upper bound, set by the
#'   researcher before the data, on the prior probability of each direction of
#'   a parameter (either side of \code{null.value}); a claim is made only if it
#'   holds under the prior most skeptical of it. When the same \code{pi.max} is
#'   applied to exchangeable parameters, it is also an upper bound on the
#'   expected proportion of parameters in the family that are not on a given
#'   side of \code{null.value}. Since null effects are among them, a judgment
#'   about how many null effects the family could contain can guide its
#'   choice. Values below 0.5 are set to 0.5, with a warning, so that no
#'   penalty is applied. Supply \strong{either} \code{pi.max} \strong{or}
#'   \code{p0}, not both.
#' @param null.value Numeric scalar or vector. The null (reference) value against
#'  which the posterior is evaluated, specified on the scale of the posterior.
#'  A single scalar applies the same null to all parameters; a vector of length \code{ncol(draws)}
#'  assigns a distinct null to each parameter. Ignored when \code{pd} is supplied directly.
#'  Defaults to \code{0}. The class of priors behind the threshold is defined
#'  for an estimation prior symmetric about the null value (see Details), which
#'  holds for flat priors and for priors centred on it.
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
#' @param order Logical. If \code{TRUE}, the returned rows are ordered by
#'   decreasing \code{pd} (strongest directional evidence first). If \code{FALSE}
#'   (the default), rows are returned in input order, so the output stays aligned
#'   with the columns of \code{draws} (or the elements of \code{pd}).
#' @param threshold Numeric scalar in \eqn{(0.5, 1)}. The threshold \eqn{c}
#'   on \emph{pd} without penalty (default \code{0.975}, a two-sided per-test
#'   \eqn{\alpha = 0.05}). The function reports the threshold implied by
#'   \eqn{c} and \code{pi.max},
#'   \eqn{c^* = c\,\bar{\pi} / ((1-\bar{\pi})(1-c) + c\,\bar{\pi})}.
#' @param alpha Optional numeric scalar in \eqn{(0, 1)}. A convenience
#'   parameterization of \code{threshold} on the error scale: when supplied, the
#'   two-sided cutoff is set to \eqn{c = 1 - \alpha/2} (so \code{alpha = 0.05}
#'   gives \code{threshold = 0.975}), \strong{overriding} any \code{threshold}
#'   value. 
#'
#' @return An object of class \code{pd_adjust} (a \code{data.frame} with a
#'   dedicated \code{\link{print.pd_adjust}} method), with one row per
#'   hypothesis. Columns: \code{pd} (the probability of direction, to be
#'   compared with \code{threshold.adj}), \code{pd.adj} (the posterior
#'   probability of the favored direction under the prior in the class most
#'   skeptical of it, to be compared with \code{threshold}), \code{pi.max},
#'   \code{m}, and \code{direction}. When
#'   \code{order = TRUE} the rows are sorted by decreasing \code{pd}. When
#'   \code{draws} are supplied, the output additionally includes
#'   \code{median.est} (posterior median per parameter) and \code{null.value}.
#'   The \code{print} method shows \eqn{m}, \eqn{\bar{\pi}}, the thresholds
#'   \eqn{c} and \eqn{c^*}, and the columns that are constant across rows
#'   (\code{null.value}, \code{direction}) in a header, followed by the per-test table; all columns remain
#'   available for programmatic access.
#'   The object also carries, as attributes, the threshold \code{threshold}
#'   (\eqn{c}), the adjusted threshold \code{threshold.adj} (\eqn{c^*}), the
#'   corresponding two-sided per-test Type I error rates \code{alpha.nominal}
#'   (\eqn{2(1-c)}) and \code{alpha.eff} (\eqn{2(1-c^*)}); retrieve them with, e.g.,
#'   \code{attr(x, "threshold.adj")}.
#'   The threshold \code{threshold.adj} is an exact function of \eqn{c} and
#'   \eqn{\bar{\pi}}, but the reported per-test \emph{rate} \code{alpha.eff} uses the identity
#'   \eqn{\alpha = 2(1 - c)}, which holds **only when the null \emph{pd} is approximately uniform on
#'   \eqn{[0.5, 1]}**.
#'
#' @references
#' Berger, J. O. (1994). An overview of robust Bayesian analysis. \emph{Test},
#' 3(1), 5–124. https://doi.org/10.1007/BF02562676
#'
#' Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A Bayesian perspective on the Bonferroni adjustment.
#' \emph{Biometrika}, 84(2), 419–427. https://doi.org/10.1093/biomet/84.2.419
#' 
#' Jeffreys, H. (1938). Significance tests when several degrees of freedom arise simultaneously. \emph{Proceedings
#' of the Royal Society of London. Series A. Mathematical and Physical Sciences}, 165(921), 161–198.
#' https://doi.org/10.1098/rspa.1938.0052
#' 
#' Storey, J. D. (2003). The positive false discovery rate: A Bayesian interpretation and the q-value. \emph{The Annals
#' of Statistics}, 31(6), 2013–2035. https://doi.org/10.1214/aos/1074290335
#'
#' Scott, J. G., & Berger, J. O. (2006). An exploration of aspects of Bayesian multiple testing.
#' \emph{Journal of Statistical Planning and Inference}, 136(7), 2144–2162.
#' https://doi.org/10.1016/j.jspi.2005.08.031
#'
#' Scott, J. G., & Berger, J. O. (2010). Bayes and empirical-Bayes multiplicity adjustment in the variable-
#' selection problem. \emph{The Annals of Statistics}, 38(5), 2587–2619. https://doi.org/10.1214/10-AOS792
#' 
#' @examples
#' # From a vector of pd values
#' pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763,
#' H5 = 0.891, H6 = 0.987)
#' pd.adjust(pd = pd_values, p0 = 0.4)
#'
#' # Equivalent call specifying pi.max directly
#' pd.adjust(pd = pd_values, pi.max = 0.4^(1/6))
#'
#' # Simulate draws
#' Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
#' mu    <- c(1, -0.1, 0.8, 0, 2, 3)
#' draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
#' colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")
#'
#' # From posterior draws
#' pd.adjust(draws = draws, pi.max = 0.7, null.value = 0)
#'
#' # Order the rows by strength of evidence
#' pd.adjust(draws = draws, pi.max = 0.7, null.value = 0, order = TRUE)
#'
#' # Mix of directional and agnostic tests with parameter-specific nulls
#' pd.adjust(draws = draws, pi.max = 0.6, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
#'           direction = c("greater", "two.sided", "greater",
#'           "two.sided", "greater", "greater"))
#'
#' @importFrom stats median
#' @export
pd.adjust <- function(pd = NULL, draws = NULL, p0 = NULL, pi.max = NULL,
                      null.value = 0, direction = NULL,
                      order = FALSE, threshold = 0.975,
                      alpha = NULL) {

  if (is.null(p0) == is.null(pi.max)) {
    stop("Specify exactly one of `p0` (largest prior probability that no claim is correct) or `pi.max` (upper bound on the prior probability of either sign of a parameter).")
  }
  if (!is.null(p0)) {
    stopifnot("`p0`: must be a single number in (0, 1)" = length(p0) == 1L && p0 > 0 && p0 < 1)
  }
  if (!is.null(pi.max)) {
    stopifnot("`pi.max`: must be a single number in (0, 1)" = length(pi.max) == 1L && pi.max > 0 && pi.max < 1)
  }
  if (!is.null(alpha)) {
    stopifnot("`alpha`: must be a single number in (0, 1)" =
                length(alpha) == 1L && is.numeric(alpha) && alpha > 0 && alpha < 1)
    threshold <- 1 - alpha / 2   # two-sided cutoff implied by the nominal rate
  }
  stopifnot("`threshold`: must be a single number in (0.5, 1)" =
              length(threshold) == 1L && is.numeric(threshold) &&
              threshold > 0.5 && threshold < 1)

  if (!is.null(pd)) {
    if (!is.numeric(pd) || any(pd < 0.5 | pd > 1, na.rm = TRUE)) {
      stop("`pd` must be numeric in [0.5, 1].")
    }
    if (!is.null(direction)) {
      warning("`direction` is ignored when `pd` is supplied: all tests are treated as direction-agnostic (`two.sided`).")
    }
    direction <- rep("two.sided", length(pd))
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

    # Per-draw indicator that each parameter lies on its claimed side.
    # For two.sided tests the claimed side is the dominant marginal direction.
    side_ind <- mapply(function(j, d) {
      if (d == "greater")   centered[, j] > 0
      else if (d == "less") centered[, j] < 0
      else if (mean(centered[, j] > 0) >= mean(centered[, j] < 0)) centered[, j] > 0
      else                  centered[, j] < 0
    }, seq_len(nc), direction)
    side_ind <- matrix(side_ind, nrow = nrow(draws), ncol = nc)

    pd <- colMeans(side_ind)            # marginal pd per hypothesis
  }

  if (is.null(pd)) stop("Either `pd` or `draws` must be provided.")

  m <- length(pd)

  # Upper bound on the prior probability of the sign opposite to the claim. If `p0` (the
  # largest prior probability that no claim in the family is correct) is
  # supplied, pi.max = p0^(1/m) for parameters independent a priori. Below 0.5 the class of priors
  # is empty, so pi.max is set to 0.5 and no penalty is applied.
  if (is.null(pi.max)) pi.max <- p0^(1 / m)
  if (pi.max < 0.5) {
    warning(sprintf("pi.max = %.3g is below 0.5 and was set to 0.5: no penalty is applied and the threshold is left unchanged.", pi.max))
    pi.max <- 0.5
  }

  # A claim is made when its posterior probability is at least `threshold` under
  # every prior in the class; the smallest posterior odds are
  # pd/(1 - pd) * (1 - pi.max)/pi.max, so the rule is pd >= threshold.adj.
  alpha.nominal <- 2 * (1 - threshold)
  threshold.adj <- (threshold * pi.max) / ((1 - pi.max) * (1 - threshold) + threshold * pi.max)
  alpha.eff <- 2 * (1 - threshold.adj)

  # Posterior probability of the favored direction under the prior in the class
  # most skeptical of it (prior probability 1 - pi.max for that direction); the
  # other side has probability 1 - pd.adj under the same prior.
  pd.adj <- pd * (1 - pi.max) / (pd * (1 - pi.max) + (1 - pd) * pi.max)

  # Row order: by decreasing evidence when `order = TRUE`; otherwise input order.
  ord <- if (isTRUE(order)) order(pd, decreasing = TRUE) else seq_len(m)

  if (from_draws) {
    out <- data.frame(
      median.est = round(apply(draws, 2, median), 4),
      null.value = null.value,
      pd         = round(pd, 4),
      pd.adj     = round(pd.adj, 4),
      pi.max     = round(rep(pi.max, m), 4),
      m          = rep(m, m),
      direction  = direction
    )
  } else {
    out <- data.frame(
      median.est = NA,
      null.value = NA,
      pd         = round(pd, 4),
      pd.adj     = round(pd.adj, 4),
      pi.max     = round(rep(pi.max, m), 4),
      m          = rep(m, m),
      direction  = direction
    )
  }

  # Preserve hypothesis identity as row names, then apply the chosen order.
  nm <- if (from_draws) colnames(draws) else names(pd)
  if (!is.null(nm) && length(nm) == m) rownames(out) <- nm

  out <- out[ord, , drop = FALSE]

  class(out) <- c("pd_adjust", "data.frame")
  attr(out, "threshold")     <- threshold
  attr(out, "threshold.adj") <- threshold.adj
  attr(out, "alpha.nominal") <- alpha.nominal
  attr(out, "alpha.eff")     <- alpha.eff
  out
}

#' Print method for \code{pd_adjust} objects
#'
#' Prints a header with the number of tests \eqn{m}, \eqn{\bar{\pi}}, the
#' unadjusted threshold \eqn{c} and the adjusted threshold \eqn{c^*} (a
#' direction is claimed when \eqn{pd \ge c^*}), followed by the per-test table
#' and a short note on \code{pd.adj}.
#' Columns that take the same value in every row (\code{null.value},
#' \code{direction}) are reported once in the header. The full set of columns remains available
#' for programmatic access via the usual \code{$} and \code{[}.
#'
#' @param x A \code{pd_adjust} object returned by \code{\link{pd.adjust}}.
#' @param digits Number of digits used when printing. Defaults to 4.
#' @param ... Further arguments passed to \code{print.data.frame}.
#'
#' @return \code{x}, invisibly.
#' 
#' @keywords internal
#' @export
print.pd_adjust <- function(x, digits = 4, ...) {
  df  <- as.data.frame(x)
  m   <- as.integer(df$m[1])
  pi.max <- df$pi.max[1]
  cat("Prior-odds threshold for directional claims\n")
  cat(sprintf("  m = %d, pi.max = %.4g\n", m, pi.max))

  cadj <- attr(x, "threshold.adj")
  if (!is.null(cadj)) {
    cat(sprintf("  c = %.4g, c* = %.4g\n", attr(x, "threshold"), cadj))
  }

  show <- df[, setdiff(names(df), c("pi.max", "m")), drop = FALSE]
  # drop columns that are entirely NA (e.g. median.est / null.value for pd input)
  keep <- !vapply(show, function(col) all(is.na(col)), logical(1))
  show <- show[, keep, drop = FALSE]
  ch <- const_cols(show, c("null.value", "direction"))
  if (!is.null(ch$text)) cat("  ", ch$text, "\n", sep = "")
  cat("\n")
  print(show[, setdiff(names(show), ch$drop), drop = FALSE], digits = digits, ...)
  if ("pd.adj" %in% names(show)) {
    cat("\npd.adj: Pr(favored direction | y) under the prior in the class most skeptical\n")
    cat("        of it; 1 - pd.adj is the probability of the other side under that prior.\n")
    cat("        Below 0.5, the priors in the class disagree on the more probable direction.\n")
    if (!is.null(cadj))
      cat(sprintf("Claim a direction when pd.adj >= %.4g, or equivalently pd >= %.4g.\n",
                  attr(x, "threshold"), cadj))
  }

  invisible(x)
}

# Columns that take the same value in every row are reported once, in the
# printed header, so that the table stays narrow.
const_cols <- function(df, cols) {
  cols <- intersect(cols, names(df))
  same <- vapply(cols, function(cl) {
    v <- df[[cl]]
    !anyNA(v) && length(unique(v)) == 1L
  }, logical(1))
  drop <- cols[same]
  text <- if (length(drop)) {
    paste(sprintf("%s = %s", sub(".", " ", drop, fixed = TRUE),
                  vapply(drop, function(cl) format(df[[cl]][1]), character(1))),
          collapse = ", ")
  }
  list(drop = drop, text = text)
}

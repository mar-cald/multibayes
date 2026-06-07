#' Prior-odds adjustment for Probability of Direction \emph{pd}
#'
#' The function accepts either a vector of pre-computed \emph{pd} values or
#' a matrix of posterior draws, from which \emph{pd} values are computed
#' internally. Both direction-agnostic and directional tests are supported:
#' the \code{direction} argument controls which formulation is applied per
#' hypothesis. The adjustment is governed by the per-test prior probability
#' that an individual effect is null, \eqn{\pi_0} (argument \code{pi0}).
#' Equivalently, one may supply the global null \eqn{\Pi_0} (argument
#' \code{p0}), the prior that all \eqn{m} hypotheses are null, from which
#' \eqn{\pi_0 = \Pi_0^{1/m}} is recovered under independence.
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
#' \strong{Per-test correction and joint statement.} By default the function
#' returns the per-test correction only (the marginal \code{pd} and its adjusted
#' value \code{pd.adj}). Setting \code{joint = TRUE} \emph{additionally} returns
#' \code{joint_cum} and orders the rows by decreasing \code{pd}, so that
#' \code{joint_cum} is the cumulative joint probability that the \eqn{k} strongest
#' claims all hold in their claimed direction. When \code{draws} are supplied this joint is computed
#' directly from the draws by intersecting the per-draw directional events, so it
#' reflects the posterior dependence among parameters; when only \code{pd} is
#' supplied it is the running product of the \emph{pd} values, valid only under
#' independence. The claimed direction for a \code{"two.sided"} test is its
#' dominant posterior side, and the last entry of \code{joint_cum} is the joint
#' probability that all \eqn{m} claims hold simultaneously.
#'
#' @param pd Numeric vector of \emph{pd} values. For direction-agnostic tests,
#'   values must be in \eqn{[0.5, 1]}. For directional tests, values are raw
#'   one-sided probabilities in \eqn{[0, 1]}. Ignored if \code{draws} is
#'   supplied.
#' @param draws Optional matrix or data frame of posterior draws (columns = parameters).
#'   If provided, \emph{pd} values are computed automatically from the draws
#'   according to \code{direction} and \code{null.value}.
#' @param p0 Numeric scalar in \eqn{(0, 1)}. The prior probability that
#'   \strong{all} hypotheses are null simultaneously (the global null,
#'   \eqn{\Pi_0}). Supply \strong{either} \code{pi0} \strong{or} \code{p0},
#'   not both. When \code{p0} is given, the per-test prior is derived under
#'   independence as \eqn{\pi_0 = \Pi_0^{1/m}}.
#' @param pi0 Numeric scalar in \eqn{(0, 1)}. The per-test prior probability that
#'   an individual effect is null, applied to every test. Supply \strong{either} \code{pi0} \strong{or} \code{p0}, not both. When
#'   \code{pi0} is given, the implied global null is \eqn{\Pi_0 = \pi_0^{m}}
#'   (valid under independence) and is returned in \code{p0}.
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
#' @param order Logical. If \code{TRUE}, the returned rows are ordered by
#'   decreasing \code{pd} (strongest directional evidence first). If \code{FALSE}
#'   (the default), rows are returned in input order, so the output stays aligned
#'   with the columns of \code{draws} (or the elements of \code{pd}); use this for
#'   programmatic work that pairs \code{pd.adj} with external vectors. Ordering is
#'   enabled automatically when \code{joint = TRUE}.
#' @param joint Logical. If \code{TRUE}, the output \emph{additionally} includes
#'   the cumulative joint probability \code{joint_cum} (see Details). Because that
#'   cumulative joint is only interpretable when accumulated from the strongest
#'   claim down, \code{joint = TRUE} orders the rows by decreasing \code{pd}
#'   (i.e. it implies \code{order = TRUE}). Defaults to \code{FALSE}, in which case
#'   only the per-test correction is returned.
#'
#' @return An object of class \code{pd_adjust} (a \code{data.frame} with a
#'   dedicated \code{\link{print.pd_adjust}} method), with one row per
#'   hypothesis. Columns: \code{pd} (values used in the adjustment),
#'   \code{pd.adj} (adjusted values), \code{pi0} (per-test null prior
#'   \eqn{\pi_0}), and \code{m} (number of tests). The global null
#'   \eqn{\Pi_0 = \pi_0^m} is \strong{not} returned: it is exact only under
#'   independence, and reporting it alongside per-test results can mislead when
#'   the tests are dependent; the per-test prior \eqn{\pi_0} is the robust
#'   quantity (the expected proportion of nulls under exchangeability). When
#'   \code{joint = TRUE} a \code{joint_cum} column is added; when
#'   \code{order = TRUE} the rows are sorted by decreasing \code{pd}. For
#'   direction-agnostic tests, both \code{pd} and \code{pd.adj} are bounded in
#'   \eqn{[0.5, 1]}; for directional tests, both are on \eqn{[0, 1]}, with values
#'   below \eqn{0.5} indicating that the data (and the adjustment) favoured the
#'   opposite direction. When \code{draws} are supplied, the output additionally
#'   includes \code{median.est} (posterior median per parameter), \code{null.value},
#'   and \code{direction}. The \code{print} method summarises the constant
#'   quantities (\code{pi0}, \code{m}) in a header and displays the per-test
#'   table; the columns themselves remain available for programmatic access.
#'
#' @examples
#' # From a vector of pd values
#' pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763,
#' H5 = 0.891, H6 = 0.987)
#' pd.adjust(pd = pd_values, p0 = 0.4)
#'
#' # Equivalent call specifying the per-test prior pi0 directly
#' pd.adjust(pd = pd_values, pi0 = 0.4^(1/6))
#'
#' # Simulate draws
#' Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
#' mu    <- c(1, -0.1, 0.8, 0, 2, 3)
#' draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
#' colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")
#'
#' # Per-test correction only (default)
#' pd.adjust(draws = draws, p0 = 0.4, null.value = 0)
#'
#' # Also return the cumulative joint statement, ordered by evidence
#' pd.adjust(draws = draws, p0 = 0.4, null.value = 0, order = TRUE, joint = TRUE)
#'
#' # Mix of directional and agnostic tests with parameter-specific nulls
#' pd.adjust(draws = draws, p0 = 0.2, null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
#'           direction = c("greater", "two.sided", "greater",
#'           "two.sided", "greater", "greater"))
#'
#' @importFrom stats median
#' @export
pd.adjust <- function(pd = NULL, draws = NULL, p0 = NULL, pi0 = NULL,
                      null.value = 0, direction = NULL,
                      order = FALSE, joint = FALSE) {

  if (is.null(p0) == is.null(pi0)) {
    stop("Specify exactly one of `p0` (global-null prior Pr(H0)) or `pi0` (per-test null prior).")
  }
  if (!is.null(p0)) {
    stopifnot("`p0`: must be a single number in (0, 1)" = length(p0) == 1L && p0 > 0 && p0 < 1)
  }
  if (!is.null(pi0)) {
    stopifnot("`pi0`: must be a single number in (0, 1)" = length(pi0) == 1L && pi0 > 0 && pi0 < 1)
  }

  if (!is.null(pd)) {
    if (!is.numeric(pd) || any(pd < 0.5 | pd > 1, na.rm = TRUE)) {
      stop("`pd` must be numeric in [0.5, 1].")
    }
    if (!is.null(direction)) {
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

  # Per-test null prior (prior_H0). If `pi0` is supplied, it is used directly: under
  # exchangeability it is the expected proportion of nulls. If `p0` (the global null)
  # is supplied, the per-test prior is derived under independence as pi0 = p0^(1/m).
  # Conversely, a supplied `pi0` implies the global null p0 = pi0^m, also under independence.
  if (!is.null(pi0)) {
    prior_H0 <- pi0
    p0 <- pi0^m
  } else {
    prior_H0 <- p0^(1 / m)
    pi0 <- prior_H0
  }

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

  # Row order: by decreasing evidence when `order = TRUE`. `joint = TRUE` also
  # forces this order, since the cumulative joint is only interpretable when
  # accumulated strongest-first; otherwise rows stay in input order.
  ord <- if (isTRUE(order) || isTRUE(joint)) order(pd, decreasing = TRUE) else seq_len(m)

  # Optional cumulative joint probability that the first k claims (in row order) hold.
  if (isTRUE(joint)) {
    if (from_draws) {
      # exact, from the draws: respects posterior dependence among parameters
      joint_cum <- numeric(m)
      running   <- rep(TRUE, nrow(side_ind))
      for (k in seq_len(m)) {
        running      <- running & side_ind[, ord[k]]
        joint_cum[k] <- mean(running)
      }
    } else {
      # only marginal pd available: running product (valid under independence)
      joint_cum <- cumprod(pd[ord])
    }
  }

  if (from_draws) {
    out <- data.frame(
      median.est = round(apply(draws, 2, median), 4),
      null.value = null.value,
      pd         = round(pd, 4),
      pd.adj     = round(pd.adj, 4),
      pi0        = round(rep(pi0, m), 4),
      m          = rep(m, m),
      direction  = direction
    )
  } else {
    out <- data.frame(
      median.est = NA,
      null.value = NA,
      pd         = round(pd, 4),
      pd.adj     = round(pd.adj, 4),
      pi0        = round(rep(pi0, m), 4),
      m          = rep(m, m),
      direction  = rep("two.sided", m)
    )
  }

  # Preserve hypothesis identity as row names, then apply the chosen order.
  nm <- if (from_draws) colnames(draws) else names(pd)
  if (!is.null(nm) && length(nm) == m) rownames(out) <- nm

  out <- out[ord, , drop = FALSE]
  if (isTRUE(joint)) out$joint_cum <- round(joint_cum, 4)   # aligned with the rows shown

  class(out) <- c("pd_adjust", "data.frame")
  out
}

#' Print method for \code{pd_adjust} objects
#'
#' Summarises the constant quantities (\code{pi0}, \code{m}) in a header and
#' prints the per-test table. The full set of columns remains available for
#' programmatic access via the usual \code{$} and \code{[}.
#'
#' @param x A \code{pd_adjust} object returned by \code{\link{pd.adjust}}.
#' @param digits Number of significant digits used when printing. Defaults to 4.
#' @param ... Further arguments passed to \code{print.data.frame}.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.pd_adjust <- function(x, digits = 4, ...) {
  df <- as.data.frame(x)
  cat("Prior-odds adjusted probability of direction\n")
  cat(sprintf("Tests: %d  | pi0 = %.4g\n\n",
              as.integer(df$m[1]), df$pi0[1]))

  show <- df[, setdiff(names(df), c("pi0", "m")), drop = FALSE]
  # drop columns that are entirely NA (e.g. mean.est / null.value for pd input)
  keep <- !vapply(show, function(col) all(is.na(col)), logical(1))
  print(show[, keep, drop = FALSE], digits = digits, ...)

  if (!is.null(df$joint_cum)) {
    cat("\njoint_cum: cumulative joint Pr(first k claims hold), in the order shown.\n")
  }
  invisible(x)
}

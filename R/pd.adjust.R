#' Prior-odds adjustment for Probability of Direction \emph{pd}
#'
#' The function accepts either a vector of pre-computed \emph{pd} values or
#' a matrix of posterior draws, from which \emph{pd} values are computed
#' internally. Both direction-agnostic and directional tests are supported:
#' the \code{direction} argument controls which formulation is applied per
#' hypothesis. The adjustment is governed by the per-test prior probability
#' that an individual effect is null, \eqn{\pi_0} (argument \code{pi0}).
#' Equivalently, one may supply the global null \eqn{\Pr(H_0)} (argument
#' \code{p0}), the prior that all \eqn{m} hypotheses are null, from which
#' \eqn{\pi_0 = \Pr(H_0)^{1/m}} is recovered under independence.
#'
#' @details
#' The adjustment follows from Bayes' theorem. Given a per-hypothesis prior
#' \eqn{\pi_0} and its complement \eqn{\pi_1 = 1 - \pi_0},
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
#' \code{pd.adjust} returns the per-test correction only (the marginal \code{pd}
#' and its adjusted value \code{pd.adj}). Family-wise summaries across hypotheses,
#' the cumulative joint statement and simultaneous credible intervals, are
#' provided separately by \code{\link{joint}}.
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
#'   \eqn{\Pr(H_0)}). Supply \strong{either} \code{pi0} \strong{or} \code{p0},
#'   not both. When \code{p0} is given, the per-test prior is derived under
#'   independence as \eqn{\pi_0 = \Pr(H_0)^{1/m}}.
#' @param pi0 Numeric scalar in \eqn{(0, 1)}. The per-test prior probability that
#'   an individual effect is null, applied to every test. Supply \strong{either} \code{pi0} \strong{or} \code{p0}, not both.
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
#'   with the columns of \code{draws} (or the elements of \code{pd}).
#' @param threshold Numeric scalar in \eqn{(0.5, 1)}. The decision cutoff applied
#'   to \emph{pd} (default \code{0.975}, a two-sided per-test
#'   \eqn{\alpha = 0.05}). It does \strong{not} change \code{pd.adj}; it is used
#'   only to report the equivalent \strong{adjusted threshold on the raw}
#'   \emph{pd}, since rejecting \eqn{pd_{adj} > c} is identical to rejecting raw
#'   \eqn{pd > c^*} with
#'   \eqn{c^* = c\,\pi_0 / ((1-\pi_0)(1-c) + c\,\pi_0)}. 
#' @param alpha Optional numeric scalar in \eqn{(0, 1)}. A convenience
#'   parameterization of \code{threshold} on the error scale: when supplied, the
#'   two-sided cutoff is set to \eqn{c = 1 - \alpha/2} (so \code{alpha = 0.05}
#'   gives \code{threshold = 0.975}), \strong{overriding} any \code{threshold}
#'   value. 
#' @param fwer Logical (default \code{FALSE}). If \code{TRUE}, the output also
#'   reports a \strong{family-wise} cutoff on the raw \emph{pd} that controls the
#'   prior-averaged FWER at the nominal level over the \eqn{m} tests, given
#'   \eqn{\pi_0}. Standard corrections (e.g. Bonferroni) implicitly assume every
#'   test is null (\eqn{\pi_0 = 1}); here the affordable per-test rate is
#'   \eqn{[1 - (1 - \alpha)^{1/m}] / \pi_0}, relaxing the standard correction by
#'   \eqn{1/\pi_0} and reducing to it as \eqn{\pi_0 \to 1}. This is a
#'   prior-averaged (Bayesian) FWER under the assumed \eqn{\pi_0}, \strong{not}
#'   strong frequentist control.
#'
#' @return An object of class \code{pd_adjust} (a \code{data.frame} with a
#'   dedicated \code{\link{print.pd_adjust}} method), with one row per
#'   hypothesis. Columns: \code{pd} (values used in the adjustment),
#'   \code{pd.adj} (adjusted values), \code{pi0} (per-test null prior
#'   \eqn{\pi_0}), and \code{m} (number of tests). When
#'   \code{order = TRUE} the rows are sorted by decreasing \code{pd}. For
#'   direction-agnostic tests, both \code{pd} and \code{pd.adj} are bounded in
#'   \eqn{[0.5, 1]}; for directional tests, both are on \eqn{[0, 1]}, with values
#'   below \eqn{0.5} indicating that the data (and the adjustment) favoured the
#'   opposite direction. When \code{draws} are supplied, the output additionally
#'   includes \code{median.est} (posterior median per parameter), \code{null.value},
#'   and \code{direction}. The \code{print} method summarises the constant
#'   quantities (\code{pi0}, \code{m}) in a header and displays the per-test
#'   table; the columns themselves remain available for programmatic access.
#'   The object also carries, as attributes, the decision cutoff \code{threshold}
#'   (\eqn{c}), the equivalent adjusted cutoff on the raw \emph{pd}
#'   \code{threshold.adj} (\eqn{c^*}), and the corresponding two-sided per-test
#'   Type I rates \code{alpha.nominal} (\eqn{2(1-c)}) and \code{alpha.eff}
#'   (\eqn{2(1-c^*)}); retrieve them with, e.g., \code{attr(x, "threshold.adj")}.
#'   The cutoffs \code{threshold.adj} and \code{threshold.fwer} are exact
#'   transformations of the decision rule, but the reported per-test \emph{rates}
#'   (\code{alpha.eff}, and the \code{fwer} calibration) use the identity
#'   \eqn{\alpha = 2(1 - c)}, which holds **only for a diffuse within-model prior** or
#'   **large \eqn{n}** (when the null \emph{pd} is approximately uniform on
#'   \eqn{[0.5, 1]}). Under an informative prior the true per-test rate is
#'   smaller, so these rates are conservative and \code{threshold.fwer} controls
#'   the FWER at or below the stated level.
#'
#' @references
#' Westfall, P. H., Johnson, W. O., & Utts, J. M. (1997). A bayesian perspective on the bonferroni adjustment.
#' \emph{Biometrika}, 84(2), 419–427.
#' 
#' Jeffreys, H. (1938). Significance tests when several degrees of freedom arise simultaneously. \emph{Proceedings
#' of the Royal Society of London. Series A. Mathematical and Physical Sciences}, 165(921), 161–198.
#' https://doi.org/10.1098/rspa.1938.0052
#' 
#' Storey, J. D. (2003). The positive false discovery rate: A bayesian interpretation and the q-value. \emph{The Annals
#' of Statistics}, 31(6), 2013–2035. https://doi.org/10.1214/aos/1074290335
#'
#' Scott, J. G., & Berger, J. O. (2006). An exploration of aspects of bayesian multiple testing. 
#' \emph{Journal of Statistical Planning and Inference}, 136(7), 2144–2162.
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
#' # Order the rows by strength of evidence
#' pd.adjust(draws = draws, p0 = 0.4, null.value = 0, order = TRUE)
#'
#' # Family-wise summaries are obtained separately with joint()
#' joint(draws, null.value = 0)               # cumulative joint statement
#' joint(draws, interval = TRUE)              # simultaneous credible intervals
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
                      order = FALSE, threshold = 0.975,
                      alpha = NULL, fwer = FALSE) {

  if (is.null(p0) == is.null(pi0)) {
    stop("Specify exactly one of `p0` (global-null prior Pr(H0)) or `pi0` (per-test null prior).")
  }
  if (!is.null(p0)) {
    stopifnot("`p0`: must be a single number in (0, 1)" = length(p0) == 1L && p0 > 0 && p0 < 1)
  }
  if (!is.null(pi0)) {
    stopifnot("`pi0`: must be a single number in (0, 1)" = length(pi0) == 1L && pi0 > 0 && pi0 < 1)
  }
  if (!is.null(alpha)) {
    stopifnot("`alpha`: must be a single number in (0, 1)" =
                length(alpha) == 1L && is.numeric(alpha) && alpha > 0 && alpha < 1)
    threshold <- 1 - alpha / 2   # two-sided cutoff implied by the nominal rate
  }
  stopifnot("`threshold`: must be a single number in (0.5, 1)" =
              length(threshold) == 1L && is.numeric(threshold) &&
              threshold > 0.5 && threshold < 1)
  stopifnot("`fwer`: must be TRUE or FALSE" = length(fwer) == 1L && is.logical(fwer))

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

  # Model-free decision threshold implied by the prior odds. Rejecting
  # `pd.adj > threshold` is identical to rejecting the raw `pd > threshold.adj`,
  # so threshold.adj is the equivalent (stricter) cutoff to apply to the
  # unadjusted pd. It is a pure transformation of the decision rule and never
  # uses the shape of the posterior, hence exact for any draws (normal, t, ...).
  # The two-sided per-test Type I rate at a cutoff `cut` is 2 * (1 - cut).
  alpha.nominal <- 2 * (1 - threshold)
  if (prior_H0 > 0.5) {
    threshold.adj <- (threshold * prior_H0) /
      ((1 - prior_H0) * (1 - threshold) + threshold * prior_H0)
  } else {
    threshold.adj <- threshold   # no adjustment applied (non-conservative prior)
  }
  alpha.eff <- 2 * (1 - threshold.adj)

  # Optional family-wise (FWER) cutoff. Standard corrections (e.g. Bonferroni)
  # control the FWER under the worst case that every test is null (pi0 = 1). If a
  # fraction pi0 of the m tests are null, the prior-averaged FWER of a per-test
  # rate a is 1 - (1 - pi0 * a)^m; setting this to the family target alpha gives
  # the affordable per-test rate a = (1 - (1 - alpha)^(1/m)) / pi0 and the
  # corresponding (two-sided) raw-pd cutoff below. This relaxes the standard
  # correction by 1/pi0 and reduces to it as pi0 -> 1. It is a prior-averaged
  # (Bayesian) FWER, not strong frequentist control.
  if (isTRUE(fwer)) {
    alpha.fwer    <- alpha.nominal                       # family-wise target
    apt_fw        <- (1 - (1 - alpha.fwer)^(1 / m)) / prior_H0
    threshold.fwer <- max(0.5, 1 - apt_fw / 2)
    alpha.fwer.pt <- 2 * (1 - threshold.fwer)            # implied per-test rate
  } else {
    alpha.fwer <- threshold.fwer <- alpha.fwer.pt <- NA_real_
  }

  # Row order: by decreasing evidence when `order = TRUE`; otherwise input order.
  ord <- if (isTRUE(order)) order(pd, decreasing = TRUE) else seq_len(m)

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

  class(out) <- c("pd_adjust", "data.frame")
  attr(out, "threshold")     <- threshold
  attr(out, "threshold.adj") <- threshold.adj
  attr(out, "alpha.nominal") <- alpha.nominal
  attr(out, "alpha.eff")     <- alpha.eff
  attr(out, "threshold.fwer") <- threshold.fwer
  attr(out, "alpha.fwer")     <- alpha.fwer
  attr(out, "alpha.fwer.pt")  <- alpha.fwer.pt
  out
}

#' Print method for \code{pd_adjust} objects
#'
#' Summarises the constant quantities (\code{pi0}, \code{m}) in a header and
#' prints the per-test table. The full set of columns remains available for
#' programmatic access via the usual \code{$} and \code{[}.
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
  df <- as.data.frame(x)
  cat("Prior-odds adjusted probability of direction\n")
  cat(sprintf("Tests: %d  |  pi0 = %.4g\n",
              as.integer(df$m[1]), df$pi0[1]))

  # Wrap header lines at the console width so they never overflow (e.g. in PDF
  # output); plain cat() would print a single long line that runs off the page.
  hdr <- function(s) cat(strwrap(s, width = getOption("width", 80L)), "", sep = "\n")

  cadj <- attr(x, "threshold.adj")
  if (!is.null(cadj)) {
    hdr(sprintf("Decision: raw pd > %.4g (adjusted from %.4g); effective per-test alpha %.3g (from %.3g)",
                cadj, attr(x, "threshold"), attr(x, "alpha.eff"), attr(x, "alpha.nominal")))
  }
  cfw <- attr(x, "threshold.fwer")
  if (!is.null(cfw) && !is.na(cfw)) {
    hdr(sprintf("Family-wise: raw pd > %.4g (FWER %.3g over m = %d under the assumed pi0; per-test alpha %.3g)",
                cfw, attr(x, "alpha.fwer"), as.integer(df$m[1]), attr(x, "alpha.fwer.pt")))
  }
  cat("\n")

  show <- df[, setdiff(names(df), c("pi0", "m")), drop = FALSE]
  # drop columns that are entirely NA (e.g. mean.est / null.value for pd input)
  keep <- !vapply(show, function(col) all(is.na(col)), logical(1))
  print(show[, keep, drop = FALSE], digits = digits, ...)

  invisible(x)
}

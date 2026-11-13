#' Prior-odds adjustment for the probability of direction
#'
#' Decides which directional claims to make in a family of tests, from the
#' probability of direction (\emph{pd}) of each parameter. Two separate choices
#' enter the decision. The threshold \eqn{c} (\code{threshold}) is the posterior
#' probability required to make a claim. The bound \eqn{\bar{\pi}}
#' (\code{pi.max}) is the largest prior probability allowed for either
#' direction of a parameter, and sets how skeptical the prior may be about a
#' claim. For each test the function returns \code{pd} and \code{pd.adj}, the
#' posterior probability of the favored direction when that direction has
#' prior probability \eqn{1 - \bar{\pi}}. A direction is claimed when
#' \code{pd.adj} is at least \eqn{c}.
#'
#' @details
#' \strong{What is computed.} \emph{pd} comes from the fitted model, whose
#' prior (the estimation prior) must be continuous and symmetric about
#' \code{null.value}. Flat priors and symmetric priors centred on
#' \code{null.value}, such as normal or Student-t priors, qualify. Each direction then has prior probability 0.5, and
#' \eqn{pd/(1 - pd)} is the Bayes factor of the favored direction against the
#' other.
#'
#' With several tests, the prior odds of each claim are lowered (Jeffreys,
#' 1938; Westfall et al., 1997). The direction of a claim is not known before
#' the data, so the function considers every prior obtained by reweighting the
#' two sides of the estimation prior so that each direction has prior
#' probability between \eqn{1 - \bar{\pi}} and \eqn{\bar{\pi}}. A claim must
#' hold under all of these priors, and it is enough to check the one most
#' skeptical of the favored direction, under which
#' \deqn{pd_{adj} = \frac{pd\,(1-\bar{\pi})}{pd\,(1-\bar{\pi}) + (1-pd)\,\bar{\pi}}.}
#' The model is not refitted, and the posterior used for estimation is not
#' changed.
#'
#' \strong{The threshold on pd.} The threshold \eqn{c} is not changed by the
#' adjustment. Because \eqn{pd_{adj}} increases with \emph{pd},
#' \eqn{pd_{adj} \ge c} holds exactly when \emph{pd} reaches
#' \deqn{c^* = \frac{c\,\bar{\pi}}{(1-c)(1-\bar{\pi}) + c\,\bar{\pi}},}
#' the value of \emph{pd} required by \eqn{c} and \eqn{\bar{\pi}} together. The
#' adjustment is thus equivalent to raising the threshold on \emph{pd} from
#' \eqn{c} to \eqn{c^*}, reported as \code{threshold.adj}.
#'
#' \strong{Choosing pi.max.} \code{pi.max} says how skeptical you are of a
#' directional claim before the data. With \code{pi.max = 0.5} no adjustment is
#' made. With \code{pi.max = 0.8}, a claim needs a Bayes factor four times
#' larger. The same \code{pi.max} is used for every test, so the parameters
#' should be exchangeable in direction, that is, no parameter is judged more
#' likely than another to be positive or negative. They can be on different
#' scales, because the sign of a parameter does not depend on its unit. The
#' proportion of effects in the family that could be null is a reference
#' value. If up to 70\% of the effects could be null, at most 30\% lie in each
#' direction, so each direction should be allowed a prior probability as low
#' as 0.3, which gives \code{pi.max = 0.7}. A larger value gives
#' fewer false claims at a cost in power. The priors considered are
#' continuous, so a null effect has probability zero under all of them, and
#' \code{pd.adj} is the posterior probability of a sign. To test whether an
#' effect is exactly null, use a prior that puts probability on the null
#' value, such as a spike-and-slab prior. A narrower prior in the fitted model
#' is no substitute, because a symmetric prior gives each direction prior
#' probability 0.5 however narrow it is. Alternatively, \code{p0} gives the largest
#' prior probability that no claim in the family is correct. For parameters
#' independent a priori, \eqn{\bar{\pi} = \Pr(H_0)^{1/m}}, which grows with the
#' number of tests \eqn{m}.
#'
#' \strong{Direction.} With \code{"two.sided"}, \emph{pd} is the larger of
#' \eqn{\Pr(\theta > \theta_\text{null} \mid y)} and
#' \eqn{\Pr(\theta < \theta_\text{null} \mid y)}, between 0.5 and 1. With
#' \code{"greater"} or \code{"less"}, \emph{pd} is the posterior probability of
#' the predicted side, between 0 and 1, and \code{pd.adj} is computed under the
#' prior most skeptical of that side. The mode can differ from one test to
#' another.
#'
#' @param pd Numeric vector of \emph{pd} values in \eqn{[0.5, 1]}. Every test
#'   is then treated as \code{"two.sided"} and \code{direction} is ignored.
#'   Ignored if \code{draws} is supplied.
#' @param draws Matrix or data frame of posterior draws, one column per
#'   parameter. \emph{pd} is computed from the draws according to
#'   \code{direction} and \code{null.value}.
#' @param p0 Number in \eqn{(0, 1)}. The largest prior probability that no
#'   claim in the family is correct, \eqn{\Pr(H_0)}, from which
#'   \eqn{\bar{\pi} = \Pr(H_0)^{1/m}}. Supply either \code{p0} or
#'   \code{pi.max}, not both.
#' @param pi.max Number in \eqn{(0, 1)}. The largest prior probability allowed
#'   for either direction of a parameter, chosen before the data. See Details
#'   for how to choose it. Values below 0.5 are set to 0.5, with a warning, so
#'   that no adjustment is made. Supply either \code{pi.max} or \code{p0}, not
#'   both.
#' @param null.value Number or vector, the null value of each parameter on the
#'   scale of the draws. A single number applies to all parameters. Ignored
#'   when \code{pd} is supplied. Defaults to \code{0}.
#' @param direction Character vector with elements \code{"two.sided"},
#'   \code{"greater"}, or \code{"less"}, one per parameter, or a single value
#'   for all. \code{NULL} (the default) means \code{"two.sided"} for all.
#' @param order Logical. If \code{TRUE}, rows are sorted by decreasing
#'   \code{pd}. If \code{FALSE} (the default), rows keep the order of the
#'   input.
#' @param threshold Number in \eqn{(0.5, 1)}. The threshold \eqn{c}, the
#'   posterior probability required to make a claim. Defaults to \code{0.975}.
#' @param alpha Optional number in \eqn{(0, 1)}. If supplied, the threshold is
#'   set to \eqn{c = 1 - \alpha/2} and \code{threshold} is ignored.
#'
#' @return An object of class \code{pd_adjust}, a data frame with one row per
#'   test and the columns \code{pd}, \code{pd.adj}, \code{pi.max}, \code{m}
#'   (the number of tests), and \code{direction}. With \code{draws}, it also
#'   has \code{median.est} (the posterior median) and \code{null.value}.
#'
#'   The thresholds are stored as attributes, \code{threshold} (\eqn{c}) and
#'   \code{threshold.adj} (\eqn{c^*}), together with \code{alpha.nominal}
#'   (\eqn{2(1-c)}) and \code{alpha.eff} (\eqn{2(1-c^*)}). The last two are
#'   per-test error rates only if \emph{pd} is uniform on \eqn{[0.5, 1]} for a
#'   null effect, which holds approximately under diffuse priors.
#'
#'   See \code{\link{print.pd_adjust}} and \code{\link{plot.pd_adjust}}.
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
#' # The same, giving pi.max directly
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
#' # Rows sorted by decreasing pd
#' pd.adjust(draws = draws, pi.max = 0.7, null.value = 0, order = TRUE)
#'
#' # Different directions and null values across parameters
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
    stop("Specify exactly one of `p0` (largest prior probability that no claim is correct) or `pi.max` (largest prior probability allowed for either direction of a parameter).")
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
    warning(sprintf("pi.max = %.3g is below 0.5 and was set to 0.5, so no adjustment is made.", pi.max))
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
#' Prints the number of tests \eqn{m}, \eqn{\bar{\pi}}, the threshold \eqn{c}
#' and the value \eqn{c^*} that \emph{pd} must reach, then one row per test
#' with \code{pd} and \code{pd.adj}. Columns with the same value in every row
#' (\code{null.value}, \code{direction}) are shown once, above the table.
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
  cat("Prior-odds adjustment for directional claims\n")
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
    cat("\npd.adj: pd when the favored direction has prior probability 1 - pi.max.\n")
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

# multibayes 0.2.0

## Breaking changes

* `pd.adjust()` reports the threshold `c*` that `pd` must reach
  (`attr(x, "threshold.adj")`): a claim is made only if its posterior
  probability is at least `c` under every prior in a class that lowers the prior
  odds of either direction, which is equivalent to `pd >= c*`.
* `pd.adj` has a new meaning: it is the posterior probability of the claimed
  direction under the prior in the class most skeptical of it, so that the
  other side has probability `1 - pd.adj` under the same prior. It is no longer
  floored at 0.5, and `pd.adj >= c` is the same rule as `pd >= c*`.
* `joint()` now returns objects of class `joint_directional` (default) or
  `joint_interval` (`interval = TRUE`). Both are still data frames.

## Other changes

* Directions are claimed when `pd.adj >= c` (equivalently `pd >= c*`), the
  analogue of `p <= alpha`; the printed rule and the documentation say so.
* The simultaneous intervals of `joint(interval = TRUE)` are described by their
  joint posterior probability (`prob`), not as coverage.
* `pi0` is an upper bound, set by the researcher, on the prior probability
  of either sign of a parameter; `p0` is the largest prior
  probability that no claim in the family is correct.
* Values of `pi0` below 0.5 (including those implied by `p0`) are set to 0.5,
  with a warning, so that no penalty is applied; `pi0 = 0.5` is accepted
  without a warning.
* Clearer printed output: the header of `pd.adjust()` shows only `m`, `pi0`,
  and the thresholds `c` and `c*`; `null.value` and `direction`, when they are
  the same for every test, are reported once in the header rather than as
  columns, so that tables fit on one line. `joint()` results print in the same
  way.
* Clearer warning when `direction` is supplied together with `pd`.

# multibayes 0.1.0

* First version.

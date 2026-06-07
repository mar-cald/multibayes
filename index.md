# multibayes

**multibayes** provides tools for Bayesian multiple testing, currently
centred on prior-odds adjustment of the Probability of Direction (*pd*).
Both direction-agnostic and directional tests are supported, and can be
mixed freely within a single call.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mar-cald/multibayes")
```

## Overview

| Function | Input | Correction type |
|----|----|----|
| [`pd.adjust()`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md) | Posterior draws or *pd* vector | Prior-odds adjustment for *pd* (agnostic or directional) |

------------------------------------------------------------------------

## `pd.adjust()`: Prior-Odds Adjustment for *pd*

The Probability of Direction (*pd*) is often used as a decision rule: a
hypothesis is accepted when *pd* exceeds a threshold. When this decision
is made repeatedly across several hypotheses, the global null
probability can be much lower than intended.
[`pd.adjust()`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md)
addresses this problem by making explicit the prior probability of the
global null.

The primary input is the per-test prior $`\pi_0`$, the probability that
an individual effect is null. Under exchangeability this is the expected
proportion of null effects in the family. Equivalently, you may specify
the global null, *p0* ($`\Pr(H_0)`$), the prior probability that **all**
$`m`$ hypotheses are null, from which the per-test prior is recovered
under independence:

``` math
\Pr(H_0) = \pi_0^{m}, \qquad \pi_0 = \Pr(H_0)^{1/m}.
```

Each *pd* is then reweighted by Bayes’ theorem:

``` math
pd_{\text{adj}} = \frac{pd \cdot \pi_1}{pd \cdot \pi_1 + (1 - pd) \cdot \pi_0}
```

Because the prior is conservative ($`\pi_0 > \pi_1`$), the adjustment
always shrinks *pd* toward its lower bound.

The function supports two testing modes, which can be mixed across
hypotheses within the same call:

- **Direction-agnostic** (`direction = "two.sided"`): *pd* =
  $`\max\big(\Pr(\hat\theta > \theta_\text{null}),\, \Pr(\hat\theta < \theta_\text{null})\big)`$,
  bounded in $`[0.5, 1]`$ by construction; $`pd_\text{adj}`$ is also
  floored at $`0.5`$.
- **Directional** (`direction = "greater"` or `"less"`): *pd* is the raw
  one-sided posterior probability on the predicted side, on $`[0, 1]`$.
  Values below $`0.5`$ indicate that the posterior is concentrated
  opposite to the predicted direction; the adjustment will further
  shrink such values, reflecting the combined weight of the data and the
  conservative prior against the hypothesis.

### Usage

``` r

library(multibayes)

# From a vector of pd values (direction-agnostic)
pd_values <- c(H1 = 0.999, H2 = 0.946, H3 = 0.813, H4 = 0.763, H5 = 0.891, H6 = 0.987)

# Specify the per-test null prior pi0 directly (primary interface)
pd.adjust(pd = pd_values, pi0 = 0.65)

# ...or specify the global null p = Pr(all null); pi0 is then p^(1/m)
pd.adjust(pd = pd_values, p0 = 0.4)

# Simulate posterior draws
Sigma <- matrix(0, nrow = 6, ncol = 6); diag(Sigma) <- 1
mu <- c(1, -0.1, 0.8, 0, 2, 3)
draws <- MASS::mvrnorm(n = 4000, mu = mu, Sigma = Sigma)
colnames(draws) <- c("H1", "H2", "H3", "H4", "H5", "H6")

# From posterior draws (per-test correction only)
pd.adjust(draws = draws, p0 = 0.4, null.value = 0)

# Also return the cumulative joint statement (this orders rows by evidence)
pd.adjust(draws = draws, p0 = 0.4, null.value = 0, joint = TRUE)

# Mix of directional and agnostic tests with parameter-specific nulls
pd.adjust(
  draws = draws,
  pi0 = 0.4,
  null.value = c(0.2, 0, 0.2, 0, 0.5, 0.5),
  direction  = c("greater", "two.sided", "greater", "two.sided", "greater", "greater")
)
```

### Output

[`pd.adjust()`](https://mar-cald.github.io/multibayes/reference/pd.adjust.md)
returns a `pd_adjust` object — a `data.frame` with a tidy `print` method
that summarises the constant quantities (`pi0`, `m`) in a header and
shows the per-test table. With `draws`, it has one row per hypothesis:

| Column | Description |
|----|----|
| `median.est` | Posterior median per parameter (the sign of `median.est − null.value` matches `pd`’s dominant direction) |
| `null.value` | Null reference value used |
| `pd` | *pd* used in the adjustment; in $`[0.5, 1]`$ for agnostic tests, $`[0, 1]`$ for directional tests |
| `pd.adj` | Adjusted *pd* after prior-odds correction; same bounds as `pd` |
| `pi0` | Per-test null prior $`\pi_0`$ — the expected proportion of nulls (robust under exchangeability) |
| `m` | Family size used |
| `direction` | Testing mode: `"greater"`, `"less"`, or `"two.sided"` |
| `joint_cum` | *(only if `joint = TRUE`)* cumulative joint probability that the first *k* hypotheses (in row order) all hold in their claimed direction — computed from the draws (dependence-aware) or, for `pd` input, as the running product of `pd` (independence) |

The global null $`\Pi_0 = \pi_0^m`$ is **not** reported as a column: it
holds only under independence and can mislead when tests are dependent,
so the per-test prior $`\pi_0`$ (`pi0`) is reported instead. You may
still *supply* the prior as `p0`; it is converted to `pi0` internally.

By default only the **per-test correction** is returned. Set
**`joint = TRUE`** to *additionally* return the **joint statement**
(`joint_cum`); this also orders the rows by decreasing `pd`, since the
cumulative joint is only interpretable strongest-first.
(**`order = TRUE`** on its own just sorts the rows by `pd`, without
adding the joint.) When a `pd` vector is supplied directly,
`median.est`, `null.value`, and `direction` are not returned. The header
is cosmetic: every quantity remains available as a column for
programmatic use (e.g. `result$pd.adj`, `result$joint_cum`).

### Choosing `pi0` (or `p0`), `null.value`, and `direction`

- **`pi0`** (primary) is your prior probability that any single effect
  is null — equivalently, the expected proportion of null effects in the
  family. Alternatively, **`p0`** encodes the prior belief that **all**
  tested hypotheses are simultaneously null, with
  $`\Pr(H_0) = \pi_0^m`$. Supply one of the two.
- **`null.value`** sets the null reference value for each parameter. A
  scalar applies the same null to all parameters; a vector allows a
  different null per parameter.
- **`direction`** selects the testing mode per hypothesis (`"greater"`,
  `"less"`, `"two.sided"`). A scalar is recycled across all hypotheses;
  a mixed vector applies different modes within the same call. For
  directional tests, values of `pd` and `pd.adj` below `0.5` indicate
  that the data contradicted the predicted direction; the adjustment
  amplifies this evidence under the conservative prior.

------------------------------------------------------------------------

## Citation

If you use **multibayes** in published research, please cite:

> Calderan, M., Gambarota, F., Toffalini, E., & Altoè, G. (2026).
> Adjusting the Probability of Direction for Multiple Testing.

Pre-print: -
[pdf](https://mar-cald.github.io/multibayes/paper/manuscript.pdf) -
[html](https://mar-cald.github.io/multibayes/paper/manuscipt.md)

------------------------------------------------------------------------

## License

GNU General Public License, Version 3

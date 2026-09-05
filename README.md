# ginteff

<!-- badges: start -->
[![R-CMD-check](https://github.com/jespernwulff/ginteff/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jespernwulff/ginteff/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

R port of the Stata
[`ginteff`](https://doi.org/10.1177/1536867X231175253) command
(Marius Radean, *Stata Journal* 23(2), 2023). Computes two- and
three-way interaction effects -- via partial derivatives or first
differences -- for fitted regression models, with delta-method
standard errors.

Built as a thin wrapper around
[`marginaleffects`](https://marginaleffects.com): a `2^k`
counterfactual hypercube of interacted-variable vertices is built,
average predictions per vertex come from
`marginaleffects::avg_predictions()`, and a signed linear combination
is then formed manually using the per-vertex variance--covariance.
That last step is what lets arbitrary `vcov =` specifications
(`"HC3"`, `~cluster`, a sandwich matrix, ...) propagate through to
the final interaction-effect SE.

## Install

```r
# install.packages("remotes")
remotes::install_github("jespernwulff/ginteff")
```

`marginaleffects` (>= 0.21.0) is the only required dependency; the
current release (1.0.0) is supported and recommended. Optional:
`ggplot2`, `MASS`, `nnet`, `sandwich` for plotting / ordered &
multinomial models / robust SEs.

## Quick start

```r
library(ginteff)

m <- glm(am ~ hp * wt * disp, data = mtcars, family = binomial)

# Three-way partial-derivative interaction
ginteff(m, dydxs = c("hp", "wt", "disp"))

# Two-way first difference, +10 hp and +1 wt
ginteff(m, firstdiff = c("hp", "wt"), nunit = c(hp = 10, wt = 1))

# Robust SEs
ginteff(m, dydxs = c("hp", "wt"), vcov = "HC3")

# Per-observation interaction effects (one row per observation, one
# column per factor-level configuration). Useful for plotting the
# distribution of effects across the sample. The default
# obseff = FALSE returns only the AIE summary.
gi <- ginteff(m, dydxs = c("hp", "wt"), obseff = TRUE)
head(gi$obseff)
```

A worked introduction lives in the package vignette:

```r
vignette("ginteff-intro", package = "ginteff")
```

The vignette covers `dydxs` vs `firstdiff`, factor & continuous
interactions, robust / cluster vcovs, per-observation effects, and
the Stata-to-R argument map.

## Argument map (Stata → R)

| Stata                  | R                                                |
|------------------------|--------------------------------------------------|
| `dydxs(x1 x2)`         | `dydxs = c("x1", "x2")`                          |
| `firstdiff(x1 x2)`     | `firstdiff = c("x1", "x2")`                      |
| `fd(...)`              | `fd = ...`                                       |
| `nunit((1) x1 x2)`     | `nunit = c(x1 = 1, x2 = 1)`                      |
| `at((mean) x3)`        | `at = list(x3 = "mean")`                         |
| `atdxs((mean) x1)`     | `atdxs = list(x1 = "mean")`                      |
| `obseff(stub)`         | `obseff = TRUE` (returned in `$obseff`)          |
| `level(95)`            | `level = 0.95`                                   |
| `intequation(#1)`      | `eqn = "outcome_label"` or integer index         |
| `predict(p11)`         | `type = "..."` (model-specific)                  |
| `vce(unconditional)`   | `vcov = ~cluster`, `"HC3"`, etc.                 |
| `noweights`            | `weights = FALSE`                                |
| `post`                 | not implemented (use `as.data.frame()`)          |

## Engines

`ginteff()` chooses an engine per model:

- **Analytic** for `lm`, `glm`, `svyglm`. Computes per-vertex
  predictions and the gradient w.r.t. β analytically (via
  `family$linkinv` and `family$mu.eta`), then forms the sandwich
  `(J' c)' V (J' c)` exactly. AIE *and* SE match Stata to printed
  precision in every case we've tested.
- **`marginaleffects` fallback** for `polr`, `multinom`, and any
  other class `marginaleffects` supports. AIEs match Stata to 4-5
  sig figs; SEs match exactly on factor # factor and, with
  marginaleffects >= 1.0.0, within ~0.05% on continuous # factor
  (1-3% with older marginaleffects, whose numerical Jacobian is
  noisier). Tested against marginaleffects 0.32.0 and 1.0.0.

## Side-by-side verification

The development repository contains parallel R / Stata harnesses that
re-run the `ginteff` Stata help-file examples (NHANES II) and a set
of synthetic-data examples through both implementations. AIEs match
to 4-5 sig figs across the board; SE behavior is summarised above.
The harness scripts and logs are not shipped with the installed
package -- they live under `verification/` in the source repo.

## Citation

If you use this package, please cite the original Stata article:

> Radean, M. (2023). ginteff: A Stata command to compute interaction
> effects from generalized linear models.
> *Stata Journal*, 23(2). <https://doi.org/10.1177/1536867X231175253>

## License

MIT © 2026 Jesper Wulff (R port).
The original Stata `ginteff` command is © Marius Radean.

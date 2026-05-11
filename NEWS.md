# ginteff 0.1.1

## Bug fixes

* `ginteff()` no longer fails with `object '<var>' not found` when the
  fitted model's formula contains an inline transformation such as
  `factor(group)`, `poly(x, 2)`, `bs(x)`, `ns(x)`, `scale(x)`, or
  `log(x)`. The internal data extractor now prefers the raw fitting
  frame (via `model$data` or by re-evaluating the original `data =`
  argument) over `model.frame(model)`, which had columns named after
  the formula expression rather than the raw variable, breaking
  `predict()`'s `predvars` re-evaluation on the counterfactual grid.
  Reported by an external user; fix verified across `factor()`,
  `poly()`, `I()`, dydxs, firstdiff, and obseff code paths
  (#regression-test T9). (#1)

## Documentation

* README quick-start now demonstrates `obseff = TRUE` and explains
  that the per-observation effects matrix is returned in `gi$obseff`.
  The default (`obseff = FALSE`) returns only the AIE summary.

# ginteff 0.1.0

* Initial GitHub release.
* Two- and three-way interaction effects via partial derivatives or
  first differences, with delta-method standard errors.
* Analytic engine for `lm` / `glm` / `svyglm` (Stata-equivalent SEs).
* `marginaleffects` fallback for `polr` / `multinom` / other classes.
* Bundled `simdata.csv` and an introductory vignette.
* R-CMD-check matrix passing on macOS, Windows, Ubuntu × R-release /
  devel / oldrel-1.

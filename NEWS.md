# ginteff 0.2.1

## marginaleffects 1.0.0 compatibility

* **Factor x factor interactions on the fallback engine work again with
  marginaleffects 1.0.0.** That release skips re-attaching non-model
  columns to `avg_predictions()` output when a multi-outcome model
  (`polr`, `multinom`, ...) meets factor-level padding of `newdata` --
  which every factor x factor grid triggers, because each vertex holds
  the interacted factors at a single level. The `by` grouping was then
  silently ignored and `ginteff()` stopped with
  `argument 1 is not a vector`. The fallback engine now supplies the
  absent factor levels itself (a handful of extra rows in their own,
  discarded group), so marginaleffects never pads. Results are
  unchanged on marginaleffects 0.32.0 and match the Stata reference
  values on 1.0.0.

* **`vcov = "stata"` means HC1 on both engines and on every
  marginaleffects version.** marginaleffects 1.0.0 changed its own
  `"stata"` alias from HC2 to HC1 (Stata's `vce(robust)`); `ginteff()`
  now maps the string itself before calling marginaleffects.
  `"classical"` / `"iid"` / `"constant"` are accepted on the fallback
  engine as well.

* **User-supplied vcov matrices are aligned by coefficient name on the
  fallback engine.** marginaleffects orders `polr` parameters
  thresholds-first, unlike `vcov(model)`; before 1.0.0 it matched a
  full named matrix positionally and returned a wildly wrong SE (0.0219
  instead of 0.00018 in the new TM9 test). 1.0.0 aligns by name and
  requires the matrix to cover *all* parameters it reports (for `polr`
  including the thresholds); `ginteff()` now aligns by name on every
  version and warns about unnamed matrices on older ones.

* **Per-variable derivative steps on the fallback engine with
  marginaleffects >= 1.0.0.** The 0.2.0 scaling fix was limited to the
  analytic engine because the fallback's numerical Jacobian was too
  noisy for it. marginaleffects 1.0.0's Jacobian is accurate enough --
  forcing a `glm` through the fallback reproduces the analytic engine's
  SE to within 1% for two- and three-way cross-partials, versus 20-40%
  off with 0.32.0 -- so the fallback now uses the same per-variable
  steps; older marginaleffects versions keep the shared `max(sd)` step.
  As a by-product, fallback SEs on continuous x factor `polr` /
  `multinom` effects now land within 0.05% of Stata's (1-3% before).

* `vcov = FALSE` stops with a clear message instead of a marginaleffects
  error about a Jacobian mismatch.

* New tests TM9-TM12 cover the four points above. The suite (73
  expectations) passes on marginaleffects 0.32.0 and 1.0.0 (R 4.2.1).

# ginteff 0.2.0

## Bug fixes

* **Per-variable derivative steps in `dydxs`.** The numerical
  cross-partial previously used a single shared step `h` scaled by the
  *largest* standard deviation among the interacted continuous
  variables. When those variables lived on very different scales, the
  small-scale variable was perturbed by many of its own standard
  deviations and the local derivative was destroyed (in a test case
  with sd ratio 1e5 the AIE came out ~4,800x too small). On the
  analytic engine (`lm`/`glm`/`svyglm`) every continuous variable now
  gets its own step proportional to its own sd. The marginaleffects
  fallback keeps the legacy shared step, which is tuned to that
  engine's numerical-Jacobian noise floor. `h` now also accepts a
  named vector for explicit per-variable steps.

* **`svyglm` design weights are now applied.** The AIE for a survey
  GLM was previously an unweighted average over the sample unless
  weights were passed manually. The default `weights = NULL` now uses
  the model's own estimation weights (svyglm design weights, lm/glm
  prior weights, polr/multinom case weights), matching Stata's
  `ginteff`, which respects estimation weights unless `noweights` is
  given. Pass `weights = FALSE` for Stata's `noweights` behaviour.

* **`weights = "column"` works on the analytic engine.** The
  documented column-name form of `weights` previously stopped with
  "only numeric `weights` are supported" for `lm`/`glm`/`svyglm`.

* **Consistent complete-case handling across vertices.** The analytic
  engine filtered each counterfactual vertex by the *baseline*
  vertex's complete-case mask; if a perturbed vertex had a different
  NA pattern, rows misaligned silently (and `obseff` could error).
  All vertices now share the common complete-case mask.

## New features

* **Full joint vcov for multi-row results.** `vcov()` on a `ginteff`
  object from the analytic engine now returns the complete
  delta-method variance-covariance matrix across result rows,
  including the off-diagonal covariances between factor contrasts.
  This enables tests of the difference between two interaction
  effects (the Stata workflow `ginteff, post` + `nlcom`):
  `a' vcov(g) a` with a contrast vector `a`. The fallback engine
  still returns a diagonal matrix.

* **`ginteffplot(output = NULL)` draws all result rows** as stacked
  panels sharing a common x-axis -- the analogue of Stata's
  `ginteffplot, xcommon(*)`.

* **Warning on ill-conditioned model vcov.** When the fallback engine
  inherits a coefficient variance-covariance matrix with non-finite
  entries or a reciprocal condition number below 1e-14 (e.g. a
  `polr` Hessian with poorly scaled interacted regressors), `ginteff()`
  now warns instead of silently propagating unreliable standard
  errors.

* **`vcov` string aliases.** `"robust"` (HC3), `"stata"` (HC1), and
  `"classical"`/`"iid"`/`"constant"` are accepted on the analytic
  engine, mirroring marginaleffects' conventions.

* **Stata-style print layout.** `print()` now mirrors Stata's
  `ginteff` output -- the point estimate, its standard error, and the
  confidence interval as a single bracketed column -- so results no
  longer wrap across two blocks in consoles and rendered documents.
  `summary()` adds the z and p columns; `as.data.frame()` is unchanged.

# ginteff 0.1.2

## Bug fixes

* `ginteff()` now applies the model's `na.action` to the raw fitting
  frame before computing interaction effects. Previously, after the
  v0.1.1 fix made `.ginteff_get_data()` return the raw frame, models
  whose fit had silently dropped rows (e.g. panel data with a lagged
  predictor that's `NA` in the first period per panel) errored with
  `factor <X> has new levels <Y>`. The dropped rows held factor
  levels that never made it into `xlev`, and `model.frame(...,
  xlev = xlev)` then refused to coerce them. The fix subsets the
  raw frame by `na.action(model)` so the AIE is computed over
  exactly the rows that contributed to the fit (and `g$n` matches
  `nobs(model)`). Reported by the same external user as the v0.1.1
  bug; new regression test T10. (#2)

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

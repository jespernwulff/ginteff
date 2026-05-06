#' ginteff: Two- and Three-Way Interaction Effects
#'
#' R port of the Stata `ginteff` command (Radean, *Stata Journal* 23(2),
#' 2023). Computes two- and three-way interaction effects -- via partial
#' derivatives or first differences -- for fitted regression models,
#' with delta-method standard errors. Built on top of the
#' \pkg{marginaleffects} package
#' (<https://marginaleffects.com>).
#'
#' The main entry point is [ginteff()]. Visualisation helper is
#' [ginteffplot()] (also dispatched via `plot()`).
#'
#' @section Method:
#' For each requested interaction effect, [ginteff()] builds a
#' counterfactual hypercube of `2^k` "vertices" over the interacted
#' variables, calls [marginaleffects::avg_predictions()] to get one
#' average prediction per vertex with the user's chosen `vcov`, and
#' then forms the signed difference-of-differences linear combination
#' manually (using the per-vertex variance-covariance from
#' [stats::vcov()] on the predictions object). For partial-derivative
#' effects the result is divided by `(2h)^k` analytically; for first
#' differences the divisor is 1.
#'
#' Bypassing [marginaleffects::hypotheses()] lets us (a) apply the
#' cross-partial scaling exactly without numeric-derivative blow-up at
#' small `h`, and (b) pass through arbitrary `vcov =` specifications
#' (`"HC3"`, `~cluster`, sandwich matrices) to the final SE.
#'
#' @docType package
#' @name ginteff-package
#' @aliases ginteff-package
#' @keywords internal
"_PACKAGE"

## Suppress R CMD check NOTEs about ggplot2 NSE references in
## ginteffplot.R. The names "x" and "y" are column references inside
## aes(), not actual global variables.
utils::globalVariables(c("x", "y"))

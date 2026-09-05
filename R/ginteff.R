## ginteff.R
## R port of Stata's ginteff command (Radean, Stata Journal 23(2), 2023).
##
## Computes two- and three-way interaction effects -- via partial
## derivatives or first differences -- with delta-method standard errors,
## by orchestrating the marginaleffects package
## (https://marginaleffects.com).
##
## The job split is:
##   * ginteff() builds the counterfactual hypercube of "vertices" needed
##     to express the interaction effect as a signed linear combination
##     of average predictions.  For lm / glm / svyglm the per-vertex
##     averages and their gradients w.r.t. the coefficients are computed
##     analytically here (exact delta method).  For every other class the
##     per-vertex averages and their joint variance-covariance come from
##     marginaleffects::avg_predictions(by = ".__vertex"), and the signed
##     combination is then formed by hand so the cross-partial scaling
##     1/(2h)^k and any user vcov override propagate exactly.
##   * For per-observation effects (obseff = TRUE), point estimates are
##     formed from a single predictions() call without the delta method.
##
## Because the heavy lifting for the fallback classes -- predict adapters
## for hundreds of model classes, weights, robust/cluster SEs,
## multi-equation models, etc. -- happens inside marginaleffects, the
## supported model space is whatever marginaleffects supports.
##
## marginaleffects compatibility notes (tested on 0.32.0 and 1.0.0):
##   * 1.0.0 drops the `by` column from avg_predictions() output for
##     multi-outcome models (polr, multinom, ...) whenever it pads
##     `newdata` for absent factor levels -- exactly what a factor x
##     factor grid looks like.  .ginteff_complete_levels() sidesteps the
##     padding by supplying the absent levels itself.
##   * 1.0.0 changed vcov = "stata" from HC2 to HC1 and now aligns
##     user-supplied vcov matrices by coefficient name (older versions
##     matched positionally).  .ginteff_me_vcov_arg() makes both
##     behaviours version-independent.
##   * 1.0.0's numerical Jacobian is accurate enough that the fallback
##     can use per-variable derivative steps like the analytic engine;
##     older versions keep the shared max(sd) step.

# =====================================================================
# Public entry point
# =====================================================================

#' Two- and Three-Way Interaction Effects (ginteff)
#'
#' @param model      A fitted model object that marginaleffects can handle.
#' @param dydxs      Character vector (length 1-3) of variables whose effect
#'                   is computed via the partial derivative. Factor variables
#'                   in `dydxs` produce the discrete change from the base
#'                   level (one result per non-base level).
#' @param firstdiff  Character vector of variables whose effect is computed
#'                   via the first difference. Must be numeric.
#' @param fd         Alias for `firstdiff`.
#' @param atdxs      Named list giving values to fix continuous variables in
#'                   `dydxs`. Each element may be a number or one of
#'                   "mean", "median", "min", "max", "zero", "asobserved"
#'                   (default), or a numeric vector matching the data length.
#' @param at         Named list of fixed values for non-interacted covariates;
#'                   same syntax as `atdxs`.
#' @param nunit      Named numeric vector of unit increases for `firstdiff`
#'                   variables (default: 1 each). May also be a single scalar.
#' @param obseff     Logical. If TRUE, per-observation interaction effects
#'                   (point estimates only) are returned in `$obseff`.
#' @param type       Prediction scale, passed to marginaleffects::predictions
#'                   (e.g. "response", "link"). Defaults to the model's
#'                   marginaleffects default.
#' @param vcov       Variance-covariance specification: `TRUE` (the
#'                   model's own), a sandwich string (`"HC0"`...`"HC5"`,
#'                   `"robust"` = HC3, `"stata"` = HC1,
#'                   `"classical"`/`"iid"` = model vcov), a one-sided
#'                   `~cluster` formula, or a coefficient vcov matrix.
#'                   Strings and formulas are resolved through the sandwich
#'                   package on the analytic path and passed to
#'                   marginaleffects on the fallback path, with the same
#'                   alias meanings on both engines (marginaleffects itself
#'                   mapped `"stata"` to HC2 before version 1.0.0). A matrix
#'                   should carry coefficient names and cover every
#'                   parameter marginaleffects reports for the model -- for
#'                   `polr` that includes the threshold parameters -- and is
#'                   aligned by name. `vcov = FALSE` is not supported.
#' @param weights    Observation weights for the average. The default
#'                   (`NULL`) uses the model's own estimation weights when
#'                   it has any -- `svyglm` design weights, `lm`/`glm` prior
#'                   weights, `polr`/`multinom` case weights -- matching
#'                   Stata's `ginteff`, which respects estimation weights
#'                   unless `noweights` is given. Pass `FALSE` to force an
#'                   unweighted average (Stata's `noweights`), a numeric
#'                   vector of length `nrow(data)`, or the name of a column
#'                   in the data.
#' @param level      Confidence level (default 0.95).
#' @param data       Optional data frame; defaults to the model's data.
#' @param h          Step size(s) for the numerical cross-partial in `dydxs`
#'                   (continuous variables only). Default: each variable
#'                   gets its own step proportional to its standard
#'                   deviation (1e-4 sd for one continuous variable, 0.01 sd
#'                   for two, 0.05 sd for three) on the analytic engine and,
#'                   with marginaleffects >= 1.0.0, on the fallback engine
#'                   too. With older marginaleffects versions the fallback
#'                   uses a single shared step scaled by the largest sd
#'                   (tuned to stay above that engine's numerical-Jacobian
#'                   noise floor). Supply a scalar to use one step for all
#'                   variables, or a named vector for per-variable steps.
#' @param eqn        Outcome / equation selector for multi-outcome models
#'                   (e.g. polr, multinom). Passed to marginaleffects via
#'                   the `type` argument or by post-filtering, depending on
#'                   the model class.
#' @param ...        Extra arguments forwarded to
#'                   marginaleffects::avg_predictions().
#'
#' @return An object of class `"ginteff"`. Print/summary/coef/vcov/confint
#'         methods are provided.
#'
#' @examples
#' \dontrun{
#'   m <- glm(am ~ hp * wt * disp + factor(cyl), data = mtcars,
#'            family = binomial)
#'
#'   # Three-way partial-derivative interaction
#'   ginteff(m, dydxs = c("hp", "wt", "disp"))
#'
#'   # Two-way first difference, +10 hp and +1 wt
#'   ginteff(m, firstdiff = c("hp", "wt"), nunit = c(hp = 10, wt = 1))
#'
#'   # Mixed: cyl factor x hp continuous
#'   ginteff(m, dydxs = c("cyl", "hp"))
#'
#'   # Robust SEs
#'   ginteff(m, dydxs = c("hp", "wt"), vcov = "HC3")
#' }
#' @export
ginteff <- function(model,
                    dydxs     = NULL,
                    firstdiff = NULL,
                    fd        = NULL,
                    atdxs     = NULL,
                    at        = NULL,
                    nunit     = NULL,
                    obseff    = FALSE,
                    type      = NULL,
                    vcov      = TRUE,
                    weights   = NULL,
                    level     = 0.95,
                    data      = NULL,
                    h         = NULL,
                    eqn       = NULL,
                    ...) {

    cl <- match.call()
    .ginteff_check_marginaleffects()

    if (isFALSE(vcov))
        stop("`vcov = FALSE` is not supported: ginteff() always computes ",
             "delta-method standard errors. Use `vcov = TRUE` (the ",
             "default) or another vcov specification.", call. = FALSE)

    ## Default `type` to match Stata's `margins`/`ginteff` behaviour
    ## (per-observation inverse-link, then average). For glm/svyglm/lm
    ## that's "response"; for polr/multinom marginaleffects names the
    ## response-scale predictor "probs". Setting type explicitly also
    ## keeps avg_predictions and predictions consistent so
    ## mean(obseff) == AIE.
    if (is.null(type)) {
        type <- if (inherits(model, c("polr", "multinom", "clm", "clmm")))
                    "probs"
                else "response"
    }

    ## ---- Resolve dydxs/firstdiff aliases ---------------------------------
    if (!is.null(fd)) {
        if (!is.null(firstdiff))
            stop("Specify only one of `fd` or `firstdiff`.")
        firstdiff <- fd
    }
    if (is.null(dydxs) && is.null(firstdiff))
        stop("One of `dydxs` or `firstdiff` is required.")

    dydxs_v <- as.character(dydxs)
    fd_v    <- as.character(firstdiff)
    all_v   <- c(dydxs_v, fd_v)

    if (length(all_v) < 2L)
        stop("Specify at least 2 interacted variables in total.")
    if (length(all_v) > 3L)
        stop("Maximum of 3 interacted variables allowed.")
    if (anyDuplicated(all_v))
        stop("Variables in `dydxs` and `firstdiff` must be distinct.")

    ## ---- Data ------------------------------------------------------------
    if (is.null(data)) data <- .ginteff_get_data(model)
    miss <- setdiff(all_v, names(data))
    if (length(miss))
        stop("Variable(s) not found in data: ",
             paste(miss, collapse = ", "))

    ## ---- Identify factors and validate -----------------------------------
    is_fac <- vapply(all_v, function(v) {
        x <- data[[v]]
        is.factor(x) || is.character(x) || is.logical(x)
    }, logical(1))
    names(is_fac) <- all_v
    if (any(is_fac[fd_v]))
        stop("Factor / discrete variables are not allowed in `firstdiff`.")
    if (length(atdxs) && any(is_fac[intersect(names(atdxs), dydxs_v)]))
        stop("Factor variables in `dydxs` cannot be fixed via `atdxs`.")

    ## ---- Apply at()/atdxs() to build the evaluation data -----------------
    base_data <- .ginteff_apply_at(data, at,
                                   skip = all_v, name = "at")
    base_data <- .ginteff_apply_at(base_data, atdxs,
                                   only = dydxs_v, name = "atdxs",
                                   allow_factor = FALSE,
                                   is_factor = is_fac)

    ## ---- nunit -----------------------------------------------------------
    units <- .ginteff_resolve_units(fd_v, nunit)
    if (length(units) && any(units == 0)) {
        warning("One or more `nunit` values are zero; ",
                "the interaction effect is then trivially 0.")
    }

    ## ---- Engine ------------------------------------------------------------
    ## Analytic for lm / glm / svyglm (exact Stata-equivalent SEs);
    ## marginaleffects for everything else (polr, multinom, ...). The
    ## analytic path also handles `vcov` strings ("HC3" etc) and ~cluster
    ## formulas via the sandwich package.
    engine <- if (inherits(model, c("glm", "lm", "svyglm")) &&
                  !inherits(model, c("polr", "multinom")))
                  "analytic"
              else
                  "marginaleffects"

    if (engine == "analytic" && !is.null(eqn))
        warning("`eqn` is ignored for single-outcome models (",
                paste(class(model), collapse = "/"), ").", call. = FALSE)

    ## Fallback models inherit whatever vcov the fitter produced; surface
    ## a warning when that matrix is itself unusable (e.g. an
    ## ill-conditioned polr/multinom Hessian from poorly scaled
    ## interacted regressors) so a bad SE is not mistaken for a good one.
    if (engine == "marginaleffects")
        .ginteff_check_model_vcov(model)

    ## ---- Weights ------------------------------------------------------------
    ## NULL (default): use the model's own estimation weights when present
    ## (svyglm design weights, lm/glm prior weights, polr/multinom case
    ## weights) -- Stata's ginteff behaviour. FALSE: Stata's `noweights`.
    if (isFALSE(weights)) {
        weights <- NULL
    } else if (is.null(weights)) {
        weights <- .ginteff_default_weights(model, nrow(base_data))
    } else if (is.character(weights) && length(weights) == 1L) {
        if (!weights %in% names(base_data))
            stop("`weights` column `", weights, "` not found in data.")
        weights <- as.numeric(base_data[[weights]])
    }
    if (is.numeric(weights) && length(weights) != nrow(base_data))
        stop("`weights` length (", length(weights),
             ") must match the evaluation data (", nrow(base_data), ").")

    ## ---- h for numerical derivatives -------------------------------------
    ## Default magnitude scales with the order of differentiation k_cont:
    ##   k=1:  1e-4 * sd(x)    k=2:  0.01 * sd(x)    k=3:  0.05 * sd(x)
    ##
    ## On the analytic engine each continuous variable gets its own step
    ## proportional to its own sd. A single shared step based on max(sd)
    ## perturbs a small-scale variable by many of its own sds when the
    ## interacted variables live on very different scales, which destroys
    ## the local derivative (the pre-0.2.0 behaviour).
    ##
    ## On the marginaleffects fallback the step policy depends on the
    ## marginaleffects version.  Before 1.0.0 that engine computed the
    ## prediction Jacobian numerically with an internal step ~1e-7--1e-8,
    ## and c'V_ap c relies on differences across vertices that must stay
    ## above that noise floor (the noise is amplified by the (2h)^(-2k)
    ## divisor); the shared max(sd) step was tuned for it, and even so
    ## the k=3 SE sat ~10-20% off the exact answer. From 1.0.0 the
    ## Jacobian is accurate enough (forced-fallback SEs within ~1% of
    ## the analytic engine for k=2 and k=3, either step policy) that the
    ## per-variable steps are used there as well. Bias from finite
    ## differencing is O(h^2); at h=0.05*sd the AIE drift on a smooth
    ## response is <0.01%.  With an old marginaleffects, prefer
    ## `firstdiff = ` for 3-way effects (its variance is exact).
    cont_dydxs <- dydxs_v[!is_fac[dydxs_v]]
    h <- .ginteff_resolve_h(h, cont_dydxs, data, engine)

    ## ---- Enumerate factor-level configurations ---------------------------
    fac_dx     <- dydxs_v[is_fac[dydxs_v]]
    fac_levels <- lapply(fac_dx, function(v) {
        x <- data[[v]]
        if (is.factor(x))      lv <- levels(x)
        else if (is.logical(x)) lv <- c(FALSE, TRUE)
        else                   lv <- sort(unique(x))
        list(base = lv[1], non_base = lv[-1])
    })
    names(fac_levels) <- fac_dx

    if (length(fac_dx) == 0L) {
        configs <- list(list())
    } else {
        grid_lvl <- expand.grid(lapply(fac_levels, `[[`, "non_base"),
                                stringsAsFactors = FALSE,
                                KEEP.OUT.ATTRS  = FALSE)
        configs <- lapply(seq_len(nrow(grid_lvl)),
                          function(i) as.list(grid_lvl[i, , drop = FALSE]))
    }

    ## ---- Loop over configs, computing the AIE for each --------------------
    ##
    ## For each config we build a stacked counterfactual data frame with
    ## one copy per vertex of the 2^k hypercube over the interacted vars.
    ## avg_predictions() returns one row per vertex; hypotheses() then
    ## reduces those rows to a single signed linear combination.

    nv     <- length(all_v)
    combos <- as.matrix(expand.grid(rep(list(c(0L, 1L)), nv),
                                    KEEP.OUT.ATTRS = FALSE))
    colnames(combos) <- all_v
    signs  <- (-1L)^(nv - rowSums(combos))

    res_list <- vector("list", length(configs))
    obs_mat  <- if (isTRUE(obseff))
                    matrix(NA_real_, nrow(base_data), length(configs))
                else NULL

    for (k in seq_along(configs)) {
        config <- configs[[k]]

        cf <- .ginteff_build_grid(base_data,
                                  combos    = combos,
                                  dydxs_v   = dydxs_v,
                                  fd_v      = fd_v,
                                  is_fac    = is_fac,
                                  fac_levels = fac_levels,
                                  config    = config,
                                  units     = units,
                                  h         = h)

        divisor <- if (length(cont_dydxs)) prod(2 * h[cont_dydxs]) else 1
        scale   <- 1 / divisor

        if (engine == "analytic") {
            res <- .ginteff_compute_analytic(
                model    = model,
                cf       = cf,
                signs    = signs,
                scale    = scale,
                vcov_arg = vcov,
                type     = type,
                weights  = weights,
                level    = level,
                want_obs = isTRUE(obseff))
        } else {
            res <- .ginteff_compute_via_me(
                model    = model,
                cf       = cf,
                signs    = signs,
                scale    = scale,
                vcov_arg = vcov,
                type     = type,
                weights  = weights,
                level    = level,
                eqn      = eqn,
                base_n   = nrow(base_data),
                ...)
        }

        res_list[[k]] <- list(
            label    = .ginteff_config_label(config, dydxs_v, fd_v,
                                             fac_levels),
            estimate = res$estimate,
            se       = res$se,
            z        = res$z,
            p        = res$p,
            lower    = res$lower,
            upper    = res$upper,
            ap       = res$ap,
            grad     = res$grad,
            V_coef   = res$V_coef
        )

        if (isTRUE(obseff)) {
            if (!is.null(res$obs)) {
                ## analytic path returned per-obs values directly.
                ## Pad with NA for rows dropped during NA filtering so the
                ## column has length nrow(base_data).
                obs_mat[, k] <- .ginteff_pad_obs(res$obs,
                                                 res$keep,
                                                 nrow(base_data))
            } else {
                obs_mat[, k] <- .ginteff_obs_effects(
                    model      = model,
                    base_data  = base_data,
                    combos     = combos,
                    signs      = signs,
                    dydxs_v    = dydxs_v,
                    fd_v       = fd_v,
                    is_fac     = is_fac,
                    fac_levels = fac_levels,
                    config     = config,
                    units      = units,
                    h          = h,
                    divisor    = divisor,
                    type       = type,
                    eqn        = eqn,
                    ...)
            }
        }
    }

    ## ---- Pack results ----------------------------------------------------
    labels <- vapply(res_list, `[[`, character(1), "label")
    aie    <- vapply(res_list, `[[`, numeric(1),  "estimate")
    se     <- vapply(res_list, `[[`, numeric(1),  "se")
    z      <- vapply(res_list, `[[`, numeric(1),  "z")
    p      <- vapply(res_list, `[[`, numeric(1),  "p")
    lo     <- vapply(res_list, `[[`, numeric(1),  "lower")
    hi     <- vapply(res_list, `[[`, numeric(1),  "upper")

    method <- if (length(fd_v) == 0L) "dydxs"
              else if (length(dydxs_v) == 0L) "firstdiff"
              else "mixed"

    ## Variance-covariance of the AIE vector.  On the analytic engine the
    ## per-config gradients w.r.t. the model coefficients are available,
    ## so the full joint vcov G V G' is formed -- including the
    ## off-diagonal covariances between configs (e.g. the two non-base
    ## contrasts of a multi-level factor), which is what a test of the
    ## difference between two interaction effects needs.  The
    ## marginaleffects fallback keeps the diagonal-only matrix.
    grads <- lapply(res_list, `[[`, "grad")
    if (length(grads) && !any(vapply(grads, is.null, logical(1)))) {
        G     <- do.call(rbind, grads)
        aie_V <- G %*% res_list[[1L]]$V_coef %*% t(G)
    } else {
        aie_V <- diag(se^2, nrow = length(se))
    }
    rownames(aie_V) <- colnames(aie_V) <- labels

    out <- list(
        call       = cl,
        aie        = stats::setNames(aie, labels),
        se         = stats::setNames(se,  labels),
        z          = stats::setNames(z,   labels),
        p          = stats::setNames(p,   labels),
        ci         = `dimnames<-`(cbind(lo, hi),
                                  list(labels, c("lower", "upper"))),
        vcov_aie   = aie_V,
        level      = level,
        method     = method,
        n          = nrow(base_data),
        vars       = all_v,
        dydxs      = dydxs_v,
        firstdiff  = fd_v,
        nunit      = units,
        obseff     = if (isTRUE(obseff))
                         `colnames<-`(obs_mat, labels)
                     else NULL,
        labels     = labels,
        is_factor  = is_fac,
        fac_levels = fac_levels,
        h          = h,
        per_vertex = lapply(res_list, `[[`, "ap")
    )
    class(out) <- "ginteff"
    out
}


# =====================================================================
# Internal: analytical engine for lm / glm / svyglm
#
# Computes per-vertex average prediction and its delta-method gradient
# w.r.t. beta, then forms the signed linear combination + variance
# without going through marginaleffects (whose numerical Jacobian for
# avg_predictions has ~1e-10 noise that catastrophically cancels in the
# 3-way cross-difference and produces a ~50% SE underestimate).
# =====================================================================

.ginteff_compute_analytic <- function(model, cf, signs, scale,
                                      vcov_arg, type, weights,
                                      level, want_obs = FALSE) {

    ## inverse link + dmu/deta. lm has no $family; treat as identity.
    fam <- if (inherits(model, c("glm", "svyglm")))
               stats::family(model)
           else
               list(linkinv = function(x) x,
                    mu.eta  = function(x) rep_len(1, length(x)))
    ## type = "link" / "xb": skip the link.
    if (identical(type, "link") || identical(type, "xb")) {
        fam <- list(linkinv = function(x) x,
                    mu.eta  = function(x) rep_len(1, length(x)))
    }

    beta <- stats::coef(model)
    ## drop aliased / NA coefficients (consistent with marginaleffects)
    valid <- !is.na(beta)
    beta  <- beta[valid]
    p     <- length(beta)

    V <- .ginteff_resolve_vcov(vcov_arg, model)
    if (!isTRUE(all.equal(dim(V), c(p, p)))) {
        ## sandwich vcovs come back at full coef length; subset.
        if (!is.null(rownames(V)) && all(names(beta) %in% rownames(V))) {
            V <- V[names(beta), names(beta), drop = FALSE]
        } else {
            stop("vcov dimension does not match coef vector after ",
                 "dropping NA coefficients.")
        }
    }

    rhs_form <- stats::delete.response(stats::terms(stats::formula(model)))
    xlev     <- model$xlevels
    contr    <- model$contrasts

    vertices <- sort(unique(cf[[".__vertex"]]))
    nvert    <- length(vertices)
    avgP     <- numeric(nvert)
    Jt       <- matrix(0, nvert, p, dimnames = list(NULL, names(beta)))

    ## Build all per-vertex model frames first and take the *common*
    ## complete-case mask.  Using vertex 1's mask for every vertex (the
    ## pre-0.2.0 behaviour) silently misaligned rows -- and broke the
    ## obseff matrix -- whenever a perturbed vertex had a different NA
    ## pattern than the baseline.
    mfs <- vector("list", nvert)
    for (k in seq_len(nvert)) {
        nd_v <- cf[cf[[".__vertex"]] == vertices[k], , drop = FALSE]
        mfs[[k]] <- stats::model.frame(rhs_form, nd_v, xlev = xlev,
                                       na.action = stats::na.pass)
    }
    keep <- Reduce(`&`, lapply(mfs, stats::complete.cases))
    n    <- sum(keep)
    if (n == 0L) stop("No rows are complete cases in every vertex.")

    w <- if (is.numeric(weights)) weights[keep] else NULL
    wsum <- if (!is.null(w)) sum(w) else n

    obs_mat <- if (isTRUE(want_obs)) matrix(0, n, nvert) else NULL

    for (k in seq_len(nvert)) {
        mf <- mfs[[k]][keep, , drop = FALSE]

        X <- stats::model.matrix(rhs_form, mf, contrasts.arg = contr)
        ## Align columns to beta's order; missing cols => 0; extra cols
        ## => dropped.
        Xa <- matrix(0, n, p, dimnames = list(NULL, names(beta)))
        common <- intersect(colnames(X), names(beta))
        if (length(common) == 0L)
            stop("model.matrix produced no columns matching coef names.")
        Xa[, common] <- X[, common]

        eta <- as.numeric(Xa %*% beta)
        mu  <- fam$linkinv(eta)
        me  <- fam$mu.eta(eta)

        if (is.null(w)) {
            avgP[k] <- mean(mu)
            Jt[k, ] <- as.numeric(crossprod(Xa, me)) / n
        } else {
            avgP[k] <- sum(w * mu) / wsum
            Jt[k, ] <- as.numeric(crossprod(Xa, w * me)) / wsum
        }

        if (isTRUE(want_obs)) obs_mat[, k] <- mu
    }

    L     <- as.numeric(sum(signs * avgP)) * scale
    g     <- as.numeric(crossprod(Jt, signs)) * scale
    names(g) <- names(beta)
    var_L <- as.numeric(t(g) %*% V %*% g)
    var_L <- max(var_L, 0)
    seL   <- sqrt(var_L)
    z_q   <- stats::qnorm(1 - (1 - level) / 2)

    obs_per <- if (isTRUE(want_obs))
                   as.numeric(obs_mat %*% signs) * scale
               else NULL

    list(
        estimate = L,
        se       = seL,
        z        = if (seL > 0) L / seL else NA_real_,
        p        = if (seL > 0) 2 * stats::pnorm(-abs(L / seL)) else NA_real_,
        lower    = L - z_q * seL,
        upper    = L + z_q * seL,
        ap       = NULL,
        obs      = obs_per,
        keep     = keep,
        grad     = g,
        V_coef   = V
    )
}

# =====================================================================
# Internal: marginaleffects-based engine (for non-lm/glm classes)
# =====================================================================

.ginteff_compute_via_me <- function(model, cf, signs, scale, vcov_arg,
                                    type, weights, level, eqn,
                                    base_n, ...) {

    nvert   <- length(unique(cf[[".__vertex"]]))
    wts_rep <- .ginteff_resolve_wts(weights, NULL, nvert = nvert,
                                    base_n = base_n)

    ## Supply absent factor levels ourselves (tagged .__vertex = 0) so
    ## marginaleffects never pads `newdata`; see .ginteff_complete_levels().
    comp <- .ginteff_complete_levels(cf)
    cf   <- comp$data
    if (comp$n_extra > 0L && !is.null(wts_rep))
        wts_rep <- c(wts_rep, rep(1, comp$n_extra))

    ap_args <- list(
        model,
        newdata    = cf,
        by         = ".__vertex",
        vcov       = .ginteff_me_vcov_arg(vcov_arg, model),
        type       = type,
        conf_level = level,
        ...)
    if (!is.null(wts_rep)) ap_args$wts <- wts_rep
    ap <- do.call(marginaleffects::avg_predictions, ap_args)

    if (is.null(ap[[".__vertex"]]))
        stop("marginaleffects::avg_predictions() did not return the ",
             "'.__vertex' grouping column, so the per-vertex averages ",
             "cannot be recovered (marginaleffects ",
             as.character(utils::packageVersion("marginaleffects")), ", ",
             paste(class(model), collapse = "/"), " model).",
             call. = FALSE)

    ## Per-vertex vcov *before* any subsetting: `[` on a predictions
    ## object does not reliably carry the attribute along.
    V_full <- stats::vcov(ap)
    if (!is.matrix(V_full) || !all(dim(V_full) == nrow(ap)))
        stop("Could not recover the per-vertex variance-covariance ",
             "matrix from marginaleffects::avg_predictions().",
             call. = FALSE)

    idx  <- which(ap[[".__vertex"]] %in% seq_len(nvert))
    keep <- .ginteff_filter_eqn(ap, eqn)
    if (!is.null(keep)) idx <- intersect(idx, keep)
    if (length(idx) != nvert)
        stop("Expected ", nvert, " per-vertex averages from ",
             "marginaleffects::avg_predictions() but found ", length(idx),
             ".", call. = FALSE)
    ord <- idx[order(ap[[".__vertex"]][idx])]

    est   <- as.numeric(ap$estimate[ord])
    V_ap  <- V_full[ord, ord, drop = FALSE]

    c_vec  <- as.numeric(signs)
    L      <- as.numeric(crossprod(c_vec, est)) * scale
    varL   <- as.numeric(crossprod(c_vec, V_ap %*% c_vec)) * scale^2
    varL   <- max(varL, 0)
    seL    <- sqrt(varL)
    zq     <- stats::qnorm(1 - (1 - level) / 2)

    list(
        estimate = L,
        se       = seL,
        z        = if (seL > 0) L / seL else NA_real_,
        p        = if (seL > 0) 2 * stats::pnorm(-abs(L / seL)) else NA_real_,
        lower    = L - zq * seL,
        upper    = L + zq * seL,
        ap       = ap[ord, , drop = FALSE],
        obs      = NULL,
        keep     = NULL
    )
}

# =====================================================================
# Internal: complete absent factor levels in the stacked grid.
#
# marginaleffects pads `newdata` with one row per level of every factor
# column that does not show all of its levels, predicts on the padded
# frame, and drops the padding afterwards.  In marginaleffects 1.0.0
# the step that re-attaches non-model columns (our `.__vertex`) to the
# predictions is skipped when the model returns several rows per
# observation (polr, multinom, ...) *and* padding took place -- the
# padding rows share one rowid, which trips the duplicate check in
# merge_original_data().  A factor x factor grid always triggers this
# because each vertex holds the interacted factors at a single level.
#
# Supplying the absent levels ourselves -- one copy of the first row per
# missing level, tagged `.__vertex = 0` -- means nothing is missing from
# marginaleffects' point of view, so it does not pad.  The extra rows
# form their own `by` group, which the caller discards; they cannot
# touch the real per-vertex averages or their covariances.  The 100-level
# ceiling mirrors marginaleffects' own (it skips padding above it).
# =====================================================================

.ginteff_complete_levels <- function(cf) {
    extra <- list()
    total <- 0L
    for (v in names(cf)) {
        x <- cf[[v]]
        if (!is.factor(x)) next
        lv   <- levels(x)
        miss <- setdiff(lv, as.character(unique(x)))
        if (!length(miss)) next
        total <- total + length(lv)
        rows  <- cf[rep(1L, length(miss)), , drop = FALSE]
        rows[[v]] <- factor(miss, levels = lv, ordered = is.ordered(x))
        extra[[v]] <- rows
    }
    if (!length(extra) || total > 100L)
        return(list(data = cf, n_extra = 0L))
    extra <- do.call(rbind, extra)
    extra[[".__vertex"]] <- 0L
    out <- rbind(cf, extra)
    rownames(out) <- NULL
    list(data = out, n_extra = nrow(extra))
}

# =====================================================================
# Internal: normalise a `vcov` argument for marginaleffects.
#
# * "stata" -> "HC1" so the alias means the same thing on both engines
#   and across marginaleffects versions (marginaleffects mapped it to
#   HC2 before 1.0.0).  "classical"/"iid"/"constant" -> TRUE.
# * A named matrix is reordered to marginaleffects' own coefficient
#   order (get_coef()), which for polr puts the thresholds first and
#   therefore differs from rownames(vcov(model)).  marginaleffects >= 1.0.0
#   aligns by name itself; older versions matched positionally and
#   silently produced wrong standard errors for such a matrix.
# =====================================================================

.ginteff_me_vcov_arg <- function(vcov_arg, model) {
    if (is.character(vcov_arg) && length(vcov_arg) == 1L) {
        key <- tolower(vcov_arg)
        if (key == "stata") return("HC1")
        if (key %in% c("classical", "constant", "iid")) return(TRUE)
        return(vcov_arg)
    }
    if (is.matrix(vcov_arg)) {
        nm <- tryCatch(names(marginaleffects::get_coef(model)),
                       error = function(e) NULL)
        rn <- rownames(vcov_arg)
        cn <- colnames(vcov_arg)
        if (is.null(rn) || is.null(cn)) {
            if (!.ginteff_me_version_ge("1.0.0"))
                warning("`vcov` is an unnamed matrix and will be matched ",
                        "positionally against marginaleffects' coefficient ",
                        "order",
                        if (!is.null(nm))
                            paste0(" (", paste(nm, collapse = ", "), ")"),
                        ". Supply dimnames to have it aligned by name.",
                        call. = FALSE)
        } else if (!is.null(nm) && all(nm %in% rn) && all(nm %in% cn)) {
            vcov_arg <- vcov_arg[nm, nm, drop = FALSE]
        }
    }
    vcov_arg
}

# =====================================================================
# Internal: installed marginaleffects version >= `v`?
# =====================================================================

.ginteff_me_version_ge <- function(v) {
    ver <- tryCatch(utils::packageVersion("marginaleffects"),
                    error = function(e) NULL)
    !is.null(ver) && ver >= v
}

# =====================================================================
# Internal: resolve a `vcov` argument to a numeric matrix.
# =====================================================================

.ginteff_resolve_vcov <- function(vcov_arg, model) {
    if (is.matrix(vcov_arg)) return(vcov_arg)
    if (is.null(vcov_arg) || isTRUE(vcov_arg))
        return(stats::vcov(model))
    if (is.character(vcov_arg)) {
        ## Aliases follow marginaleffects' conventions so the same string
        ## means the same thing on both engines.
        key <- tolower(vcov_arg)
        if (key %in% c("classical", "constant", "iid"))
            return(stats::vcov(model))
        if (!requireNamespace("sandwich", quietly = TRUE))
            stop("vcov = '", vcov_arg, "' requires the 'sandwich' package.")
        if (key == "robust") return(sandwich::vcovHC(model))
        if (key == "stata")  return(sandwich::vcovHC(model, type = "HC1"))
        return(sandwich::vcovHC(model, type = vcov_arg))
    }
    if (inherits(vcov_arg, "formula")) {
        if (!requireNamespace("sandwich", quietly = TRUE))
            stop("vcov = formula requires the 'sandwich' package.")
        return(sandwich::vcovCL(model, cluster = vcov_arg))
    }
    stop("Unsupported `vcov` specification.")
}

# =====================================================================
# Internal: per-variable h resolution for dydxs cross-partials.
#
# Returns a named numeric vector over the continuous dydxs variables.
# Defaults scale each variable's step by its own sd on the analytic
# engine and on the marginaleffects fallback from marginaleffects 1.0.0;
# older marginaleffects versions keep the legacy shared max(sd) step
# (tuned to their numerical-Jacobian noise floor).
# =====================================================================

.ginteff_resolve_h <- function(h, cont_dydxs, data, engine) {
    k_cont <- length(cont_dydxs)
    if (k_cont == 0L)
        return(if (is.null(h)) 1e-4 else as.numeric(h)[1L])

    if (!is.null(h)) {
        if (length(h) == 1L && is.null(names(h)))
            return(stats::setNames(rep(as.numeric(h), k_cont), cont_dydxs))
        if (is.null(names(h)))
            stop("`h` must be a single number or a named vector.")
        bad <- setdiff(names(h), cont_dydxs)
        if (length(bad))
            stop("`h` names not among the continuous dydxs variables: ",
                 paste(bad, collapse = ", "))
        out <- .ginteff_resolve_h(NULL, cont_dydxs, data, engine)
        out[names(h)] <- as.numeric(h)
        return(out)
    }

    sds <- vapply(cont_dydxs,
                  function(v) {
                      x <- as.numeric(data[[v]])
                      s <- stats::sd(x, na.rm = TRUE)
                      if (!is.finite(s) || s == 0) 1 else s
                  }, numeric(1))
    c_k <- switch(min(k_cont, 3L),
                  `1` = 1e-4,
                  `2` = 0.01,
                  `3` = 0.05)
    if (identical(engine, "analytic") || .ginteff_me_version_ge("1.0.0"))
        stats::setNames(c_k * sds, cont_dydxs)
    else
        stats::setNames(rep(c_k * max(sds), k_cont), cont_dydxs)
}

# =====================================================================
# Internal: the model's own estimation weights (Stata parity).
#
# Returns a numeric weights vector aligned with the evaluation data, or
# NULL when the model is unweighted / weights cannot be recovered /
# the length does not match.  Uniform weights collapse to NULL (the
# unweighted fast path gives the same average).
# =====================================================================

.ginteff_default_weights <- function(model, n_expected) {
    w <- tryCatch(stats::weights(model), error = function(e) NULL)
    if (is.null(w) || !length(w))
        w <- tryCatch(stats::model.weights(stats::model.frame(model)),
                      error = function(e) NULL)
    if (is.null(w) || !is.numeric(w) || !length(w)) return(NULL)
    if (anyNA(w)) w <- w[!is.na(w)]
    if (length(w) != n_expected) return(NULL)
    if (isTRUE(all(w == w[1L]))) return(NULL)
    as.numeric(w)
}

# =====================================================================
# Internal: warn when a fallback model's own vcov is unusable.
# =====================================================================

.ginteff_check_model_vcov <- function(model) {
    V <- tryCatch(stats::vcov(model), error = function(e) NULL)
    if (is.null(V) || !is.matrix(V)) return(invisible(NULL))
    dV <- diag(V)
    bad_finite <- !all(is.finite(dV))
    rc <- if (bad_finite) NA_real_
          else tryCatch(rcond(V), error = function(e) NA_real_)
    if (bad_finite || (is.finite(rc) && rc < 1e-14)) {
        warning("The model's coefficient variance-covariance matrix ",
                "appears ill-conditioned",
                if (bad_finite) " (non-finite entries)"
                else sprintf(" (reciprocal condition number %.1e)", rc),
                ". Interaction-effect standard errors inherit it and may ",
                "be unreliable. Consider rescaling the interacted ",
                "variables or refitting with a better-conditioned ",
                "implementation (e.g. ordinal::clm for ordered models).",
                call. = FALSE)
    }
    invisible(NULL)
}

# =====================================================================
# Internal: pad a per-obs vector back to nrow(base_data) length, filling
# NA for rows that were dropped during model-frame NA filtering.
# =====================================================================

.ginteff_pad_obs <- function(obs_per, keep, n_full) {
    if (is.null(keep) || length(obs_per) == n_full) return(obs_per)
    out <- rep(NA_real_, n_full)
    out[keep] <- obs_per
    out
}

# =====================================================================
# Internal: build the counterfactual grid (one config)
# =====================================================================

.ginteff_build_grid <- function(base_data, combos, dydxs_v, fd_v,
                                is_fac, fac_levels, config, units, h) {

    n     <- nrow(base_data)
    nvert <- nrow(combos)

    out <- vector("list", nvert)
    for (i in seq_len(nvert)) {
        nd <- base_data
        for (v in c(dydxs_v, fd_v)) {
            hi <- combos[i, v] == 1L
            x  <- nd[[v]]
            if (v %in% fd_v) {
                ## first difference: low = obs, high = obs + nunit
                if (hi) nd[[v]] <- x + units[[v]]
            } else if (is_fac[[v]]) {
                ## factor in dydxs: low = base, high = config[[v]]
                target <- if (hi) config[[v]] else fac_levels[[v]]$base
                if (is.factor(x)) {
                    nd[[v]] <- factor(rep(as.character(target), n),
                                      levels = levels(x))
                } else if (is.logical(x)) {
                    nd[[v]] <- rep(as.logical(target), n)
                } else {
                    nd[[v]] <- rep(target, n)
                }
            } else {
                ## continuous in dydxs: central diff around obs, using
                ## that variable's own step
                h_v <- if (!is.null(names(h)) && v %in% names(h)) h[[v]]
                       else h[[1L]]
                nd[[v]] <- x + (if (hi) h_v else -h_v)
            }
        }
        nd[[".__vertex"]] <- i
        out[[i]] <- nd
    }

    do.call(rbind, out)
}


# =====================================================================
# Internal: filter avg_predictions output to a single outcome category
# (for polr/multinom). Returns row indices to keep, or NULL if no
# filtering needed.
# =====================================================================

.ginteff_filter_eqn <- function(ap, eqn) {
    if (!"group" %in% names(ap)) return(NULL)
    grps <- as.character(ap$group)
    uniq <- unique(grps)
    if (length(uniq) == 1L) return(NULL)
    if (is.null(eqn))
        stop("Model produces predictions for multiple outcomes (",
             paste(uniq, collapse = ", "),
             "); specify which one with `eqn = `.")
    target <- if (is.numeric(eqn)) uniq[as.integer(eqn)] else as.character(eqn)
    if (!target %in% uniq)
        stop("`eqn = ", target, "` not in the model's outcome groups: ",
             paste(uniq, collapse = ", "))
    which(grps == target)
}

# =====================================================================
# Internal: per-observation IE values (point estimates only)
# =====================================================================
#
# For obseff = TRUE we want a length-n vector of IE_i values.  We could
# call predictions() once on the stacked grid, but to also handle the
# `wts` argument cleanly we just reshape the result manually.

.ginteff_obs_effects <- function(model, base_data, combos, signs,
                                 dydxs_v, fd_v, is_fac, fac_levels,
                                 config, units, h, divisor,
                                 type, eqn = NULL, ...) {
    cf <- .ginteff_build_grid(base_data,
                              combos     = combos,
                              dydxs_v    = dydxs_v,
                              fd_v       = fd_v,
                              is_fac     = is_fac,
                              fac_levels = fac_levels,
                              config     = config,
                              units      = units,
                              h          = h)

    p <- marginaleffects::predictions(model, newdata = cf,
                                      type = type, vcov = FALSE, ...)
    if ("group" %in% names(p)) {
        uniq <- unique(as.character(p$group))
        if (length(uniq) > 1L) {
            if (is.null(eqn))
                stop("Multi-outcome model: specify `eqn = `.")
            target <- if (is.numeric(eqn)) uniq[as.integer(eqn)]
                      else as.character(eqn)
            p <- p[as.character(p$group) == target, , drop = FALSE]
        }
    }
    ## The reshape below assumes rows arrive in cf (vertex-major) order;
    ## restore it explicitly in case predictions() sorted by group.
    if ("rowid" %in% names(p)) p <- p[order(p$rowid), , drop = FALSE]
    n     <- nrow(base_data)
    nvert <- nrow(combos)
    ## Each vertex contributes n rows in order; reshape to n x nvert.
    M <- matrix(p$estimate, nrow = n, ncol = nvert, byrow = FALSE)
    as.numeric((M %*% signs) / divisor)
}


# =====================================================================
# Internal: data extraction
#
# Prefer the *raw* fitting data (with untransformed columns) over
# `model.frame(model)` (which contains post-transformation columns named
# after the formula expression, e.g. `factor(group)` instead of `group`).
# Using the raw frame matters whenever the user's formula contains an
# inline transform -- factor(), poly(), I(x^2), bs(), ns(), scale(),
# log(), etc. With a post-transformation frame, downstream `predict()`
# calls re-evaluate `predvars` against our counterfactual newdata and
# fail with `object '<raw_var>' not found`, and any user reference to
# `<raw_var>` in `dydxs`/`firstdiff` also breaks because the column has
# been renamed.
#
# Resolution order:
#   1. `model$data` -- some fitters (lme4, glmmTMB, ...) cache the raw
#      fitting frame here.
#   2. Re-evaluate the `data =` argument from `model$call` in the
#      formula's environment -- works for lm/glm/polr/multinom/svyglm
#      whenever the original data frame is still in scope.
#   3. Fall back to `model.frame(model)` -- post-transformation frame
#      that works only for plain (no inline transform) formulas.
# =====================================================================

.ginteff_get_data <- function(model) {

    raw <- NULL

    ## (1) Some fitters store the raw frame on the object.
    d <- model$data
    if (!is.null(d) && is.data.frame(d) && nrow(d) > 0L)
        raw <- as.data.frame(d)

    ## (2) Re-evaluate the original `data =` argument.
    if (is.null(raw)) {
        d <- tryCatch(eval(stats::getCall(model)$data,
                           envir = environment(stats::formula(model))),
                      error = function(e) NULL)
        if (!is.null(d) && is.data.frame(d) && nrow(d) > 0L)
            raw <- as.data.frame(d)
    }

    if (!is.null(raw)) {
        ## Apply the same na.action that lm/glm/polr/multinom did
        ## during fitting.  Otherwise rows the model never saw stay
        ## in the frame and (a) `model.matrix(..., xlev = xlev)` errors
        ## with "factor X has new levels Y" if the dropped rows held
        ## levels that didn't survive into `xlev`, and (b) the AIE is
        ## averaged over a different sample than the fit.  Common
        ## trigger: panel data with lagged predictors -- the first
        ## period per panel has NA on the lag, glm drops those rows,
        ## but the corresponding factor(year) levels never make it
        ## into `xlev`.
        omit <- stats::na.action(model)
        if (!is.null(omit)) {
            idx <- as.integer(omit)
            if (length(idx) && all(idx >= 1L & idx <= nrow(raw)))
                raw <- raw[-idx, , drop = FALSE]
        }
        return(raw)
    }

    ## (3) Last resort: post-transformation model frame.  This will not
    ## carry raw columns through inline transforms, so users with
    ## `factor(group)` / `poly(x,2)` / `I(x^2)` in the formula can hit
    ## "object 'group' not found" downstream.  Fix: pass `data = ` to
    ## ginteff() explicitly, or re-fit with the transform pre-evaluated
    ## as a column.
    d <- tryCatch(stats::model.frame(model), error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0L)
        return(as.data.frame(d))

    stop("Could not extract data from `model`. Pass `data = ...`.")
}

# =====================================================================
# Internal: at()/atdxs() application
# =====================================================================

.ginteff_apply_at <- function(data, spec, skip = character(),
                              only = NULL, name = "at",
                              allow_factor = TRUE,
                              is_factor = NULL) {
    if (is.null(spec)) return(data)
    if (!is.list(spec) || is.null(names(spec)))
        stop("`", name, "` must be a named list.")
    for (v in names(spec)) {
        if (v %in% skip)
            stop("`", name, "` cannot set interacted variable `", v, "`.")
        if (!is.null(only) && !v %in% only)
            stop("`", name, "` may only set variables in `dydxs`. ",
                 "Offending: `", v, "`.")
        if (!v %in% names(data))
            stop("Variable `", v, "` not found in data.")
        x <- data[[v]]
        is_fac_v <- if (!is.null(is_factor)) isTRUE(is_factor[[v]])
                    else is.factor(x) || is.character(x) || is.logical(x)
        if (is_fac_v && !allow_factor)
            stop("`", name, "` cannot fix factor variable `", v, "`.")

        new <- .ginteff_resolve_at_value(x, spec[[v]], v, name)
        if (length(new) == 1L) {
            if (is.factor(x))
                data[[v]] <- factor(rep(as.character(new), nrow(data)),
                                    levels = levels(x))
            else
                data[[v]] <- rep(new, nrow(data))
        } else if (length(new) == nrow(data)) {
            if (is.factor(x))
                data[[v]] <- factor(as.character(new), levels = levels(x))
            else
                data[[v]] <- new
        } else {
            stop("`", name, "` for `", v, "` has length ", length(new),
                 "; expected 1 or ", nrow(data), ".")
        }
    }
    data
}

.ginteff_resolve_at_value <- function(x, val, v, name) {
    if (is.numeric(val) || is.logical(val)) return(val)
    if (is.factor(val)) return(as.character(val))
    if (is.character(val) && length(val) == 1L &&
        val %in% c("mean", "median", "min", "max", "zero",
                   "asobserved", "base")) {
        return(switch(val,
                      mean       = mean(as.numeric(x), na.rm = TRUE),
                      median     = stats::median(as.numeric(x), na.rm = TRUE),
                      min        = min(as.numeric(x), na.rm = TRUE),
                      max        = max(as.numeric(x), na.rm = TRUE),
                      zero       = 0,
                      asobserved = x,
                      base       = if (is.factor(x)) levels(x)[1]
                                   else if (is.logical(x)) FALSE
                                   else sort(unique(x))[1]))
    }
    if (is.character(val)) return(val)
    stop("Unrecognised value in `", name, "$", v, "`.")
}

# =====================================================================
# Internal: nunit resolution
# =====================================================================

.ginteff_resolve_units <- function(fd_v, nunit) {
    if (length(fd_v) == 0L) return(numeric(0))
    out <- stats::setNames(rep(1, length(fd_v)), fd_v)
    if (is.null(nunit)) return(out)

    if (is.numeric(nunit) && length(nunit) == 1L && is.null(names(nunit))) {
        out[] <- nunit
        return(out)
    }
    if (is.null(names(nunit)))
        stop("`nunit` must be named, or a single scalar applied to all.")
    miss <- setdiff(names(nunit), fd_v)
    if (length(miss))
        stop("`nunit` names not in firstdiff: ",
             paste(miss, collapse = ", "))
    out[names(nunit)] <- as.numeric(nunit)
    out
}

# =====================================================================
# Internal: weights resolution -- pass to marginaleffects::wts argument
# =====================================================================

.ginteff_resolve_wts <- function(weights, data, nvert = 1L, base_n = NULL) {
    if (is.null(weights)) return(NULL)
    if (is.character(weights) && length(weights) == 1L) return(weights)
    if (is.numeric(weights)) {
        n_target <- if (!is.null(base_n)) base_n else nrow(data)
        if (length(weights) != n_target)
            stop("`weights` length must match nrow of evaluation data.")
        return(rep(weights, times = nvert))
    }
    weights
}

# =====================================================================
# Internal: ensure marginaleffects is installed
# =====================================================================

.ginteff_check_marginaleffects <- function() {
    if (!requireNamespace("marginaleffects", quietly = TRUE)) {
        stop("ginteff() requires the 'marginaleffects' package. ",
             "Install with: install.packages('marginaleffects')",
             call. = FALSE)
    }
}

# =====================================================================
# Internal: human-readable label for a factor-level configuration
# =====================================================================

.ginteff_config_label <- function(config, dydxs_v, fd_v, fac_levels) {
    parts <- character(0)
    for (v in dydxs_v) {
        if (v %in% names(config)) {
            base <- fac_levels[[v]]$base
            parts <- c(parts,
                       sprintf("d%s[%s>%s]", v,
                               as.character(base),
                               as.character(config[[v]])))
        } else {
            parts <- c(parts, sprintf("d%s", v))
        }
    }
    for (v in fd_v) parts <- c(parts, sprintf("FD(%s)", v))
    paste(parts, collapse = " # ")
}

# =====================================================================
# Methods: print / summary / coef / vcov / confint
# =====================================================================

## Shared formatter: one row per effect, the CI as a single bracketed
## column so the table fits a standard console / rendered code block
## without wrapping.  print() mirrors Stata's ginteff output (Statistic,
## Std. Err., CI); summary() adds the z and p columns.
.ginteff_format_table <- function(aie, se, z, p, ci, labels, level,
                                  digits = 4, zp = FALSE) {
    num <- function(v) formatC(v, digits = digits, format = "g")
    tab <- data.frame(
        AIE = num(aie),
        SE  = num(se),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    if (isTRUE(zp)) {
        tab$z       <- formatC(z, digits = 3, format = "f")
        tab$`P>|z|` <- formatC(p, digits = 3, format = "f")
    }
    tab$ci <- sprintf("[%s, %s]", num(ci[, 1]), num(ci[, 2]))
    names(tab)[names(tab) == "ci"] <- sprintf("%.0f%% CI", 100 * level)
    rownames(tab) <- labels
    tab
}

#' @export
print.ginteff <- function(x, digits = 4, ...) {
    cat("\nInteraction Effects (ginteff)\n")
    cat(strrep("-", 60), "\n", sep = "")
    cat("Method   : ", x$method, "\n", sep = "")
    cat("Variables: ", paste(x$vars, collapse = ", "), "\n", sep = "")
    if (length(x$firstdiff)) {
        u <- x$nunit
        cat("nunit    : ",
            paste(names(u), "=", u, collapse = ", "),
            "\n", sep = "")
    }
    cat("N        : ", x$n, "\n", sep = "")
    cat("\nAverage interaction effect(s):\n")
    print(.ginteff_format_table(x$aie, x$se, x$z, x$p, x$ci,
                                x$labels, x$level, digits))
    invisible(x)
}

#' @export
summary.ginteff <- function(object, ...) {
    out <- list(
        call    = object$call,
        table   = data.frame(
            AIE   = object$aie,
            SE    = object$se,
            z     = object$z,
            p     = object$p,
            lower = object$ci[, 1],
            upper = object$ci[, 2],
            row.names = object$labels,
            check.names = FALSE
        ),
        n      = object$n,
        method = object$method,
        vars   = object$vars,
        level  = object$level
    )
    class(out) <- "summary.ginteff"
    out
}

#' @export
print.summary.ginteff <- function(x, digits = 4, ...) {
    cat("\nCall:\n  "); print(x$call)
    cat("\nMethod : ", x$method, "\n", sep = "")
    cat("N      : ", x$n,      "\n", sep = "")
    cat("Level  : ", x$level,  "\n\n", sep = "")
    print(.ginteff_format_table(x$table$AIE, x$table$SE, x$table$z,
                                x$table$p,
                                cbind(x$table$lower, x$table$upper),
                                rownames(x$table), x$level, digits,
                                zp = TRUE))
    invisible(x)
}

#' @export
coef.ginteff <- function(object, ...) object$aie

#' Variance-covariance matrix of the AIE vector
#'
#' Returns the delta-method variance-covariance matrix of the average
#' interaction effects stored in `object$aie`. On the analytic engine
#' (`lm` / `glm` / `svyglm`) this is the full joint matrix -- the
#' off-diagonal entries are the covariances between rows (e.g. between
#' the non-base contrasts of a multi-level factor), so a test of the
#' difference between two interaction effects is
#' `(a' vcov(g) a)` with the appropriate contrast vector `a`. On the
#' marginaleffects fallback (`polr`, `multinom`, ...) the matrix is
#' diagonal. Useful for chaining with downstream tools (e.g.
#' \pkg{multcomp}, \pkg{car}::deltaMethod).
#'
#' @param object a `ginteff` object
#' @param ... ignored
#' @export
vcov.ginteff <- function(object, ...) object$vcov_aie

#' @export
confint.ginteff <- function(object, parm, level = 0.95, ...) {
    z <- stats::qnorm(1 - (1 - level) / 2)
    ci <- cbind(lower = object$aie - z * object$se,
                upper = object$aie + z * object$se)
    rownames(ci) <- object$labels
    if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
    ci
}

#' Convert a ginteff object to a tidy data frame
#' @param x a `ginteff` object
#' @param ... ignored
#' @param row.names ignored
#' @param optional ignored
#' @export
as.data.frame.ginteff <- function(x, row.names = NULL, optional = FALSE, ...) {
    data.frame(
        term     = x$labels,
        estimate = x$aie,
        std.error = x$se,
        statistic = x$z,
        p.value  = x$p,
        conf.low = x$ci[, 1],
        conf.high = x$ci[, 2],
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

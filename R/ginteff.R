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
##     of average predictions, calls marginaleffects::avg_predictions(),
##     and then asks marginaleffects::hypotheses() to evaluate that
##     combination (which automatically carries the delta-method SE
##     through any vcov override the user supplied).
##   * For per-observation effects (obseff = TRUE), point estimates are
##     formed from a single predictions() call without the delta method.
##
## Because all heavy lifting -- factor handling, predict adapters for
## hundreds of model classes, weights, robust/cluster SEs, multi-equation
## models, etc. -- happens inside marginaleffects, this file is short
## and the supported model space is whatever marginaleffects supports.

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
#' @param vcov       Variance-covariance specification, passed straight to
#'                   marginaleffects (e.g. TRUE, "HC3", ~cluster, a matrix).
#' @param weights    Observation weights for the average, passed as `wts`.
#'                   Either a numeric vector, the name of a column in the
#'                   data, or NULL.
#' @param level      Confidence level (default 0.95).
#' @param data       Optional data frame; defaults to the model's data.
#' @param h          Step size for numerical differentiation in `dydxs`
#'                   (continuous vars). Default sqrt(.Machine$double.eps)
#'                   times a heuristic scale.
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

    ## ---- h for numerical derivatives -------------------------------------
    ## Default scales with the order of differentiation k_cont:
    ##   k=1 (single continuous dydxs var):  1e-4 * max sd(x)
    ##   k=2 (typical 2-way cross-partial):  0.01 * max sd(x)
    ##   k=3 (3-way cross-partial):           0.05 * max sd(x)
    ##
    ## The larger step at higher k is necessary because marginaleffects
    ## computes the prediction Jacobian numerically with an internal step
    ## ~1e-7--1e-8; c'V_ap c relies on small differences across vertices
    ## that fall below that noise floor when h is small, and the noise is
    ## amplified by (2h)^(-2k). Bias from finite differencing is O(h^2);
    ## at h=0.05*sd the AIE drift on a smooth response is <0.01%.
    ##
    ## A residual quirk: for k=3 the SE in the stable plateau tends to be
    ## ~10% lower than Stata's analytic answer. This is a finite-difference
    ## artifact in the c'V_ap c term. If you need analytical-quality SEs
    ## for 3-way effects, prefer the FD path (`firstdiff = `) -- the FD
    ## variance is exact (no h-divisor amplifies the noise).
    cont_dydxs <- dydxs_v[!is_fac[dydxs_v]]
    if (is.null(h)) {
        k_cont <- length(cont_dydxs)
        if (k_cont) {
            sds <- vapply(cont_dydxs,
                          function(v) {
                              x <- as.numeric(data[[v]])
                              s <- stats::sd(x, na.rm = TRUE)
                              if (!is.finite(s) || s == 0) 1 else s
                          }, numeric(1))
            scale_sd <- max(sds)
            h <- switch(min(k_cont, 3L),
                        `1` = 1e-4 * scale_sd,
                        `2` = 0.01 * scale_sd,
                        `3` = 0.05 * scale_sd)
        } else {
            h <- 1e-4
        }
    }

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

    ## Pick engine: analytic for lm / glm / svyglm (exact Stata-equivalent
    ## SEs); marginaleffects for everything else (polr, multinom, ...).
    ## The analytic path also handles `vcov` strings ("HC3" etc) and
    ## ~cluster formulas via the sandwich package.
    engine <- if (inherits(model, c("glm", "lm", "svyglm")) &&
                  !inherits(model, c("polr", "multinom")))
                  "analytic"
              else
                  "marginaleffects"

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

        cont_dx_count <- sum(!is_fac[dydxs_v])
        divisor       <- if (cont_dx_count > 0L) (2 * h)^cont_dx_count else 1
        scale         <- 1 / divisor

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
            ap       = res$ap
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

    ## Diagonal vcov of the AIE vector.  The off-diagonal cross-config
    ## covariances would require joint Jacobians of all configs against
    ## the model coefficients; for the typical case (one config = pure
    ## continuous, or one config per non-base factor level evaluated
    ## independently) the diagonal is what users reach for via vcov().
    aie_V <- diag(se^2, nrow = length(se))
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

    ## Determine which rows of base_data survive NA filtering.  We use
    ## the first vertex (no perturbation that introduces NAs) as the
    ## reference; FD vertices that add a unit shouldn't change which
    ## rows have NAs in the model variables.
    rows_kept <- NULL
    obs_mat   <- if (isTRUE(want_obs)) matrix(NA_real_, 0, nvert) else NULL

    for (k in seq_len(nvert)) {
        nd_v <- cf[cf[[".__vertex"]] == vertices[k], , drop = FALSE]
        mf <- stats::model.frame(rhs_form, nd_v, xlev = xlev,
                                 na.action = stats::na.pass)
        keep <- stats::complete.cases(mf)
        if (k == 1L) rows_kept <- keep   # for obseff padding later
        mf  <- mf[keep, , drop = FALSE]
        n   <- nrow(mf)
        if (n == 0L) stop("No complete-case rows in vertex ", k, ".")

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

        if (is.null(weights)) {
            avgP[k] <- mean(mu)
            Jt[k, ] <- as.numeric(crossprod(Xa, me)) / n
        } else {
            w <- if (is.numeric(weights)) weights[which(rows_kept)] else NULL
            if (is.null(w))
                stop("Analytic engine: only numeric `weights` are supported.")
            wsum <- sum(w)
            avgP[k] <- sum(w * mu) / wsum
            Jt[k, ] <- as.numeric(crossprod(Xa, w * me)) / wsum
        }

        if (isTRUE(want_obs)) {
            if (k == 1L)
                obs_mat <- matrix(0, n, nvert)
            obs_mat[, k] <- mu
        }
    }

    L     <- as.numeric(sum(signs * avgP)) * scale
    g     <- as.numeric(crossprod(Jt, signs)) * scale
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
        keep     = rows_kept
    )
}

# =====================================================================
# Internal: marginaleffects-based engine (for non-lm/glm classes)
# =====================================================================

.ginteff_compute_via_me <- function(model, cf, signs, scale, vcov_arg,
                                    type, weights, level, eqn,
                                    base_n, ...) {

    wts_rep <- .ginteff_resolve_wts(weights, NULL,
                                    nvert = length(unique(cf[[".__vertex"]])),
                                    base_n = base_n)
    ap_args <- list(
        model,
        newdata    = cf,
        by         = ".__vertex",
        vcov       = vcov_arg,
        type       = type,
        conf_level = level,
        ...)
    if (!is.null(wts_rep)) ap_args$wts <- wts_rep
    ap <- do.call(marginaleffects::avg_predictions, ap_args)

    keep <- .ginteff_filter_eqn(ap, eqn)
    if (!is.null(keep)) ap <- ap[keep, , drop = FALSE]

    ord    <- order(ap$.__vertex)
    est    <- as.numeric(ap$estimate[ord])
    V_full <- stats::vcov(ap)
    if (!is.null(keep) && nrow(V_full) != length(ord))
        V_full <- V_full[keep, keep, drop = FALSE]
    V_ap   <- V_full[ord, ord, drop = FALSE]

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
        ap       = ap,
        obs      = NULL,
        keep     = NULL
    )
}

# =====================================================================
# Internal: resolve a `vcov` argument to a numeric matrix.
# =====================================================================

.ginteff_resolve_vcov <- function(vcov_arg, model) {
    if (is.matrix(vcov_arg)) return(vcov_arg)
    if (is.null(vcov_arg) || isTRUE(vcov_arg))
        return(stats::vcov(model))
    if (is.character(vcov_arg)) {
        if (!requireNamespace("sandwich", quietly = TRUE))
            stop("vcov = '", vcov_arg, "' requires the 'sandwich' package.")
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
                ## continuous in dydxs: central diff around obs
                nd[[v]] <- x + (if (hi) h else -h)
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
    n     <- nrow(base_data)
    nvert <- nrow(combos)
    ## Each vertex contributes n rows in order; reshape to n x nvert.
    M <- matrix(p$estimate, nrow = n, ncol = nvert, byrow = FALSE)
    as.numeric((M %*% signs) / divisor)
}


# =====================================================================
# Internal: data extraction
# =====================================================================

.ginteff_get_data <- function(model) {
    d <- tryCatch(stats::model.frame(model), error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0L) {
        d <- tryCatch(eval(stats::getCall(model)$data,
                           envir = environment(stats::formula(model))),
                      error = function(e) NULL)
    }
    if (is.null(d))
        stop("Could not extract data from `model`. Pass `data = ...`.")
    as.data.frame(d)
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
    tab <- data.frame(
        AIE   = x$aie,
        SE    = x$se,
        z     = x$z,
        p     = x$p,
        lower = x$ci[, 1],
        upper = x$ci[, 2],
        check.names = FALSE
    )
    rownames(tab) <- x$labels
    ## Use significant digits so very small AIEs (e.g. 1e-6) still print
    ## with non-zero figures, while moderate ones still align nicely.
    print(format(tab, digits = digits, scientific = FALSE,
                 zero.print = "0"))
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
    print(format(x$table, digits = digits, scientific = FALSE,
                 zero.print = "0"))
    invisible(x)
}

#' @export
coef.ginteff <- function(object, ...) object$aie

#' Variance-covariance matrix of the AIE vector
#'
#' Returns the delta-method variance-covariance matrix of the average
#' interaction effects stored in `object$aie`. Useful for chaining with
#' downstream tools (e.g. \pkg{multcomp}, \pkg{car}::deltaMethod).
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

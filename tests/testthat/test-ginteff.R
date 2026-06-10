## Smoke tests for ginteff().
##
## All tests reuse a shared linear-model fixture so the analytical AIE
## values are exact polynomials in the coefficients (no approximation
## bias).  Two more tests touch the binomial GLM path and the obseff
## machinery.

# Fixture ---------------------------------------------------------------------

make_fixture <- function() {
    set.seed(7)
    n  <- 200
    x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
    y  <- 0.5 + 0.3*x1 + 0.2*x2 - 0.1*x3 +
          0.4*x1*x2 + 0.1*x1*x2*x3 + rnorm(n, sd = 0.5)
    dat <- data.frame(y, x1, x2, x3)
    list(dat = dat,
         m   = lm(y ~ x1 * x2 * x3, data = dat),
         x3  = x3,
         n   = n)
}

# T1: 2-way FD matches analytical AIE = beta_{12} + beta_{123} * mean(x3)
test_that("2-way firstdiff matches the analytical AIE on a linear model", {
    fx <- make_fixture()
    g  <- ginteff(fx$m, firstdiff = c("x1", "x2"),
                  nunit = c(x1 = 1, x2 = 1))
    b  <- coef(fx$m)
    aie_hand <- mean(b["x1:x2"] + b["x1:x2:x3"] * fx$x3)
    expect_equal(unname(g$aie[[1]]), unname(aie_hand), tolerance = 1e-6)
})

# T2: 2-way dydxs (cross-partial) reproduces the FD answer exactly for
# a linear model (where the cross-partial *is* the second mixed coef).
test_that("2-way dydxs reproduces firstdiff on a linear model", {
    fx <- make_fixture()
    g_fd <- ginteff(fx$m, firstdiff = c("x1", "x2"),
                    nunit = c(x1 = 1, x2 = 1))
    g_dx <- ginteff(fx$m, dydxs = c("x1", "x2"))
    expect_equal(unname(g_dx$aie[[1]]),
                 unname(g_fd$aie[[1]]),
                 tolerance = 1e-3)
})

# T3: 3-way FD matches the analytical AIE = beta_{123}
test_that("3-way firstdiff matches the analytical 3-way AIE", {
    fx <- make_fixture()
    g  <- ginteff(fx$m, firstdiff = c("x1", "x2", "x3"))
    b  <- coef(fx$m)
    expect_equal(unname(g$aie[[1]]),
                 unname(b["x1:x2:x3"]),
                 tolerance = 1e-6)
})

# T4: HC3 vcov flows through and changes the SE
test_that("vcov = 'HC3' propagates to the interaction SE", {
    skip_if_not_installed("sandwich")
    fx <- make_fixture()
    g_cl  <- ginteff(fx$m, firstdiff = c("x1", "x2"))
    g_hc3 <- ginteff(fx$m, firstdiff = c("x1", "x2"), vcov = "HC3")
    ## Same point estimate, different SE
    expect_equal(unname(g_hc3$aie[[1]]), unname(g_cl$aie[[1]]),
                 tolerance = 1e-8)
    expect_false(isTRUE(all.equal(unname(g_hc3$se[[1]]),
                                  unname(g_cl$se[[1]]),
                                  tolerance = 1e-3)))
})

# T5: at() and atdxs() are applied without error
test_that("at() and atdxs() are accepted and applied", {
    fx <- make_fixture()
    g  <- ginteff(fx$m, dydxs = c("x1", "x2"),
                  at    = list(x3 = 0),
                  atdxs = list(x1 = "mean"))
    expect_s3_class(g, "ginteff")
    expect_length(g$aie, 1L)
    expect_true(is.finite(g$se))
})

# T6: obseff produces a per-observation matrix whose mean equals AIE
test_that("obseff = TRUE returns a per-obs matrix with mean == AIE", {
    fx <- make_fixture()
    g  <- ginteff(fx$m, firstdiff = c("x1", "x2"), obseff = TRUE)
    expect_equal(dim(g$obseff), c(fx$n, 1L))
    expect_equal(mean(g$obseff[, 1]), unname(g$aie[[1]]), tolerance = 1e-6)
})

# T7: factor x continuous yields one row per non-base factor level, all
# with finite SEs
test_that("factor x continuous interaction emits one row per non-base level", {
    set.seed(11)
    nf <- 400
    fx  <- factor(sample(c("a", "b", "c"), nf, replace = TRUE))
    xc  <- rnorm(nf)
    eta <- (-0.2 + (fx == "b") * 0.6 + (fx == "c") * -0.3) +
           (0.4 + (fx == "b") * -0.5 + (fx == "c") * 0.2) * xc
    yf  <- rbinom(nf, 1, plogis(eta))
    df  <- data.frame(yf, fx, xc)
    m   <- glm(yf ~ fx * xc, data = df, family = binomial)
    g   <- ginteff(m, dydxs = c("fx", "xc"))
    expect_length(g$aie, 2L)
    expect_true(all(is.finite(g$se)))
    expect_true(all(g$se < 1))
})

# T8: tidy as.data.frame output has the expected column set
test_that("as.data.frame.ginteff returns a tidy frame", {
    fx <- make_fixture()
    g  <- ginteff(fx$m, firstdiff = c("x1", "x2", "x3"))
    df <- as.data.frame(g)
    expect_named(df, c("term", "estimate", "std.error", "statistic",
                       "p.value", "conf.low", "conf.high"))
    expect_equal(nrow(df), 1L)
})

# T9: regression test for the inline-transform bug.  Before the fix in
# .ginteff_get_data, ginteff() pulled `model.frame(model)` as the data
# source -- which has columns named after the formula expression
# (`factor(group)`, `poly(x1,2)`) rather than the raw variable.  When the
# downstream predict() re-evaluated `predvars` against the counterfactual
# grid, R couldn't find the raw column and errored with
#   "object 'group' not found"
# (or the dydxs lookup itself failed: "Variable(s) not found in data: x1").
# After the fix, `.ginteff_get_data` prefers `model$data` / the original
# `data =` argument, which still contain the raw columns.
#
# The test checks that AIE / SE from a model with inline `factor()` match
# the AIE / SE from the same model fit with a pre-factored column.
test_that("inline factor() / poly() in formula doesn't break ginteff()", {
    set.seed(1)
    n  <- 500
    df <- data.frame(
        y     = rbinom(n, 1, 0.3),
        x1    = rnorm(n),
        x2    = rnorm(n),
        group = sample(1:5, n, replace = TRUE)
    )
    df$group_f <- factor(df$group)

    ## ----- inline factor(group) -----
    m_inline <- glm(y ~ x1 * x2 + factor(group),
                    family = binomial(link = "probit"), data = df)
    m_pre    <- glm(y ~ x1 * x2 + group_f,
                    family = binomial(link = "probit"), data = df)

    g_inline <- ginteff(m_inline, dydxs = c("x1", "x2"))
    g_pre    <- ginteff(m_pre,    dydxs = c("x1", "x2"))
    expect_equal(unname(g_inline$aie), unname(g_pre$aie), tolerance = 1e-10)
    expect_equal(unname(g_inline$se),  unname(g_pre$se),  tolerance = 1e-10)

    ## ----- inline poly(x1, 2) -----
    m_poly <- glm(y ~ poly(x1, 2) * x2 + group_f,
                  family = binomial(link = "probit"), data = df)
    g_poly <- ginteff(m_poly, dydxs = c("x1", "x2"))
    expect_true(is.finite(unname(g_poly$aie)))
    expect_true(is.finite(unname(g_poly$se)))

    ## ----- firstdiff path with inline factor -----
    g_fd <- ginteff(m_inline, firstdiff = c("x1", "x2"),
                    nunit = c(x1 = 1, x2 = 1))
    expect_true(is.finite(unname(g_fd$aie)))
    expect_true(is.finite(unname(g_fd$se)))

    ## ----- obseff path with inline factor -----
    g_obs <- ginteff(m_inline, dydxs = c("x1", "x2"), obseff = TRUE)
    expect_equal(dim(g_obs$obseff), c(n, 1L))
    expect_equal(mean(g_obs$obseff[, 1]), unname(g_obs$aie), tolerance = 1e-8)
})

# T10: panel data with NA-dropped rows + inline factor() in the formula.
# Regression test for the v0.1.1 follow-up:  glm() applied na.action
# silently dropped first-period-per-panel rows (where the lagged predictor
# is NA), so xlev was incomplete for the dropped rows' factor levels.
# Before the fix, .ginteff_get_data() returned the full raw frame and
# the analytic engine then tripped on
#   "factor factor(year) has new levels 1900".
# After the fix, the same na.action filter is applied to the raw data
# so the AIE is computed over exactly the rows the model was fit on
# (and ginteff()'s reported N matches glm's nobs()).
test_that("inline factor() + NA-dropped rows: ginteff respects na.action", {
    set.seed(42)
    n_per <- 30; n_grp <- 10
    df <- expand.grid(year = 1900:(1900 + n_per - 1), id = 1:n_grp)
    df$x  <- rnorm(nrow(df))
    df$z  <- rnorm(nrow(df))
    df$y  <- rbinom(nrow(df), 1, 0.3)
    df    <- df[order(df$id, df$year), ]
    df$x_lag <- ave(df$x, df$id, FUN = function(v) c(NA, head(v, -1)))

    m <- glm(y ~ x_lag * z + factor(year),
             family = binomial(link = "probit"), data = df)
    ## Sanity: glm dropped 10 rows (year=1900, x_lag is NA), and 1900
    ## is missing from xlev as a result.
    expect_true(length(stats::na.action(m)) == 10L)
    expect_false("1900" %in% m$xlevels[["factor(year)"]])

    g <- ginteff(m, dydxs = c("x_lag", "z"))
    expect_true(is.finite(unname(g$aie)))
    expect_true(is.finite(unname(g$se)))
    ## ginteff's N should equal glm's surviving-sample size, not nrow(df).
    expect_equal(g$n, stats::nobs(m))

    ## Should equal the result on the pre-filtered data.
    df_clean <- df[!is.na(df$x_lag), ]
    m_clean  <- glm(y ~ x_lag * z + factor(year),
                    family = binomial(link = "probit"), data = df_clean)
    g_clean  <- ginteff(m_clean, dydxs = c("x_lag", "z"))
    expect_equal(unname(g$aie), unname(g_clean$aie), tolerance = 1e-10)
    expect_equal(unname(g$se),  unname(g_clean$se),  tolerance = 1e-10)
})

# T11: per-variable h -- regression test for the shared-h scaling bug.
# Pre-0.2.0, dydxs used one step h = c_k * max(sd) for ALL continuous
# variables; with sd(xa)/sd(xb) = 1e5 the small variable was perturbed
# by ~1000 of its own sds and the AIE came out orders of magnitude too
# small. With per-variable steps the dydxs estimate must agree with a
# small-unit first-difference quotient (an h-free reference).
test_that("dydxs handles continuous variables on wildly different scales", {
    set.seed(3)
    n  <- 4000
    xa <- rnorm(n, sd = 1000)
    xb <- rnorm(n, sd = 0.01)
    eta <- -0.5 + 0.0008 * xa + 30 * xb + 0.012 * xa * xb
    y  <- rbinom(n, 1, plogis(eta))
    d  <- data.frame(y, xa, xb)
    m  <- suppressWarnings(glm(y ~ xa * xb, data = d, family = binomial))

    g_dx <- ginteff(m, dydxs = c("xa", "xb"))
    ## h-free reference: tiny first differences, scaled to a rate
    g_fd <- ginteff(m, firstdiff = c("xa", "xb"),
                    nunit = c(xa = 1, xb = 1e-5))
    ref  <- unname(g_fd$aie) / (1 * 1e-5)
    expect_equal(unname(g_dx$aie), ref, tolerance = 1e-3)

    ## named-vector h is accepted and equivalent at the same steps
    g_h <- ginteff(m, dydxs = c("xa", "xb"),
                   h = c(xa = 10, xb = 1e-4))
    expect_equal(unname(g_h$aie), ref, tolerance = 1e-3)
})

# T12: svyglm -- the AIE must be the design-weighted average, and the
# default weights must equal passing the design weights explicitly.
test_that("svyglm uses design weights by default", {
    skip_if_not_installed("survey")
    set.seed(4)
    n  <- 1500
    z  <- rnorm(n); f <- factor(rbinom(n, 1, 0.5))
    pw <- runif(n, 0.5, 4)
    y  <- rbinom(n, 1, plogis(-0.3 + 0.5 * z - 0.4 * (f == "1") +
                              0.6 * z * (f == "1")))
    d  <- data.frame(y, z, f, pw)
    des <- survey::svydesign(ids = ~1, weights = ~pw, data = d)
    m   <- survey::svyglm(y ~ z * f, design = des,
                          family = stats::quasibinomial())

    g_def <- ginteff(m, dydxs = c("z", "f"), obseff = TRUE)
    ## point estimate == design-weighted mean of the per-obs effects
    expect_equal(unname(g_def$aie),
                 stats::weighted.mean(g_def$obseff[, 1],
                                      stats::weights(m)),
                 tolerance = 1e-8)
    ## and == passing the design weights explicitly
    g_exp <- ginteff(m, dydxs = c("z", "f"),
                     weights = as.numeric(stats::weights(m)))
    expect_equal(unname(g_def$aie), unname(g_exp$aie), tolerance = 1e-10)
    expect_equal(unname(g_def$se),  unname(g_exp$se),  tolerance = 1e-10)
    ## weights = FALSE forces the unweighted average (Stata noweights)
    g_now <- ginteff(m, dydxs = c("z", "f"), weights = FALSE,
                     obseff = TRUE)
    expect_equal(unname(g_now$aie), mean(g_now$obseff[, 1]),
                 tolerance = 1e-8)
    expect_false(isTRUE(all.equal(unname(g_def$aie), unname(g_now$aie),
                                  tolerance = 1e-6)))
})

# T13: weights as a column name (documented API) works on the analytic
# engine and matches the same weights passed as a numeric vector.
test_that("weights = 'column' works on the analytic engine", {
    set.seed(5)
    n <- 600
    d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                    w  = runif(n, 0.2, 3))
    d$y <- rbinom(n, 1, plogis(0.2 + 0.4 * d$x1 - 0.3 * d$x2 +
                               0.5 * d$x1 * d$x2))
    m <- glm(y ~ x1 * x2, data = d, family = binomial)
    g_name <- ginteff(m, dydxs = c("x1", "x2"), weights = "w")
    g_num  <- ginteff(m, dydxs = c("x1", "x2"), weights = d$w)
    expect_equal(unname(g_name$aie), unname(g_num$aie), tolerance = 1e-10)
    expect_equal(unname(g_name$se),  unname(g_num$se),  tolerance = 1e-10)
})

# T14: full joint vcov across factor-contrast rows (analytic engine).
# diag must reproduce se^2; off-diagonals must be non-zero (the two
# contrasts share the base-level cells); and the difference test
# a' V a must be positive and smaller than the no-covariance value
# when the covariance is positive.
test_that("vcov.ginteff returns the full joint matrix on the analytic engine", {
    set.seed(11)
    n  <- 800
    fx <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
    f2 <- factor(rbinom(n, 1, 0.5))
    eta <- -0.2 + 0.5 * (fx == "b") - 0.3 * (fx == "c") +
           0.4 * (f2 == "1") +
           0.6 * (fx == "b") * (f2 == "1") - 0.5 * (fx == "c") * (f2 == "1")
    y  <- rbinom(n, 1, plogis(eta))
    d  <- data.frame(y, fx, f2)
    m  <- glm(y ~ fx * f2, data = d, family = binomial)
    g  <- ginteff(m, dydxs = c("fx", "f2"))

    V <- vcov(g)
    expect_equal(dim(V), c(2L, 2L))
    expect_equal(unname(diag(V)), unname(g$se^2), tolerance = 1e-10)
    expect_true(abs(V[1, 2]) > 0)
    expect_equal(V[1, 2], V[2, 1], tolerance = 1e-12)
    ## difference between the two contrasts: var(a'theta), a = (1, -1)
    a <- c(1, -1)
    v_diff <- as.numeric(t(a) %*% V %*% a)
    expect_true(v_diff > 0)
    expect_equal(v_diff,
                 unname(g$se[1]^2 + g$se[2]^2 - 2 * V[1, 2]),
                 tolerance = 1e-12)
})

# T15: mixed factor x continuous x continuous three-way dydxs equals the
# factor-contrast of the two-way continuous cross-partial computed at
# the factor's two levels via at().  This is the textbook definition of
# the mixed three-way effect; it is also the branch where the original
# Stata command pairs its internal derivative variables incorrectly
# (its output equals the factor contrast of an *own* second derivative
# instead -- verified against Stata's own margins).  Guard the R
# behaviour against regressions.
test_that("mixed factor x cont x cont dydxs equals the at()-difference identity", {
    set.seed(21)
    n  <- 2000
    f  <- factor(rbinom(n, 1, 0.5), c(0, 1), c("lo", "hi"))
    x1 <- rnorm(n, sd = 5)     # deliberately different scales
    x2 <- rnorm(n, sd = 0.5)
    eta <- -0.3 + 0.4 * (f == "hi") + 0.05 * x1 + 0.6 * x2 +
           0.04 * x1 * x2 + 0.03 * (f == "hi") * x1 +
           0.5 * (f == "hi") * x2 - 0.06 * (f == "hi") * x1 * x2
    y  <- rbinom(n, 1, pnorm(eta))
    d  <- data.frame(y, f, x1, x2)
    m  <- glm(y ~ f * x1 * x2, data = d, family = binomial("probit"))

    g3 <- ginteff(m, dydxs = c("f", "x1", "x2"))
    g_hi <- ginteff(m, dydxs = c("x1", "x2"), at = list(f = "hi"))
    g_lo <- ginteff(m, dydxs = c("x1", "x2"), at = list(f = "lo"))
    expect_equal(unname(g3$aie),
                 unname(g_hi$aie - g_lo$aie),
                 tolerance = 1e-8)
})

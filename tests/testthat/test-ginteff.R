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

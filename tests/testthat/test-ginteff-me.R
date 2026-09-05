## Tests for the *marginaleffects* fallback engine in ginteff().
##
## The analytic engine in test-ginteff.R covers lm/glm/svyglm.  These
## tests cover the path used for everything else (polr, multinom,
## probit-via-polr, ...).  Each test pins the AIE / SE produced by the
## fallback to the value Stata's `ginteff` produces on the *same model
## fit* on the *same data* (simdata.csv, bundled in inst/extdata).
## The reference numbers are taken from
## verification/me_fallback_examples.log -- see also the diff table at
## the bottom of verification/me_fallback_examples.R.
##
## Two kinds of expectations:
##   * AIEs are pinned tight (1e-5 relative tolerance) -- they match
##     Stata to 4-5 sig figs in every case.
##   * SEs are pinned with a wider tolerance per case, reflecting the
##     known finite-difference / Hessian-conditioning effects
##     documented in CLAUDE.md (2026-05-06 entry).  With marginaleffects
##     >= 1.0.0 the continuous x factor SEs land within 0.05% of Stata;
##     the wider tolerances keep the suite green on older versions.

skip_unless_me <- function() {
    skip_if_not_installed("marginaleffects")
}

load_fixture <- function() {
    path <- system.file("extdata", "simdata.csv", package = "ginteff")
    if (!nzchar(path)) skip("simdata.csv not bundled with ginteff")
    dat <- utils::read.csv(path)
    dat$female_f <- factor(dat$female, levels = c(0L, 1L),
                           labels = c("male", "female"))
    dat$race_f   <- factor(dat$race,   levels = c(1L, 2L, 3L),
                           labels = c("white", "black", "other"))
    dat$health_o <- ordered(dat$health, levels = 1:5)
    dat$health_f <- factor(dat$health, levels = 1:5)
    dat
}

# ---------------------------------------------------------------------------
# TM1: ologit (polr) -- continuous x factor.  Reference: Stata
#   ologit health c.age##i.female i.race
#   ginteff, dydxs(age female) predict(outcome(#2))
#     AIE = -0.00028479,  SE = 0.00018072
# ---------------------------------------------------------------------------
test_that("polr (logit), age x female: AIE matches Stata; SE within 5%", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ age * female_f + race_f,
                    data = dat, Hess = TRUE)
    g <- ginteff(m, dydxs = c("age", "female_f"), eqn = 2L)
    expect_equal(unname(g$aie),  -0.00028479, tolerance = 1e-4)
    ## SE drift is the marginaleffects FD-cross-partial step (~3% high
    ## on marginaleffects < 1.0.0; 0.04% from 1.0.0)
    expect_equal(unname(g$se),    0.00018072, tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# TM2: ologit (polr) -- factor x factor.  Reference: Stata
#   ologit health i.female##i.race age
#   ginteff, dydxs(female race) predict(outcome(#2))
#     AIE [1#2] =  0.00085646, SE = 0.00943743
#     AIE [1#3] = -0.00480777, SE = 0.00976232
# Two-row result; SE comparison uses the relative tolerance directly
# because factor#factor has NO finite-difference step at all => SEs
# should match Stata to ~6 decimals.
# ---------------------------------------------------------------------------
test_that("polr (logit), female x race: factor#factor SEs match Stata exactly", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ female_f * race_f + age,
                    data = dat, Hess = TRUE)
    g <- ginteff(m, dydxs = c("female_f", "race_f"), eqn = 2L)
    expect_equal(unname(g$aie), c( 0.00085646, -0.00480777), tolerance = 1e-4)
    expect_equal(unname(g$se),  c( 0.00943743,  0.00976232), tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# TM3: oprobit (polr method = probit).  Reference: Stata
#   oprobit health c.age##i.female i.race
#   ginteff, dydxs(age female) predict(outcome(#2))
#     AIE = -0.00040625,  SE = 0.00017371
# ---------------------------------------------------------------------------
test_that("polr (probit), age x female: AIE matches Stata; SE within 5%", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ age * female_f + race_f,
                    data = dat, method = "probit", Hess = TRUE)
    g <- ginteff(m, dydxs = c("age", "female_f"), eqn = 2L)
    expect_equal(unname(g$aie), -0.00040625, tolerance = 1e-4)
    expect_equal(unname(g$se),   0.00017371, tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# TM4: mlogit (multinom) -- continuous x factor.  Reference: Stata
#   mlogit health c.age##i.female i.race
#   ginteff, dydxs(age female) predict(outcome(#3))
#     AIE = -0.00281553, SE = 0.00099417
# ---------------------------------------------------------------------------
test_that("multinom, age x female (eqn=3): AIE matches Stata; SE within 3%", {
    skip_unless_me()
    skip_if_not_installed("nnet")
    dat <- load_fixture()
    m <- nnet::multinom(health_f ~ age * female_f + race_f,
                        data = dat, trace = FALSE)
    g <- ginteff(m, dydxs = c("age", "female_f"), eqn = 3L)
    expect_equal(unname(g$aie), -0.00281553, tolerance = 1e-4)
    expect_equal(unname(g$se),   0.00099417, tolerance = 0.03)
})

# ---------------------------------------------------------------------------
# TM5: mlogit (multinom) -- factor x factor.  Reference: Stata
#   mlogit health i.female##i.race age
#   ginteff, dydxs(female race) predict(outcome(#3))
#     AIE [1#2] = -0.05500156, SE = 0.04561435
#     AIE [1#3] = -0.08747977, SE = 0.04938842
# AIE has a slightly wider gap here (0.04% vs 0.001% for polr) because
# multinom's BFGS converges to a slightly different optimum on this
# data; the 0.5% AIE tolerance accommodates that.  SEs still match to
# ~4 decimals.
# ---------------------------------------------------------------------------
test_that("multinom, female x race (eqn=3): factor#factor SEs match Stata", {
    skip_unless_me()
    skip_if_not_installed("nnet")
    dat <- load_fixture()
    m <- nnet::multinom(health_f ~ female_f * race_f + age,
                        data = dat, trace = FALSE)
    g <- ginteff(m, dydxs = c("female_f", "race_f"), eqn = 3L)
    expect_equal(unname(g$aie), c(-0.05500156, -0.08747977), tolerance = 5e-3)
    expect_equal(unname(g$se),  c( 0.04561435,  0.04938842), tolerance = 1e-3)
})

# ---------------------------------------------------------------------------
# TM6: `eqn` actually selects the right outcome.  We pull the per-vertex
# AvgP column for two different eqn values and check the AIE comes out
# correspondingly different (and that the AvgP for eqn=k matches what
# marginaleffects::predictions returns for group=k directly).
# ---------------------------------------------------------------------------
test_that("eqn = k actually filters the marginaleffects predictions to group k", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ age * female_f + race_f,
                    data = dat, Hess = TRUE)
    g2 <- ginteff(m, dydxs = c("age", "female_f"), eqn = 2L)
    g3 <- ginteff(m, dydxs = c("age", "female_f"), eqn = 3L)
    expect_false(isTRUE(all.equal(unname(g2$aie), unname(g3$aie))))
    ## With 5 outcome levels, the per-eqn AIEs should sum (over k) to
    ## roughly zero -- changing X reallocates probability mass between
    ## levels but the total is conserved.
    aies <- vapply(1:5, function(k)
                       unname(ginteff(m, dydxs = c("age","female_f"), eqn = k)$aie),
                   numeric(1))
    expect_equal(sum(aies), 0, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# TM7: obseff = TRUE through the marginaleffects fallback returns a
# per-observation matrix whose mean (over rows) equals the AIE point
# estimate.  Same invariant the analytic engine guarantees in T6.
# ---------------------------------------------------------------------------
test_that("obseff = TRUE through the fallback: mean(obseff) == AIE", {
    skip_unless_me()
    skip_if_not_installed("nnet")
    dat <- load_fixture()
    m <- nnet::multinom(health_f ~ age * female_f + race_f,
                        data = dat, trace = FALSE)
    g <- ginteff(m, dydxs = c("age", "female_f"), eqn = 3L,
                 obseff = TRUE)
    expect_equal(dim(g$obseff), c(nrow(dat), 1L))
    expect_equal(mean(g$obseff[, 1]), unname(g$aie), tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# TM8: AIE for FD path through the marginaleffects fallback matches
# Stata even when the model's vcov is degraded.  Reference (Stata):
#   ologit health c.age##c.height i.female i.race weight
#   ginteff, firstdiff(age height) nunit((1) age height) predict(outcome(#2))
#     AIE = 9.520e-06   (SE diverges -- see verification log)
# We deliberately do NOT pin the SE here: MASS::polr's vcov on this
# fit is itself broken (see CLAUDE.md 2026-05-06), so ginteff faithfully
# propagates a wrong SE.  Pinning the AIE confirms the FD code path
# through the fallback is sound.
# ---------------------------------------------------------------------------
test_that("polr 2-way FD via fallback: AIE matches Stata (SE not pinned)", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ age * height + female_f + race_f + weight,
                    data = dat, Hess = TRUE)
    ## Since 0.2.0 ginteff() warns that this fit's vcov is unusable
    ## (NaN diagonal entries from the ill-conditioned BFGS Hessian)
    ## rather than silently propagating it.
    expect_warning(
        g <- ginteff(m,
                     firstdiff = c("age", "height"),
                     nunit     = c(age = 1, height = 1),
                     eqn       = 2L),
        "ill-conditioned")
    expect_equal(unname(g$aie), 9.520e-06, tolerance = 5e-3)
    expect_true(is.finite(unname(g$se)))
    expect_true(unname(g$se) > 0)
})

# ---------------------------------------------------------------------------
# TM9: a user-supplied *named* coefficient vcov matrix on the fallback must
# reproduce the default (vcov = TRUE) result.  polr's vcov(m) lists the
# slopes before the thresholds while marginaleffects' get_coef() puts the
# thresholds first; ginteff aligns the matrix by name before handing it
# over (marginaleffects < 1.0.0 matched positionally and returned a wildly
# wrong SE for exactly this matrix).
# ---------------------------------------------------------------------------
test_that("named vcov matrix on the fallback matches vcov = TRUE", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ age * female_f + race_f,
                    data = dat, Hess = TRUE)
    g0 <- ginteff(m, dydxs = c("age", "female_f"), eqn = 2L)
    g1 <- ginteff(m, dydxs = c("age", "female_f"), eqn = 2L, vcov = vcov(m))
    expect_equal(unname(g1$aie), unname(g0$aie), tolerance = 1e-10)
    expect_equal(unname(g1$se),  unname(g0$se),  tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# TM10: factor x factor through the fallback with obseff = TRUE.  This is
# the grid on which marginaleffects 1.0.0 drops the `.__vertex` grouping
# column (factor-level padding + multi-outcome model); ginteff completes
# the absent levels itself so no padding happens.  The per-observation
# matrix must average to the AIE and carry one column per non-base race
# level, and the AIE / SE must still hit the Stata pins from TM2.
# ---------------------------------------------------------------------------
test_that("factor x factor via the fallback with obseff: mean(obseff) == AIE", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ female_f * race_f + age,
                    data = dat, Hess = TRUE)
    g <- ginteff(m, dydxs = c("female_f", "race_f"), eqn = 2L, obseff = TRUE)
    expect_equal(dim(g$obseff), c(nrow(dat), 2L))
    expect_equal(unname(colMeans(g$obseff)), unname(g$aie), tolerance = 1e-8)
    expect_equal(unname(g$aie), c( 0.00085646, -0.00480777), tolerance = 1e-4)
    expect_equal(unname(g$se),  c( 0.00943743,  0.00976232), tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# TM11: the level-completion rows are inert.  An unused factor column with
# an absent level forces completion on a grid that otherwise needs none;
# the result must be identical to the plain run.
# ---------------------------------------------------------------------------
test_that("level-completion rows do not alter per-vertex averages", {
    skip_unless_me()
    skip_if_not_installed("MASS")
    dat <- load_fixture()
    m <- MASS::polr(health_o ~ age * female_f + race_f,
                    data = dat, Hess = TRUE)
    g0 <- ginteff(m, dydxs = c("age", "female_f"), eqn = 2L)
    dat2 <- dat
    dat2$junk <- factor("a", levels = c("a", "b"))
    g1 <- ginteff(m, dydxs = c("age", "female_f"), eqn = 2L, data = dat2)
    expect_equal(unname(g1$aie), unname(g0$aie), tolerance = 1e-12)
    expect_equal(unname(g1$se),  unname(g0$se),  tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# TM12: the derivative-step policy on the fallback follows the installed
# marginaleffects version -- per-variable steps (each variable's own sd)
# from 1.0.0, the shared max(sd) step before that.
# ---------------------------------------------------------------------------
test_that("fallback step policy follows the marginaleffects version", {
    skip_unless_me()
    dat <- load_fixture()
    h   <- .ginteff_resolve_h(NULL, c("age", "height"), dat, "marginaleffects")
    sds <- c(age = stats::sd(dat$age), height = stats::sd(dat$height))
    if (utils::packageVersion("marginaleffects") >= "1.0.0")
        expect_equal(h, 0.01 * sds)
    else
        expect_equal(unname(h), rep(0.01 * max(sds), 2L))
})

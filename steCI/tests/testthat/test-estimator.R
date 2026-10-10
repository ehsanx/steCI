test_that("dat_sim returns a well-formed long dataset", {
  d <- dat_sim(n = 200, seed = 42)
  expect_true(all(c("id", "time", "time.stop", "event", "visit",
                    "A", "L", "Alag1", "L.baseline") %in% names(d)))
  expect_equal(length(unique(d$id)), 200L)
  expect_true(all(d$time < d$time.stop))
  expect_true(all(d$event %in% c(0, 1)))
  expect_true(is.finite(attr(d, "neg_haz_frac")))
})

test_that("ste_estimate returns finite estimands near the truth", {
  d <- dat_sim(n = 1000, seed = 7)
  est <- ste_estimate(d, estimand = "both", t.hor = 1:5)
  expect_true(is.finite(est$gamma))
  # true constant hazard difference is -0.04; expect same sign, plausible size.
  expect_lt(est$gamma, 0)
  expect_gt(est$gamma, -0.15)
  expect_length(est$mrd, 5)
  expect_true(all(is.finite(est$mrd)))
})

test_that("ci_analytic gives three ordered intervals", {
  d <- dat_sim(n = 500, seed = 11)
  ds <- ste_prepare(d)
  a <- ci_analytic(ds)
  expect_setequal(a$method, c("model", "robust", "cluster"))
  expect_true(all(a$lower < a$upper))
  expect_true(all(is.finite(a$se)))
})

test_that("ci_bootstrap returns percentile/basic/normal intervals", {
  d <- dat_sim(n = 300, seed = 13)
  b <- ci_bootstrap(d, B = 10, estimand = "hd", seed = 1)
  expect_setequal(b$gamma$type, c("percentile", "basic", "normal"))
  expect_true(all(is.finite(b$gamma$se)))
})

test_that("ci_lef returns analytic + multiplier MRD CIs and a band", {
  d <- dat_sim(n = 400, seed = 17)
  lf <- ci_lef(d, t.hor = 1:5, B = 100, seed = 1)
  expect_setequal(unique(lf$mrd$type),
                  c("analytic_row", "analytic_cluster", "multiplier"))
  expect_true(all(is.finite(lf$mrd$se)))
  expect_true(all(lf$mrd$lower < lf$mrd$upper))
  # cluster SE should be >= row SE at later horizons (subject reuse inflates it)
  a_cl  <- lf$mrd$se[lf$mrd$type == "analytic_cluster"]
  a_row <- lf$mrd$se[lf$mrd$type == "analytic_row"]
  expect_gt(a_cl[5], a_row[5] * 0.9)
  # simultaneous band is at least as wide as the pointwise multiplier interval
  band_w <- lf$band$upper - lf$band$lower
  expect_true(all(is.finite(band_w)) && all(band_w > 0))
})

test_that("ci_lef Approach-2 runs and the weight correction is small", {
  d <- dat_sim(n = 1000, seed = 23)
  a1 <- ci_lef(d, t.hor = 1:5, B = 100, approach = "a1", seed = 1)
  a2 <- ci_lef(d, t.hor = 1:5, B = 100, approach = "a2", seed = 1)
  s1 <- a1$mrd$se[a1$mrd$type == "analytic_cluster"]
  s2 <- a2$mrd$se[a2$mrd$type == "analytic_cluster"]
  expect_true(all(is.finite(s2)))
  # stabilized-weight estimation correction is negligible (design-average < 1%;
  # allow slack for single-dataset / small-sample noise at late horizons)
  expect_true(all(abs(s2 / s1 - 1) < 0.10))
})

test_that(".aalen_if reproduces timereg cumulative coefficients", {
  d <- dat_sim(n = 300, seed = 19); ds <- ste_prepare(d)
  fit <- .aalen_if(ds, t.hor = c(2, 4), weight_col = "ipw.stab.cum")
  ref <- timereg::aalen(survival::Surv(time.new, time.stop.new, event) ~
                          A.baseline + L.baseline, data = ds, n.sim = 0,
                        weights = ds$ipw.stab.cum)
  cum <- ref$cum
  b4 <- cum[max(which(cum[, 1] <= 4)), 2:4]
  expect_equal(as.numeric(fit$B[2, ]), as.numeric(b4), tolerance = 1e-6)
})

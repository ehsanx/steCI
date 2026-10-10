# Tie handling (ties = "timereg" / "simultaneous") and exact grid times (time.grid).

# Continuous-time data on a 1/10 sub-visit grid: stop times rounded UP to the grid, so the
# trial-relative stop times of different subjects coincide (ties) -- and, without time.grid,
# coincide only up to the last bit when reached from different trial origins.
grid_data <- function(n, seed, g = 10) {
  d <- dat_sim(n = n, seed = seed)
  d$time.stop <- d$time + ceiling((d$time.stop - d$time) * g) / g
  d
}

test_that(".semiaalen reproduces timereg on tie-free data (estimate and all three variances)", {
  d <- dat_sim(n = 300, seed = 19); ds <- ste_prepare(d)
  expect_false(any(duplicated(ds$time.stop.new[ds$event == 1])))
  w <- ds$ipw.stab.cum
  f0 <- timereg::aalen(survival::Surv(time.new, time.stop.new, event) ~
                         const(A.baseline) + L.baseline,
                       data = ds, n.sim = 0, weights = w, robust = 1)
  f1 <- timereg::aalen(survival::Surv(time.new, time.stop.new, event) ~
                         const(A.baseline) + L.baseline,
                       data = ds, n.sim = 0, weights = w, robust = 1, clusters = ds$id)
  s <- .semiaalen(ds, clusters = ds$id)
  expect_equal(s$gamma, unname(f0$gamma[1, 1]), tolerance = 1e-12)
  expect_equal(s$var_model[1, 1], unname(f0$var.gamma[1, 1]), tolerance = 1e-12)
  # 1640 rows > max.clust = 1000: timereg bins the row-level clusters; so must we
  expect_gt(nrow(ds), 1000)
  expect_equal(s$var_robust[1, 1], unname(f0$robvar.gamma[1, 1]), tolerance = 1e-12)
  expect_equal(s$var_cluster[1, 1], unname(f1$robvar.gamma[1, 1]), tolerance = 1e-12)
})

test_that("ties = 'simultaneous' equals ties = 'timereg' on tie-free data", {
  d <- dat_sim(n = 300, seed = 7)
  eT <- ste_estimate(d, estimand = "both", ties = "timereg")
  eS <- ste_estimate(d, estimand = "both", ties = "simultaneous")
  expect_equal(eS$gamma, eT$gamma, tolerance = 1e-12)
  expect_equal(eS$mrd, eT$mrd, tolerance = 1e-12)
  aT <- ci_analytic(attr(eT, "dat.seq"), ties = "timereg")
  aS <- ci_analytic(attr(eT, "dat.seq"), ties = "simultaneous")
  expect_equal(aS$se, aT$se, tolerance = 1e-12)
})

test_that("time.grid removes last-bit near-ties and changes nothing else", {
  d <- grid_data(300, 29)
  ds0 <- ste_prepare(d); ds1 <- ste_prepare(d, time.grid = 10)
  near <- function(ds) { et <- sort(unique(ds$time.stop.new[ds$event == 1])); sum(diff(et) < 1e-9) }
  expect_gt(near(ds0), 0)
  expect_equal(near(ds1), 0)
  expect_lte(max(abs(ds1$time.stop.new - ds0$time.stop.new)), 1e-12)
  keep <- setdiff(names(ds0), "time.stop.new")
  expect_identical(ds1[keep], ds0[keep])
  expect_identical(ste_prepare(d, time.grid = NULL), ds0)   # default: untouched
  expect_error(ste_prepare(dat_sim(n = 100, seed = 1), time.grid = 10), "not on a 1/10 grid")
})

test_that("ties = 'simultaneous' is deterministic on tied data and equals the M6 fit", {
  d <- grid_data(300, 31)
  set.seed(1); s0 <- .Random.seed
  e <- ste_estimate(d, estimand = "both", ties = "simultaneous", time.grid = 10)
  a <- ci_analytic(attr(e, "dat.seq"), ties = "simultaneous")
  expect_identical(.Random.seed, s0)                 # no jitter, no random numbers
  expect_gt(sum(duplicated(attr(e, "dat.seq")$time.stop.new[attr(e, "dat.seq")$event == 1])), 0)
  lf <- ci_lef(d, B = 20, seed = 1, time.grid = 10)
  expect_equal(e$mrd, lf$point$mrd, tolerance = 1e-12)
  expect_true(all(is.finite(a$se)))
})

test_that("resampling intervals refuse a point estimate from the other tie convention", {
  d <- grid_data(150, 37)
  eT <- ste_estimate(d, estimand = "hd", ties = "timereg", time.grid = 10)
  expect_error(ci_bootstrap(d, B = 2, estimand = "hd", ties = "simultaneous", time.grid = 10,
                            point = eT), "ties")
  eS <- ste_estimate(d, estimand = "hd", ties = "simultaneous", time.grid = 10)
  b <- ci_bootstrap(d, B = 5, estimand = "hd", ties = "simultaneous", time.grid = 10,
                    point = eS, seed = 1)
  expect_identical(b$gamma$estimate[1], eS$gamma)    # centred on exactly the supplied estimate
  expect_identical(b$ties, "simultaneous")
})

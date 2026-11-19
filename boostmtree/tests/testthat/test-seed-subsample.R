# Regression test for a fixed control$seed reusing one subsample for every tree.
#
# Passing the same seed to rfsrc() at every boosting iteration drew the same
# subjects out of bag each time. Under cv.flag = TRUE the others were never out
# of bag, so their cross-validated mu stayed at the starting value: a flat line
# at the mean response of the other subjects.
#
# gamma.i.list records a subject's step only at iterations where that subject
# was out of bag, so its non-NULL entries are the out-of-bag set.

seed.fixture <- function() {
  set.seed(3)
  simLong(n = 30, n.time = 4, model = 2)$data.list
}

seed.fit <- function(d, seed, M = 6) {
  boostmtree(
    x = d$features, tm = d$time, id = d$id, y = d$y,
    M = M, cv.flag = TRUE, verbose = FALSE,
    control = boostmtree.control(seed = seed)
  )
}

oob.sets <- function(fit) {
  lapply(fit$gamma.i.list[[1]], function(gamma.m) {
    which(!vapply(gamma.m, is.null, logical(1)))
  })
}

test_that("a fixed seed draws a different subsample at each iteration", {
  d <- seed.fixture()
  sets <- oob.sets(seed.fit(d, seed = 11))

  expect_length(sets, 6L)
  expect_true(all(lengths(sets) > 0L))
  # Measured: one distinct set across all 6 iterations before the fix.
  expect_gt(length(unique(sets)), 1L)
})

test_that("a fixed seed lets every subject be out of bag over a long run", {
  d <- seed.fixture()
  fit <- seed.fit(d, seed = 11, M = 40)

  # Measured before the fix: the same 11 of 30 subjects at every iteration.
  expect_equal(sort(unique(unlist(oob.sets(fit)))), seq_len(30L))
})

test_that("the same seed reproduces the fit exactly", {
  d <- seed.fixture()
  fit.a <- seed.fit(d, seed = 11)
  fit.b <- seed.fit(d, seed = 11)

  expect_identical(fit.a$mu, fit.b$mu)
  expect_identical(oob.sets(fit.a), oob.sets(fit.b))
  expect_false(identical(oob.sets(seed.fit(d, seed = 12)), oob.sets(fit.a)))
})

test_that("a seeded fit neither reads nor moves the caller's RNG state", {
  d <- seed.fixture()
  set.seed(1)
  before <- .Random.seed
  fit.a <- seed.fit(d, seed = 11)
  expect_identical(.Random.seed, before)

  set.seed(2)
  fit.b <- seed.fit(d, seed = 11)
  expect_identical(fit.a$mu, fit.b$mu)
})

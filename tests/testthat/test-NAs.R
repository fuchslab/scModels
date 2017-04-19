
test_that("NA parameters in density function", {
  expect_true(is.na(dmpb(NA, 5, 3, 20)))
  expect_true(is.na(dmpb(1, NA, 3, 20)))
  expect_true(is.na(dmpb(1, 5, NA, 20)))
  expect_true(is.na(dmpb(1, 5, 3, NA)))
})

test_that("NA parameters in distribution function", {
  expect_true(is.na(pmpb(NA, 5, 3, 20)))
  expect_true(is.na(pmpb(2, NA, 3, 20)))
  expect_true(is.na(pmpb(2, 5, NA, 20)))
  expect_true(is.na(pmpb(2, 5, 3, NA)))
})

test_that("NA parameters in quantile function", {
  expect_true(is.na(qmpb(NA, 5, 3, 20)))
  expect_true(is.na(qmpb(0.2, NA, 3, 20)))
  expect_true(is.na(qmpb(0.2, 5, NA, 20)))
  expect_true(is.na(qmpb(0.2, 5, 3, NA)))
})

test_that("NA parameters in RNG function", {
  expect_error(rmpb(NA, 5, 3, 20))
  expect_warning(expect_true(is.na(rmpb(1, NA, 3, 20))))
  expect_warning(expect_true(is.na(rmpb(1, 5, NA, 20))))
  expect_warning(expect_true(is.na(rmpb(1, 5, 3, NA))))
})

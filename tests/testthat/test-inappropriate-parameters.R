
test_that("Wrong parameters in density function", {
  expect_warning(expect_true(is.nan(dmpb(0, -1, 3, 20))))
  expect_warning(expect_true(is.nan(dmpb(1, -1, 3, 20))))
  expect_warning(expect_true(is.nan(dmpb(1, 5, -1, 20))))
  expect_warning(expect_true(is.nan(dmpb(1, 5, 3, -1))))
})

test_that("Wrong parameters in distribution function", {
  expect_warning(expect_true(is.nan(pmpb(2, -1, 3, 20))))
  expect_warning(expect_true(is.nan(pmpb(2, 5, -1, 20))))
  expect_warning(expect_true(is.nan(pmpb(2, 5, 3, -1))))
})

test_that("Wrong parameters in quantile function", {
  expect_warning(expect_true(is.nan(qmpb(-0.5, 5, 3, 20))))
  expect_warning(expect_true(is.nan(qmpb(0.2, -1, 3, 20))))
  expect_warning(expect_true(is.nan(qmpb(0.2, 5, -1, 20))))
  # expect_true(is.nan(qmpb(0.2, 5, 3, -1)))
})

# test_that("Wrong parameters in RNG function", {
#   expect_error(rmpb(NA, 5, 3, 20))
#   expect_warning(expect_true(is.na(rmpb(1, NA, 3, 20))))
#   expect_warning(expect_true(is.na(rmpb(1, 5, NA, 20))))
#   expect_warning(expect_true(is.na(rmpb(1, 5, 3, NA))))
# })

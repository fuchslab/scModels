
test_that("Discrete probabilities sum to unity", {
  expect_equal(sum(dmpb(0:100, 5, 3, 20)), 1)
})

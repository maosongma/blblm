test_that("multiplication works", {
  m = 3
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  co <- coef(fit)
  # since m is the number of subdata + the intercept
  expect_equal(length(co), m+1)
})

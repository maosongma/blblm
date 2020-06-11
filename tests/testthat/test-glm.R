test_that("multiplication works", {
  fit <- blblm_glm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  expect_s3_class(fit, "blblm")
})

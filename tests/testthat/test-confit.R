test_that("confint works", {
  x=blblm(mpg~wt*hp,data=mtcars,m=3,B=100)
  expect_equal(class(confint(x, c("wt","hp"))[1,1]),'numeric')
})
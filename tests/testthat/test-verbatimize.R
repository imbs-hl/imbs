context("Check verbatimize")

x <- c('foo', 'bar', 'buz')
x_pasted <- '\\verb=foo=, \\verb=bar=, \\verb=buz'

test_that("verbatimze", {
  expect_that(verbatimze(x) == x_pasted, is_true())
  expect_that(length(verbatimize(x)), equals(1))
  expect_that(verbatimize(x), is_a('character'))
})

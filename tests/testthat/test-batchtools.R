
context("batchtools")

test_that("load_or_create_registry", {
  
  reg_dir <- tempdir()
  
  expect_class(load_or_create_registry(reg_dir), "Registry")
  
})
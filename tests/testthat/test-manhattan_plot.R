context("Check manhattan plot")

test_data <- data.frame(chromosome = sample(1:22, 1e6, TRUE), 
                        position = sample(1:23467548, 1e6, TRUE), 
                        pvalue = runif(1e6))

test_that("manhattan plot", {
  
  p <- manhattan_plot(test_data, 
                      header = list(chr = "chromosome", 
                                    pos = "position", 
                                    pval = "pvalue"), 
                      thin = TRUE, 
                      thin.param = list(logp.thin.thresh = 2.5, 
                                        logp.thin.rate=0.05, 
                                        buffer = 10000))
  
  expect_class(p, "ggplot")
  
})


test_that("manhattan plot input checks", {
  
  header <- list(chr = "chromosome", 
                 pos = "position", 
                 pval = "pvalue")
  thin.param <- list(logp.thin.thresh = 5,
                     logp.thin.rate = 0.1,
                     buffer = 5L)
  
  expect_error(manhattan_plot(test_data[, 1:2]))
  expect_error(manhattan_plot(test_data))
  expect_error(manhattan_plot(test_data, header = list()))
  expect_error(manhattan_plot(test_data, header = list(foo = "bar", fooo = "baar", foooo = 7)))
  expect_error(manhattan_plot(test_data, header = list(foo = "bar", fooo = "baar", foooo = "baaar")))
  expect_error(manhattan_plot(test_data, threshold = c(8.4, 289), header = header))
  expect_error(manhattan_plot(test_data, threshold = 8, header = header))
  expect_error(manhattan_plot(test_data, threshold = 5e-8, header = header, thin = 3))
  expect_error(manhattan_plot(test_data, threshold = 5e-8, header = header, thin = TRUE, thin.param = list()))
  expect_error(manhattan_plot(test_data, threshold = 5e-8, header = header, thin = TRUE, thin.param = list(foo = "bar", fooo = "baar", foooo = 7)))
  
})

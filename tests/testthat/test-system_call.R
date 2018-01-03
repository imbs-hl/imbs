
context("system_call")

test_that("system_call", {
  
  dummy_dir <- "ainSxnElGCdoU4Pz3CpMP3IsR56zQmIQRIpriHMCrWMPGfk4XuNke8iMgA8ftDopozDotAfMjF8Z0rY0h9qiSvE9mDsJ6dHwLlii"
  
  expect_error(system_call("cd", dummy_dir))
  expect_error(system_call(1))
  
  expect_character(system_call("ls"))
  
})
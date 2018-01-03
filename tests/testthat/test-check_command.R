
context("checkCommand")

test_that("checkCommand", {
  
  dummy_command <- "Kb1VL4bg1TaJ92cYqwNiz5nvZy96VfiFueDV9wjcU8gP54N0uJXZ0t1gSwewcLrhBTI47KIl74rzmoxAWxAACzjr8x0IiHQdRDND"
  command <- "cd"
  
  expect_identical(check_command, checkCommand)
  expect_identical(test_command, testCommand)
  expect_identical(assert_command, assertCommand)
  
  expect_true(check_command(command))
  expect_true(test_command(command))
  expect_equal(assert_command(command), command)
  
  expect_string(check_command(dummy_command))
  expect_equal(check_command(dummy_command), sprintf("%s: command not found", dummy_command))
  expect_equal(check_command(1), "Must be of type 'string', not 'double'")
  
  expect_false(test_command(1))
  expect_false(test_command(dummy_command))
  
  expect_error(assert_command(1))
  expect_error(assert_command(dummy_command))
  
})
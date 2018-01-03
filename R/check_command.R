
#' Check if command exists
#'
#' @param x  [\code{string}]\cr
#'           Command to check.
#'
#' @return If the check is successful, i.e. the system command is found, returns \code{TRUE}. If the check is not successful, returns a string with the error message.
#' @export
#'
#' @import BBmisc checkmate
#'
check_command <- function(x) {
  
  cs <- checkmate::checkString(x, na.ok = FALSE, min.chars = 1, null.ok = FALSE)
  
  if(checkmate::testString(cs)) {
    return(cs)
  }
  
  BBmisc::suppressAll({
    try_cmd <- try(system2(command = x,
                           stderr = TRUE,
                           stdout = TRUE),
                   silent = TRUE)
  })
  
  if (BBmisc::is.error(try_cmd)) {
    sprintf("%s: command not found", x)
  } else {
    TRUE
  }
  
}

#' @export
#' @rdname check_command
#' 
#' @details See \link[checkmate]{assert} for details.
assert_command <- checkmate::makeAssertionFunction(check_command)

#' @export
#' @rdname check_command
test_command <- checkmate::makeTestFunction(check_command)

#' @export
#' @rdname check_command
checkCommand <- check_command

#' @export
#' @rdname check_command
assertCommand <- assert_command

#' @export
#' @rdname check_command
testCommand <- test_command

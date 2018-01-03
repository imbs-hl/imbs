
#' Wrapper for \code{system2}
#'
#' @param bin  [\code{string}]\cr
#'             The system command to be invoked as character string.
#' @param ...  [any]\cr
#'             Additional arguments passed to \code{system2}.
#'
#' @details Output of '\code{stdout}' and '\code{stderr}' is captured in a character vector and checked for the exit status. If the status is not \code{0} an error is thrown. Otherwise the captured output is returned as a string.
#'
#' @return The captured output as a string.
#' @export
#' 
#' @import checkmate
#'
system_call <- function(bin, ...) {
  
  assertions <- checkmate::makeAssertCollection()
  
  assertCommand(bin, add = collection)
  
  checkmate::reportAssertions(assertions)
  
  sys_out <- system2(command = bin,
                     stdout = TRUE,
                     stderr = TRUE, ...)
  
  if (!is.null(attr(sys_out, "status"))) {
    stop(sprintf("%s exited with status %d.",
                 basename(bin),
                 attr(sys_out, "status")))
  }
  
  return(sys_out)
  
}

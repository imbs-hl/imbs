
#' Load or create a \code{batchtools} registry
#'
#' @param file.dir       [\code{string}]\cr
#'                       File path where the registry is stored.
#' @param work.dir       [\code{string}]\cr
#'                       File path where the working dir of the registry is set to.
#' @param writeable      [\code{flag}]\cr
#'                       If a registry is loaded, make it writeable?
#' @param overwrite      [\code{flag}]\cr
#'                       Should an existing registry be removed and recreated?
#' @param ...            [any]\cr
#'                       Further arguments passed to \code{\link[batchtools]{makeRegistry}} and \code{\link[batchtools]{loadRegistry}}.
#'
#' @return A \code{\link[batchtools]{Registry}}.
#' @export
#'
#' @import checkmate batchtools
#'
load_or_create_registry <- function(file.dir, work.dir = getwd(), 
                                    writeable = TRUE, 
                                    overwrite = FALSE, ...) {
  
  assertions <- makeAssertCollection()
  
  checkmate::assert_directory(dirname(file.dir), add = assertions)
  checkmate::assert_directory(work.dir, add = assertions)
  checkmate::assert_flag(writeable, add = assertions)
  checkmate::assert_flag(overwrite, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  dots <- list(...)
  
  tryCatch({
    load_args <- names(formals(batchtools::loadRegistry))
    reg <- do.call(batchtools::loadRegistry,
                   args = c(file.dir = file.dir,
                            work.dir = work.dir,
                            writeable = writeable,
                            dots[intersect(names(dots), load_args)]))
    if (overwrite) {
      batchtools::removeRegistry(reg = reg)
    } else {
      return(reg)
    }
  },
  error = function(e) {
    message(e$message)
    message("This deletes all files in ", file.dir, ". Proceeding in 5 seconds...")
    Sys.sleep(5)
    message("Recursively removing files in ", file.dir, "...")
    unlink(file.dir, recursive = TRUE, force = TRUE)
  })
  
  make_args <- names(formals(batchtools::makeRegistry))
  reg <- do.call(batchtools::makeRegistry,
                 args = c(file.dir = file.dir,
                          work.dir = work.dir,
                          dots[intersect(names(dots), make_args)]))
  
  return(reg)
  
}

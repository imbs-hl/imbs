
##' Paste strings using named elements
##'
##' @title Named paste
##' @param ... one or more R objects, to be converted to character vectors.
##' @param sep a character string to separate the terms.
##' @return Named vector of concatenated values.
##' @examples 
##' named_paste("a", c(b = "b", c = "c"))
##' @aliases named_paste0
##' @export
named_paste <- function(..., sep = " ") {
  args <- list(...)
  result <- paste(..., sep = sep)
  long_idx <- which(sapply(args, length) > 1)
  names(result) <- names(args[[long_idx[1]]])
  
  if (length(long_idx) > 1) {
    warning("More than one argument with length >1 found, using names of the first argument.")
  }
  result
}

##' @rdname named_paste
named_paste0 <- function(..., sep = "") {
  named_paste(..., sep = sep)
}

##' Paste file names using named elements
##'
##' @title Named file path
##' @param ... character vectors.
##' @param sep the path separator to use.
##' @return Named vector of concatenated file names.
##' @examples 
##' named_file_path("home", 
##'                 c(wright = "myWork", schillert = "mywork"), 
##'                 "beratungen")
##' @export
named_file_path <- function(..., sep = .Platform$file.sep) {
  named_paste(..., sep = sep)
}

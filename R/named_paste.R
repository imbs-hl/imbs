
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
  ## Check for multiple part with more than one element
  n_el <- sapply(args, length)
  long_idx <- which(n_el > 1)
    if (length(long_idx) > 1) {
    ## Compare names of multiple element parts
    nms_raw <- lapply(args[long_idx], names)
    nms <- nms_raw[!sapply(nms_raw, is.null)]
    uni_nms <- unique(nms)

    if (length(uni_nms) > 1) {
      print(nms)
      stop('Names do not match')
    } else {
      names(result) <- unlist(uni_nms)
      return(result)
    }
    names(result)
    return(result)
  } else {
    if (length(long_idx) == 1){
      names(result) <- names(unlist(args[long_idx]))  
      return(result)
    } else {
      warning('No names specified')
    }
  }
   result
}



##' @rdname named_paste
##' @export
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
##'                 c(homer = "projects", marge = "work"), 
##'                 "doo")
##' @export
named_file_path <- function(..., sep = .Platform$file.sep) {
  named_paste(..., sep = sep)
}

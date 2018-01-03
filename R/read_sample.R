
##' Reads SNPTEST sample files, ignores the second line though
##'
##' @title Read SNPTEST sample file
##' @param file Name of the sample file
##' @return data.frame containing the data
##' @export 

read_sample <- function(file){
  header <- read.table(file, as.is = TRUE, nrows = 1)
  body <-  read.table(file, as.is = TRUE, skip = 2,
                      header = FALSE)
  names(body) <- header
  body[body == -9] <- NA
  body
}
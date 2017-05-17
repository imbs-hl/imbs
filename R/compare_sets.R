
#' Compare two vectors with simple set theory
#' 
#' @param a [vector(n)] First set of elements
#' @param b [vector(m)] Second set of elements
#' @param label [char(2)] Alternative labels for the sets. If NULL, names of 
#'   \code{a} and \code{b} are taken.
#' @return List with names of sets (\code{label}), size of sets (\code{n}) and
#'   actual content (\code{sets}).
#' @export
#' 
#' @examples
#' 
#' x <- 1:10
#' y <- 4:12
#' compare_sets(a = x, b = y)

compare_sets <- function(a, b, label = NULL){
  sets <- list(a_and_b = sort(intersect(a, b)),
               a_not_b = sort(setdiff(a, b)),
               b_not_a = sort(setdiff(b, a)))
  
  ns <- sapply(sets, length)
  
  if (is.null(label)) {
    lab_a <- deparse(substitute(a))
    lab_b <- deparse(substitute(b))
  } else {
    lab_a <- label[1]
    lab_b <- label[2]
  }
  
  labels <- c(sprintf('%s AND %s', lab_a, lab_b),
              sprintf('%s NOT %s', lab_a, lab_b),
              sprintf('%s NOT %s', lab_b, lab_a))
  names(labels) <- c('a_and_b', 'a_not_b', 'b_not_a')
  
  message(sprintf('Number of %s: %d\n', 
                  labels, ns))
  
  list(label = labels,
       n = ns,
       sets = sets)
  
}

## x <- 1:10
## y <- 4:12
## cset <- compare_sets(a = x, b = y)

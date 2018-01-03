#' Combine rows with the same content
#'
#' @param df           [\code{data.frame}]\cr
#'                     \code{data.frame} object of which to combine rows with 
#'                     same content to multirows.
#' @param cols         [\code{integer}]\cr
#'                     Vector of column indices to apply multirow combination to.
#' @param cmidrule     [\code{flag}]\cr
#'                     Should cmidrules inserted between groups of combinations?
#'
#' @return A \code{data.frame} with multirow(s) inserted.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- data.frame(FIRST = rep("A", 28), 
#'                  SECOND = c(rep("A", 4), rep ("B", 6), rep("A", 8), rep("B", 10)), 
#'                  THIRD = rep(c(LETTERS[0:13 %% 24 + 1]), each = 2), 
#'                  FORTH = LETTERS[0:27 %% 24 + 1])
#'
#' print(
#'   xtable(
#'     do.multirow(df, cols = 1:3)), 
#'   booktabs = TRUE, 
#'   sanitize.text.function = identity, 
#'   include.rownames = FALSE
#' )
#' }
do.multirow <- function(df, cols=1:ncol(df), cmidrule = TRUE) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assertDataFrame(df, add = assertions)
  checkmate::assertInteger(cols, add = assertions)
  checkmate::assertFlag(cmidrule, add = assertions)
  
  checkmate::reportAssertions()
  
  df <- as.data.frame(df)
  for (c in cols) {
    runs <- rle(as.character(df[, c]))
    if (any(runs$lengths > 1)) {
      tmp <- rep("", nrow(df))
      if (cmidrule) {
        
        r <- c(1, 1 + head(cumsum(runs$lengths), -1))
        
        adjust <- sapply(mapply(seq, from = r, to = data.table::shift(r, type = "lead", fill = nrow(df) + 1) - 1, SIMPLIFY = FALSE),
                         function(rows) {
                           if (length(cols[cols > c]) == 0) {
                             "0pt"
                           } else {
                             num_cmids <- as.integer(max(
                               sapply(cols[cols > c],
                                      function(x) {
                                        rl <- rle(as.character(df[rows, x]))
                                        length(which(rl$lengths > 1)) - 1
                                      }
                               )
                             ))
                             sprintf("-\\dimexpr%d\\cmidrulewidth+%f\\aboverulesep+%f\\belowrulesep\\relax",
                                     num_cmids, ceiling(num_cmids/2)*0.925, floor(num_cmids/2)*0.925
                             )
                           }
                         })
        
        tmp[c(1, 1 + head(cumsum(runs$lengths), -1))] <- paste("\\multirow{", runs$lengths, "}{*}[", adjust, "]{", df[c(1, 1 + head(cumsum(runs$lengths), -1)), c], "}", sep = "")
      } else {
        tmp[c(1, 1 + head(cumsum(runs$lengths), -1))] <- paste("\\multirow{", runs$lengths, "}{*}{", df[c(1, 1 + head(cumsum(runs$lengths), -1)), c], "}", sep = "")
      }
      df[, c] <- tmp
      if (cmidrule) {
        tmp <- df[, 1]
        tmp[c(1 + head(cumsum(runs$lengths), -1))] <- paste("\\cmidrule(lr){", c, "-", ncol(df), "}", df[c(1 + head(cumsum(runs$lengths), -1)), 1], sep = "")
        df[, 1] <- tmp
      }
    }
  }
  return(df)
}

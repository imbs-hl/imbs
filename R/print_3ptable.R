
##' print xtable objekts with table notes using threeparttable
##'
##' @title print_3ptable
##' @param ... everything you would give to print.xtable 
##' @param notes explanation of the table footnotes, see example
##' @return a well formated Latex-table
##' @examples 
##' \dontrun{
##' notes <- paste('\\begin{tablenotes}',
##'                '\\item [1] Das ist der Parameterschätzer $\\beta$',
##'                '\\item [2] Das ist der p-Wert',
##'                '\\end{tablenotes}',
##'                sep = '\n')
##'                
##' print_3ptable(notes = notes, xtable(data2, caption = "Regression Coefficients",
##'                                     label = "tab:3", digits = c(0, 4, 4, 2, 2),
##'                                     display = c('s', 'f','f','f','e'), ),
##'               include.rownames = TRUE)
##'               }
##' @export
##' 
print_3ptable <- function(notes, ...){
  tmp_table <- print.xtable(...,
                            latex.environments = c(getOption("xtable.latex.environments", c("center")), 'threeparttable'),
                            sanitize.text.function=function(x){x},
                            print.results = FALSE)
  final_table <- sub('\\end{tabular}', paste0('\\end{tabular}', '\n', notes),
                     tmp_table, fixed=TRUE)
  cat(final_table)
}

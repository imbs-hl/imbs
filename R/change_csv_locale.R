
#' Convert csv into csv2 files and vice versa. 
#' Can cope with vectors of input files
#' 
#'@title Changes delimiters in csv files
#'
#' @param in_files Input csv file(s) [character]
#' @param out_files Output csv file(s) [character]
#' @param from Language/locale of source file must be \code{'en'} or \code{'de'}
#'
#' @return Nothing, invoked for side effect of writing new files
#' @import utils
#' @export
change_csv_locale <- function(in_files, out_files, from = c('de', 'en')){
  sep <- c(de = ';', en = ',')[from]  
  dec <- c(de = ',', en = '.')[from]  
  
  dat_list <- lapply(in_files, utils::read.table, 
                     header = TRUE, as.is = TRUE,
                     dec = dec, sep = sep)
  mapply(utils::write.table, x = dat_list, file = out_files,
         MoreArgs = list(row.names = FALSE, col.names = TRUE,
         dec = dec, sep = sep))
  
}


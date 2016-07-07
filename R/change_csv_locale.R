
#' Convert csv into csv2 files and vice versa. 
#' Can cope with vectors of input files
#' 
#'@title Changes delimiters in csv files
#'
#' @param in_file Input csv file(s) [character]
#' @param out_file Output csv file(s) [character]
#' @param from Language/locale of source file must be \code{'en'} or \code{'de'}
#'
#' @return Nothing, invoked for side effect of writing new files
#' @export
change_csv_locale <- function(in_file, out_file, from = c('de', 'en')){
  sep <- c(de = ';', en = ',')[from]  
  dec <- c(de = ',', en = '.')[from]  
  
  dat_list <- lapply(in_files, read.table, 
                     header = TRUE, as.is = TRUE,
                     dec = dec, sep = sep)
  lapply(dat_list, write.table, 
         row.names = FALSE, col.names = TRUE,
         dec = dec, sep = sep)
  
}


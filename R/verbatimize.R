##' Encapsulate elements of vector in LaTeX verbatim commands
##'
##' @title Verbatimize string
##' @param x Character vector to be encapsulated
##' @param sep String; How pasted elements should be separated, Default is ", ".
##' @return String of pasted values
##' @export
verbatimize <- function(x, sep = ', '){
    ## Determine possible character to use as verb border
    borders <- c('=', '+', ';', '-', '_', letters)

    chars <- unlist(lapply(x, strsplit, split = ""))
    border_cands <- setdiff(borders, chars)

    if(length(border_cands) < 1){
        stop('No candidate for separating verbatim chunks found!')
    } else {
        border <- border_cands[1]
    }
    
    verb <- paste('\\verb', border, x, border, sep = "")
    paste(verb, collapse = sep)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_csv_locale.R
\name{change_csv_locale}
\alias{change_csv_locale}
\title{Changes delimiters in csv files}
\usage{
change_csv_locale(in_files, out_files, from = c("de", "en"))
}
\arguments{
\item{in_files}{Input csv file(s) [character]}

\item{out_files}{Output csv file(s) [character]}

\item{from}{Language/locale of source file must be \code{'en'} or \code{'de'}}
}
\value{
Nothing, invoked for side effect of writing new files
}
\description{
Convert csv into csv2 files and vice versa. 
Can cope with vectors of input files
}

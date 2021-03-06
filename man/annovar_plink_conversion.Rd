% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annovar.R
\name{annovar_plink_conversion}
\alias{annovar_plink_conversion}
\title{Convert PLINK bim file to ANNOVAR input}
\usage{
annovar_plink_conversion(bim.file, output.prefix, chr)
}
\arguments{
\item{bim.file}{[\code{string}]\cr
The PLINK \code{bim} file.}

\item{output.prefix}{[\code{string}]\cr
The basename of the output files.}

\item{chr}{[\code{numeric}]\cr
Numeric vector of chromosomes to consider.}
}
\value{
Nothing.
}
\description{
Convert PLINK bim file to ANNOVAR input
}
\details{
Uses \code{awk} to transform information from \code{bim} file to valid input to ANNOVAR.
}

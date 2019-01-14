
#' Convert PLINK bim file to ANNOVAR input
#'
#' @param bim.file            [\code{string}]\cr
#'                            The PLINK \code{bim} file.
#' @param output.prefix       [\code{string}]\cr
#'                            The basename of the output files.
#' @param chr                 [\code{numeric}]\cr
#'                            Numeric vector of chromosomes to consider.
#'
#' @details Uses \code{awk} to transform information from \code{bim} file to valid input to ANNOVAR.
#'
#' @return Nothing.
#' 
#' @export
#'
annovar_plink_conversion <- function(bim.file, output.prefix, chr) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(bim.file, add = assertions)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  if (!missing(chr)) {
    checkmate::assert_numeric(chr, lower = 0, upper = 26, finite = TRUE, any.missing = FALSE, all.missing = FALSE, null.ok = FALSE, add = assertions)
    chr_pattern <- sprintf("$1~/^(%s)$/", paste(chr, collapse = "|"))
  } else {
    chr_pattern <- ""
  }
  
  checkmate::reportAssertions(assertions)
  
  
  system_call(
    bin = "awk",
    args = c(
      sprintf("'%s{l=length($5); end=$4+(l-1); if ($5==\"*\") {$5=\"-\"}; if ($6==\"*\") {$6=\"-\"}; print $1,$4,end,$5,$6}'", chr_pattern),
      bim.file,
      ">", sprintf("%s.avinput", output.prefix)
    )
  )
  
}

#' Annotate variant regions
#'
#' @param avinput.file        [\code{string}]\cr
#'                            The file name of the ANNOVAR input file.
#' @param output.prefix       [\code{string}]\cr
#'                            The basename of the output files.
#' @param db.path             [\code{string}]\cr
#'                            The path to the database files. Default is \code{/imbs/external_data/annotation_and_references/annovar/humandb}.
#' @param build               [\code{string}]\cr
#'                            Human genome build code. Default is 'hg19'.
#' @param exec                [\code{string}]\cr
#'                            Path to \code{annotate_variation} script.
#'                            
#' @details Use ANNOVAR to find if variants hit exons or hit intergenic regions, or hit introns, or hit a non-coding RNA genes.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
annovar_annotate_region <- function(avinput.file, output.prefix, db.path = "/imbs/external_data/annotation_and_references/annovar/humandb", build = "hg19", exec = "annotate_variation.pl") {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(avinput.file, add = assertions)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_directory(db.path, add = assertions)
  
  checkmate::assert_choice(build, c("hg18", "hg19", "hg38"), add = assertions)
  
  assert_command(exec, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(
      sprintf("-out %s", output.prefix),
      sprintf("-build %s", build),
      avinput.file,
      db.path
    )
  )
  
}
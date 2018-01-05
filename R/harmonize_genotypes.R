
#' Genotype harmonization
#' 
#' @description Harmonization of genotype data stored using different file formats with different and potentially unknown strands. Linkage disequilibrium (LD) patterns are used to determine the correct strand GC and AT SNPs. This is a simple wrapper for \code{GenotypeHarmonizer}.
#' 
#' @details \url{https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer}
#'
#' @param input           [\code{string}]\cr
#'                        The base path of the data to align. The extensions are determined based on the input data type.
#' @param ref             [\code{string}]\cr
#'                        The base path of the reference data used for alignment. The extensions are determined based on the input data type. If not specified the input data is simply converted to the specified output type.
#' @param output          [\code{string}]\cr
#'                        The base path of the output data.
#' @param inputType       [\code{string}]\cr
#'                        The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path.
#' @param refType         [\code{string}]\cr
#'                        The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path.
#' @param outputType      [\code{string}]\cr
#'                        The output data type. Defaults to --inputType.
#' @param ...             [\code{character}]\cr
#'                        Other options passed to \code{GenotypeHarmonizer}.
#' @param exec            [\code{string}]\cr
#'                        Path of \code{GenotypeHarmonizer} executable.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
harmonize_genotypes <- function(input, ref, output,
                                inputType, refType, outputType,
                                ...,
                                exec = "GenotypeHarmonizer") {
  
  inputTypes <- list(
    PED_MAP = c("ped", "map"),
    VCF = c("vcf.gz", "vcf.gz.tbi"),
    PLINK_BED = c("bed", "bim", "fam"),
    SHAPEIT2 = c("haps", "sample"),
    GEN = c("gen", "sample")
  )
  
  outputTypes <- list(
    PED_MAP = c("ped", "map"),
    PLINK_BED = c("bed", "bim", "fam"),
    SHAPEIT2 = c("haps", "sample"),
    GEN = c("gen", "sample")
  )
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assertTRUE(
    any(
      sapply(inputTypes, 
             function(extensions) {
               checkmate::testFile(sprintf("%s.%s", input, extensions))
             }
      )
    ), add = assertions
  )
  
  if (!missing(ref)) {
  checkmate::assertTRUE(
    any(
      sapply(inputTypes, 
             function(extensions) {
               checkmate::testFile(sprintf("%s.%s", ref, extensions))
             }
      )
    ), add = assertions
  )
    ref <- sprintf("--ref %s", ref)
  } else {
    ref <- ""
  }
  
  checkmate::assertPathForOutput(dirname(output), add = assertions)
  
  if (!missing(inputType)) {
    checkmate::assertChoice(inputType, names(inputTypes))
    inputType <- sprintf("--inputType %s", inputType)
  } else {
    inputType <- ""
  }
  if (!missing(refType)) {
    checkmate::assertChoice(refType, names(inputTypes))
    refType <- sprintf("--inputType %s", refType)
  } else {
    refType <- ""
  }
  if (!missing(outputType)) {
    checkmate::assertChoice(outputType, names(outputTypes))
    outputType <- sprintf("--inputType %s", outputType)
  } else {
    outputType <- ""
  }
  
  assertCommand(exec, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(exec,
              args = c("--input", input,
                       inputType,
                       ref,
                       refType,
                       "--output", output,
                       outputType,
                       ...))
  
}
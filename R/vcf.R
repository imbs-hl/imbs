
#' Normalize VCF file to bi-allelic variants
#' 
#' Converting VCF files to plink format has never been easier. However, there are a few issues related to some intrinsic limitations of the plink format. The first is related to the fact that variants in a plink file are bi-allelic only, while variants in a VCF file can be multi-allelic. The second is related to an intrinsic limitation of plink which makes indel definitions ambiguous. Here is an example: is the following variant an insertion or a deletion compared to the GRCh37 reference?
#' 
#' 20 31022441 A AG
#' 
#' There is no way to tell, as the plink format does not record this information.
#' 
#' Keeping this in mind, we are going to split mulit-allelic variants into bi-allelic ones, left-normalize indels, and assign unique idetifiers.
#'
#' @param vcf.file           [\code{string}]\cr
#'                           The input VCF file path.
#' @param ref.file           [\code{string}]\cr
#'                           A human reference genome \code{fasta} file to normalize indels against.
#' @param output.file        [\code{string}]\cr
#'                           The output VCF file path.
#' @param bcftools.exec      [\code{string}]\cr
#'                           Path of bcftools executable. 
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by bcftools
#'                           Default is determined by SLURM environment variables and at least 1. 
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
vcf_normalization <- function(vcf.file, ref.file, output.file,
                              bcftools.exec = "bcftools",
                              num.threads) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, add = assertions)
  checkmate::assert_file(ref.file, add = assertions)
  checkmate::assert_directory(dirname(output.file), add = assertions)
  
  assert_command(bcftools.exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  # Run bcftools
  system_call(
    bin = bcftools.exec,
    args = c(
      "norm", "-Ou", "-m", "-any", vcf.file, # split multi-allelic alleles
      "|",
      bcftools.exec, "norm", "-Ou", "-f", ref.file, # normalize indels
      "|",
      bcftools.exec, "annotate", "--threads", num.threads, "-Oz", "-o", output.file, "-x", "ID", "-I", "+'%CHROM:%POS:%REF:%ALT'" # assign unique identifier
    )
  )
  
}
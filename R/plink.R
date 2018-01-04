
#' Write SNP or sample IDs to file
#' 
#' @param ids Vector of SNP or Sample IDs
#' @param file File name to write to
#' @param update Logical, Default is FALSE. Skip writing file if it exists
#'   
#' @param sample Logical. Should be TRUE for samples (the default) and FALSE for
#'   SNPs
#'   
#' @return NULL
#' @import utils
#' @export
write_plink_ids <- function(ids, file, update, sample = TRUE) {
  if (sample) {
    out <- cbind(ids, ids)
  } else {
    out <- ids
  }
  if (update | !file.exists(file)) {
    utils::write.table(
      out,
      file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
}

#' Subset samples and variants using PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param remove             [\code{string}]\cr
#'                           File name of a list of sample FIDs and IIDs to remove.
#' @param keep               [\code{string}]\cr
#'                           File name of a list of sample FIDs and IIDs to keep.
#' @param exclude            [\code{string}]\cr
#'                           File name of a list of variant names to remove.
#' @param extract            [\code{string}]\cr
#'                           File name of a list of variant names to keep.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by SLURM environment variables and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_subset <- function(bfile, output.prefix, remove, keep, exclude, extract, ...,
                         exec = "plink",
                         num.threads = max(1, as.integer(Sys.getenv("SLURM_NPROCS")), na.rm = TRUE),
                         memory = max(5000, as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE)) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
  checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
  checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  if (missing(remove)) {
    remove <- ""
  } else {
    checkmate::assert_file(remove, add = assertions)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(keep)) {
    keep <- ""
  } else {
    checkmate::assert_file(keep, add = assertions)
    keep <- sprintf("--keep %s", keep)
  }
  
  if (missing(exclude)) {
    exclude <- ""
  } else {
    checkmate::assert_file(exclude, add = assertions)
    exclude <- sprintf("--exclude %s", exclude)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    checkmate::assert_file(extract, add = assertions)
    extract <- sprintf("--extract %s", extract)
  }
  
  assert_command(exec, add = assertions)
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c("--bfile", bfile,
             "--threads", num.threads,
             "--memory", memory,
             remove,
             keep,
             extract,
             exclude,
             "--make-bed",
             "--allow-extra-chr",
             "--out", output.prefix, ...)
  )
  
}

#' Convert VCF files to PLINK binary files
#'
#' @param vcf.file           [\code{string}]\cr
#'                           The VCF file path.
#' @param ref.file           [\code{string}]\cr
#'                           A human reference genome \code{fasta} file to normalize indels against.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param build              [\code{string}]\cr
#'                           Human genome build code to split the X chromosome into pseudo-autosomal region and pure X. Default is 'b37' (alias 'hg19').
#' @param bcftools.exec      [\code{string}]\cr
#'                           Path of bcftools executable.
#' @param plink.exec         [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by SLURM environment variables and at least 5000.
#'
#' @details Based on the best practices on \url{http://apol1.blogspot.de/2014/11/best-practice-for-converting-vcf-files.html}. The procedure will take the VCF file, strip the variant IDs, split multi-allelic sites into bi-allelic sites, assign names to make sure indels will not become ambiguous, and finally convert to PLINK format. See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_conversion <- function(vcf.file, ref.file, output.prefix, ..., 
                             build = "b37",
                             bcftools.exec = "bcftools", plink.exec = "plink",
                             num.threads = max(1, as.integer(Sys.getenv("SLURM_NPROCS")), na.rm = TRUE),
                             memory = max(5000, as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE)) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, add = assertions)
  checkmate::assert_file(ref.file, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  checkmate::assertChoice(build, c("b36", "hg18", "b37", "hg19", "b38", "hg38"), add = assertions)
  assert_command(bcftools.exec, add = assertions)
  assert_command(plink.exec, add = assertions)
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  # Run PLINK
  system_call(
    bin = bcftools.exec,
    args = c("norm", "-Ou", "-m", "-any", vcf.file,
             "|",
             bcftools.exec, "norm", "-Ou", "-f", ref.file,
             "|",
             bcftools.exec, "annotate", "-Ob", "-x", "ID", "-I", "+'%CHROM:%POS:%REF:%ALT'",
             "|",
             plink.exec, "--bcf", "/dev/stdin",
             "--keep-allele-order",
             "--vcf-idspace-to", "_",
             "--const-fid",
             "--allow-extra-chr", "0",
             "--split-x", build, "no-fail",
             "--threads", num.threads,
             "--memory", memory,
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
}

#' Merge two PLINK datasets
#'
#' @param old.prefix       [\code{string}]\cr
#'                         The basename of the first binary PLINK file set.
#' @param new.prefix       [\code{string}]\cr
#'                         The basename of the second binary PLINK file set.
#' @param merge.mode       [\code{int}]\cr
#'                         Merge mode.
#' @param output.prefix    [\code{string}]\cr
#'                         The basename of the output binary PLINK file set.
#' @param ...              [\code{character}]\cr
#'                         Additional arguments passed to PLINK.
#' @param exec             [\code{string}]\cr
#'                         Path of PLINK executable.
#' @param num.threads      [\code{int}]\cr
#'                         Number of CPUs usable by PLINK.
#'                         Default is determined by SLURM environment variables and at least 1.
#' @param memory           [\code{int}]\cr
#'                         Memory for PLINK in Mb.
#'                         Default is determined by SLURM environment variables and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_merge <- function(old.prefix, new.prefix, merge.mode,
                        output.prefix,
                        ...,
                        exec = "plink",
                        num.threads = max(1, as.integer(Sys.getenv("SLURM_NPROCS")), na.rm = TRUE),
                        memory = max(5000, as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE)) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(sprintf("%s.bed", old.prefix), add = assertions)
  checkmate::assert_file(sprintf("%s.bim", old.prefix), add = assertions)
  checkmate::assert_file(sprintf("%s.fam", old.prefix), add = assertions)
  
  checkmate::assert_file(sprintf("%s.bed", new.prefix), add = assertions)
  checkmate::assert_file(sprintf("%s.bim", new.prefix), add = assertions)
  checkmate::assert_file(sprintf("%s.fam", new.prefix), add = assertions)
  
  checkmate::assert_int(merge.mode, lower = 1, upper = 7, add = assertions)
  
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  assert_command(exec, add = assertions)
  
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  # Make a merge with PLINK
  system_call(
    bin = exec,
    args = c("--bfile ", old.prefix,
             "--bmerge ", new.prefix,
             "--threads", num.threads,
             "--memory", memory,
             "--merge-mode ", merge.mode,
             "--out", output.prefix,
             ...)
  )
  
}

#' LD pruning with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param window.size        [\code{int}]\cr
#'                           Window size in kilobase.
#' @param step.size          [\code{int}]\cr
#'                           Step size in kilobase.
#' @param threshold          [\code{number}]\cr
#'                           Pairwise \eqn{r^2} threshold.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by SLURM environment variables and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_ld_pruning <- function(bfile, output.prefix, 
                             window.size, step.size, threshold, ...,
                             exec = "plink",
                             num.threads = max(1, as.integer(Sys.getenv("SLURM_NPROCS")), na.rm = TRUE),
                             memory = max(5000, as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE)) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
  checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
  checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_int(window.size, lower = 2, add = assertions)
  checkmate::assert_int(step.size, lower = 1, add = assertions)
  checkmate::assert_number(threshold, lower = 0, upper = 1, add = assertions)
  
  assert_command(exec, add = assertions)
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c("--bfile", bfile,
             "--threads", num.threads,
             "--memory", memory,
             "--keep-allele-order",
             "--indep-pairwise", sprintf("%dkb", window.size), sprintf("%dkb",step.size), threshold,
             "--allow-extra-chr",
             "--out", output.prefix, ...)
  )
  
}

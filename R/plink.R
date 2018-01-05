
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
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_subset <- function(bfile, output.prefix, remove, keep, exclude, extract, ...,
                         bed.file = NULL, bim.file = NULL, fam.file = NULL,
                         exec = "plink",
                         num.threads,
                         memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
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
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  }  
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000,     
                 min(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000,           
                     num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE),                 
                 na.rm = TRUE)  
  }  
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             remove,
             keep,
             extract,
             exclude,
             "--make-bed",
             "--keep-allele-order",
             "--allow-extra-chr", "0",
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
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details Based on the best practices on \url{http://apol1.blogspot.de/2014/11/best-practice-for-converting-vcf-files.html}. The procedure will take the VCF file, strip the variant IDs, split multi-allelic sites into bi-allelic sites, assign names to make sure indels will not become ambiguous, and finally convert to PLINK format. See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_vcf_conversion <- function(vcf.file, ref.file, output.prefix, ..., 
                             build = "b37",
                             bcftools.exec = "bcftools", plink.exec = "plink",
                             num.threads,
                             memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  checkmate::assert_file(vcf.file, add = assertions)
  checkmate::assert_file(ref.file, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  checkmate::assertChoice(build, c("b36", "hg18", "b37", "hg19", "b38", "hg38"), add = assertions)
  
  assert_command(bcftools.exec, add = assertions)
  assert_command(plink.exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)  
  } 
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {  
    memory = max(5000,       
                 min(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000,              
                     num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE),                   na.rm = TRUE)  
  } 
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
             "--vcf-idspace-to", "_",
             "--const-fid",
             "--keep-allele-order",
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
#' @param first.prefix     [\code{string}]\cr
#'                         The basename of the first binary PLINK file set.
#' @param second.prefix    [\code{string}]\cr
#'                         The basename of the second binary PLINK file set.
#' @param merge.mode       [\code{int}]\cr
#'                         Merge mode.
#' @param output.prefix    [\code{string}]\cr
#'                         The basename of the output binary PLINK file set.
#' @param ...              [\code{character}]\cr
#'                         Additional arguments passed to PLINK.
#' @param first.bed.file   [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param first.bim.file   [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param first.fam.file   [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param second.bed.file  [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param second.bim.file  [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param second.fam.file  [\code{string}]\cr
#'                         Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec             [\code{string}]\cr
#'                         Path of PLINK executable.
#' @param num.threads      [\code{int}]\cr
#'                         Number of CPUs usable by PLINK.
#'                         Default is determined by SLURM environment variables and at least 1.
#' @param memory           [\code{int}]\cr
#'                         Memory for PLINK in Mb.
#'                         Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate
#'
plink_merge <- function(first.prefix, second.prefix, merge.mode,
                        output.prefix,
                        ...,
                        first.bed.file = NULL, first.bim.file = NULL, first.fam.file = NULL,
                        second.bed.file = NULL, second.bim.file = NULL, second.fam.file = NULL,
                        exec = "plink",
                        num.threads,
                        memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(first.prefix)) {
    checkmate::assert_file(sprintf("%s.bed", first.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", first.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", first.prefix), add = assertions)
    first <- sprintf("--bfile %s", first.prefix)
  } else {
    checkmate::assert_file(first.bed.file, add = assertions)
    checkmate::assert_file(first.bim.file, add = assertions)
    checkmate::assert_file(first.fam.file, add = assertions)
    first <- sprintf("--bed %s --bim %s --fam %s", first.bed.file, first.bim.file, first.fam.file)
  }
  
  if (!missing(second.prefix)) {
    checkmate::assert_file(sprintf("%s.bed", second.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", second.prefix), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", second.prefix), add = assertions)
    second <- sprintf("--bmerge %s", second.prefix)
  } else {
    checkmate::assert_file(second.bed.file, add = assertions)
    checkmate::assert_file(second.bim.file, add = assertions)
    checkmate::assert_file(second.fam.file, add = assertions)
    second <- sprintf("--bmerge %s %s %s", second.bed.file, second.bim.file, second.fam.file)
  }
  
  checkmate::assert_int(merge.mode, lower = 1, upper = 7, add = assertions)
  
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {   
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE) 
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {   
    memory = max(5000,           
                 min(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000,                    
                     num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE),      
                 na.rm = TRUE) 
  }  
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  # Make a merge with PLINK
  system_call(
    bin = exec,
    args = c(first,
             second,
             "--threads", num.threads,
             "--memory", memory,
             "--keep-allele-order",
             "--allow-extra-chr", "0",
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
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
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
                             bed.file = NULL, bim.file = NULL, fam.file = NULL,
                             exec = "plink",
                             num.threads,
                             memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_int(window.size, lower = 2, add = assertions)
  checkmate::assert_int(step.size, lower = 1, add = assertions)
  checkmate::assert_number(threshold, lower = 0, upper = 1, add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {  
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE) 
  }  
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) { 
    memory = max(5000,       
                 min(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000,                   
                     num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE),                
                 na.rm = TRUE)  
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             "--keep-allele-order",
             "--allow-extra-chr", "0",
             "--indep-pairwise", sprintf("%dkb", window.size), sprintf("%dkb",step.size), threshold,
             "--out", output.prefix, ...)
  )
  
}

#' Sex imputation with PLINK
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param f.values           [\code{numeric(2)}]\cr
#'                           \code{numeric} of length \code{2} with minimum F coefficients to impute females (\code{[1]}) or males (\code{[2]}).
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_CPU} and \code{num.threads * SLURM_MEM_PER_NODE} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate tools
#'
plink_sex_imputation <- function(bfile, output.prefix, 
                                 f.values, ...,
                                 bed.file = NULL, bim.file = NULL, fam.file = NULL,
                                 exec = "plink",
                                 num.threads,
                                 memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_numeric(f.values, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, all.missing = FALSE, len = 2, null.ok = FALSE, add = assertions)
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {     
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), na.rm = TRUE)   
  }   
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  if (missing(memory)) {     
    memory = max(5000,                   
                 min(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000,                       num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE),                   
                 na.rm = TRUE)   
  }   
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             "--impute-sex", f.values,
             "--keep-allele-order",
             "--allow-extra-chr", "0",
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
}

#' Merge multiple PLINK datasets
#'
#' @param bfile              [\code{string}]\cr
#'                           The basename of the binary PLINK files.
#' @param output.prefix      [\code{string}]\cr
#'                           The basename of the new binary PLINK files.
#' @param merge.list         [\code{string}]\cr
#'                           File listing other PLINK datasets to be merged.
#' @param merge.mode         [\code{int}]\cr
#'                           Merge mode.
#' @param ...                [\code{character}]\cr
#'                           Additional arguments passed to PLINK.
#' @param bed.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param bim.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param fam.file           [\code{string}]\cr
#'                           Alternative to \code{bfile} interface. Specify \code{bed}, \code{bim} and \code{fam} files individually.
#' @param exec               [\code{string}]\cr
#'                           Path of PLINK executable.
#' @param num.threads        [\code{int}]\cr
#'                           Number of CPUs usable by PLINK.
#'                           Default is determined by SLURM environment variables and at least 1.
#' @param memory             [\code{int}]\cr
#'                           Memory for PLINK in Mb.
#'                           Default is determined by minimum of SLURM environment variables \code{SLURM_MEM_PER_NODE} and \code{num.threads * SLURM_MEM_PER_CPU} and at least 5000.
#'
#' @details See PLINK manual \url{https://www.cog-genomics.org/plink/1.9/}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
#' @import checkmate tools
#'
plink_merge_list <- function(bfile, output.prefix, 
                                 merge.list, merge.mode, ...,
                                 bed.file = NULL, bim.file = NULL, fam.file = NULL,
                                 exec = "plink",
                                 num.threads,
                                 memory) {
  
  assertions <- checkmate::makeAssertCollection()
  
  if (!missing(bfile)) {
    checkmate::assert_file(sprintf("%s.bed", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.bim", bfile), add = assertions)
    checkmate::assert_file(sprintf("%s.fam", bfile), add = assertions)
    input <- sprintf("--bfile %s", bfile)
  } else {
    checkmate::assert_file(bed.file, add = assertions)
    checkmate::assert_file(bim.file, add = assertions)
    checkmate::assert_file(fam.file, add = assertions)
    input <- sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file)
  }
  
  checkmate::assert_string(output.prefix, add = assertions)
  checkmate::assert_directory(dirname(output.prefix), add = assertions)
  
  checkmate::assert_file(merge.list, add = assertions)
  
  if (!missing(merge.mode)) {
  checkmate::assert_int(merge.mode, lower = 1, upper = 7, add = assertions)
    merge.mode <- sprintf("--merge-mode %d", merge.mode)
  } else {
    merge.mode <- ""
  }
  
  assert_command(exec, add = assertions)
  
  if (missing(num.threads)) {
    num.threads <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                       na.rm = TRUE)
  }
  checkmate::assert_int(num.threads, lower = 1, add = assertions)
  
  if (missing(memory)) {
    memory = max(5000, 
                 min(as.integer(Sys.getenv("SLURM_MEM_PER_NODE")) - 1000, 
                     num.threads * as.integer(Sys.getenv("SLURM_MEM_PER_CPU")) - 1000, na.rm = TRUE), 
                 na.rm = TRUE)
  }
  checkmate::assert_int(memory, lower = 1000, add = assertions)
  
  checkmate::reportAssertions(assertions)
  
  system_call(
    bin = exec,
    args = c(input,
             "--threads", num.threads,
             "--memory", memory,
             "--merge-list", merge.list,
             merge.mode,
             "--keep-allele-order",
             "--allow-extra-chr", "0",
             "--make-bed",
             "--out", output.prefix, ...)
  )
  
}


#' Genotype harmonization
#' 
#' @description Harmonization of genotype data stored using different file formats with different and potentially unknown strands. Linkage disequilibrium (LD) patterns are used to determine the correct strand GC and AT SNPs. This is a simple wrapper for \code{GenotypeHarmonizer}.
#' 
#' @details \url{https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer}
#'
#' @param input                   [\code{string}]\cr
#'                                The base path of the data to align. The extensions are determined based on the input data type.
#' @param ref                     [\code{string}]\cr
#'                                The base path of the reference data used for alignment. The extensions are determined based on the input data type. If not specified the input data is simply converted to the specified output type.
#' @param output                  [\code{string}]\cr
#'                                The base path of the output data.
#' @param input.type              [\code{string}]\cr
#'                                The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path.
#' @param ref.type                [\code{string}]\cr
#'                                The input data type. If not defined will attempt to automatically select the first matching dataset on the specified path.
#' @param output.type             [\code{string}]\cr
#'                                The output data type. Defaults to --input.type.
#' @param input.prob              [\code{number}]\cr
#'                                The minimum posterior probability to call genotypes in the input data. Defaults to \code{0.4}.
#' @param force.chr               [\code{int} or \code{X, Y, MT}]\cr
#'                                \code{SHAPEIT2} does not output the sequence name in the first column of the haplotype file and for \code{GEN} files this can also be the case. Use this option to force the chromosome for all variants. This option is only valid in combination with \code{input.type} \code{SHAPEIT2} and \code{input.type} \code{GEN}.
#' @param call.rate.filter        [\code{number}]\cr
#'                                The minimum call rate to include variant from input data.
#' @param chr.filter              [\code{int} or \code{X, Y, MT}]\cr
#'                                Filter input data on chromosome.
#' @param hwe.filter              [\code{number}]\cr
#'                                The minimum hardy weinberg equilibrium p-value to include variant from input data.
#' @param maf.filter              [\code{number}]\cr
#'                                The minimum minor allele frequency to include variant from input data.
#' @param sample.filter.list      [\code{string}]\cr
#'                                Path to file with samples IDs to include from input data. For plink data and oxford sample files only the sample id (column 2) is used.
#' @param variant.filter.list     [\code{string}]\cr
#'                                Path to file with variant IDs to include from input data.
#' @param mach.r2.filter          [\code{number}]\cr
#'                                The minimum MACH R2 measure to include SNPs.
#' @param variant.pos.filter.list [\code{string}]\cr
#'                                Path to file with variant \code{CHR\\tPOS} or \code{CHR:POS} to include from input data.
#' @param ambiguous.snp.filter    [\code{flag}]\cr
#'                                Filter out ambiguous SNPs (A/T, C/G) SNPs.
#' @param update.id               [\code{flag}]\cr
#'                                Update the variant identifiers using the reference data. The identifiers of the output data will be the same as the reference data.
#' @param min.ld                  [\code{number}]\cr
#'                                The minimum LD (\eqn{r^2}) between the variant to align and potential supporting variants. Defaults to \code{0.3}.
#' @param min.variants            [\code{int}]\cr
#'                                The minimum number of supporting variant before before we can do an alignment. Defaults to \code{3}.
#' @param variants                [\code{int}]\cr
#'                                Number of flanking variants to consider. Defaults to \code{100}.
#' @param check.ld                [\code{flag}]\cr
#'                                Also check the LD structure of non AT and non GC variants. Variants that do not pass the check are excluded.
#' @param maf.align               [\code{number}]\cr
#'                                If there are not enough variants in LD and the minor allele frequency (MAF) of a variant <= the specified value in both study as in reference then the minor allele can be used as a backup for alignment. Defaults to \code{0}.
#' @param update.reference.allele [\code{flag}]\cr
#'                                Make sure the output data uses the same reference allele as the reference data set.
#' @param exec                    [\code{string}]\cr
#'                                Path of \code{GenotypeHarmonizer} executable. You can also give a \code{JAVA} call: \code{java -Xmx5g -jar <path/to/GenotypeHarmonizer.jar>}.
#'
#' @return Captured system output as \code{character} vector.
#' @export
#'
harmonize_genotypes <- function(input, ref, output,
                                input.type, ref.type, output.type,
                                input.prob = 0.4, force.chr, call.rate.filter = 0, chr.filter, hwe.filter = 1, maf.filter = 0, sample.filter.list, variant.filter.list, mach.r2.filter = 0, variant.pos.filter.list, ambiguous.snp.filter = FALSE,
                                update.id = FALSE, min.ld = 0.3, min.variants = 3, variants = 100, check.ld = FALSE, maf.align = 0, update.reference.allele = FALSE,
                                exec = "GenotypeHarmonizer") {
  
  input.types <- list(
    PED_MAP = c("ped", "map"),
    VCF = c("vcf.gz", "vcf.gz.tbi"),
    PLINK_BED = c("bed", "bim", "fam"),
    SHAPEIT2 = c("haps", "sample"),
    GEN = c("gen", "sample")
  )
  
  output.types <- list(
    PED_MAP = c("ped", "map"),
    PLINK_BED = c("bed", "bim", "fam"),
    SHAPEIT2 = c("haps", "sample"),
    GEN = c("gen", "sample")
  )
  
  assertions <- checkmate::makeAssertCollection()
  
  # Basic I/O arguments ----
  checkmate::assertTRUE(
    any(
      sapply(input.types, 
             function(extensions) {
               checkmate::testFile(sprintf("%s.%s", input, extensions))
             }
      )
    ), add = assertions
  )
  
  if (!missing(ref)) {
    checkmate::assertTRUE(
      any(
        sapply(input.types, 
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
  
  checkmate::assertDirectory(dirname(output), add = assertions)
  
  if (!missing(input.type)) {
    checkmate::assertChoice(input.type, names(input.types))
    input.type <- sprintf("--inputType %s", input.type)
  } else {
    input.type <- ""
  }
  if (!missing(ref.type)) {
    checkmate::assertChoice(ref.type, names(input.types))
    ref.type <- sprintf("--refType %s", ref.type)
  } else {
    ref.type <- ""
  }
  if (!missing(output.type)) {
    checkmate::assertChoice(output.type, names(output.types))
    output.type <- sprintf("--outputType %s", output.type)
  } else {
    output.type <- ""
  }
  
  # Input filtering arguments ----
  checkmate::assertNumber(input.prob, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
  input.prob <- sprintf("--inputProb %f", input.prob)
  if (input.type %in% c("SHAPEIT2", "GEN")) {
    checkmate::assertChoice(as.character(force.chr), c(1:25, "X", "Y", "MT"), add = assertions)
    force.chr <- sprintf("--forceChr %s", as.character(force.chr))
  } else {
    force.chr <- ""
  }
  checkmate::assertNumber(call.rate.filter, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
  call.rate.filter <- sprintf("--callRateFilter %f", call.rate.filter)
  if (!missing(chr.filter)) {
    checkmate::assertChoice(as.character(chr.filter), c(1:25, "X", "Y", "MT"), add = assertions)
    chr.filter <- sprintf("--chrFilter %s", as.character(chr.filter))
  } else {
    chr.filter <- ""
  }
  checkmate::assertNumber(maf.filter, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
  maf.filter <- sprintf("--mafFilter %f", maf.filter)
  if (!missing(sample.filter.list)) {
    checkmate::assertFile(sample.filter.list, add = assertions)
    sample.filter.list <- sprintf("--sampleFilterList %s", sample.filter.list)
  } else {
    sample.filter.list <- ""
  }
  if (!missing(variant.filter.list)) {
    checkmate::assertFile(variant.filter.list, add = assertions)
    variant.filter.list <- sprintf("--variantFilterList %s", variant.filter.list)
  } else {
    variant.filter.list <- ""
  }
  checkmate::assertNumber(mach.r2.filter, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
  mach.r2.filter <- sprintf("--machR2Filter %f", mach.r2.filter)
  if (!missing(variant.pos.filter.list)) {
    checkmate::assertFile(variant.pos.filter.list, add = assertions)
    variant.pos.filter.list <- sprintf("--variantPosFilterList %s", variant.pos.filter.list)
  } else {
    variant.pos.filter.list <- ""
  }
  checkmate::assertFlag(ambiguous.snp.filter, add = assertions)
  if (ambiguous.snp.filter) {
    ambiguous.snp.filter <- "--ambiguousSnpFilter"
  } else {
    ambiguous.snp.filter <- ""
  }
  
  # Strand alignment arguments ----
  checkmate::assertFlag(update.id, add = assertions)
  if (update.id) {
    update.id <- "--update-id"
  } else {
    update.id <- ""
  }
  checkmate::assertNumber(min.ld, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
  min.ld <- sprintf("--min-ld %f", min.ld)
  checkmate::assertInt(min.variants, lower = 1, null.ok = FALSE, add = assertions)
  min.variants <- sprintf("--min-variants %f", min.variants)
  checkmate::assertInt(variants, lower = 1, null.ok = FALSE, add = assertions)
  variants <- sprintf("--variants %f", variants)
  checkmate::assertFlag(check.ld, add = assertions)
  if (check.ld) {
    check.ld <- "--check-ld"
  } else {
    check.ld <- ""
  }
  checkmate::assertNumber(maf.align, lower = 0, upper = 1, finite = TRUE, null.ok = FALSE, add = assertions)
  maf.align <- sprintf("--mafAlign %f", maf.align)
  checkmate::assertFlag(update.reference.allele, add = assertions)
  if (update.reference.allele) {
    update.reference.allele <- "--update-reference-allele"
  } else {
    update.reference.allele <- ""
  }
  
  if (grepl("java", exec)) {
    exec <- unlist(strsplit(exec, " "))
    n <- length(exec)
    assertCommand(exec[1], add = assertions)
    checkmate::assertFile(exec[n])
    java_opts <- exec[2:n]
    exec <- exec[1]
  } else {
    assertCommand(exec, add = assertions)
    java_opts <- ""
  }
  
  checkmate::reportAssertions(assertions)
  
  system_call(exec,
              args = c(java_opts,
                       "--input", input,
                       input.type,
                       ref,
                       ref.type,
                       "--output", output,
                       output.type,
                       input.prob, force.chr, 
                       call.rate.filter,
                       chr.filter,
                       hwe.filter,
                       maf.filter, 
                       sample.filter.list, 
                       variant.filter.list, 
                       mach.r2.filter,
                       variant.pos.filter.list,
                       ambiguous.snp.filter,
                       update.id, 
                       min.ld , 
                       min.variants,
                       variants ,
                       check.ld ,
                       maf.align , 
                       update.reference.allele))
  
}
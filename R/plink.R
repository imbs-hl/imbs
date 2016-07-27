
#' Write SNP or sample IDs to file
#'
#' @param ids Vector of SNP or Sample IDs
#' @param file File name to write to
#' @param sample Logical. Should be TRUE for samples and FALSE for SNPs
#'
#' @return NULL
#' @export
write_plink_ids <- function(ids, file, sample = TRUE){
  if (sample) {
      out <- cbind(ids)
  } else {
      out <- ids
  }
  write.table(out, file, row.names = FALSE, col.names = FALSE,
              quote = FALSE)
  }
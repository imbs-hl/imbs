
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
#' @export
write_plink_ids <- function(ids, file, update, sample = TRUE) {
  if (sample) {
    out <- cbind(ids, ids)
  } else {
    out <- ids
  }
  if (update | !file.exists(file)) {
    write.table(
      out,
      file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
}
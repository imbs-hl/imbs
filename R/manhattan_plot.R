
#' Create a manhattan plot
#' 
#' Create a manhattan plot from a \code{data.frame}.
#'
#' @param data         [\code{\link{data.frame}}]\cr
#'                     A \code{data.table} object containing at least 
#'                     information about chromosome, position and p value.
#' @param threshold    [\code{number}]\cr
#'                     Single numeric to draw a threshold line at. On p value 
#'                     scale. Optional.
#' @param thin         [\code{logical}]\cr
#'                     Should the data be thinned?
#' @param thin.param   [\code{list}]\cr
#'                     A list with the following elements: 
#'                     \code{logp.thin.thresh} (\code{number}), 
#'                     \code{logp.thin.rate} (\code{number}) and 
#'                     \code{buffer} (\code{integer}) indicating at which rate 
#'                     points with a low log p value should be removed and how 
#'                     many base pairs at the start and end of each chromosome 
#'                     should not be excluded from thinning.
#' @param header       [\code{list}]\cr
#'                     A list with the following elements:
#'                     \code{chr}, \code{pos} and \code{pval} (all 
#'                     \code{character}) indicating the columns to use for 
#'                     chromosome, position and p value information.
#'                     
#' @return A \code{\link[ggplot2]{ggplot2}} object.
#' @export
#'
#' @import data.table ggplot2
manhattan_plot <- function(data,
                           threshold,
                           thin = FALSE,
                           thin.param = list(logp.thin.thresh = 5,
                                             logp.thin.rate = 0.1,
                                             buffer = 1e6),
                           header = list(chr = "CHR",
                                         pos = "BP",
                                         pval = "P")) {
  
  checkmate::assertDataFrame(data, min.cols = 3)
  checkmate::assertFlag(thin)
  if(thin) {
    checkmate::assertList(thin.param, len = 3, unique = TRUE, names = "strict", types = "numeric")
    checkmate::assertSubset(c("logp.thin.thresh", "logp.thin.rate", "buffer"), names(thin.param))
  }
  checkmate::assertList(header, len = 3, unique = TRUE, names = "strict", types = "character")
  checkmate::assertSubset(names(header), c("chr", "pos", "pval"))
  checkmate::assertSubset(unlist(header), names(data))
  if(!missing(threshold)) {
    checkmate::assertNumber(threshold, lower = 0, upper = 1)
  }
  
  data <- data.table::as.data.table(copy(data))
  
  pos <- header$pos
  chr <- header$chr
  pval <- header$pval
  
  # Get maximum base pair position per chromosome
  maxBP <- data.table::data.table("chr" = 1:22, key = "chr")
  data.table::setnames(maxBP, old = "chr", new = chr)
  data.table::setkeyv(maxBP, chr)
  data.table::setkeyv(data, chr)
  maxBP <- data[, .(maxBP = max(get(pos))), by = chr][maxBP]
  maxBP[is.na(maxBP), maxBP := 0]
  
  # Shift the chromosome number to add the previous maximum base pair position
  # to each chromosome
  maxBP[, CHR := get(chr)+1]
  maxBP[, eval(chr) := NULL]
  
  # Chromosome 1 doesn't need to be shifted
  maxBP[CHR == 23, c("CHR", "maxBP") := list(1, 0)]
  data.table::setkeyv(maxBP, "CHR")
  
  # Define the shift per chromosome
  maxBP[, shift := cumsum(as.numeric(maxBP))]
  maxBP[, maxBP := NULL]
  
  # Merge the genomic base pair position with GWAS results
  manhattanData <- data[, .SD, .SDcols = c(pos, chr, pval)][maxBP]
  manhattanData[, GBP := shift + get(pos)]
  manhattanData[, CHR := factor(get(chr))]
  manhattanData[, LOGP := -log10(get(pval))]
  
  # Get number of SNPs tested for Bonferroni adjusted threshold
  snpCount <- manhattanData[, .N]
  bThresh <- -log10(0.05/snpCount)
  
  if(thin) {
    manhattanData[, THINABLE := (get(pos) > (min(get(pos)) + thin.param$buffer)) &
                    (get(pos) < (max(get(pos)) - thin.param$buffer)) &
                    (LOGP < thin.param$logp.thin.thresh), by = get(chr)]
    manhattanData <-manhattanData[c(sample(which(THINABLE), sum(THINABLE, na.rm = TRUE)*thin.param$logp.thin.rate),
                                    which(!THINABLE))]
  }
  
  p <- ggplot2::ggplot(data = manhattanData, aes(x = GBP, y = -log10(get(pval)))) +
    geom_point(aes(col = CHR), size = 0.5) +
    geom_abline(slope = 0, intercept = bThresh, color = "red") +
    scale_color_manual(values = rep(c("black", "darkgray"), 15)) +
    scale_x_continuous(breaks = setkey(manhattanData[, .(midCHR = mean(GBP)), by = CHR], CHR)$midCHR,
                       labels = sort(unique(manhattanData$CHR)),
                       minor_breaks = c(manhattanData$shift, manhattanData[, max(GBP)])) +
    theme(legend.position = "none") +
    ylab(expression(-log[10](p))) +
    xlab("Chromosome")
  
  if(!missing(threshold)) {
    p <- p + geom_abline(slope = 0, intercept = -log10(threshold), color = "blue")
  } else {
    threshold <- 1
  }
  
  p <- p + ylim(c(0, max(c(manhattanData$LOGP,
                           bThresh,
                           -log10(threshold)))))
  
  return(p)
  
}
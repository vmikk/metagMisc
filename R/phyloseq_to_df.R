
#' @title Convert phyloseq object to data frame (for exporting).
#'
#' @param physeq A phyloseq-class object
#' @param addtax Add taxonomy to the result table
#' @param addtot Add total OTU abundance to the result table
#' @param abund_sort Sort OTUs by total abundance
#'
#' @return Data frame with OTU taxonomy and abundance
#' @export
#' @seealso \code{\link{phyloseq-class}}
#'
#' @examples
#' data(GlobalPatterns)
#' GlobalPatternsDF <- phyloseq_to_df(GlobalPatterns)
#' str(GlobalPatternsDF)
#'
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, abund_sort = T){

  require(phyloseq)

  if(addtax==TRUE) { res <- data.frame(OTU = taxa_names(physeq), tax_table(physeq), otu_table(physeq), stringsAsFactors = F) }
  if(addtax==FALSE){ res <- data.frame(OTU = taxa_names(physeq), otu_table(physeq), stringsAsFactors = F) }

  if(abund_sort == TRUE){
    otus <- res[, which(colnames(res) %in% sample_names(physeq))]
    res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
  }

  if(addtot == TRUE){
    res$Total <- rowSums(res[,which(colnames(res) %in% sample_names(physeq))])
  }

  rownames(res) <- NULL
  return(res)
}

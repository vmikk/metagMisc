
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
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, sorting = "abundance"){

  # require(phyloseq)

  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }

  ## Prepare data frame
  res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)

  ## Add taxonomy
  if(addtax == TRUE){
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    res <- cbind(res, taxx[match(x = res$OTU, table = rownames(taxx)), ])

    ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
    res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
  }

  ## Reorder OTUs
  if(!is.null(sorting)){

    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }

    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )

      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }

  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }

  rownames(res) <- NULL
  return(res)
}

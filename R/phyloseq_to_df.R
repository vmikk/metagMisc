
#' @title Convert phyloseq object to data frame (for exporting).
#'
#' @param physeq A phyloseq-class object
#' @param addtax Logical; add taxonomy to the result table (default, TRUE)
#' @param addtot Logical; add total OTU abundance to the resulting table (default, FALSE)
#' @param sorting Reorder OTUs according to their total abundance ("abundance", default) or taxonomy (sorting = "taxonomy"); if reordering is not required use 'sorting = NULL'
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
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){

  # require(phyloseq)

  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }

  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }

  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }

    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }

  ## Add taxonomy
  if(addtax == TRUE){
    
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)

    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]

    ## Add taxonomy table to the data
    res <- cbind(res, taxx)

    ## Add max tax rank column
    if(addmaxrank == TRUE){

      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)

      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]

    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
    
    } # end of addmaxrank
  }   # end of addtax

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

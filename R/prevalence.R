
#' @title Compute prevalence of each species.
#'
#' @param physeq A phyloseq-class object
#' @param add_tax Logical, add taxonomy to the results
#'
#' @return Data frame
#'
#' @examples
#' data(esophagus)
#' prevalence(esophagus)
#'
#' data(GlobalPatterns)
#' head( prevalence(GlobalPatterns, add_tax = T) )
#'
prevalence <- function(physeq, add_tax = TRUE){
  prevdf <- apply(X = otu_table(physeq),
                  MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})

  ## Add taxonomy and total read counts to this data.frame
  prevdf <- data.frame(Prevalence = prevdf,
                       TotalAbundance = taxa_sums(physeq))

  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(physeq))
  }
  return(prevdf)
}

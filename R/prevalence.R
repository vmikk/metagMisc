
#' @title Compute prevalence and abundance summary of each species/OTU.
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

  ## Check if taxa are rows
  trows <- taxa_are_rows(physeq)

  ## Extract OTU table
  otutab <- as.data.frame(otu_table(physeq))

  ## Transpose OTU table (species should be arranged by rows)
  if(trows == FALSE){
    otutab <- t(otutab)
  }

  ## Estimate prevalence (number of samples with OTU present)
  prevdf <- apply(X = otutab,
                  # MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
                  MARGIN = 1,
                  FUN = function(x){sum(x > 0)})

  ## Add total and average read counts per OTU
  prevdf <- data.frame(Prevalence = prevdf,
                       TotalAbundance = taxa_sums(physeq),
                       MeanAbundance = rowMeans(otutab),
                       MedianAbundance = apply(otutab, 1, median))

  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(physeq))
  }
  return(prevdf)
}

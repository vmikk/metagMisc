phyloseq_otu_to_df <- function(physeq, taxa_as_rows = TRUE){

  res <- as.data.frame(phyloseq::otu_table(physeq))
  current_orientation <- phyloseq::taxa_are_rows(physeq)

  ## Transpose if needed to match desired orientation
  if(taxa_as_rows == TRUE){
    ## Want taxa as rows
    if(current_orientation == FALSE){
      ## Currently samples as rows, need to transpose
      res <- t(res)
    }
    ## If current_orientation == TRUE, no change needed
  } else {
    ## Want samples as rows (taxa as columns)
    if(current_orientation == TRUE){
      ## Currently taxa as rows, need to transpose
      res <- t(res)
    }
    ## If current_orientation == FALSE, no change needed
  }

  return(res)
}

phyloseq_SRS <- function(physeq, Cmin, drop_zeros = FALSE, set_seed = TRUE, seed = 1, ...){

  ## Extract OTU table
  data <- as.data.frame(otu_table(physeq))

  ## Transpose data (columns should be samples and rows are the OTU counts)
  if(!taxa_are_rows(physeq)){
    data <- t(data)
  }

  ## Perform SRS-normalization
  res <- SRS(data = data, Cmin = Cmin, ...)
  rownames(res) <- rownames(data)

  ## Transpose back
  if(!taxa_are_rows(physeq)){
    res <- t(res)
  }

  ## Replace OTU table
  otu_table(physeq) <- otu_table(res, taxa_are_rows = taxa_are_rows(physeq))

  ## Remove zero-OTUs
  if(drop_zeros == TRUE){
    physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  }

  return(physeq)
}

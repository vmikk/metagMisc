
## Standardize OTU abundance (with vegan)
phyloseq_standardize_otu_abundance <- function(physeq, method="total", ...){

  ## Check the orientation of the OTU table
  trows <- taxa_are_rows(physeq)
  if(trows == TRUE){ marg <- 2 } else { marg <- 1 }

  ## Extact OTU table
  comm <- as(object = otu_table(physeq), Class = "matrix")

  ## Standardize community table
  comm_std <- decostand(comm, method, MARGIN = marg, ...)

  ## Replace old otu_table with the new one
  otu_table(physeq) <- otu_table(comm_std, taxa_are_rows = trows)

  return(physeq)
}

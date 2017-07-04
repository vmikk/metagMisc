
## Count number of gene copies for each feature (e.g., KEGG) in each sample
predict_metagenomes <- function(otu_tab, func_tab, NSTI_present = TRUE, rel_abund = TRUE){

  # extract OTUs from functional table
  funcs <- as.matrix( func_tab[match(x = otu_tab[, 1], table = func_tab[,1]), -1] )
  otus <- as.matrix( otu_tab[, -1] )

  # remove NSTI scores if present
  if(NSTI_present == TRUE){ funcs <- funcs[, -ncol(funcs)] }

  # standardize OTU count to relative abundance
  if(rel_abund == TRUE){ otus <- decostand(otus, method = "total", MARGIN = 2) }

  # Multiply gene counts for each OTU by the abundance of that OTU in each each sample
  # and sum across all OTUs
  res <- crossprod(otus, funcs)
  return(res)
}

## Count number of gene copies for each feature (e.g., KEGG) in each sample
predict_metagenomes <- function(otu_tab, func_tab, NSTI_present = TRUE, rel_abund = TRUE, add_sub_tabs = FALSE){

  # Subset OTUs for the features present in the functional table
  otu_vs_func <- otu_tab[, 1] %in% func_tab[, 1]
  if(any(!otu_vs_func)){
    cat("Warning: OTUs that are no present in the table with functional features will be remved.\n")
    otu_tab <- otu_tab[otu_vs_func, ]
  }
  otus <- as.matrix( otu_tab[, -1] )

  # Extract OTUs from functional table
  match_otus <- match(x = otu_tab[, 1], table = func_tab[,1])
  if(any(is.na(match_otus))){ match_otus <- na.omit(match_otus) }
  funcs <- as.matrix( func_tab[match_otus, -1] )

  # Remove NSTI scores if present
  if(NSTI_present == TRUE){ funcs <- funcs[, -ncol(funcs)] }

  # Standardize OTU count to relative abundance
  if(rel_abund == TRUE){ otus <- decostand(otus, method = "total", MARGIN = 2) }

  # Multiply gene counts for each OTU by the abundance of that OTU in each each sample
  # and sum across all OTUs
  res <- crossprod(otus, funcs)

  # Add subsetted tables to the results
  if(add_sub_tabs == TRUE){
    attr(res, which = "OTUs") <- otu_tab
    attr(res, which = "Funcs") <- func_tab[match_otus, ]
  }

  return(res)
}

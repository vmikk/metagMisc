
## Function to rename species with tax ranks
phyloseq_rename_with_tax <- function(physeq, taxrank = "Family"){

  ## Extract taxonomy table
  txx <- as.data.frame(tax_table(physeq), stringsAsFactors = F)

  ## Check if
  if(!taxrank %in% colnames(txx)){
    stop(taxrank, " is not in columns of tax_table.\n")
  }

  ## Extract new names from tax table
  newnames <- txx[, taxrank]

  ## Check for the uniqueness of taxa names
  if(nrow(txx) != length(unique(newnames))){
    stop("Taxonomy column should contain only unique names.\n")
  }

  ## Rename taxa
  # taxa_names(physeq) <- newnames[ match(x = taxa_names(physeq), table = rownames(txx)) ]  # reorder names
  taxa_names(physeq) <- newnames   # names should be in the same order

  return(physeq)
}


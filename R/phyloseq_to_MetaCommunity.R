
## Convert phyloseq object to MetaCommunity class from entropart package
phyloseq_to_MetaCommunity <- function(physeq, wei = NULL, ...){

  require(entropart)

  ## Extract OTU table
  otus <- as.data.frame(otu_table(physeq))
  if(taxa_are_rows(physeq) == FALSE){ otus <- t(otus) }

  ## Prepare equal community weights
  if(is.null(wei)){
    wei <- rep(1, ncol(otus))
  }
  Weights <- data.frame(Communities = colnames(otus), Weights = wei)

  ## Convert phyloseq class to MetaCommunity class
  MC <- MetaCommunity(otus, Weights)

  return(MC)
}

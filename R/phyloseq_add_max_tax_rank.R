
## Add the lowest level of taxonomic classification to the taxonomy table
phyloseq_add_max_tax_rank <- function(physeq, abbreviate = F, ...){
  # ... passed to abbreviate_taxa_names

  ## Determine the lowest level of taxonomic classification
  sp_ranks <- get_max_taxonomic_rank(physeq)                                   # [metagMisc]
  sp_ranks$RankName <- as.character(sp_ranks$RankName)

  ## Extract the names of the lowest ranks
  sp_ranks <- plyr::adply(
    .data = sp_ranks,
    .margins = 1,
    .fun = function(z){ data.frame(LowestTaxon = z[, z$RankName]) })

  ## Abbreviate taxa names and make them unique
  if(abbreviate == TRUE){
    sp_ranks$LowestTaxon <- abbreviate_taxa_names(sp_ranks$LowestTaxon, ...)   # [metagMisc]
  }

  ## Replace taxonomy table with the new one
  rownames(sp_ranks) <- sp_ranks$TaxaName
  sp_ranks$TaxaName <- NULL
  tax_table(physeq) <- tax_table(as.matrix(sp_ranks))

  return(physeq)
}

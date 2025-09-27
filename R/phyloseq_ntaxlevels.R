
## Estimate number of unique (non-NA) taxonomic levels
phyloseq_ntaxlevels <- function(physeq, add_all_samps = TRUE){

  ## Melt phyloseq data
  mm <- speedyseq::psmelt(physeq)
  setDT(mm)

  ## Convert all factors to character
  factor_cols <- names(mm)[sapply(mm, is.factor)]
  if(length(factor_cols) > 0) {
    mm[, (factor_cols) := lapply(.SD, as.character), .SDcols = factor_cols]
  }

  ## Get taxonomic ranks
  tranks <- phyloseq::rank_names(physeq)
  
  ## Replace NAs in taxonomy
  for(col in tranks) {
    mm[is.na(get(col)), (col) := "NA"]
  }

  ## Reorder columns (taxonomy first)
  other_cols <- setdiff(names(mm), tranks)
  setcolorder(mm, c(tranks, other_cols))

  ## Create combined tax ranks
  tranks_cmb <- character(length(tranks) - 1)
  for(i in seq_len(length(tranks) - 1)) {
    tranks_cmb[i] <- paste(tranks[1:(i+1)], collapse = "_")
  }
  
  ## Create lookup table for ranks
  tranks_cmbs <- data.table::data.table(
    TaxRank = tranks[-1], 
    CombinedRank = tranks_cmb
  )

  ## Create combined taxonomic columns
  for(i in seq_len(length(tranks) - 1)) {
    rank_cols <- tranks[1:(i+1)]
    combined_col <- tranks_cmb[i]
    mm[, (combined_col) := do.call(paste, c(.SD, sep = "_")), .SDcols = rank_cols]
  }

  ## For each taxonomic rank and sample, count unique non-zero, non-NA taxa
  res_list <- vector("list", nrow(tranks_cmbs))
  
  for(i in seq_len(nrow(tranks_cmbs))) {
    tax_rank <- tranks_cmbs$TaxRank[i]
    combined_rank <- tranks_cmbs$CombinedRank[i]
    
    ## Count unique taxa per sample for this rank
    temp_res <- mm[Abundance > 0 & get(tax_rank) != "NA", 
                   .(N.tax.levels = uniqueN(get(combined_rank))), 
                   by = Sample]
    temp_res[, TaxRank := tax_rank]
    res_list[[i]] <- temp_res
  }
  
  ## Combine all results
  res <- rbindlist(res_list, use.names = TRUE)
  
  ## Add results for all samples combined if requested
  if(add_all_samps == TRUE && uniqueN(mm$Sample) > 1) {
    res_all_list <- vector("list", nrow(tranks_cmbs))
    
    for(i in seq_len(nrow(tranks_cmbs))) {
      tax_rank <- tranks_cmbs$TaxRank[i]
      combined_rank <- tranks_cmbs$CombinedRank[i]
      
      ## Count unique taxa across all samples for this rank
      n_unique <- mm[Abundance > 0 & get(tax_rank) != "NA", 
                     uniqueN(get(combined_rank))]
      
      res_all_list[[i]] <- data.table(
        TaxRank = tax_rank,
        Sample = "All_samples", 
        N.tax.levels = n_unique
      )
    }
    
    res_all <- rbindlist(res_all_list, use.names = TRUE)
    res <- rbindlist(list(res, res_all), use.names = TRUE)
  }
  
  ## Reorder columns
  setcolorder(res, c("TaxRank", "Sample", "N.tax.levels"))
  
  return(res)
}

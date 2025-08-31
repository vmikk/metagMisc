
## Estimate number of unique (non-NA) taxonomic levels
phyloseq_ntaxlevels <- function(physeq, add_all_samps = TRUE){

  ## Melt phyloseq data object into large data.frame
  mm <- phyloseq::psmelt(physeq)

  ## Convert all factors to character
  i <- sapply(mm, is.factor)
  mm[i] <- lapply(mm[i], as.character)
  rm(i)

  ## Replace NAs in taxonomy
  tranks <- phyloseq::rank_names(physeq)
  for(i in tranks){
    nas <- is.na( mm[, i] )
    if(any(nas)){
      mm[, i][ which(nas) ] <- "NA"
    }
    rm(nas)
  }
  rm(i)

  ## Reorder columns (taonomy first)
  mm <- mm[, c(tranks, colnames(mm)[ which(!colnames(mm) %in% tranks)]) ]

  ## Create combined tax ranks
  tranks_ids <- list()
  for(i in 1:(length(tranks)-1)){
    tranks_ids[[i]] <- 1:(i+1)
  }
  rm(i)


  tranks_cmb <- laply(.data = tranks_ids, .fun = function(z){ paste(tranks[z], collapse = "_") })  # rank names
  names(tranks_ids) <- tranks_cmb

  tranks_cmbs <- data.frame(TaxRank = tranks[-1], CombinedRank = tranks_cmb, stringsAsFactors = F)

  tranks_ccc <- llply(.data = tranks_ids, .fun = function(z){
    data.frame(Cmb = apply(mm[,z], 1, paste, collapse= "_"), stringsAsFactors = F)
  })

  tranks_ccc <- do.call(cbind, tranks_ccc)
  colnames(tranks_ccc) <- tranks_cmb
  mm <- cbind(mm, tranks_ccc)
  rm(tranks_ccc, tranks_cmb, tranks_ids)


  ## Count number of OTUs for each sample
  count_taxlevels <- function(z, TaxRank = TaxRank){

    ## Remove zero-OTUs
    zero_otu <- z$Abundance > 0
    if(any(zero_otu)){ z <- z[zero_otu, ] }

    ## Remove NA tax levels
    na_tax <- z[, TaxRank] %in% "NA"
    if(any(na_tax)){ z <- z[!na_tax, ] }

    ## Find which column to take for the aggregation
    TaxRankCmb <- tranks_cmbs[ which(tranks_cmbs$TaxRank == TaxRank), "CombinedRank"]

    ## Count number of unique levels per taxonomic rank selected
    rez <- data.frame(N.tax.levels = length(unique(z[, TaxRankCmb])), stringsAsFactors = F)

    return(rez)
  }

  ## Count number of tax levels for each sample and each rank
  res <- mdply(
    .data = data.frame(TaxRank = tranks_cmbs$TaxRank, stringsAsFactors = F),
    .fun = function(...){ plyr::ddply(.data = mm, .variables = "Sample", .fun = count_taxlevels, ...) })  # TaxRank = TaxRank

  ## Count for all samples (if there are > 1 samples)
  if(add_all_samps == TRUE & length(unique(mm$Sample)) > 1){
    res_all <- mdply(
      .data = data.frame(TaxRank = tranks_cmbs$TaxRank, stringsAsFactors = F),
      .fun = function(...){ count_taxlevels(mm, ...) })
    res_all$Sample <- "All_samples"

    ## Combine results
    res <- rbind(res, res_all)
  }
  return(res)
}

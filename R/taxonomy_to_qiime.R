
## Prepare taxonomy in QIIME-style
taxonomy_to_qiime <- function(x, dropNA = TRUE, add_OTUID = TRUE){
  # x = data.frame or phyloseq with tax_table()
  # dropNA = remove missing tax ranks

  ## If phyloseq is provided - extract OTU table
  if(class(x) %in% c("phyloseq", "taxonomyTable")){
    tx <- as.data.frame(phyloseq::tax_table(x), stringsAsFactors = F)
  }
  ## If data frame is provided
  if(class(x) %in% "data.frame"){
    tx <- x
  }

  ## Column names abbreviations
  tax.levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax.abbr <- c("k", "p", "c", "o", "f", "g", "s")
  taxx <- data.frame(Level = c(tax.levels, tolower(tax.levels)), Abbr = tax.abbr, stringsAsFactors = F)

  ## Rename tax columns
  colnames(tx) <- taxx[match(x = colnames(tx), table = taxx$Level), "Abbr"]

  ## Function to merge a single taxonomy string
  tax_merge <- function(z, dropNA = TRUE){
    # z = data frame with single row

    ## Remove missing ranks
    if(dropNA == TRUE){
      nn <- is.na(z)
      if(any(nn)){ z <- z[, -which(nn)] }
    }

    ## Add tax levels to the rank names
    rez <- paste(colnames(z), z, sep = "__")

    ## Merge tax ranks into a single string
    rez <- paste(rez, collapse=";")
    return(rez)
  } # end of tax_merge


  ## Batch merging
  res <- plyr::aaply(.data = tx, .margins = 1, .fun = tax_merge, dropNA = dropNA, .expand = F)

  ## Add column with OTU name
  if(add_OTUID == TRUE){
    res <- data.frame(OTU_ID = rownames(tx), taxonomy = res, stringsAsFactors = F)
  }

  return(res)
}

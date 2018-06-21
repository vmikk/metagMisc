
#' @title Prepare taxonomy in QIIME-style
#' @description This function merges multiple columns with taxonomic ranks into a single string in QIIME-style.
#' @param x data.frame or phyloseq with \code{\link[phyloseq]{tax_table}} slot
#' @param dropNA Logical; if TRUE, missing tax ranks will be removed
#' @param add_OTUID Logical; if TRUE, OTU name will be added to the first column of the resulting table
#' @param custom_tax_ranks Data frame with two columns, taxonomic ranks (with the same names as columns in x) and the corresponding single letter abbreviations
#' 
#' @return character vector (if 'add_OTUID = FALSE') or data.frame with taxonomy annotations.
#' @export
#'
#' @examples
#' # Load data
#' data("GlobalPatterns")
#' 
#' # Subset data
#' GP <- metagMisc::phyloseq_filter_top_taxa(GlobalPatterns, n = 20)
#' 
#' # Prepare taxonomy in QIIME-style (result = data frame with 2 columns)
#' qtax <- taxonomy_to_qiime(GP)
#' head(qtax)
#' 
#' # Vector with OTU taxonomy
#' qvec <- taxonomy_to_qiime(GP, add_OTUID = FALSE)
#' head(qvec)
#' 
taxonomy_to_qiime <- function(x, dropNA = TRUE, add_OTUID = TRUE, custom_tax_ranks = NULL){

  ## If phyloseq is provided - extract OTU table
  if(class(x) %in% c("phyloseq", "taxonomyTable")){
    tx <- as.data.frame(phyloseq::tax_table(x), stringsAsFactors = F)
  }
  ## If data frame is provided
  if(class(x) %in% "data.frame"){
    tx <- x
  }

  ## Default column names abbreviations
  if(!is.null(custom_tax_ranks)){
    tax.levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    tax.abbr <- c("k", "p", "c", "o", "f", "g", "s")
    taxx <- data.frame(Level = c(tax.levels, tolower(tax.levels)), Abbr = tax.abbr, stringsAsFactors = F)
  } else {
  ## Custom column names abbreviations

    taxx <- custom_tax_ranks
    
    ## Convert factors to character
    taxx <- data.frame(lapply(taxx, as.character), stringsAsFactors=FALSE)
    
    ## Rename columns
    colnames(taxx) <- c("Level", "Abbr")
  }

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

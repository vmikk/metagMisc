
#' @title Split phyloseq object by taxonomic rank.
#' @description This function splits a phyloseq object by a specified taxonomic rank, returning a list of objects whose components each correspond to a taxonomic rank.
#' @param physeq Phyloseq object
#' @param TaxRank Taxonomy rank name (e.g., "Phylum")
#' @param drop_NA Logical, remove ranks with NAs
#'
#' @return List of phyloseq objects
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' # Subset data
#' GP <- subset_taxa(GlobalPatterns, Phylum %in% c("Acidobacteria", "Actinobacteria", "Firmicutes", "Verrucomicrobia"))
#' phyloseq_sep_tax(GP, TaxRank = "Phylum")
#'
phyloseq_sep_tax <- function(physeq, TaxRank = "Phylum", drop_NA = FALSE){
  ## TO DO - change from for-loop to plyr
  ##       - verify that the taxonomic names are unique across higher ranks (e.g. no same class-names within different phylums)

  require(phyloseq)
  require(plyr)

  ## Test if the taxonomic table is present
  if( is.null(tax_table(physeq)) ){
    stop("Nothing subset. No taxonomyTable in physeq.\n")
  }

  ## Test if the specified TaxRank is valid
  if( !TaxRank %in% rank_names(physeq) ){
  	stop("Specified taxonomic rank is missing in the taxonomy table of physeq.\n")
  }

  ## Extract tax table
  txtbl <- as.data.frame(tax_table(physeq), stringsAsFactors = F)

  ## Remove missing values if required
  if(drop_NA == TRUE){
    nas <- is.na(txtbl[, TaxRank])
    if(any(nas)){ txtbl <- txtbl[-which(nas), ] }
  }

  ## Repace NAs with character "NA"
  if(drop_NA == FALSE){
    nas <- is.na(txtbl[, TaxRank])
    if(any(nas)){ txtbl[which(nas), TaxRank] <- "NAs" }
  }

  # Get unique levels within TaxRank
  rnk <- as.vector( unique(txtbl[, TaxRank]) )

  ## Remove missing values if required
  # if(drop_NA == TRUE){ rnk <- na.omit(rnk) }

  ## TO DO - subset by each tax rank
  # dlply(.data=txtbl, .variables=TaxRank, .fun=function(x, ps){
  #   tax_subs <- as(x, "matrix")
  #   tax_table(ps) <- tax_table(tax_subs)
  #   return(ps)
  #   }, ps = physeq)

  ## Subset by each tax rank
  res <- list()
  for(i in 1:length(rnk)){
    tax_subs <- txtbl[which(txtbl[,TaxRank] == rnk[i]), ]
    tax_subs <- as(tax_subs, "matrix")
    ps <- physeq
    tax_table(ps) <- tax_table(tax_subs)
    res[[i]] <- ps
    rm(ps, tax_subs)
  }
  names(res) <- rnk

  return(res)
}

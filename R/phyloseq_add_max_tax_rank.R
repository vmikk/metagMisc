
#' @title Add the lowest level of taxonomic classification to the taxonomy table of phyloseq-class object.
#' @param physeq A phyloseq-class object
#' @param abbreviate Logical; if TRUE, lowest taxon name will be abbreviated (useful for plotting). Default is FALSE
#' @param ... Additional arguments may be passed to \code{\link{abbreviate_taxa_names}}
#' @return Phyloseq object with modified taxonomy table.
#' @export
#'
#' @examples
#' # Load data 
#' data("GlobalPatterns")
#' 
#' # Subset data to 5% of the most abundant taxa (just for the example)
#' GP <- phyloseq_filter_top_taxa(GlobalPatterns, perc = 5)
#'
#' # Add the lowest level of taxonomic classification to the taxonomy table
#' GP <- phyloseq_add_max_tax_rank(GP, abbreviate = T)
#'
#' # Replace taxa names with abbreviated lowest taxa ranks,
#' # e.g., this names could be placed on ordination plot
#' taxa_names(GP) <- as.data.frame(tax_table(GP))$LowestTaxon
#' tail(taxa_names(GP))
#' 
phyloseq_add_max_tax_rank <- function(physeq, abbreviate = F, ...){

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
  phyloseq::tax_table(physeq) <- phyloseq::tax_table( as.matrix(sp_ranks) )

  return(physeq)
}


#' @title Remove taxa with abundance less then a certain fraction of total abundance.
#'
#' @param physeq Phyloseq object
#' @param frac The minimum cutoff for the OTU abundance in the table. This number is a fraction, not a percent.
#' @details
#' If frac = 0.0001, this will retain all OTU's that have at least a 0.01% total abundance in the OTU table.
#' If you wanted to retain OTUs with at least 1% total abundance, you must specify, 0.01.
#'
#' @return
#' @export
#' @seealso http://qiime.org/scripts/filter_otus_from_otu_table.html
#' @examples
#' data("esophagus")
#' phyloseq_filter_taxa_tot_fraction(esophagus, frac = 0.01)
#'
phyloseq_filter_taxa_tot_fraction <- function(physeq, frac = 0.01){
  require(phyloseq)

  # Estimate total abundance of OTUs
  tot <- sum(taxa_sums(physeq))

  # Remove OTUs
  res <- filter_taxa(physeq, function(x){ ( sum(x)/tot ) > frac }, prune = TRUE)
  return(res)
}

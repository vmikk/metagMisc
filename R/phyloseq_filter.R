

#' @title Remove taxa with small mean relative abundance.
#'
#' @param physeq Phyloseq object
#' @param frac The minimum cutoff for the relative OTU abundance
#' @details This function searches for taxa with small mean relative abundance and removes them. Result will be returned with original counts in the abundance table.
#' @return Phyloseq object with a subset of taxa.
#' @export
#'
#' @examples
#' data("esophagus")
#' phyloseq_filter_taxa_rel_abund(esophagus, frac = 0.01)
#' phyloseq_filter_taxa_rel_abund(esophagus, frac = 0.1)
#'
phyloseq_filter_taxa_rel_abund <- function(physeq, frac = 1e-4){
  require(phyloseq)

  ## Transform OTU counts to relative abundance
  rel <- transform_sample_counts(physeq, function(x) x / sum(x) )

  ## Filter OTUs
  rel.subs <- filter_taxa(rel, function(x){ mean(x) > frac }, prune = FALSE)

  ## if prune = TRUE
  # tn <- taxa_names(rel.subs)              # OTUs to preserve
  # tr <- setdiff(taxa_names(physeq), tn)   # OTUs to remove

  ## Taxa to remove
  tr <- names(rel.subs)[ which(rel.subs == FALSE) ]

  ## If all taxa should be removed
  if(length(tr) == ntaxa(physeq)){
    stop("Error: all taxa will be removed with the specified 'frac' cutoff.\n")
  }

  ## If there is nothing to remove
  if(length(tr) == 0){
    res <- physeq
    cat("Warning: no taxa removed.\n")
  }

  ## Keep taxa which satisfies the truncation threshold
  if(length(tr) > 0){
    res <- prune_taxa(taxa = rel.subs, physeq)
  }

  return(res)
}




#' @title Remove taxa with abundance less then a certain fraction of total abundance.
#'
#' @param physeq Phyloseq object
#' @param frac The minimum cutoff for the OTU abundance in the table. This number is a fraction, not a percent.
#' @details
#' If frac = 0.0001, this will retain all OTU's that have at least a 0.01% total abundance in the OTU table.
#' If you wanted to retain OTUs with at least 1% total abundance, you must specify, 0.01.
#'
#' @return Phyloseq object with a subset of taxa.
#' @export
#' @seealso http://qiime.org/scripts/filter_otus_from_otu_table.html
#' @examples
#' data("esophagus")
#' phyloseq_filter_taxa_tot_fraction(esophagus, frac = 0.01)
#'
phyloseq_filter_taxa_tot_fraction <- function(physeq, frac = 0.01){
  require(phyloseq)

  ## Estimate total abundance of OTUs
  tot <- sum(taxa_sums(physeq))

  ## Remove OTUs
  res <- filter_taxa(physeq, function(x){ ( sum(x)/tot ) > frac }, prune = TRUE)
  return(res)
}

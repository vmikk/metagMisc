

#' @title Remove taxa with small mean relative abundance.
#'
#' @param physeq A phyloseq-class object
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
#' @param physeq A phyloseq-class object
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




#' @title Filter low-prevalence OTUs.
#' @description This function will remove taxa (OTUs) with low prevalence, where prevalence is the fraction of total samples in which an OTU is observed.
#' @param physeq A phyloseq-class object
#' @param prev.trh Prevalence threshold (default, 0.05 = 5\% of samples)
#' @param abund.trh Abundance threshold (default, NULL)
#' @details Abundance threshold defines if the OTU should be preserved if its abundance is larger than threshold (e.g., >= 50 reads).
#' @return  Phyloseq object with a subset of taxa.
#' @seealso \code{\link{phyloseq_prevalence_plot}}
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' GlobalPatterns  # 19216 taxa
#'
#' # OTUs that are found in at least 5% of samples
#' phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.05, abund.trh = NULL)  # 15389 taxa
#'
#' # The same, but if total OTU abundance is >= 10 reads it'll be preserved too
#' phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.05, abund.trh = 10)    # 15611 taxa
#'
phyloseq_filter_prevalence <- function(physeq, prev.trh = 0.05, abund.trh = NULL){

  ## Check for the low-prevalence species (compute the total and average prevalences of the features in each phylum)
  prevdf_smr <- function(prevdf){
    ddply(prevdf, "Phylum", function(df1){ data.frame(Average = mean(df1$Prevalence), Total = sum(df1$Prevalence))})
  }
  # prevdf_smr( prevalence(physeq) )

  ## Check the prevalence threshold
  # phyloseq_prevalence_plot(prevdf, physeq)

  ## Define prevalence threshold as % of total samples
  ## This function is located in 'phyloseq_prevalence_plot.R' file
  prevalenceThreshold <- prev.trh * nsamples(physeq)

  ## Calculate prevalence (number of samples with OTU) and OTU total abundance
  prevdf <- prevalence(physeq)

  ## Which taxa to preserve
  if(is.null(abund.trh)) { tt <- prevdf$Prevalence >= prevalenceThreshold }
  if(!is.null(abund.trh)){ tt <- (prevdf$Prevalence >= prevalenceThreshold | prevdf$TotalAbundance > abund.trh) }
  keepTaxa <- rownames(prevdf)[tt]

  ## Execute prevalence filter
  res <- prune_taxa(keepTaxa, physeq)
  return(res)
}



#' @title Filter rare OTUs based on minimum abundance threshold.
#' @description This function performs sample-wise OTU abundance trimming.
#' @param physeq A phyloseq-class object
#' @param minabund Abundance threshold (default, 10)
#' @param rm_zero_OTUs Logical, remove OTUs with zero total abundance
#' @details OTUs can be considered as rare if they comprise fewer than X (e.g., 10) sequences within a sample. This function is intented to censore OTU abundance (unsing an arbitrary threshold) on a sample-wise basis.
#' @return Phyloseq object with a filtered data.
#' @export
#'
#' @examples
#' # Load data
#' data(GlobalPatterns)
#'
#' # Trim GlobalPatterns data (19216 taxa)
#' GP <- phyloseq_filter_sample_wise_abund_trim(GlobalPatterns) # 10605 taxa
#'
#' # Compare raw and trimmed data
#' phyloseq_compare(GlobalPatterns, GP, cols = c("GlobalPatterns", "Trimmed GlobalPatterns"))
#'
phyloseq_filter_sample_wise_abund_trim <- function(physeq, minabund = 10, rm_zero_OTUs = TRUE){

  ## Censore OTU abundance
  res <- transform_sample_counts(physeq, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })

  ## Remove zero-OTUs
  if(rm_zero_OTUs == TRUE){
    res <- prune_taxa(taxa_sums(res) > 0, res)
  }
  return(res)
}

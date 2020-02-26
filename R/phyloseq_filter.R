
#' @title Remove samples from phyloseq object that have less than n taxa
#'
#' @param physeq A phyloseq-class object
#' @param mintaxa Minimum number of taxa that should be present in a sample (default, 10)
#'
#' @return Trimmed phyloseq object (All samples will have >= N taxa)
#' @export
#'
#' @examples
#' data("esophagus")
#' esophagus
#' phyloseq_richness_filter(esophagus, mintaxa = 30)
#' phyloseq_richness_filter(esophagus, mintaxa = 100)
#'
phyloseq_richness_filter <- function(physeq, mintaxa = 10){

  ## Estimate number of OTUs per sample
  sp <- phyloseq::estimate_richness(physeq, measures = "Observed")
  samples_to_keep <- rownames(sp)[ which(sp$Observed >= mintaxa) ]

  if(length(samples_to_keep) == 0){
    stop("All samples will be removed.\n")
  }

  if(length(samples_to_keep) == phyloseq::nsamples(physeq)){
    cat("All samples will be preserved\n")
    res <- physeq
  }

  if(length(samples_to_keep) < phyloseq::nsamples(physeq)){
    res <- phyloseq::prune_samples(samples = samples_to_keep, x = physeq)
  }

  return(res)
}



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

  # require(phyloseq)

  ## Transform OTU counts to relative abundance
  rel <- phyloseq::transform_sample_counts(physeq, function(x) x / sum(x) )

  ## Filter OTUs
  rel.subs <- phyloseq::filter_taxa(rel, function(x){ mean(x) > frac }, prune = FALSE)

  ## if prune = TRUE
  # tn <- taxa_names(rel.subs)              # OTUs to preserve
  # tr <- setdiff(taxa_names(physeq), tn)   # OTUs to remove

  ## Taxa to remove
  tr <- names(rel.subs)[ which(rel.subs == FALSE) ]

  ## If all taxa should be removed
  if(length(tr) == phyloseq::ntaxa(physeq)){
    stop("Error: all taxa will be removed with the specified 'frac' cutoff.\n")
  }

  ## If there is nothing to remove
  if(length(tr) == 0){
    res <- physeq
    cat("Warning: no taxa removed.\n")
  }

  ## Keep taxa which satisfies the truncation threshold
  if(length(tr) > 0){
    res <- phyloseq::prune_taxa(taxa = rel.subs, physeq)
  }

  return(res)
}




#' @title Remove taxa with abundance less then a certain fraction of total abundance.
#'
#' @param physeq A phyloseq-class object
#' @param frac The minimum cutoff for the OTU abundance in the table. This number is a fraction, not a percent.
#' @details
#' If frac = 0.0001, this will retain all OTU's that have at least a 0.01\% total abundance in the OTU table.
#' If you wanted to retain OTUs with at least 1\% total abundance, you must specify, 0.01.
#'
#' @return Phyloseq object with a subset of taxa.
#' @export
#' @seealso http://qiime.org/scripts/filter_otus_from_otu_table.html
#' @examples
#' data("esophagus")
#' phyloseq_filter_taxa_tot_fraction(esophagus, frac = 0.01)
#'
phyloseq_filter_taxa_tot_fraction <- function(physeq, frac = 0.01){

  # require(phyloseq)

  ## Estimate total abundance of OTUs
  tot <- sum(phyloseq::taxa_sums(physeq))

  ## Remove OTUs
  res <- phyloseq::filter_taxa(physeq, function(x){ ( sum(x)/tot ) > frac }, prune = TRUE)
  return(res)
}




#' @title Filter low-prevalence OTUs.
#' @description This function will remove taxa (OTUs) with low prevalence, where prevalence is the fraction of total samples in which an OTU is observed.
#' @param physeq A phyloseq-class object
#' @param prev.trh Prevalence threshold (default, 0.05 = 5\% of samples)
#' @param abund.trh Abundance threshold (default, NULL)
#' @param threshold_condition Indicates type of prevalence and abundance conditions, can be "OR" (default) or "AND"
#' @param abund.type Character string indicating which type of OTU abundance to take into account for filtering ("total", "mean", or "median")
#' @details
#' Abundance threshold defines if the OTU should be preserved if its abundance is larger than threshold (e.g., >= 50 reads).
#' Parameter "threshold_condition" indicates whether OTU should be kept if it occurs in many samples AND/OR it has high abundance.
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
#' phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")  # 15639 taxa
#'
#' # Include only taxa with more than 10 reads (on average) in at least 10% samples
#' phyloseq_filter_prevalence(GlobalPatterns, prev.trh = 0.1, abund.trh = 10, abund.type = "mean", threshold_condition = "AND")  # 4250 taxa
#'
phyloseq_filter_prevalence <- function(physeq, prev.trh = 0.05, abund.trh = NULL, threshold_condition = "OR", abund.type = "total"){

  ## Threshold validation
  if(prev.trh > 1 | prev.trh < 0){ stop("Prevalence threshold should be non-negative value in the range of [0, 1].\n") }
  if(abund.trh <= 0){ stop("Abundance threshold should be non-negative value larger 0.\n") }

  ## Check for the low-prevalence species (compute the total and average prevalences of the features in each phylum)
  prevdf_smr <- function(prevdf){
    ddply(prevdf, "Phylum", function(df1){ data.frame(Average = mean(df1$Prevalence), Total = sum(df1$Prevalence))})
  }
  # prevdf_smr( prevalence(physeq) )

  ## Check the prevalence threshold
  # phyloseq_prevalence_plot(prevdf, physeq)

  ## Define prevalence threshold as % of total samples
  ## This function is located in 'phyloseq_prevalence_plot.R' file
  prevalenceThreshold <- prev.trh * phyloseq::nsamples(physeq)

  ## Calculate prevalence (number of samples with OTU) and OTU total abundance
  prevdf <- prevalence(physeq)

  ## Get the abundance type
  if(abund.type == "total") { prevdf$AbundFilt <- prevdf$TotalAbundance }
  if(abund.type == "mean")  { prevdf$AbundFilt <- prevdf$MeanAbundance }
  if(abund.type == "median"){ prevdf$AbundFilt <- prevdf$MedianAbundance }

  ## Which taxa to preserve
  if(is.null(abund.trh)) { tt <- prevdf$Prevalence >= prevalenceThreshold }
  if(!is.null(abund.trh)){
    ## Keep OTU if it either occurs in many samples OR it has high abundance
    if(threshold_condition == "OR"){
      tt <- (prevdf$Prevalence >= prevalenceThreshold | prevdf$AbundFilt >= abund.trh)
    }

    ## Keep OTU if it occurs in many samples AND it has high abundance
    if(threshold_condition == "AND"){
      tt <- (prevdf$Prevalence >= prevalenceThreshold & prevdf$AbundFilt >= abund.trh)
    }
  }

  ## Extract names for the taxa we whant to keep
  keepTaxa <- rownames(prevdf)[tt]

  ## Execute prevalence filter
  res <- phyloseq::prune_taxa(keepTaxa, physeq)
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
phyloseq_filter_sample_wise_abund_trim <- function(physeq, minabund = 10, relabund = FALSE, rm_zero_OTUs = TRUE){

  ## Censore OTU abundance
  if(relabund == FALSE){     # trim based on absolute OTU counts

    res <- phyloseq::transform_sample_counts(physeq, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })

  } else {                   # trim based on relative abundances within sample, but return original counts

    if(!minabund > 0 & minabund <= 1){
      stop("Error: for relative abundance trimmin 'minabund' should be in (0,1] interval.\n")
    }

    ## Convert data to relative abundances
    res <- phyloseq_standardize_otu_abundance(physeq, method = "total")

    ## Remove relative abundances less than the threshold value
    res <- phyloseq::transform_sample_counts(res, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })

    ## Sample sums and data orientation
    smps <- phyloseq::sample_sums(physeq)
    if(phyloseq::taxa_are_rows(physeq) == TRUE){
      mar <- 2
    } else {
      mar <- 1
    }

    ## Convert back to counts by multiplying relative abundances by sample sums
    phyloseq::otu_table(res) <- phyloseq::otu_table(
      sweep(x = phyloseq::otu_table(res), MARGIN = mar, STATS = smps, FUN = `*`),
      taxa_are_rows = phyloseq::taxa_are_rows(physeq))
  }

  ## Remove zero-OTUs
  if(rm_zero_OTUs == TRUE){
    res <- phyloseq::prune_taxa(taxa_sums(res) > 0, res)
  }
  return(res)
}



#' @title Extract the most abundant taxa.
#' @param physeq A phyloseq-class object
#' @param perc Percentage of the most abundant taxa to retain
#' @param n Number of the most abundant taxa to retain (this argument will override perc argument)
#' @return Phyloseq object with a filtered data.
#' @export
#'
#' @examples
#'
phyloseq_filter_top_taxa <- function(physeq, perc = 10, n = NULL){

  ## Arguments validation
  if(perc <= 0 | perc > 100){ stop("Error: percentage should be in 1-100 range.\n") }

  ## Get total abundances for all taxa
  taxx <- sort(phyloseq::taxa_sums(physeq), decreasing = TRUE)

  ## Find how many taxa to preserve (if percentage is specified)
  if(is.null(n)){
    n <- phyloseq::ntaxa(physeq) * perc / 100
    n <- floor(n)
  }

  ## Extract names for the taxa that should be preserved
  keepTaxa <- names(taxx)[1:n]

  ## Extract this taxa
  physeq_pruned <- phyloseq::prune_taxa(keepTaxa, physeq)

  return(physeq_pruned)
}



#' @title Check the range of the top-taxa filtering values to determine the optimal threshold.
#' @description This function performs taxa filtering by retaining the most abundant taxa.
#' A range of abundance percentages (5 - 95\%) will be explored.
#' @param physeq A phyloseq-class object
#' @param show_plot Logical; if TRUE, shows the plot on screen
#' @return ggplot-object.
#' @export
#'
#' @examples
#'
phyloseq_filter_top_taxa_range <- function(physeq, show_plot = TRUE){
  percs <- seq(5, 95, 5)

  fr <- plyr::mlply(.data = data.frame(perc = percs), .fun = function(...){ phyloseq_filter_top_taxa(physeq, ...) })
  names(fr) <- percs

  fr_tab <- plyr::ldply(.data = fr, .fun = function(z){
    sz <- phyloseq::sample_sums(z)
    res <- data.frame(Sample = names(sz), Preserved = sz)
    return(res)
  })

  pp <- ggplot(data = fr_tab, aes(x = perc, y = Preserved, group = Sample)) +   # color = Sample
    geom_vline(xintercept=75, color="grey", linetype = "longdash") +
    geom_line() +
    geom_point() +
    labs(x = "Number of most abundant taxa retained, %", y = "Percentage of total sample abundance") +
    theme(legend.position = "none")

  if(show_plot == TRUE){ print(pp) }
  invisible(pp)
}

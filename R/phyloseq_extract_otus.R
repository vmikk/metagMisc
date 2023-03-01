
#' @title Extract common species (OTUs) between samples.
#' @description This function will subset phyloseq object to the OTUs that are present in all samples. By default all samples will be taken into account, otherwise it's possible to take a subset of samples.
#'
#' @param x A phyloseq-class object
#' @param samp_names Character vector with sample names (default, all samples)
#'
#' @return Phyloseq object with the subset of data.
#' @export
#'
#' @examples
#' data(esophagus)
#'
#' # Extract OTUs that are present in all samples
#' ps <- phyloseq_extract_shared_otus(esophagus, samp_names = sample_names(esophagus))
#' ps
#' otu_table(ps)
#'
#' # Extract shared OTUs between two samples
#' ps2 <- phyloseq_extract_shared_otus(esophagus, samp_names = c("B", "C"))
#' ps2
#' otu_table(ps2)
#'
phyloseq_extract_shared_otus <- function(x, samp_names = sample_names(x)){

  require(phyloseq)

  ## Test if the sample names are valid
  if( any(!samp_names %in% sample_names(x)) ){
    stop("Check the sample names, not all of them are present in the phyloseq object.\n")
  }

  ## Extract samples
  if(length(unique(samp_names)) != length(unique(sample_names(x)))){
    x <- prune_samples(samples = samp_names, x = x)
  }

  ## Count number of OTUs (non-zero)
  n_tot_otu <- ntaxa( prune_taxa(taxa_sums(x) > 0, x) )

  ## Subset to OTUs that are present in both samples
  x <- try( filter_taxa(x, function(z){ sum(z >= 1) == length(samp_names) }, TRUE) )

  ## Count number of shared OTUs
  if(! "try-error" %in% class(x)){
    n_shared_otu <- ntaxa(x)
  } else {
    cat("WARNING: no shared OTUs found!\n")
    n_shared_otu <- 0
  }
  
  shared_ratio <- n_shared_otu * 100 / n_tot_otu

  ## Add attributes to the results
  attr(x, which = "TotalNumberOfOTUs")  <- n_tot_otu
  attr(x, which = "NumberOfSharedOTUs") <- n_shared_otu
  attr(x, which = "SharedOTURatio") <- shared_ratio

  return(x)
}



#' @title Extract non-shared species (OTUs) between samples.
#' @description This function will subset phyloseq object to the OTUs unique for each sample. By default all samples will be taken into account, otherwise it's possible to take a subset of samples.
#'
#' @param x A phyloseq-class object
#' @param samp_names Character vector with sample names (default, all samples)
#'
#' @return Phyloseq object with the subset of data.
#' @export
#'
#' @examples
#' data(esophagus)
#'
#' # Extract unique (non-shared) OTUs for each sample in the entire dataset
#' ps <- phyloseq_extract_non_shared_otus(esophagus, samp_names = sample_names(esophagus))
#' ps
#' otu_table(ps)
#'
#' # Extract OTUs that are unique between two samples
#' ps2 <- phyloseq_extract_non_shared_otus(esophagus, samp_names = c("B", "C"))
#' ps2
#' otu_table(ps2)
#'
phyloseq_extract_non_shared_otus <- function(x, samp_names = sample_names(x)){

  require(phyloseq)

  # test if the sample names are valid
  if( any(!samp_names %in% sample_names(x)) ){
    stop("Check the sample names, not all of them are present in the phyloseq object.\n")
  }

  # extract samples
  xx <- prune_samples(samples = samp_names, x = x)

  # subset to OTUs that are present only in 1 sample
  xx <- try( filter_taxa(xx, function(z){ sum(z >= 1) == 1 }, TRUE) )

  if("try-error" %in% class(xx)){
    cat("WARNING: no shared OTUs found!\n")
  }

  return(xx)
}


#' @title Determine which OTUs appeared or disappeared in comparison with reference sample group.
#'
#' @param phys A phyloseq-class object
#' @param ref_level Sample name(s) relative to which appearance or disappearance of OTU will be measured
#'
#' @return List with OTU names that appeared, disappeared or remained.
#' @export
#'
#' @examples
#'
phyloseq_otu_appearance <- function(phys, ref_level){

  ## Remove tax & phy slots
  if(!is.null(phyloseq::phy_tree(phys, errorIfNULL = F))) { phys@phy_tree <- NULL }
  if(!is.null(phyloseq::tax_table(phys, errorIfNULL = F))){ phys@tax_table <- NULL }

  ## Split phyloseq to reference sample(s) and all other
  pref <- phyloseq::prune_samples(phyloseq::sample_names(phys) %in% ref_level, phys)
  poth <- phyloseq::prune_samples(!phyloseq::sample_names(phys) %in% ref_level, phys)

  ## Drop missing OTUs
  pref <- phyloseq::prune_taxa(taxa_sums(pref) > 0, pref)
  poth <- phyloseq::prune_taxa(taxa_sums(poth) > 0, poth)

  ## Extract OTU names
  nref <- phyloseq::taxa_names(pref)
  noth <- phyloseq::taxa_names(poth)

  res <- list()
  res$disappeared <- nref[ !nref %in% noth ]
  res$remained <- intersect(nref, noth)
  res$appeared <- noth[ !noth %in% nref ]

  return(res)
}

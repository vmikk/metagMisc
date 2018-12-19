
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

  # test if the sample names are valid
  if( any(!samp_names %in% sample_names(x)) ){
    stop("Check the sample names, not all of them are present in the phyloseq object.\n")
  }

  # extract samples
  xx <- prune_samples(samples = samp_names, x = x)

  # subset to OTUs that are present in both samples
  xx <- filter_taxa(xx, function(z){ sum(z >= 1) == length(samp_names) }, TRUE)

  return(xx)
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
  xx <- filter_taxa(xx, function(z){ sum(z >= 1) == 1 }, TRUE)

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
  if(!is.null(phy_tree(phys, errorIfNULL = F))) { phys@phy_tree <- NULL }
  if(!is.null(tax_table(phys, errorIfNULL = F))){ phys@tax_table <- NULL }

  ## Split phyloseq to reference sample(s) and all other
  pref <- prune_samples(sample_names(phys) %in% ref_level, phys)
  poth <- prune_samples(!sample_names(phys) %in% ref_level, phys)

  ## Drop missing OTUs
  pref <- prune_taxa(taxa_sums(pref) > 0, pref)
  poth <- prune_taxa(taxa_sums(poth) > 0, poth)

  ## Extract OTU names
  nref <- taxa_names(pref)
  noth <- taxa_names(poth)

  res <- list()
  res$disappeared <- nref[ !nref %in% noth ]
  res$remained <- intersect(nref, noth)
  res$appeared <- noth[ !noth %in% nref ]

  return(res)
}

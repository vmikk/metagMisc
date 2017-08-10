
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

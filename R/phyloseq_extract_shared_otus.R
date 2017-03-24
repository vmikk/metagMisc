
## Extract common species (OTUs) between samples
# by default all samples will be taken
phyloseq_extract_shared_otus <- function(x, samp_names = sample_names(x)){
  # x = phyloseq object
  # samp_names = character vector with sample names

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

## Example:
data(esophagus)
phyloseq_extract_shared_otus(esophagus, samp_names = sample_names(esophagus))
phyloseq_extract_shared_otus(esophagus, samp_names = c("B", "C"))

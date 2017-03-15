

# Multiple rarefaction (Resample an OTU table such that all samples have the same library size)
phyloseq_mult_raref <- function(x, SampSize = NULL, MinSizeTreshold = NULL, iter = 1000, ...){
  # x = physeq object
  # SampSize = rarefaction depth (number of reads to sample)
  # MinSizeTreshold = remove samples with number of reads less then this treshold
  # iter = number of rarefication iterations

  require(plyr)
  require(phyloseq)

  # Filter samples by number of reads
  if(!is.null(MinSizeTreshold)){ x <- prune_samples(sample_sums(x) >= MinSizeTreshold, x) }

  # Define rarefication depth
  if(is.null(SampSize)){ SampSize <- 0.9*min(sample_sums(x)) }

  # rarefy without replacement
  res <- mlply(
    .data = 1:iter,
    .fun = function(z){ rarefy_even_depth(x, rngseed=z, sample.size=SampSize, replace=F, verbose = F) },
    .progress = "text")

  return(res)
}


#' @title Multiple rarefaction of phyloseq-object
#' @description This function will resample an OTU table such that all samples have the same library size. Resampling will be performed multiple times.
#' @param x A phyloseq-class object
#' @param SampSize Rarefaction depth (number of reads to sample)
#' @param MinSizeTreshold Remove samples with number of reads less then this treshold
#' @param iter Number of rarefication iterations
#' @param replace Logical, whether to sample with replacement (TRUE) or without replacement (FALSE, default)
#' @param seeds Integer vector used for the reproducible random subsampling (should be of the same length as the number of iterations)
#' @param ... Additional arguments will be passed to \code{\link{rarefy_even_depth}}
#'
#' @details
#' If the sample size ('SampSize') is not specified, rarefaction will be made for the depth equeal to 0.9 * minimal observed sample size.
#' By default, sampling is performed without replacement ('replace = FALSE'), which differs from the default behaviour of \code{\link{rarefy_even_depth}}.
#'
#' @return List of rarefied phyloseq-objects.
#' @seealso \code{\link{rarefy_even_depth}}
#' @export
#'
#' @examples
#' data("esophagus")
#' eso_raref <- phyloseq_mult_raref(esophagus, iter = 10)
#'
#' # Discard samples with number of reads < 210 and perform resampling with replacement
#' eso_raref_t <- phyloseq_mult_raref(esophagus, replace = T, MinSizeTreshold = 210, SampSize = 210, iter = 10)
#' sample_sums(esophagus)
#' sample_sums(eso_raref_t[[1]])
#'
#' # Do not remove OTUs from the dataset that are no longer observed in any sample (have a count of zero in every sample)
#' phyloseq_mult_raref(esophagus, trimOTUs = F, replace = T, MinSizeTreshold = 210, SampSize = 210, iter = 10)
#'
phyloseq_mult_raref <- function(x, SampSize = NULL, MinSizeTreshold = NULL, iter = 1000, replace = F, seeds = NULL, ...){

  require(plyr)
  require(phyloseq)

  ## Sanity check for the random number generator
  if(!is.null(seeds)){
    if(length(seeds) != iter){ stop("Error: lenght of 'seeds' should be the same as the number of iterations.\n") }
    if(length(seeds) != length(unique(seeds))){ warning("Warning: Provided seeds are not unique which leads to the identical results of random sampling.\n") }
    if(!isTRUE(all(seeds == floor(seeds)))){ stop("Error: Seeds must only contain integer values.\n") }
  }

  ## Filter samples by number of reads
  if(!is.null(MinSizeTreshold)){ x <- prune_samples(sample_sums(x) >= MinSizeTreshold, x) }

  ## Define rarefication depth
  if(is.null(SampSize)){ SampSize <- 0.9*min(sample_sums(x)) }

  ## Prepare seed values
  if(is.null(seeds)){ seeds <- 1:iter }

  ## Rarefy
  res <- mlply(
    .data = seeds,
    .fun = function(z){ rarefy_even_depth(x, rngseed=z, sample.size=SampSize, replace=replace, verbose = F, ...) },
    .progress = "text")

  return(res)
}

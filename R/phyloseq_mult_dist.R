## Function to perform multiple rarefactions and average sample dissimilarity across rarefactions
phyloseq_mult_dist <- function(physeq, dissimilarity = "bray", meanfun = mean,
    SampSize = min(sample_sums(physeq)), iter = 1000, parallel = FALSE, verbose = TRUE, ...){
  # meanfun = The calculation to use for the average (mean or median)

  ## TO DO: - add data validation (e.g. unifrac & phy_tree)
  ## TO DO: - add ... for `phyloseq::distance` (now ... are passed to `phyloseq_mult_raref`)
  ## TO DO: - add averaging type - mean/median:
  ##   meanfun <- match.fun(meanfun)  # Which function to use for averaging
  ## TO DO: - add taxa subsets (e.g. for a particular phylum)

  ## Multiple rarefaction
  if(verbose == TRUE){ cat("..Multiple rarefaction\n") }
  phys_raref <- phyloseq_mult_raref(physeq, SampSize = SampSize, iter = iter, multithread = parallel, ...)


  ## Estimate dissimilarities
  if(verbose == TRUE){ cat("..Dissimilarity estimation\n") }

  if(length(dissimilarity) == 1){  # single dissimilarity coefficient
    res <- mult_dissim(phys_raref, method = dissimilarity, average = T)
  } else {                         # multiple dissimilarity coefficients (e.g., dissimilarity = c("bray", "unifrac"))
    res <- plyr::mlply(
        .data = data.frame(method = dissimilarity, stringsAsFactors = F),
        .fun = function(...){ mult_dissim(phys_raref, average = T, ...) })
    names(res) <- dissimilarity
  }

  ## Add rarefaction attributes to the results
  attributes(res) <- c(attributes(res), attributes(phys_raref)[c("RarefactionDepth", "RarefactionReplacement")])
  return(res)
}
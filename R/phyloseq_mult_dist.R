
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


#' @title Average multiple distance matrices.
#' @description This function can be used in order to avergage beta diversity over multiple rarefaction iterations.
#' @param dlist List of distance matrices (class 'dist' or 'matrix')
#'
#' @return Distance matrix of the same class as input matrices ('dist' or 'matrix') with averaged values.
#' @export
#'
#' @examples
#' library(plyr)
#' # Generate dummy data (list with 3 5x5 matrices)
#' ddd <- rlply(.n = 3, .expr = function(){ as.dist( matrix(sample(1:1000, size = 25, replace = T), ncol = 5)) })
#' ddd
#'
#' # Average matrices
#' mult_dist_average(ddd)
#'
mult_dist_average <- function(dlist){
  # dlist = list of `dist` or `matrix` objects

  ## Input = list of matrices
  if(class(dlist[[1]]) %in% "matrix"){
    res <- Reduce("+", dlist) / length(dlist)
  }

  ## Input = list of dist
  if(class(dlist[[1]]) %in% "dist"){
    res <- Reduce("+", plyr::llply(.data = dlist, .fun = as.matrix)) / length(dlist)
    res <- as.dist(res)
  }

  return(res)
}


#' @title Compute beta diversity for each rarefaction iteration.
#'
#' @param x List of phyloseq objects (result of \code{\link{phyloseq_mult_raref}})
#' @param method A character string with the name of supported dissimilarity index (see \code{\link{distanceMethodList}})
#' @param average Logical; if TRUE, dissimilarity averaged over rarefication iterations will be returned; if FALSE, list of dissimilarity matrices will be returned.
#' @param ... Additional arguments will be passed to \code{\link[phyloseq]{distance}}
#'
#' @return List of 'dist'-matrices (if average = FALSE) or a single 'dist' (if average = TRUE).
#' @export
#'
#' @examples
#' # Load data
#' data(esophagus)
#' sample_sums(esophagus)  # samples has different number of sequences
#'
#' # Perform multiple rarefaction (sample 200 sequences from each sample, repeat the procedure 100 times)
#' esor <- phyloseq_mult_raref(esophagus, SampSize = 200, iter = 100)
#' sample_sums(esor[[1]])  # rarefied data
#'
#' # Estimate sample dissimilarity independently for each iteration
#' eso_dis <- mult_dissim(esor, method = "unifrac", average = F)
#' eso_dis[[1]]   # unweighted UniFrac distances for the first rarefaction iteration
#'
#' # Average sample dissimilarities over all rarefaction iterations
#' eso_dis_avg <- mult_dissim(esor, method = "unifrac", average = T)
#' eso_dis_avg    # mean unweighted UniFrac distances
#'
mult_dissim <- function(x, method = "bray", average = T, ...){

  # require(phyloseq)
  # require(plyr)

  physeq_dissim <- plyr::llply(
    .data = x,
    .fun = function(z, ...){ phyloseq::distance(physeq = z, type = "samples", ...) },
    method = method,
    .progress = "text")

  if(average == TRUE){
    physeq_dissim_avg <- mult_dist_average(dlist = physeq_dissim)
    return(physeq_dissim_avg)
  } else {
    return(physeq_dissim)
  }
}


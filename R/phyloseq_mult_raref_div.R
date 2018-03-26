
#' @title Average diversity estimates across multiple rarefaction iterations.
#' @description This function performs multiple rarefaction and estimates average diversity (e.g., Simpson or Shannon indices) for each sample.
#' @param physeq A phyloseq-class object
#' @param SampSize Rarefaction depth (number of reads to sample)
#' @param iter Number of rarefication iterations
#' @param divindex Alpha-diversity measures to estimate (for the supported indices see \code{\link[phyloseq]{estimate_richness}})
#' @param parallel Logical; if TRUE, attempts to run the function on multiple cores
#' @param verbose Logical; if TRUE, progress messages from the function will be printed
#' @param ... Additional argument may be passed to \code{\link{phyloseq_mult_raref}}
#'
#' @return A data.frame of the richness estimates.
#' @export
#'
#' @examples
phyloseq_mult_raref_div <- function(physeq, SampSize = min(sample_sums(physeq)), iter = 1000,
  divindex = c("Observed", "Shannon"), parallel = FALSE, verbose = TRUE, ...){

  ## Progress indicator
  if(verbose == TRUE){
    progr <- "text"
  } else {
    progr <- "none"
  }

  ## Multiple rarefaction
  if(verbose == TRUE){ cat("..Multiple rarefaction\n") }
  phys_raref <- phyloseq_mult_raref(physeq, SampSize = SampSize, iter = iter, multithread = parallel, ...)

  ## Estimate richness and diversity
  phys_diversity <- function(x, indices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")){
    res <- phyloseq::estimate_richness(x, split = T, measures = indices)
    res <- dfRowName(res, name = "Sample")
    return(res)
  }

  if(verbose == TRUE){ cat("..Diversity esimation\n") }
  DIV <- plyr::ldply(
    .data = phys_raref,
    .fun = phys_diversity,
    .id = "Iteration",
    indices = divindex,
    .progress = progr)

  ## bug in ldply? column name is "X1" instead of "Iteration"
  colnames(DIV)[1] <- "Iteration"

  ## Reshape data to long-format
  DIV.resh <- reshape2::melt(
    data = DIV,
    id.vars = c("Iteration", "Sample"),
    variable.name = "Index",
    value.name = "Value")

  ## Average diversity estimates over all iterations
  if(verbose == TRUE){ cat("..Averaging diversity estimates\n") }
  DIV.means <- plyr::ddply(
    .data = DIV.resh,
    .variables = c("Sample", "Index"),
    .fun = function(x) { ggplot2::mean_cl_boot(x$Value) })  # mean_cl_boot = ggplot2-wrapper for Hmisc

  colnames(DIV.means) <- c("Sample", "Index", "Estimate", "CI.lower", "CI.upper")

  return(DIV.means)
}

#' @title Inspect random samples (and taxa) in a phyloseq object
#'
#' @description Inspect a random subset of samples from a
#'   \code{\link[phyloseq]{phyloseq-class}} object. Optionally, also
#'   subsample taxa (OTUs/ASVs) for compact viewing.
#'
#' @param x A \code{\link[phyloseq]{phyloseq-class}} object
#' @param n Number of samples to return
#' @param n_otu Optional number of taxa (OTUs/ASVs) to retain; set to \code{NULL}
#'   to keep all taxa
#' @return A phyloseq object pruned to \code{n} samples and up to \code{n_otu} taxa.
#' @seealso \code{\link[phyloseq]{prune_samples}}, \code{\link[phyloseq]{prune_taxa}},
#'   \code{\link[phyloseq]{nsamples}}, \code{\link[phyloseq]{ntaxa}}
#' @export
#' @author Inspired by \code{\link[car]{some}} in the car-package by John Fox.
#' @examples
#' library(phyloseq)
#' data(enterotype)
#' 
#' phyloseq_some(enterotype)
#' 
#'
phyloseq_some <- function(x, n = 10, n_otu = 10){

  ## Input data = data.frame or matrix
  if(any(class(x) %in% c("data.frame", "matrix"))){
    if(nrow(x) > n){
      show <- sample(x = 1:nrow(x), size = n)
      res <- x[show, ]
    } else {
      warning("Cannot take a sample larger than the population.\n")
      res <- x
    }
  }

  ## Input data = phyloseq object
  if("phyloseq" %in% class(x)){

    ## How many samples?    
    nsamps <- phyloseq::nsamples(x)

    ## Subset samples
    if(nsamps > n){
      samps <- sample(x = 1:nsamps, size = n)
      res <- phyloseq::prune_samples(1:nsamps %in% samps, x)
    } else {
      warning("Number of samples in phyloseq <= n \n")
      res <- x
    }
 
    ## How many taxa?
    ntax <- phyloseq::ntaxa(res)

    ## Subset OTUs
    if(!is.null(n_otu) & ntax > n_otu){
      otus <- sample(x = 1:ntax, size = n_otu)
      res <- phyloseq::prune_taxa(1:ntax %in% otus, res)
    } else {
      warning("Number of OTUs in phyloseq <= n_otu \n")
    }

  } # end of phyloseq

  return(res)
}

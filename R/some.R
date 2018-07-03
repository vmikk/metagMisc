#' @title Check random lines in the data (including phyloseq objects)
#'
#' @description Returns several random lines from the data.frame, matrix or phyloseq.
#' Usage is quite similar to the \code{head} or \code{tail} functions.
#'
#' @param x Coluld be a "data.frame", "matrix", or "phyloseq" class
#' @param n Number of records (or samples in phyloseq) that will be returned
#' @param n_otu Number of OTUs to show (if the input is phyloseq-object)
#' @return Part of the data.frame or matrix with n lines or a subset of phyloseq-object.
#' @author Adapted from \code{some} in the car-package by John Fox.
#' @examples
#' # Load phyloseq data
#' data(enterotype)
#' some(enterotype)
#' 
#' # If input is data.frame or matrix
#' x <- matrix(rnorm(500), ncol=5)
#' some(x)
#'
some <- function(x, n = 10, n_otu = 10){

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

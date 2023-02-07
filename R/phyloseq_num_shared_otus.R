
#' @title Estimate the number of shared and non-shared species (OTUs/ESvs) between all pairwise comparisons of samples.
#' @description
#' This function counts the number of OTUs that are common (shared) and
#' distinct (non-shared) between a pair of samples. All pairwise sample
#' combinations are considered, and the output is stored as a matrix.
#'
#' @param physeq A phyloseq-class object
#' @param use_Matrix Logical, use sparse matrix class to perform the analysis (default, TRUE)
#'
#' @export
#' @examples
#' data(esophagus)
#'
#' # Count number of shared and non-shared OTUs
#' ps <- phyloseq_num_shared_otus(esophagus)
#' ps
#'
#' # The same, but using base-R matrix class (not recommended for large data)
#' ps <- phyloseq_num_shared_otus(esophagus, use_Matrix = FALSE)
#' ps
#'
phyloseq_num_shared_otus <- function(physeq, use_Matrix = TRUE){

  require(phyloseq)
  if(use_Matrix == TRUE){ require(Matrix) }

  ## Extract OTU table
  x <- as( phyloseq::otu_table(physeq), "matrix" )

  ## Convert to presence-absence
  x[ x > 0 ] <- 1

  ## Convert to the sparse matrix
  if(use_Matrix == TRUE){
    # x <- as(x, "Matrix")
    x <- as(x, "sparseMatrix")
  }

  ## Estimate the number of shared OTUs
  if(taxa_are_rows(physeq) == TRUE){
    otus_shared <- t(x) %*% x
  } else {
    otus_shared <- x %*% t(x)
  }

  ## Estimate the number of non-shared OTUs
  ## Asymmetric (A-B, B-A)
  otus_nonshared_asymm <-  abs(
    sweep(otus_shared, MARGIN = 2, diag(otus_shared))   # subtract the diagonal
    )

  ## Total number of non-shared OTUs (symmetric)
  otus_nonshared <- otus_nonshared_asymm + t(otus_nonshared_asymm)

  ## Prepare results
  res <- list(
    shared               = otus_shared,
    nonshared_total      = otus_nonshared,
    nonshared_asymmetric = otus_nonshared_asymm
    )
  return(res)
}


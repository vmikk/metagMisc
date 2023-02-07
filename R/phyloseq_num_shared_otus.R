
#' @title Estimate the number of shared and non-shared species (OTUs/ESvs) between all pairwise comparisons of samples.
#' @description
#' This function counts the number of OTUs that are common (shared) and
#' distinct (non-shared) between a pair of samples. All pairwise sample
#' combinations are considered, and the output is stored as a matrix.
#'
#' @param physeq A phyloseq-class object
#' @param use_Matrix Logical, use sparse matrix class to perform the analysis (default, TRUE)
#'
#' @details
#' This function uses the sparse matrix class implemented in the package
#' \code{\link[Matrix]{Matrix}} to improve computational speed and
#' decrease memory footprint.
#'
#' @return
#' The function's output is a list containing three matrices
#' (rows and columns correspond to samples):
#'
#' \itemize{
##'  \item{\strong{"shared"}}{ - The number of shared OTUs between the samples}
##'
##'  \item{\strong{"nonshared_total"}}{ - The number of non-shared OTUs between the samples.
##'   This matrix contains the total number of non-shared OTUs between the two samples.
##'   The matrix is symmetric, meaning that the values in the matrix are equal on both sides of the diagonal.}
##'
##'  \item{\strong{"nonshared_asymmetric"}}{ - The number of non-shared OTUs between the samples.
##'    The matrix is asymmetric, with the upper and lower triangular parts of
##'    the matrix representing the number of unique OTUs present only in one of the two samples.}
##' }
#'
#' @export
#' @seealso
#' \code{\link[metagMisc]{phyloseq_extract_shared_otus}},
#' \code{\link[metagMisc]{phyloseq_extract_non_shared_otus}}
#'
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
#' # To obtain a list of shared or non-shared OTUs between a pair of samples, use
#' phyloseq_extract_shared_otus(x = esophagus, samp_names = c("B", "C"))
#' phyloseq_extract_non_shared_otus(x = esophagus, samp_names = c("B", "C"))
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


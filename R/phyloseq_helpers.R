

#' Extract OTU table as data.frame with controlled rows and columns orientation
#'
#' @description 
#' Extracts the OTU abundance table from a phyloseq object and converts it to 
#' a data.frame with control over whether taxa should be arranged as rows or columns.
#'
#' @param physeq A \code{phyloseq} object containing an OTU table (or \code{otu_table} object)
#' @param taxa_as_rows Logical or NULL. 
#'   If \code{NULL} (default), the orientation of the OTU table will not be changed. 
#'   If \code{TRUE}, taxa will be arranged as rows and samples as columns in the output data.frame. 
#'   If \code{FALSE}, samples will be rows and taxa will be columns.
#'
#' @return A data.frame containing OTU abundances with the specified orientation.
#'   Row and column names are preserved from the original phyloseq object.
#'
#' @details
#' This function handles the orientation automatically by checking the internal 
#' structure of the phyloseq object and transposing the data if necessary to match 
#' the requested output format. This is particularly useful when preparing data 
#' for functions that expect specific data orientations (e.g., vegan functions 
#' typically expect samples as rows).
#'
#' @importFrom phyloseq otu_table taxa_are_rows
#'
#' @export
#'
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' 
#' # Extract OTU table with taxa as rows (default)
#' otu_taxa_rows <- phyloseq_otu_to_df(GlobalPatterns)
#' dim(otu_taxa_rows)
#' 
#' # Extract OTU table with samples as rows (for vegan compatibility)
#' otu_samp_rows <- phyloseq_otu_to_df(GlobalPatterns, taxa_as_rows = FALSE)
#' dim(otu_samp_rows)
#'
phyloseq_otu_to_df <- function(physeq, taxa_as_rows = NULL){

  res <- as.data.frame(phyloseq::otu_table(physeq))
  current_orientation <- phyloseq::taxa_are_rows(physeq)

  ## If taxa_as_rows is NULL, use the default orientation of the OTU table
  if(is.null(taxa_as_rows)){
    taxa_as_rows <- current_orientation
  }

  ## Validation
  if(!is.logical(taxa_as_rows)){
    stop("taxa_as_rows must be TRUE or FALSE\n")
  }

  ## Transpose if needed to match desired orientation

  ## Want taxa as rows
  if(taxa_as_rows == TRUE){

    ## Currently samples as rows, need to transpose
    if(current_orientation == FALSE){
      res <- t(res)
    }
    ## If current_orientation == TRUE, no change needed

  ## Want samples as rows (taxa as columns)
  } else {

    ## Currently taxa as rows, need to transpose
    if(current_orientation == TRUE){
      res <- t(res)
    }
    ## If current_orientation == FALSE, no change needed
  }

  return(res)
}

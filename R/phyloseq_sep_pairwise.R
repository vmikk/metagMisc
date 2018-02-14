
#' @title Split phyloseq-class object for pairwise comparisons by sample-level variable.
#' @param phys A phyloseq-class object
#' @param group Variable name (contained in \code{\link{sample_data}})
#' @param comparis Character "all" for all possible pairwise comparisons (default), or a matrix with comparison names to use (similar to the output of \code{\link{combn}} with m=2)
#' @param drop_zeroes Logical, indicating weather OTUs with zero abundance or samples with zero total abundance should be removed
#' @return List with phyloseq objects (each with only two groups based on selected variable).
#' @export
#'
#' @examples
#'
phyloseq_sep_pairwise <- function(phys, group, comparis = "all", drop_zeroes = TRUE){

  ## Split by groups
  phgr <- phyloseq_sep_variable(phys, variable = group, drop_zeroes = drop_zeroes)

  ## All possible comparisons
  if(is.character(comparis) & comparis == "all"){
    compars <- combn(x = names(phgr), m = 2)
  } else if(is.matrix(comparis)){
  ## Use the provided comparisons (e.g., not all combinations, or another order of pairs)
    ## Data validation
    if(nrow(comparis) != 2){
      stop("Matrix of pairwise comparisons should have two rows and be similar to the output of 'combn(..., m=2)' function.\n")
    }
    if(any(!unique(as.vector(compars)) %in% names(phgr))){
      stop("Some sample groups provided in matrix of pairwise combinations 'comparis' are missing in phyloseq object 'phys'.\n")
    }

    compars <- comparis

  } else {
    stop("'comparis' should be a character (e.g., 'all') or a matrix similar to the output of 'combn(..., m=2)' function.\n")
  }

  ## Iterate over combinations and merge phyloseqs
  res <- list()
  for(i in 1:ncol(compars)){
    res[[i]] <- phyloseq::merge_phyloseq(
      phgr[[compars[1, i]]],
      phgr[[compars[2, i]]]
      )
  }
  names(res) <- aaply(.data = compars, .margins = 2, .fun = paste, collapse="__")

  ## Remove zero OTUs
  if(drop_zeroes == TRUE){
    res <- llply(.data = res, .fun = function(z){ prune_taxa(taxa_sums(z) > 0, z) })
  }

  return(res)
}

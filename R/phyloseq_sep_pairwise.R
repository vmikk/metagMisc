
#' @title Split phyloseq-class object for pairwise comparisons by sample-level variable.
#' @param phys A phyloseq-class object
#' @param group Variable name (contained in \code{\link[phyloseq]{sample_data}})
#' @param comparis Character "all" for all possible pairwise comparisons (default), or a matrix with comparison names to use (similar to the output of \code{\link{combn}} with m=2)
#' @param drop_zeroes Logical, indicating weather OTUs with zero abundance or samples with zero total abundance should be removed
#' @return List with phyloseq objects (each with only two groups based on selected variable).
#' @export
#'
#' @examples
#' data(enterotype)
#' 
#' # Split enterotype data for pairwise comparisons by sequencing technology
#' ent <- phyloseq_sep_pairwise(enterotype, group = "SeqTech", drop_zeroes = TRUE)
#' ent
#' 
#' # Perform differential abundance testing for each of the pairwise comparisons
#' 
phyloseq_sep_pairwise <- function(phys, group, comparis = "all", drop_zeroes = TRUE){

  ## Split by groups
  phgr <- phyloseq_sep_variable(phys, variable = group, drop_zeroes = drop_zeroes)

  ## All possible comparisons
  if(is.character(comparis) & comparis == "all"){
    compars <- combn(x = names(phgr), m = 2)

    ## Issue warning if there will be too many comparisons
    if(ncol(compars) > 100){
      warning("There will be too many possible pairwise comparisons (> 100). Please verify that it is intentional.\n")
    }

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
  names(res) <- plyr::aaply(.data = compars, .margins = 2, .fun = paste, collapse="__")

  ## Remove zero OTUs
  if(drop_zeroes == TRUE){
    res <- plyr::llply(
      .data = res,
      .fun = function(z){ phyloseq::prune_taxa(phyloseq::taxa_sums(z) > 0, z) })
  }

  return(res)
}


## Split phyloseq object for pairwise comparisons (by groupping variable)
phyloseq_sep_pairwise <- function(phys, group = "FactoryZone", comparis = "hardcoded", drop_zeroes = TRUE){
  # comparis = "all" (for all possible pairwise comparisons) or "hardcoded"

  ## Split by groups
  phgr <- phyloseq_sep_variable(phys, variable = group, drop_zeroes = FALSE)

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

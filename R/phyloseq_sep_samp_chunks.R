
#' @title Separate phyloseq-class object into a number of chunks by sample.
#' @description This function splits a phyloseq object by samples, returning a list of phyloseqs with N samples each.
#' @param physeq A phyloseq-class object
#' @param nchunks A number of samples per chunks
#' @param drop_zeroes Logical, indicating weather OTUs with zero abundance or samples with zero total abundance should be removed
#' @return List with phyloseq objects.
#' @export
#'
#' @examples
#' data("enterotype")
#' phyloseq_sep_samp_chunks(enterotype, nchunks = 5)
#'
phyloseq_sep_samp_chunks <- function(physeq, nchunks = NULL, drop_zeroes = T){
    
  ## Check if nchunks is specified
  if(is.null(nchunks)){
    stop("Error: number of samples is not specified.\n")
  }

  ## If there is just a single sample in phyloseq
  if(phyloseq::nsamples(physeq) == 1){
    cat("Warning: there is only one sample in the resulting list.\n")
    res <- list()
    res[[1]] <- physeq
  } else {
  ## If there are multiple samples

  ## Extract sample names
  smp <- phyloseq::sample_names(physeq)

  ## Split samples into N groups
  chunks <- chunk(x = smp, n = nchunks)
  # chunks <- chunk(x = 1:nsamples(physeq), n = nchunks)

  ## Extract samples
  res <- plyr::llply(.data = chunks, .fun = function(samples){ 
    phyloseq::prune_samples(sample_names(physeq) %in% samples, physeq)
    })

  }

  ## Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }

  return(res)
}


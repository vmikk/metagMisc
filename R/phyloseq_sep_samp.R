
#' @title Separate phyloseq-class object by sample.
#' @description This function splits a phyloseq object by sample, returning a list of objects whose components each correspond to a separate sample.
#' @param physeq A phyloseq-class object.
#' @param drop_zeroes Logical, indicating weather OTUs with zero abundance or samples with zero total abundance should be removed.
#' @return List with phyloseq objects.
#' @export
#'
#' @examples
#' data("esophagus")
#' phyloseq_sep_samp(esophagus)
#'
phyloseq_sep_samp <- function(physeq, drop_zeroes = T){
    require(phyloseq)
    require(plyr)

    if(nsamples(physeq) == 1){
        cat("Warning: there is only one sample in the resulting list.\n")
    }

    ## Check if there are samples with zero-totals
    ssum <- sample_sums(physeq)
    if(any(ssum == 0)){
      if(drop_zeroes == FALSE){
        stop("Error: there are samples with zero total abundance. Use 'drop_zeroes = T' to remove them.\n")
      }
      ## Remove samples with missing data
      if(drop_zeroes == TRUE){
        snames <- paste(names(ssum)[which(ssum == 0)], collapse = ", ")
        cat("Warning: samples with zero total abundance were removed.\n(Removed samples: ", snames, ")\n", sep = "")
        physeq <- prune_samples(ssum > 0, physeq)
      }
    }

    ## Extract sample names
    smp <- data.frame(samples = sample_names(physeq), stringsAsFactors = F)

    ## Extract samples
    res <- mlply(.data = smp, .fun = function(samples){ prune_samples(samples, x = physeq) })
    names(res) <- smp$samples

    ## Remove taxa with zero abundance
    if(drop_zeroes == TRUE){
        res <- llply(.data = res, .fun = function(x){ prune_taxa(taxa_sums(x) > 0, x) })
    }

    return(res)
}

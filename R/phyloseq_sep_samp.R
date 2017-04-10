

## Separate by sample
# This function splits a phyloseq object by sample, returning a list of objects whose components each correspond to a separate sample
phyloseq_sep_samp <- function(physeq, drop_zeroes = T){
    require(phyloseq)
    require(plyr)

    if(nsamples(physeq) == 1){
        cat("Warning: there is only one sample in the resulting list.\n")
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
# phyloseq_sep_samp(physeq)

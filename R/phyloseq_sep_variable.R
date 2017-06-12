
#' @title Split phyloseq-class object by sample-level variable.
#' @description This function splits a phyloseq object by sample meta-data, returning a list of objects whose components each correspond to a group of samples.
#' @param physeq A phyloseq-class object
#' @param variable Variable name (contained in \code{\link{sample_data}})
#' @param drop_zeroes Logical, indicating weather OTUs with zero abundance withing a group of samples should be removed
#' @return List with phyloseq objects.
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' phyloseq_sep_variable(GlobalPatterns, variable = "SampleType")
#'    
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
    require(phyloseq)
    require(plyr)

    ## Check the input
    if(is.null(sample_data(physeq, errorIfNULL = T))){
        stop("Sample data is missing in the phyloseq-object.\n")
    }

    if(nsamples(physeq) == 1){
        cat("Warning: there is only one sample in the resulting list.\n")
    }

    ## Extract sample meta-data
    smp <- data.frame(
        SID = sample_names(physeq),
        as(object = sample_data(physeq), Class = "data.frame"),
        stringsAsFactors = F)

    ## Exatract sample names by the specified variable
    svv <- dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })

    ## Extract samples by groupping variable
    res <- llply(.data = svv, .fun = function(z){ prune_samples(z, x = physeq) })
    
    ## Remove taxa with zero abundance
    if(drop_zeroes == TRUE){
        res <- llply(.data = res, .fun = function(x){ prune_taxa(taxa_sums(x) > 0, x) })
    }

    return(res)
}

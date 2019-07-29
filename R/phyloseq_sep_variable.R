
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
#'
#' # Split data by sample type (e.g., Soil, Ocean, etc.)
#' phyloseq_sep_variable(GlobalPatterns, variable = "SampleType")
#'
#' # Do not remove OTUs with total zero abundance within each sample type
#' phyloseq_sep_variable(GlobalPatterns, variable = "SampleType", drop_zeroes = FALSE)
#'
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
    
    # require(phyloseq)
    # require(plyr)

    ## Check the input
    if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
        stop("Sample data is missing in the phyloseq-object.\n")
    }

    ## Extract sample meta-data
    mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")

    if(!variable %in% colnames(mtd)){
        stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
    }

    if(class(mtd[, variable]) %in% c("integer", "numeric") ){
        if( length( unique(mtd[, variable]) ) > 5){
          stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
        } else {
          warning("Groupping variable is numeric and it was coerced to factor.\n")
          mtd[, variable] <- factor(mtd[, variable])
        }
    }

    if(length(table(mtd[, variable])) == 1){
        cat("Warning: there is only one group of samples in the resulting list.\n")
    }

    ## Add sample IDs to the meta-data
    smp <- data.frame(
        SID = phyloseq::sample_names(physeq),
        mtd,
        stringsAsFactors = F)

    ## Extract sample names by the specified variable
    svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })

    ## Extract samples by groupping variable
    res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })

    ## Remove taxa with zero abundance
    if(drop_zeroes == TRUE){
        res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
    }

    return(res)
}

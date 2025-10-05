#' @title Combine samples by summing OTU abundance
#' @description Sums OTU abundances of replicates or sample groups into single samples.
#' This function can either combine all samples into one single sample or 
#' combine samples within groups defined by a sample metadata variable.
#'
#' @param physeq A phyloseq-class object
#' @param variable Character string; name of the sample metadata variable to use for grouping. 
#'   If NULL (default), all samples are combined into a single sample.
#'
#' @return A phyloseq object with combined samples. Original sample metadata is removed
#'   when combining all samples, but taxonomy table and phylogenetic tree (if present) are preserved.
#'
#' @export
#' @seealso \code{\link{phyloseq_merge_samples}}, \code{\link{phyloseq_sep_variable}}
#' @examples
#' library(phyloseq)
#' data("GlobalPatterns")
#' 
#' # Combine all samples into a single sample
#' gp_total <- phyloseq_combine_samples(GlobalPatterns)
#' 
#' # Combine samples by SampleType
#' gp_by_type <- phyloseq_combine_samples(GlobalPatterns, variable = "SampleType")
#'
phyloseq_combine_samples <- function(physeq, variable = NULL){

  ## Function to combine all samples into a single one
  comb_phys <- function(phys){
    ## Remove sample metadata
    if(!is.null(phyloseq::sample_data(phys, errorIfNULL = FALSE))){ phys@sam_data <- NULL }
    
    ## Estimate OTU total abundance 
    tss <- phyloseq::taxa_sums(phys)

    ## Replace OTU table
    ott <- matrix(data = tss, ncol = 1, dimnames = list(names(tss), "Total"))
    phyloseq::otu_table(phys) <- phyloseq::otu_table(ott, taxa_are_rows = TRUE)
    
    return(phys)
  }

  ### All samples into a single one (no variables specified)
  if(is.null(variable)){
    res <- comb_phys(physeq)
  } else {
  ### Combine samples by group

    ## Extract slots from phyloseq object (later we'll return them)
    ## and remove them to save RAM
    if(!is.null(phyloseq::tax_table(physeq, errorIfNULL=F))){
      taxx <- phyloseq::tax_table(physeq)
      taxpresent <- TRUE
      physeq@tax_table <- NULL
    } else {
      taxpresent <- FALSE   # no tax_table in phyloseq
    }

    if(!is.null(phyloseq::phy_tree(physeq, errorIfNULL=F))){
      phyy <- phyloseq::phy_tree(physeq)
      phypresent <- TRUE
      physeq@phy_tree <- NULL
    } else {
      phypresent <- FALSE   # no phy_tree in phyloseq
    }


    ## Split samples by group
    pps <- phyloseq_sep_variable(physeq, variable = variable, drop_zeroes = FALSE)

    ## Combine OTUs by sample group
    ppc <- plyr::llply(.data = pps, .fun = comb_phys)

    ## Combine phyloseq objects into a single object
    ## using for-loop and merge_phyloseq  -- very slow for large number of groups
    # res <- ppc[[1]]
    # for(i in 2:length(ppc)){
    #   res <- phyloseq::merge_phyloseq(res, ppc[[i]])
    # }
    # rm(i)

    ## Extract OTU tables
    ppc <- plyr::llply(.data = ppc, .fun = function(x){ phyloseq_otu_to_df(x) })

    ## Merge OTU tables and rename samples
    ppc <- do.call(cbind, ppc)
    colnames(ppc) <- names(pps)

    ## Create phyloseq object with merged data
    res <- phyloseq::phyloseq(
              phyloseq::otu_table(ppc, taxa_are_rows = TRUE))

    ## Recover phyloseq slots
    if(taxpresent == TRUE){ res <- phyloseq::merge_phyloseq(res, taxx) }
    if(phypresent == TRUE){ res <- phyloseq::merge_phyloseq(res, phyy) }
  }

  return(res)
}

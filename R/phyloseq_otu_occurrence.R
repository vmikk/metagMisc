
#' @title Estimate species occurrence (within groups of samples)
#'
#' @param physeq A phyloseq-class object
#' @param variable Character string defining a variable name of the sample groups (this variable should be present in \code{\link[phyloseq]{sample_data}}) or NULL (no sample groups)
#' @param taxa_frequency Logical; if TRUE (default), relative frequency of species occurence within a sample group will be returned; if FALSE, number of samples within each sample group will be returned for each taxa
#' @param drop_zeroes Logical; if TRUE, taxa with total zero abundance will be removed
#' @param justdf Logical; if TRUE, return only a data frame with taxa occurrences; if FALSE (default), modified phyloseq object will be returned
#' @param long Logical; if TRUE, data frame with taxa occurrences will be returned in long format (with a single column defining the sample group); if FALSE (default), species occurrences will be returned in a wide format (sample groups as columns)
#'
#' @return A phyloseq-class object (if justdf = FALSE) or a data frame (justdf = TRUE).
#' @export
#'
#' @examples
#'
phyloseq_otu_occurrence <- function(physeq, variable = NULL,
  taxa_frequency = TRUE, drop_zeroes = FALSE, justdf = FALSE, long = FALSE){

  require(plyr)
  require(reshape2)

  ## Function to collapse samples into occurrences for a single sample group
  single_group_occurrence <- function(phys, rel = FALSE){

    ## Transform OTU abundances into presence-absence form
    sp_count <- apply(X = otu_table(phys),
        MARGIN = ifelse(taxa_are_rows(phys), yes = 1, no = 2),
        FUN = function(x){sum(x > 0)})

    ## Absolute occurrence (e.g., number of samples with the species)
    if(rel == FALSE){
      rez <- data.frame(Taxa = names(sp_count), Occurrence = sp_count, stringsAsFactors = F)
    }

    ## Relative occurrence (frequency)
    if(rel == TRUE){
      rez <- data.frame(Taxa = names(sp_count), Occurrence = sp_count / nsamples(phys), stringsAsFactors = F)
    }

    rownames(rez) <- NULL
    return(rez)
  }

  ## If no sample groups are specified
  if(is.null(variable)){
    res <- single_group_occurrence(physeq, rel = taxa_frequency)

  ## If multiple groups are specified
  } else {

    ## Split phyloseq by group
    pg <- phyloseq_sep_variable(physeq, variable, drop_zeroes = FALSE)

    ## Count species occurrence within each group
    resl <- ldply(.data = pg, .fun = single_group_occurrence, rel = taxa_frequency, .id = "SampleGroup")

    ## Reshape species occurrences into a wide format (samples as columns)
    res <- dcast(data = resl, Taxa ~ SampleGroup, value.var = "Occurrence")

  }

  ## Remove species with zero abundance (only from the wide-format data)
  if(drop_zeroes == TRUE){
    spsums <- rowSums(res[,-1]) > 0
    if(any(!spsums)){
      res <- res[spsums, ]
    }
  }

  ## Return only species occurrences
  if(justdf == TRUE){

    ## Return single-group data
    if(is.null(variable)){ return(res) }

    ## Return multiple-group data
    if(long == FALSE & !is.null(variable)){
      return(res)                    # return data in wide format
    } else {
      colnames(resl)[1] <- variable  # rename sample group variable in long data
      return(resl)                   # return data in long format
    }

  }

  ## Return transformed phyloseq object
  if(justdf == FALSE){

    ## Remove sample metadata
    if(!is.null(sample_data(physeq, errorIfNULL = T))){
      physeq@sam_data <- NULL
    }

    ## Prepare new OTU table
    rownames(res) <- res$Taxa
    res$Taxa <- NULL

    ## Replace OTU table with the new one
    otu_table(physeq) <- otu_table(res, taxa_are_rows = T)

    return(physeq)
  }
}


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
#' # Load data
#' data("GlobalPatterns")
#'
#' ## Return phyloseq-object with sample groups instead of samples
#' # With absolute counts (e.g., number of samples with the species)
#' gpa <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = F)
#' gpa
#' head( otu_table(gpa) )
#'
#' # With relative frequency of species occurence within a sample group
#' gpr <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = T)
#' gpr
#' head( otu_table(gpr) )
#'
#'
#' ## Return just a data frame with taxa occurrences
#' # In wide format (sample groups as columns)
#' gpw <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = T, drop_zeroes = F, justdf = T, long = F)
#' head(gpw)
#'
#' # In long format (with a single column with sample type)
#' gpl <- phyloseq_otu_occurrence(GlobalPatterns, variable = "SampleType", taxa_frequency = T, drop_zeroes = F, justdf = T, long = T)
#' head(gpl)
#'
phyloseq_otu_occurrence <- function(physeq, variable = NULL,
  taxa_frequency = "percentage", drop_zeroes = FALSE, justdf = FALSE, long = FALSE){

  # taxa_frequency = "count"      - just counts OTU occurence (number of samples)
  # taxa_frequency = "percentage" - proportion of sampling units that contains the OTU (sample units in which OTU occurred / total number of sample units examined)
  # taxa_frequency = "relfreq"    - relative frequency of OTU (frequency of each OTU / sum of the frequency of all OTUs)

  ## Function to collapse samples into occurrences for a single sample group
  single_group_occurrence <- function(phys, rel = "count"){

    ## Transform OTU abundances into presence-absence form
    sp_count <- apply(X = phyloseq::otu_table(phys),
        MARGIN = ifelse(phyloseq::taxa_are_rows(phys), yes = 1, no = 2),
        FUN = function(x){sum(x > 0)})

    ## Absolute occurrence (e.g., number of samples with the species)
    if(rel == "count"){
      rez <- data.frame(Taxa = names(sp_count), Occurrence = sp_count, stringsAsFactors = F)
    } else {

      ## Frequency of occurrence
      rez <- data.frame(Taxa = names(sp_count), Occurrence = sp_count / nsamples(phys), stringsAsFactors = F)

      ## Relative frequency (OTU frequency / sum of the frequency of all OTUs)
      if(rel == "relfreq"){
        rez$Occurrence <- rez$Occurrence / sum(rez$Occurrence)
      }
    }

    rownames(rez) <- NULL
    return(rez)
  } # end of single_group_occurrence

  ## If no sample groups are specified
  if(is.null(variable)){
    res <- single_group_occurrence(physeq, rel = taxa_frequency)

  ## If multiple groups are specified
  } else {

    ## Split phyloseq by group
    pg <- phyloseq_sep_variable(physeq, variable, drop_zeroes = FALSE)

    ## Count species occurrence within each group
    resl <- plyr::ldply(.data = pg, .fun = single_group_occurrence, rel = taxa_frequency, .id = "SampleGroup")

    ## Reshape species occurrences into a wide format (samples as columns)
    setDT(resl)
    res <- dcast(data = resl, Taxa ~ SampleGroup, value.var = "Occurrence")
    setDF(res)

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
    if(!is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))){
      physeq@sam_data <- NULL
    }

    ## Prepare new OTU table
    rownames(res) <- res$Taxa
    res$Taxa <- NULL

    ## Replace OTU table with the new one
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(res, taxa_are_rows = T)

    return(physeq)
  }
}

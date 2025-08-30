
#' @title Perform rarefaction and average relative OTU abundance
#'
#' @param physeq A phyloseq-class object
#' @param SampSize Rarefaction depth (number of reads to sample)
#' @param iter Number of rarefication iterations
#' @param parallel Logical or integer; if TRUE, attempts to run the function on multiple cores; integer defines the number of cores to use (if it is set to TRUE, all cores will be used)
#' @param verbose Logical; if TRUE, informational messages will be shown on screen
#' @param ... Additional arguments will be passed to \code{\link[metagMisc]{phyloseq_mult_raref}}
#'
#' @return phyloseq object with relative OTU abundance averaged across multiple rarefactions.
#' @export
#'
#' @examples
#' # Load data
#' data(GlobalPatterns)
#' 
#' # Subset data for faster processing
#' GP <- prune_taxa(taxa_sums(GlobalPatterns) >= 100, GlobalPatterns)
#' 
#' # Perform multiple rarefaction and average relative OTU abundances across rarefaction iterations
#' res <- phyloseq_mult_raref_avg(GP)
#' res
#' 
phyloseq_mult_raref_avg <- function(physeq, SampSize = NULL, iter = 1000, parallel = FALSE, verbose = TRUE, ...){

  # require(compositions)
  # require(plyr)
  # require(reshape2)

  ## Progress indicator
  if(verbose == TRUE){
    progress <- "text"
  } else {
    progress <- "none"
  }

  ## Transpose data if necessary
  if(taxa_are_rows(physeq) == FALSE){
    cat("Warning: transposing input data (OTUs should be arrenged by rows, and samples by columns)\n")
    physeq <- t(physeq)
  }

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

  if(!is.null(phyloseq::sample_data(physeq, errorIfNULL=F))){
    samm <- phyloseq::sample_data(physeq)
    sampresent <- TRUE
    physeq@sam_data <- NULL
  } else {
    sampresent <- FALSE   # no sample_data in phyloseq
  }


  ## Rarefy (do not remove zero-OTUs)
  if(verbose == TRUE){ cat("..Multiple rarefaction\n") }
  phys_raref <- phyloseq_mult_raref(physeq, SampSize = SampSize, iter = iter, multithread = parallel, trimOTUs = F, ...)

  ## Rename rarefied samples (add raref ID)
  if(verbose == TRUE){ cat("..Sample renaming\n") }
  samps_names <- vector()
  for(i in 1:length(phys_raref)){
    tmp_name <- paste(phyloseq::sample_names(phys_raref[[i]]), i, sep="__")
    samps_names <- c(samps_names, tmp_name)
    phyloseq::sample_names(phys_raref[[i]]) <- tmp_name
    rm(tmp_name)
  }
  rm(i)

  ### Combine rarefied OTU tables into a single table
  if(verbose == TRUE){ cat("..Rarefied data merging\n") }
  ## using for-loop and merge_phyloseq  -- very slow for large number of iterations
  # phys <- phys_raref[[1]]
  # for(i in 2:length(phys_raref)){
  #   phys <- phyloseq::merge_phyloseq(phys, phys_raref[[i]])
  # }
  # rm(i)

  ## Extract OTU tables
  phys_tabs <- plyr::llply(.data = phys_raref, .fun = function(x){
    phyloseq_otu_to_df(x)
  })

  ## Merge OTU tables
  # phys_tabs <- do.call("cbind", phys_tabs)  # loses rownames
  phys_tabs <- do.call(cbind, phys_tabs)      # why does it alters column names??
  colnames(phys_tabs) <- samps_names

  ## Create metadata
  metad <- data.frame(
    SampleID = gsub(pattern = "__[0-9]+$", replacement = "", x = samps_names),  # remove rarefaction ID
    # SampleRarefName = samps_names,
    stringsAsFactors = FALSE
    )
  rownames(metad) <- samps_names

  ## Create phyloseq object with merged data
  phys <- phyloseq::phyloseq(
            phyloseq::otu_table(phys_tabs, taxa_are_rows = TRUE),
            phyloseq::sample_data(metad))

  ## Extract rarefaction attributes and remove rarefied data to save RAM
  ## Add rarefaction parameters as attributes to the phyloseq object
  RarefactionDepth <- attr(phys_raref, which = "RarefactionDepth")
  RarefactionReplacement <- attr(phys_raref, which = "RarefactionReplacement")
  rm(phys_raref, phys_tabs)

  ## Split by sample
  if(verbose == TRUE){ cat("..Splitting by sample\n") }
  smps <- phyloseq_sep_variable(phys, variable = "SampleID", drop_zeroes = T)

  ## Average relative OTU abundances within each sample across rarefaction iterations
  if(verbose == TRUE){ cat("..OTU abundance averaging within rarefaction iterations\n") }

  ## OTU averaging function shortcut
  OTU_avg <- function(z){

    ## Function source code is in 'phyloseq_average.R'
    otuavg <- OTU_average(z, avg_type = "arithmetic", result = "matrix", verbose = FALSE)
              # zeroimp = FALSE, meth = "CZM"

    ## Add column name with OTU IDs
    otuavg <- dfRowName(x = otuavg, name = "OTU")
    return(otuavg)  # rows = OTUs
  }

  ## Batch averaging by sample
  smps_avg <- plyr::ldply(
    .data = smps,
    .fun = OTU_avg,
    .id = "SampleID",
    .progress = progress)

  ## Rename samples (for llply instead of ldply)
  # for(i in 1:length(smps_avg)){
  #   colnames(smps_avg[[i]]) <- names(smps)[i]
  # }
  # rm(i)

  ## Reshape data (rows = OTUs, columns = samples)
  smps_avg_wide <- reshape2::dcast(data = smps_avg, formula = OTU ~ SampleID, fill = 0, value.var = "Average")
  rownames(smps_avg_wide) <- smps_avg_wide$OTU
  smps_avg_wide$OTU <- NULL

  ## Re-create phyloseq object
  if(verbose == TRUE){ cat("..Re-create phyloseq object\n") }
  res <- phyloseq::phyloseq(
            phyloseq::otu_table(smps_avg_wide, taxa_are_rows = TRUE)
            )

  ## Recover phyloseq slots
  if(taxpresent == TRUE){ res <- phyloseq::merge_phyloseq(res, taxx) }
  if(phypresent == TRUE){ res <- phyloseq::merge_phyloseq(res, phyy) }
  if(sampresent == TRUE){ res <- phyloseq::merge_phyloseq(res, samm) }

  ## Add rarefaction parameters as attributes to the phyloseq object
  attr(res, which = "RarefactionDepth") <- RarefactionDepth
  attr(res, which = "RarefactionReplacement") <- RarefactionReplacement

  return(res)
}

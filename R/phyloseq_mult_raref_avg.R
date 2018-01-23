
## Perform rarefaction and average relative OTU abundance
phyloseq_mult_raref_avg <- function(z, SampSize = NULL, iter = 1000, parallel = FALSE, ...){

  # require(compositions)
  # require(plyr)
  # require(reshape2)

  ## Extract slots from phyloseq object (later we'll return them)
  ## and remove them to save RAM
  if(!is.null(phyloseq::tax_table(physeq, errorIfNULL=F))){
    taxx <- phyloseq::tax_table(physeq)
    taxpresent <- TRUE
    physeq@tax_table <- NULL
  } else {
    taxpresent <- FALSE   # no tax_table in pheloseq
  }

  if(!is.null(phyloseq::phy_tree(physeq, errorIfNULL=F))){
    phyy <- phyloseq::phy_tree(physeq)
    phypresent <- TRUE
    physeq@phy_tree <- NULL
  } else {
    phypresent <- FALSE   # no phy_tree in pheloseq
  }

  if(!is.null(phyloseq::sample_data(physeq, errorIfNULL=F))){
    samm <- phyloseq::sample_data(physeq)
    sampresent <- TRUE
    physeq@sam_data <- NULL
  } else {
    sampresent <- FALSE   # no sample_data in pheloseq
  }

  ## Rarefy (do not remove zero-OTUs)
  cat("..Multiple rarefaction\n")
  phys_raref <- phyloseq_mult_raref(z, SampSize = SampSize, iter = iter, multithread = parallel, trimOTUs = F, ...)

  ## Rename rarefied samples (add raref ID)
  cat("..Sample renaming\n")
  for(i in 1:length(phys_raref)){
    phyloseq::sample_names(phys_raref[[i]]) <- paste(phyloseq::sample_names(phys_raref[[i]]), i, sep="_")
  }
  rm(i)

  ## Combine rarefied OTU tables into a single table
  cat("..Rarefied data merging\n")
  phys <- phys_raref[[1]]
  for(i in 2:length(phys_raref)){
    phys <- phyloseq::merge_phyloseq(phys, phys_raref[[i]])
  }
  rm(i)

  ## Split by sample
  cat("..Splitting by sample\n")
  smps <- phyloseq_sep_variable(phys, variable = "MiSeqSample", drop_zeroes = T)

  ## Average relative OTU abundances within each sample across rarefaction iterations
  cat("..OTU abundance averaging within rarefaction iterations\n")

  ## Function to average OTU table
  OTU_average <- function(x){
    # x = phyloseq

    ## Convert to relative OTU abundance
    xx <- phyloseq::transform_sample_counts(x, function(OTU) OTU/sum(OTU) )

    ## Extract OTU table
    otus <- as.data.frame(phyloseq::otu_table(xx))

    ## Transpose OTU abundance table (samples must be ROWS from this step!)
    if(phyloseq::taxa_are_rows(x) == TRUE){
      otus <- t(otus)
    }

    ## Transform the data using the the centred log-ratio (Aitchison compositions)
    otucomp <- suppressWarnings( compositions::acomp(otus) )

    ## Average proportions
    # TO DO: add possibilty to specify a robust estimator ('robust = TRUE')
    otuavg <- suppressWarnings( compositions::mean.acomp(otucomp) )
    otuavg <- as.matrix(otuavg)    # it will be transposed here
    colnames(otuavg) <- "Average"  # rename average proporion column

    ## Add column name with OTU IDs
    otuavg <- dfRowName(x = otuavg, name = "OTU")

    return(otuavg)  # rows = OTUs
  }

  ## Batch averaging by sample
  smps_avg <- plyr::ldply(.data = smps, .fun = OTU_average, .id = "SampleID", .progress = "text")

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
  cat("..Re-create phyloseq object\n")
  res <- phyloseq::phyloseq(
            phyloseq::otu_table(smps_avg_wide, taxa_are_rows = TRUE),
            phyloseq::tax_table(z),
            phyloseq::sample_data(z))

  ## Recover phyloseq slots
  if(taxpresent == TRUE){ res <- phyloseq::merge_phyloseq(res, taxx) }
  if(phypresent == TRUE){ res <- phyloseq::merge_phyloseq(res, phyy) }
  if(sampresent == TRUE){ res <- phyloseq::merge_phyloseq(res, samm) }
  
  ## Add rarefaction parameters as attributes to the phyloseq object
  attr(res, which = "RarefactionDepth") <- attr(phys_raref, which = "RarefactionDepth")
  attr(res, which = "RarefactionReplacement") <- attr(phys_raref, which = "RarefactionReplacement")

  return(res)
}

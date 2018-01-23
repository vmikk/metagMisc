
## Perform rarefaction and average relative OTU abundance
phyloseq_mult_raref_avg <- function(physeq, SampSize = NULL, iter = 1000, parallel = FALSE, ...){

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
  phys_raref <- phyloseq_mult_raref(physeq, SampSize = SampSize, iter = iter, multithread = parallel, trimOTUs = F, ...)

  ## Rename rarefied samples (add raref ID)
  cat("..Sample renaming\n")
  samps_names <- vector()
  for(i in 1:length(phys_raref)){
    tmp_name <- paste(phyloseq::sample_names(phys_raref[[i]]), i, sep="__")
    samps_names <- c(samps_names, tmp_name)
    phyloseq::sample_names(phys_raref[[i]]) <- tmp_name
    rm(tmp_name)
  }
  rm(i)

  ### Combine rarefied OTU tables into a single table
  cat("..Rarefied data merging\n")
  ## using for-loop and merge_phyloseq  -- very slow for large number of iterations
  # phys <- phys_raref[[1]]
  # for(i in 2:length(phys_raref)){
  #   phys <- phyloseq::merge_phyloseq(phys, phys_raref[[i]])
  # }
  # rm(i)

  ## Extract OTU tables
  phys_tabs <- plyr::llply(.data = phys_raref, .fun = function(x){
    as.data.frame(phyloseq::otu_table(x))
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
  cat("..Splitting by sample\n")
  smps <- phyloseq_sep_variable(phys, variable = "SampleID", drop_zeroes = T)

  ## Average relative OTU abundances within each sample across rarefaction iterations
  cat("..OTU abundance averaging within rarefaction iterations\n")
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


## Function to average OTU table
OTU_average <- function(xx){
  # xx = phyloseq
  
  ## Remove OTUs with zero abundance
  # xx <- phyloseq::prune_taxa(phyloseq::taxa_sums(xx) > 0, xx)   # zeroes are already removed

  ## Convert to relative OTU abundance
  xx <- phyloseq::transform_sample_counts(xx, function(OTU) OTU/sum(OTU) )

  ## Extract OTU table
  otus <- as.data.frame(phyloseq::otu_table(xx))

  ## Transpose OTU abundance table (samples must be ROWS from this step!)
  if(phyloseq::taxa_are_rows(xx) == TRUE){
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

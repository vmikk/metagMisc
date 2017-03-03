## Functions to transform phyloseq OTU tables

## Cumulative sum scaling (CSS) normalization
# median scaling factor across samples will be used as default
# Paulson et al. Nature Methods 10, 1200–1202 (2013) doi:10.1038/nmeth.2658
# http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html
css_scaling <- function(physeq, norm = TRUE, log = TRUE, ...){
  require(metagenomeSeq)

  MGS <- phyloseq_to_metagenomeSeq(physeq)

  # get the normalized count matrix
  otu_norm <- MRcounts(MGS, norm = norm, log = log, ...)
  # exportMat(datt.norm, file = "tmp.txt")    # norm = TRUE, log = TRUE
  ## save sample statistics (sample scaling factor, quantile value, number of identified features and library size):
  # exportStats(dattM.sf, file = "tmp_stats.txt")

  # Substitue raw abundance to the css-normalized data
  physeq.tr <- physeq
  otu_table(physeq.tr) <- otu_table(otu_norm, taxa_are_rows = T)
  return(physeq.tr)
}


## Variance stabilizing transformation (VST)
# For downstream analysis use sample covariate information (blind = F) !!!!!
vst_blind_scaling <- function(physeq, dropneg = F, dropmissing = T, ...){
  ## dropneg = replace negative transformed values with 0

  require(DESeq2)

  dsc <- phyloseq_to_deseq2(physeq, design = formula(~ 1))
  # otu_norm <- varianceStabilizingTransformation(dsc, blind = T, ...)
  # otu_norm <- assay(otu_norm)

  dsc <- estimateSizeFactors(dsc)
  dsc <- estimateDispersions(dsc)
  otu_norm <- getVarianceStabilizedData(dsc)
  # Negative values probably correspond to “less than one count”

  # Set to zero all values less than zero
  if(dropneg == TRUE){
    otu_norm[otu_norm < 0.0] <- 0.0
  }

  # Substitue raw abundance to the variance stabilized data
  physeq.tr <- physeq
  otu_table(physeq.tr) <- otu_table(otu_norm, taxa_are_rows = T)

  # Remove missing OTUs
  if(dropneg == TRUE & dropmissing == TRUE){
    physeq.tr <- prune_taxa(taxa_sums(physeq.tr) > 0, physeq.tr)
  }

  return(physeq.tr)
}


## Regularized-log transformation
# rlog is preferable to the vst if the size factors vary widely
rlog_blind_scaling <- function(physeq, dropneg = F, dropmissing = T, ...){
  require(DESeq2)

  dsc <- phyloseq_to_deseq2(physeq, design = formula(~ 1))
  otu_norm <- rlogTransformation(dsc, blind = T, ...)
  otu_norm <- assay(otu_norm)

  # Set to zero all values less than zero
  if(dropneg == TRUE){
    otu_norm[otu_norm < 0.0] <- 0.0
  }

  # Substitue raw abundance to the rlog transformed data
  physeq.tr <- physeq
  otu_table(physeq.tr) <- otu_table(otu_norm, taxa_are_rows = T)

  # Remove missing OTUs
  if(dropneg == TRUE & dropmissing == TRUE){
    physeq.tr <- prune_taxa(taxa_sums(physeq.tr) > 0, physeq.tr)
  }

  return(physeq.tr)
}


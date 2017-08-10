## Functions to transform phyloseq OTU tables

## TO DO:
# - add warning about negative values
# - add function description
# - add documentation and links (e.g., to DESeq2 manual)
# - make single wrapper function for all methods ??


#' @title Cumulative sum scaling (CSS) normalization of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param norm Logical, return normalized counts
#' @param log Logical, apply a logarithmic transform (log2) to the normalized count data
#' @param ... Additional arguments will be passed to \code{\link{MRcounts}}
#'
#' @details Median scaling factor across samples will be used as default.
#' @return Phyloseq object with transformed counts in OTU table.
#' @export
#' @references
#' Paulson et al. Nature Methods 10, 1200–1202 (2013) doi:10.1038/nmeth.2658. http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html
#' @seealso \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{rlog}}, \code{\link{varianceStabilizingTransformation}}
#' @examples
#'
phyloseq_transform_css <- function(physeq, norm = TRUE, log = TRUE, ...){
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


#' @title Variance stabilizing transformation (VST) of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param dropneg Logical, replace negative transformed values with 0
#' @param dropmissing Logical, remove missing data
#@param ... Additional arguments will be passed to \code{\link{varianceStabilizingTransformation}
#'
#' @details For downstream analysis it could be better to use sample covariate information (blind = FALSE in \code{\link{varianceStabilizingTransformation}}).
#' @return Phyloseq object with transformed counts in OTU table.
#' @export
#' @seealso \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{phyloseq_transform_css}}, \code{\link{varianceStabilizingTransformation}}, \code{\link{rlog}}
#' @examples
#'
#'
phyloseq_transform_vst_blind <- function(physeq, dropneg = F, dropmissing = T, ...){

  require(DESeq2)

  # Add dummy sample data (phyloseq_to_deseq2 doesn't work without sample_data)
  if(is.null( sample_data(physeq, errorIfNULL = F) )){
    smpdat_nul <- TRUE
    smpdat <- data.frame(TMP = rep(1, times = nsamples(physeq)))
    rownames(smpdat) <- sample_names(physeq)
    sample_data(physeq) <- smpdat
  }

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

  # Remove dummy sample data if present
  if(smpdat_nul == TRUE){
    pp@sam_data <- NULL
  }

  return(physeq.tr)
}


#' @title Regularized-log (rlog) transformation of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param dropneg Logical, replace negative transformed values with 0
#' @param dropmissing Logical, remove missing data
#' @param ... Additional arguments will be passed to \code{\link{rlogTransformation}}
#'
#' @details rlog transformation (this function) is preferable to the vst (\code{\link{phyloseq_transform_vst_blind}}) if the size factors vary widely.
#' @return Phyloseq object with transformed counts in OTU table.
#' @export
#' @seealso \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_css}}, \code{\link{rlog}}, \code{\link{varianceStabilizingTransformation}}
#' @examples
#'
#'
phyloseq_transform_rlog_blind <- function(physeq, dropneg = F, dropmissing = T, ...){
  require(DESeq2)

  # Add dummy sample data (phyloseq_to_deseq2 doesn't work without sample_data)
  if(is.null( sample_data(physeq, errorIfNULL = F) )){
    smpdat_nul <- TRUE
    smpdat <- data.frame(TMP = rep(1, times = nsamples(physeq)))
    rownames(smpdat) <- sample_names(physeq)
    sample_data(physeq) <- smpdat
  }

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

  # Remove dummy sample data if present
  if(smpdat_nul == TRUE){
    pp@sam_data <- NULL
  }

  return(physeq.tr)
}



#' @title Log-transformation of OTU abundance table.
#'
#' @param physeq A phyloseq-class object
#' @param ... Additional arguments (e.g., "logbase" for the logarithm base) will be passed to \code{\link{decostand}}
#'
#' @details
#' Logarithmic transformation as suggested by Anderson et al. (2006) will be applied to the OTU table in phyloseq-object. First, non-integer data will be divided by smallest positive value. Second, log(x) + 1 for x > 0.
#' Default value of the logarithm base ("logbase" parameter) is 2.
#' This function is a wrapper to \code{\link{decostand}} in vegan.
#' @return Phyloseq object with log-transformed counts in OTU table.
#' @export
#' @seealso \code{\link{decostand}}
#' @references  Anderson MJ, Ellingsen KE and McArdle BH. (2006) Multivariate dispersion as a measure of beta diversity. Ecology Letters 9, P. 683–693.
#'
#' @examples
#' data(enterotype)
#' physeq_transform_anderson_log(enterotype)
#'
physeq_transform_anderson_log <- function(physeq, ...){
  require(vegan)
  otus <- as.data.frame( otu_table(physeq) )
  otu_log <- decostand(otus, method = "log", ...)
  rownames(otu_log) <- taxa_names(physeq)
  otu_table(physeq) <- otu_table(otu_log, taxa_are_rows = T)
  return(physeq)
}

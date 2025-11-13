#' Normalize phyloseq object using Scaling with Ranked Subsampling (SRS)
#'
#' @description 
#' Applies Scaling with Ranked Subsampling (SRS) normalization to a phyloseq object. 
#' SRS is a superior alternative to rarefying that preserves the original community 
#' structure with greater reproducibility and minimal subsampling error. Unlike rarefying, 
#' SRS maintains relative frequencies of taxa and provides deterministic results.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq-class}} object containing OTU/ASV abundance data
#' @param Cmin Numeric. The target number of counts per sample after normalization (minimum library size)
#' @param drop_zeros Logical, If \code{TRUE}, remove taxa with zero abundance after normalization. Default is \code{FALSE}
#' @param set_seed Logical, If \code{TRUE}, set random seed for reproducibility. Default is \code{TRUE}
#' @param seed Numeric. Random seed value when \code{set_seed = TRUE}. Default is 1
#' @param ... Additional arguments passed to the underlying SRS function
#'
#' @details
#' SRS (Scaling with Ranked Subsampling) is a two-step normalization procedure described 
#' by Beule & Karlovsky (2020):
#' 
#' \enumerate{
#'   \item \strong{Scaling step}: All OTU/ASV counts are divided by a scaling factor so 
#'         that the sum of scaled counts equals the target \code{Cmin}. Relative frequencies 
#'         remain unchanged.
#'   \item \strong{Ranked subsampling step}: Non-integer count values are converted to 
#'         integers using an algorithm that minimizes subsampling error while maintaining 
#'         the target total count.
#' }
#' 
#' SRS provides several advantages over traditional rarefying:
#' \itemize{
#'   \item Higher reproducibility (often deterministic results)
#'   \item Better preservation of alpha diversity measures
#'   \item Minimal variance in relative abundances of taxa
#'   \item Improved preservation of community structure for beta diversity analysis
#' }
#'
#' @return A \code{\link[phyloseq]{phyloseq-class}} object with SRS-normalized OTU/ASV table. 
#'   All samples will have exactly \code{Cmin} total counts after normalization.
#'
#' @references
#' Beule L, Karlovsky P. (2020) Improved normalization of species count data in ecology by scaling with ranked subsampling (SRS): application to microbial communities // PeerJ 8:e9593, DOI:10.7717/peerj.9593
#'
#' @seealso \code{\link{phyloseq_SRS_lowcount}} for SRS normalization with preservation of low-count samples
#' @importFrom phyloseq otu_table taxa_are_rows prune_taxa taxa_sums
#' @export
#'
#' @examples
#' library(phyloseq)
#' library(SRS)
#' 
#' # Load example data
#' data(esophagus)
#' 
#' # Check original library sizes
#' sample_sums(esophagus)
#' 
#' # Normalize to 100 reads per sample
#' ps_norm <- phyloseq_SRS(esophagus, Cmin = 100)
#' 
#' # Verify normalization
#' sample_sums(ps_norm)  # All should equal 100
#' 
#' # Remove zero-abundance taxa after normalization
#' ps_norm_clean <- phyloseq_SRS(esophagus, Cmin = 100, drop_zeros = TRUE)
#' ntaxa(ps_norm_clean)  # without zero-abundance taxa
#' ntaxa(ps_norm)        # with zero-abundance taxa, original number of taxa preserved
#' 
phyloseq_SRS <- function(physeq, Cmin, drop_zeros = FALSE, set_seed = TRUE, seed = 1, ...){

  ## Extract OTU table
  ## and transpose data (columns should be samples and rows are the OTU counts)
  data <- phyloseq_otu_to_df(physeq, taxa_as_rows = TRUE)

  ## Perform SRS-normalization
  res <- SRS(data = data, Cmin = Cmin, ...)
  rownames(res) <- rownames(data)

  ## Transpose back
  if(!taxa_are_rows(physeq)){
    res <- t(res)
  }

  ## Replace OTU table
  otu_table(physeq) <- otu_table(res, taxa_are_rows = taxa_are_rows(physeq))

  ## Remove zero-OTUs
  if(drop_zeros == TRUE){
    physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  }

  return(physeq)
}


#' SRS normalization with preservation of low-count samples
#'
#' @description 
#' Applies SRS normalization only to samples with sufficient read counts while preserving 
#' low-count samples unchanged. This approach prevents the loss of potentially valuable 
#' samples that have sequencing depths below the normalization threshold.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq-class}} object containing OTU/ASV abundance data
#' @param treshold Numeric. The minimum sequencing depth threshold. Samples with counts 
#'   below this value are preserved unchanged, while samples at or above this threshold 
#'   are normalized to exactly this value using SRS. Default is 10,000
#'
#' @details
#' This function implements a mixed normalization strategy:
#' 
#' \enumerate{
#'   \item \strong{Sample splitting}: Samples are divided into two groups based on the 
#'         \code{treshold} parameter
#'   \item \strong{Selective normalization}: Only high-count samples (\code{>= treshold}) 
#'         are normalized to exactly \code{treshold} counts using SRS
#'   \item \strong{Preservation}: Low-count samples (\code{< treshold}) remain unchanged 
#'         to avoid information loss
#'   \item \strong{Merging}: Both groups are recombined into a single phyloseq object
#' }
#' 
#' This approach is particularly useful when:
#' \itemize{
#'   \item The dataset contains samples with highly variable sequencing depths
#'   \item Some samples have low sequencing depth but contain unique or important taxa
#'   \item Complete removal of low-coverage samples would result in significant data loss
#'   \item Downstream analyses can handle mixed normalization states
#' }
#' 
#' \strong{Note}: The resulting dataset will have samples with different total counts 
#' (low-count samples retain their original depths, high-count samples are normalized 
#' to \code{treshold}). Consider whether this is appropriate for your specific analysis.
#'
#' @return A \code{\link[phyloseq]{phyloseq-class}} object where high-count samples 
#'   have been SRS-normalized to the threshold value, while low-count samples retain 
#'   their original sequencing depths.
#'
#' @seealso \code{\link{phyloseq_SRS}} for standard SRS normalization of all samples
#'
#' @importFrom phyloseq prune_samples sample_sums merge_phyloseq
#' @export
#'
#' @examples
#' # Load required packages
#' library(phyloseq)
#' library(SRS)
#' 
#' # Load example data
#' data(esophagus)
#' 
#' # Check distribution of library sizes
#' summary(sample_sums(esophagus))
#' 
#' # Apply mixed normalization (threshold = 220)
#' ps <- phyloseq_SRS_lowcount(esophagus, treshold = 220)
#' sample_sums(ps)
#'
phyloseq_SRS_lowcount <- function(physeq, treshold = 10000){
  
  ## Split data into low- and high- count parts
  pl <- prune_samples(sample_sums(physeq) < treshold, physeq)
  ph <- prune_samples(sample_sums(physeq) >= treshold, physeq)

  ## SRS-normalize high-count samples
  ps <- phyloseq_SRS(ph, Cmin = treshold, drop_zeros = FALSE)

  ## Merge phyloseq objects
  res <- merge_phyloseq(pl, ps)
  return(res)
}

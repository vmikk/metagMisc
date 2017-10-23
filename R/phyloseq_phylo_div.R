
#' @title Estimate phylogenetic diversity (PD), mean pairwise distance (MPD), and mean nearest taxon distance (MNTD).
#' @description Currently only non-abundance-weighted estimates are implemented.
#' @param physeq A phyloseq-class object (phylogenetic tree is required)
#' @param measures Character vector with diversity indices names ("PD", "MPD", "MNTD", "SES.PD", "SES.MPD", "SES.MNTD")
#' @param ... Additional arguments (standardize, null.model, abundance.weights) may be passed to \code{\link[PhyloMeasures]{pd.query}}
#' @details Implementation of the phylogenetic diversity measures is based on the \code{\link[PhyloMeasures]{PhyloMeasures}} package which is much faster than corresponding functions from \code{\link[picante]{picante}} package.
#'
#' MNTD index sometimes is also reffered as MNND (mean nearest neighbour distance).
#'
#' Diversity values could be standardized to the mean and expectation of the metric (SES, standardized effect size) by passing additional argument ('standardize = T') to the function or directly calling esimation of this measures ('measures = c("SES.PD", "SES.MPD", "SES.MNTD")').
#' Standardization leads to some of the most commonly used measures, including NRI (net relatedness index, the standardized version of MPD), NTI (nearest taxon index, derived from MNTD) and PDI (phylogenetic diversity index, derived from PD).
#' These metrics describe how diffrent the observed species assemblages is from other possible assemblages of the same size.
#' With this standardization, values of 0 are consistent with random phylogenetic structure, while phylogenetic clustering is associated with negative values, and phylogenetic overdispersion with positive values.
#' The default null model preserves only the species richness of samples (abundance is not maintained), which leads to the taxon shuff approach ("taxa.labels" null model from picante).
#' Note that standardization is performed by evaluating analytical expressions for the mean and variance of the basic measures, rather than by the slow and inexact randomization techniques (as in \code{\link[metagMisc]{phyloseq_phylo_ses}}).
#' @return A data.frame of the diversity estimates.
#' @export
#' @seealso
#' \code{\link[metagMisc]{phyloseq_phylo_ses}},
#' \code{\link[PhyloMeasures]{mpd.query}}, \code{\link[PhyloMeasures]{mntd.query}}, and \code{\link[PhyloMeasures]{pd.query}} from \code{\link[PhyloMeasures]{PhyloMeasures}} package;
#' \code{\link[picante]{mpd}}, \code{\link[picante]{mntd}}, and \code{\link[picante]{pd}} from \code{\link[picante]{picante}} package;
#' \code{\link[phyloseq]{estimate_richness}}, \code{\link[phyloseq]{phy_tree}}
#' @references
#' Tsirogiannis C., Sandel B. (2016) PhyloMeasures: a package for computing phylogenetic biodiversity measures and their statistical moments. Ecography, 39, 709–714. doi:10.1111/ecog.01814
#'
#' Webb C., Ackerly D., McPeek M., and Donoghue M. (2002) Phylogenies and community ecology. Annual Review of Ecology and Systematics, 33, 475-505.
#'
#' Faith D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#' @examples
#' # Load data
#' data(esophagus)
#'
#' # Estimate phylogenetic diversity (PD) and mean pairwise phylogenetic distance (MPD)
#' phyloseq_phylo_div(esophagus, measures = c("PD", "MPD"))
#'
#' # Standardized effect size of phylogenetic diversity (PDI)
#' phyloseq_phylo_div(esophagus, measures = "PD", standardize = T)
#' phyloseq_phylo_div(esophagus, measures = "SES.PD")  # the same
#'
phyloseq_phylo_div <- function(physeq, measures=c("PD", "MPD", "MNTD", "SES.PD", "SES.MPD", "SES.MNTD"), ...){

  require(PhyloMeasures)

  ## Prepare community matrix = matrix with binary (0/1) values, where each row represents a tip set
  ## Each column name in the matrix must match a tip label on the input tree

  ## Check the orientation of the OTU table
  trows <- taxa_are_rows(physeq)

  ## Extact OTU table
  comm <- as(object = otu_table(physeq), Class = "matrix")
  if(trows == TRUE){ comm <- t(comm) }

  ## Scale OTU abundance to presence/absence
  comm <- ifelse(comm > 0, 1, 0)

  ## Calculate diversity metrics for each community
  res <- vector("list")  # initialize results

  ## Basic metrics
  if("PD" %in% measures){
    res <- c(res, list(PD = pd.query(tree = phy_tree(physeq), matrix = comm, ...) ))
  }
  if("MPD" %in% measures){
    res <- c(res, list(MPD = mpd.query(tree = phy_tree(physeq), matrix = comm, ...) ))
  }
  if("MNTD" %in% measures){
    res <- c(res, list(MNTD = mntd.query(tree = phy_tree(physeq), matrix = comm, ...) ))
  }

  ## Standardized effect sizes
  if("SES.PD" %in% measures){
    res <- c(res, list(SES.PD = pd.query(tree = phy_tree(physeq), matrix = comm, standardize = T, ...) ))
  }
  if("SES.MPD" %in% measures){
    res <- c(res, list(SES.MPD = mpd.query(tree = phy_tree(physeq), matrix = comm, standardize = T, ...) ))
  }
  if("SES.MNTD" %in% measures){
    res <- c(res, list(SES.MNTD = mntd.query(tree = phy_tree(physeq), matrix = comm, standardize = T, ...) ))
  }


  ## Combine results
  res <- do.call("cbind", res)
  rownames(res) <- sample_names(physeq)
  res <- as.data.frame(res)

  return(res)
}

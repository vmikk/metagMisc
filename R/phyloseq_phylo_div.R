
#' @title Estimate phylogenetic diversity (PD), mean pairwise distance (MPD), and mean nearest taxon distance (MNTD).
#' @description Currently only non-abundance-weighted estimates are implemented.
#' @param physeq A phyloseq-class object (phylogenetic tree is required)
#' @param measures Character vector with diversity indices names ("PD", "MPD", "MNTD")
#' @param ... Additional arguments (standardize, null.model, abundance.weights) may be passed to \code{\link[PhyloMeasures]{pd.query}}
#' @details Implementation of the phylogenetic diversity measures is based on the \code{\link[PhyloMeasures]{PhyloMeasures}} package which is much faster than corresponding functions from \code{\link[picante]{picante}} package.
#'
#' MNTD index sometimes is also reffered as MNND (mean nearest neighbour distance).
#' @return A data.frame of the diversity estimates.
#' @export
#' @seealso
#' \code{\link[PhyloMeasures]{mpd.query}}, \code{\link[PhyloMeasures]{mntd.query}}, and \code{\link[PhyloMeasures]{pd.query}} from \code{\link[PhyloMeasures]{PhyloMeasures}} package;
#' \code{\link[picante]{mpd}}, \code{\link[picante]{mntd}}, and \code{\link[picante]{pd}} from \code{\link[picante]{picante}} package;
#' \code{\link[phyloseq]{estimate_richness}}, \code{\link[phyloseq]{phy_tree}}
#' @references
#' Faith D.P. (1992) Conservation evaluation and phylogenetic diversity. Biological Conservation, 61, 1-10.
#'
#' Webb C., Ackerly D., McPeek M., and Donoghue M. (2002) Phylogenies and community ecology. Annual Review of Ecology and Systematics, 33, 475-505.
#' @examples
#'
phyloseq_phylo_div <- function(physeq, measures=c("PD", "MPD", "MNTD"), ...){
  # measures = vector with diversity indices names
  # ... =

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

  if("PD" %in% measures){
    res <- c(res, list(PD = pd.query(tree = phy_tree(physeq), matrix = comm, ...) ))
  }
  if("MPD" %in% measures){
    res <- c(res, list(PD = mpd.query(tree = phy_tree(physeq), matrix = comm, ...) ))
  }
  if("MNTD" %in% measures){
    res <- c(res, list(PD = mntd.query(tree = phy_tree(physeq), matrix = comm, ...) ))
  }

  ## Combine results
  res <- do.call("cbind", res)
  rownames(res) <- sample_names(physeq)
  res <- as.data.frame(res)

  return(res)
}


## Estimate phylogenetic diversity (PD), mean pairwise distance (MPD), and mean nearest taxon distance (MNTD)
phyloseq_phylo_div <- function(physeq, measures=c("PD", "MPD", "MNTD"), ...){
  # measures = vector with diversity indices names
  # ... = other arguments (standardize, null.model, abundance.weights) may be passed to PhyloMeasures::pd.query

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

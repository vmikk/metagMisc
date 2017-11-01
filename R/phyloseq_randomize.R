
#' @title Randomize abundance table and phylogeny in phyloseq objects (for null model testing and simulation).
#'
#' @param physeq A phyloseq-class object
#' @param null_model Character string defining the null model (for the description of supported models see \code{\link[picante]{randomizeMatrix}})
#' @param verbose Logical; if TRUE, additional information messages will be displayed
#' @param ... Additional arguments may be passed to \code{\link[picante]{randomizeMatrix}}
#' @details Currently only null models from picante package are implemented.
#' @return A phyloseq-class object with randomized abundance table and/or phylogeny.
#' @export
#' @seealso \code{\link[picante]{randomizeMatrix}}, \code{\link[vegan]{commsim}}
#' @examples
#'
phyloseq_randomize <- function(physeq, null_model = "phylogeny.pool", verbose = T, ...){
  # c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")

  require(phyloseq)
  require(picante)
  require(vegan)

  ## Picante models
  pm <- c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")

  ## Vegan models
  # Binary null models
  v1 <- c("r00", "r0", "r0_old", "r1", "r2", "c0", "swap", "tswap", "curveball", "quasiswap", "backtracking", "backtrack")

  # Quantitative Models for Counts with Fixed Marginal Sums
  v2 <- c("r2dtable", "quasiswap_count")

  # Quantitative swap models
  v3 <- c("swap_count", "abuswap_r", "abuswap_c")

  # Quantitative Swap and Shuffle
  v4 <- c("swsh_samp", "swsh_both", "swsh_samp_r", "swsh_samp_c", "swsh_both_r", "swsh_both_c")

  # Quantitative Shuffle Methods
  v5 <- c("r00_ind", "r0_ind", "c0_ind", "r00_samp", "r0_samp", "c0_samp", "r00_both", "r0_both", "c0_both")

  # All vegan models
  vv <- c(v1, v2, v3, v4, v5)

  ## Print implementation details
  if(verbose == TRUE){
    if(null_model %in% pm){
      cat("Randomization null model is based on implementation from the 'picante' package.\n")
    }
    if(null_model %in% vv){
      cat("Randomization null model is based on implementation from the 'vegan' package.\n")
    }
    cat("Please cite it in the publications.\n")
  }

  ## TO DO - vegan models
  if(null_model %in% vv){
    stop("Vegan null models are currently not yet implemented.\n")
  }

  ### Extract phyloseq components
  ## OTU table
  comm <- as.data.frame(otu_table(physeq))
  if(!taxa_are_rows(physeq)){
    comm <- t(comm)
  }

  ## Phylo tree
  tree_null <- is.null(phy_tree(physeq, errorIfNULL=F))
  if(!tree_null){
    phy <- phy_tree(physeq)
  }


  if(null_model == "taxa.labels"){
    if(tree_null){ # No phylogeny -> shuffle taxa names
      rownames(comm) <- sample(taxa_names(physeq))
    } else {       # Phylogeny present -> shuffle tip names
      phy <- picante::tipShuffle(phy)
    }
  } # end of "taxa.labels"

  ## Randomize community data matrices with picante::randomizeMatrix
  if(null_model %in% c("frequency", "richness", "independentswap", "trialswap")){

    if(null_model == "sample.pool"){ null_model <- "richness" }  # this is the same models?
    # https://github.com/skembel/picante/blob/649edc7938b878429914c617e22c67198a8c189a/R/phylodiversity.R#L179

      comm <- t( randomizeMatrix(t(comm), null.model = null_model, ...) )
  }
  if(null_model == "phylogeny.pool"){
    if(tree_null){ stop("Error: phylogeny tree is not available; therefore 'phylogeny.pool' model is not applicable.\n") }
    comm <- t( randomizeMatrix(t(comm), null.model = "richness", ...) )
    phy <- picante::tipShuffle(phy)
  }

  ## Replace phyloseq slots with the randomized ones
  if(!taxa_are_rows(physeq)){  # transpose OTU table back
    comm <- t(comm)
  }
  otu_table(physeq) <- otu_table(comm, taxa_are_rows = TRUE)

  if(!tree_null){
    phy_tree(physeq) <- phy
  }

  return(physeq)
}

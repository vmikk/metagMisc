
#' @title Convert phyloseq data to entropart MetaCommunity object.
#' @description The phyloseq data is converted to the MetaCommunity object (entropart package),
#' which can then be used for measurement and partitioning of diversity, including
#' species-neutral, phylogenetic and functional diversity metrics.
#' @param physeq A phyloseq-class object
#' @param wei A vector of positive numbers equal to community weights (could be NULL for equal weights)
#'
#' @return A \code{\link[entropart]{MetaCommunity}} class object
#' @seealso \code{\link[entropart]{MetaCommunity}}
#' @export
#'
#' @examples
#'
phyloseq_to_MetaCommunity <- function(physeq, wei = NULL){

  # require(entropart)

  ## Extract OTU table
  otus <- as.data.frame(phyloseq::otu_table(physeq))
  if(phyloseq::taxa_are_rows(physeq) == FALSE){ otus <- t(otus) }

  ## Prepare equal community weights
  if(is.null(wei)){
    wei <- rep(1, ncol(otus))
  }
  Weights <- data.frame(Communities = colnames(otus), Weights = wei)

  ## Convert phyloseq class to MetaCommunity class
  MC <- entropart::MetaCommunity(otus, Weights)

  return(MC)
}

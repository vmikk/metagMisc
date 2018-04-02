
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
#' # Load data
#' data("esophagus")
#' 
#' library(entropart)
#' 
#' # Convert phyloseq to MetaCommunity class
#' eso <- phyloseq_to_MetaCommunity(esophagus)
#' summary(eso)
#' 
#' # Estimate diversity (Shannon diversity, q=1)
#' ad <- AlphaDiversity(eso, q = 1, Correction = "None")
#' summary(ad)
#' plot(ad)
#' 
#' # Diversity partitioning into alpha and beta components
#' dp <- DivPart(q = 1, eso, Correction = "None")
#' summary(dp)
#' plot(dp)
#' 
#' # Diversity profile
#' dpr <- DivProfile(q.seq = seq(0, 2, 0.1), eso, Correction = "None")
#' summary(dpr)
#' plot(dpr)
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

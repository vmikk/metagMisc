
#' Convert between Bray-Curtis and Jaccard dissimilarity measures
#'
#' @description 
#' Converts dissimilarity values between Bray-Curtis and Jaccard measures 
#' using a mathematical transformation.
#'
#' @param x Numeric vector or distance matrix containing dissimilarity values (0-1 range)
#' @param to Character string specifying conversion direction. Either "bray" (convert 
#'   from Jaccard to Bray-Curtis) or "jaccard" (convert from Bray-Curtis to Jaccard)
#'
#' @details
#' The conversion formulas are based on the mathematical relationship between these 
#' dissimilarity measures:
#' \itemize{
#'   \item Jaccard to Bray-Curtis: \code{x / (2 - x)}
#'   \item Bray-Curtis to Jaccard: \code{(2 * x) / (1 + x)}
#' }
#' 
#' Both measures range from 0 (identical) to 1 (completely dissimilar), but they 
#' weight differences differently. Bray-Curtis gives more weight to abundant species, 
#' while Jaccard treats all species equally.
#'
#' @return Numeric vector or matrix of the same structure as input, containing 
#'   converted dissimilarity values
#'
#' @examples
#' # Convert Jaccard to Bray-Curtis
#' jaccard_vals <- c(0.2, 0.5, 0.8)
#' bray_vals <- jaccard_bray(jaccard_vals, to = "bray")
#' 
#' # Convert back to Jaccard
#' converted_back <- jaccard_bray(bray_vals, to = "jaccard")
#' 
#' # Works with distance matrices too
#' \dontrun{
#' library(vegan)
#' data(dune)
#' jac_dist <- vegdist(dune, method = "jaccard")
#' bray_dist <- jaccard_bray(jac_dist, to = "bray")
#' }
#'
#' @export
jaccard_bray <- function(x, to = "bray"){

  ## Convert Bray-Curtis to Jaccard dissimilarity
  if(to == "jaccard"){
    res <- (2*x)/(1+x)
  }

  ## Convert Jaccard to Bray-Curtis dissimilarity 
  if(to == "bray"){
    res <- x/(2-x)
  }

  return(res)
}

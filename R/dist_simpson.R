
#' Calculate Simpson's dissimilarity coefficient
#'
#' @description 
#' Computes pairwise Simpson's dissimilarity coefficient between samples using either 
#' abundance-weighted or presence-absence data.
#'
#' @param x A numeric matrix or data frame containing community data where rows are samples 
#'   and columns are taxa (species/OTUs)
#' @param abundance Logical, if \code{TRUE} (default), compute the abundance-
#'   weighted version. If \code{FALSE}, compute the binary (presence–absence) version.
#'
#' @details
#' 
#' \strong{Presence–absence} (\code{abundance = FALSE}):
#' Let \eqn{a} be the number of shared taxa and \eqn{b} and \eqn{c} the numbers
#' unique to each sample. Simpson's dissimilarity (Lennon et al., 2001) is
#' \deqn{β_{sim} = \frac{\min(b,,c)}{a + \min(b,,c)}\ ,}
#' which equals 0 for identical lists and 1 when no taxa are shared. Implemented
#' via \code{vegan::betadiver(x > 0, method = "sim")}, i.e. any positive value
#' is treated as presence.
#'
#' \strong{Abundance-weighted} (\code{abundance = TRUE}):
#' Rows are standardized to relative abundances (row sums = 1). For a pair of
#' samples, let \eqn{U} and \eqn{V} be the sums of relative abundances (in each
#' sample) over taxa shared by the pair. Following Chao et al. (2006), the
#' Simpson \emph{similarity} is
#' \deqn{S = \frac{UV}{,UV + \min{U(1 - V),; V(1 - U)},},}
#' and the reported \emph{dissimilarity} is \eqn{1 - S}. This reduces to the
#' binary form when abundances are equal across taxa.
#'
#' \strong{Data handling and edge cases:}
#' \itemize{
#' \item Pairs with no shared taxa (\eqn{U = 0} or \eqn{V = 0}) return
#' dissimilarity 1.
#' \item Any pair involving an all-zero row returns dissimilarity 1.
#' }
#'
#' \emph{Estimator note:} The abundance-weighted calculation is a simple
#' plug-in estimator and does not include unseen-shared-species bias correction
#' (Chao et al., 2006); in undersampled communities this can inflate
#' dissimilarity.
#'
#' @return A distance object of class \code{dist} containing pairwise Simpson 
#'   dissimilarity values. The object includes an \code{abundance_weighted} 
#'   attribute indicating the calculation method used.
#'
#' @references
#' 
#' Lennon JJ, Koleff P, Greenwood JJD, Gaston KJ (2001) The geographical structure of British bird distributions: 
#' Diversity, spatial turnover and scale // Journal of Animal Ecology, V. 70, P. 966–979, 
#' DOI:10.1046/j.0021-8790.2001.00563.x
#' 
#' Chao A, Chazdon RL, Colwell RK, Shen T-J (2006) Abundance-based similarity indices and their estimation 
#' when there are unseen species in samples // Biometrics, V. 62 (2), P. 361–371, 
#' DOI:10.1111/j.1541-0420.2005.00489.x
#'
#' @examples
#' library(vegan)
#' data(BCI)
#' 
#' # Abundance-weighted Simpson dissimilarity (default)
#' dist_abund <- dist_simpson(BCI)
#' 
#' # Presence-absence Simpson dissimilarity  
#' dist_pa <- dist_simpson(BCI, abundance = FALSE)
#' 
#' @importFrom vegan decostand betadiver
#' @export
#' 
dist_simpson <- function(x, abundance = TRUE){

  ## Function for Simpson dissimilarity for a pair of samples
  abund_simps_pair <- function(x){
    # x = data.frame or matrix with two rows (samples), relative species abundances

    ## Remove double zeros and subset to shared species
    xx <- x[, colSums(x > 0) == 2, drop = FALSE]

    if(ncol(xx) > 0){

    ## Estimate sum of total [relative] abundances of shared species for each community
    U <- rowSums(xx[1,, drop = FALSE])
    V <- rowSums(xx[2,, drop = FALSE])

    ## Estimate Simpson's / Lennon's dissimilarity coefficient (see Chao et al., 2006 Biometrics; Table 2)
    UV <- U*V
    dd <- UV / (UV + min(U - UV, V - UV))
    dd <- 1 - dd
    } else {
      # No shared species
      dd <- 1
    }
    return(dd)
  }

  ## Convert to relative abundances
  x <- decostand(x, method = "total", MARGIN = 1)

  ## Abundance-based Simpson's dissimilarity
  if(abundance == TRUE){

    ## Initialize distance matrix
    dd <- matrix(data = NA, nrow = nrow(x), ncol = nrow(x),
                dimnames = list(rownames(x), rownames(x)))

    ## Estimate pairwise dissimilarities
    cmb <- combn(x = 1:nrow(x), m = 2)

    for(i in 1:ncol(cmb)){
      s1 <- cmb[1, i]
      s2 <- cmb[2, i]
      ds <- abund_simps_pair(x = x[c(s1, s2), ])
      dd[s1, s2] <- ds
      dd[s2, s1] <- ds
      rm(ds)
    }

    diag(dd) <- 0
    dd <- as.dist(dd)
    attr(x = dd, which = "abundance_weighted") <- TRUE

  ## Presence-absence Simpson's dissimilarity
  } else {

    dd <- vegan::betadiver(x = x, method = "sim")
    attr(x = dd, which = "abundance_weighted") <- FALSE
  }

  return(dd)
}



## Function to estimate abundance-based version of Simpson's dissimilarity coefficient
dist_simpson_abund <- function(x){
  # x = otu table, taxa = columns

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
  return(dd)
}
# e.g. dist_simpson_abund(BCI)

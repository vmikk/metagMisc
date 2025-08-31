gaussian_kernel <- function(d, bandwidth = 1, invert = TRUE) {

  # Validate inputs
  if(! "dist" %in% class(d)) {
    stop("Input must be an object of class 'dist'.\n")
  }
  if(bandwidth <= 0) {
    stop("Bandwidth must be positive.\n")
  }
  
  # Compute Gaussian kernel
  res <- exp(- (d^2) / (2 * bandwidth^2))
  
  ## Convert back to dissimilarity
  if(invert == TRUE){
    res <- 1 - res
  }
  
  return(res)
}


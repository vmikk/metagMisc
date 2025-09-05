
#' Gaussian kernel transformation
#'
#' @param d A distance matrix (class `dist`).
#' @param bandwidth A positive number representing the bandwidth of the kernel
#'   (default is 1). The bandwidth controls the spread of the kernel; a smaller
#'   bandwidth leads to a kernel that drops off more quickly.
#' @param invert Logical (default is TRUE); if TRUE, converts the Gaussian kernel into a measure of dissimilarity.
#'
#' @details
#' The Gaussian kernel transformation is used to convert distances into a
#' measure of similarity or dissimilarity. It gives more weight to points that
#' are close to each other and rapidly decreases this weight as the distance
#' increases. This feature is particularly useful in ecological and genetic
#' data, where closer relationships (e.g., closely related species or similar
#' genetic sequences) are often more important than distant ones. The
#' bandwidth parameter allows for control over the sensitivity to distance. 
#' A smaller bandwidth makes the kernel more sensitive to small distances,
#' emphasizing local structure, while a larger bandwidth smooths out local
#' variations, focusing on global structure.
#'
#' @return An object of class `dist`, with values in the range between 0 and 1.
#' These values can be interpreted as the probability of similarity 
#' (if `invert = FALSE`) or dissimilarity (if `invert = TRUE`).
#' 
#' @export
#'
#' @examples
#' d <- dist(1:10)
#' dk <- gaussian_kernel(d)
#' 
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


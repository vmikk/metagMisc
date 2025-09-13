
#' Gaussian kernel transformation
#'
#' @param d A distance matrix (class `dist`).
#' @param bandwidth A positive number representing the bandwidth of the kernel
#'   (default is 1). The bandwidth controls the spread of the kernel; a smaller
#'   bandwidth leads to a kernel that drops off more quickly.
#' @param invert Logical (default is TRUE); if TRUE, converts the Gaussian kernel into a measure of dissimilarity.
#' @param normalize_to_max Logical (default is FALSE); if TRUE, normalizes the
#'   transformed values to the [0, 1] range using `max_distance` as the maximum
#'   possible distance for the original metric.
#' @param max_distance A positive number (default is 1) indicating the maximum
#'   possible dissimilarity in the original metric. Used only when
#'   `normalize_to_max = TRUE`.
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
#' When `normalize_to_max = TRUE`, the output is rescaled to ensure the full
#' [0, 1] range is attainable for the selected `bandwidth` and chosen
#' `max_distance` (denoted as D_max):
#' - If `invert = TRUE` (dissimilarity): delta(d) = 1 - exp(-d^2/(2*sigma^2)), and
#'   delta*(d) = delta(d) / delta(D_max).
#' - If `invert = FALSE` (similarity): s(d) = exp(-d^2/(2*sigma^2)), and
#'   s*(d) = {s(d) - s(D_max)} / {1 - s(D_max)}.
#'
#' This preserves monotonicity and ensures standardized scaling across datasets
#' when a common D_max is used (e.g., 1 for Bray–Curtis or Jaccard dissimilarities).
#'
#' @return An object of class `dist`, with values in the range between 0 and 1.
#' These values can be interpreted as the probability of similarity 
#' (if `invert = FALSE`) or dissimilarity (if `invert = TRUE`).
#' 
#' @export
#'
#' @examples
#' library(ggplot2)
#' 
#' d  <- dist(1:100) / 100
#' dk <- gaussian_kernel(d)
#  shepard_plot(dis = d, ord = dk) + labs(x = "Original distance", y = "Transformed distance")
#' 
#' dk <- rbind(
#' data.frame(d = as.vector(d), Bandwidth = "1", dk = as.vector(gaussian_kernel(d, bandwidth = 1))),
#' data.frame(d = as.vector(d), Bandwidth = "1.5", dk = as.vector(gaussian_kernel(d, bandwidth = 1.5))),
#' data.frame(d = as.vector(d), Bandwidth = "2", dk = as.vector(gaussian_kernel(d, bandwidth = 2))),
#' data.frame(d = as.vector(d), Bandwidth = "4", dk = as.vector(gaussian_kernel(d, bandwidth = 4))) )
#' 
#' ggplot(dk, aes(x = d, y = dk, color = Bandwidth)) +
#    geom_line() +
#'   labs(x = "Original distance", y = "Transformed distance") + theme_classic()
#'
#' # Normalize by the maximum possible distance
#' dkn <- rbind(
#' data.frame(d = as.vector(d), Bandwidth = "1", dk = as.vector(gaussian_kernel(d, bandwidth = 1, normalize_to_max = TRUE))),
#' data.frame(d = as.vector(d), Bandwidth = "2", dk = as.vector(gaussian_kernel(d, bandwidth = 2, normalize_to_max = TRUE))))
#' 
#' ggplot(dkn, aes(x = d, y = dk, color = Bandwidth)) +
#'   annotate("segment", x=-Inf, xend=Inf,y=-Inf, yend=Inf, color = "grey70", linetype = "dashed") +
#'   geom_line() +
#'   labs(x = "Original distance", y = "Transformed distance (normalized)") + theme_classic()
#' 
gaussian_kernel <- function(d, bandwidth = 1, invert = TRUE,
                            normalize_to_max = FALSE, max_distance = 1) {

  # Validate inputs
  if(! "dist" %in% class(d)) {
    stop("Input must be an object of class 'dist'.\n")
  }
  if(bandwidth <= 0) {
    stop("Bandwidth must be positive.\n")
  }
  if(normalize_to_max) {
    if(!is.numeric(max_distance) || length(max_distance) != 1 || !is.finite(max_distance) || max_distance <= 0) {
      stop("When 'normalize_to_max = TRUE', 'max_distance' must be a single positive, finite number.\n")
    }
    if(max(as.vector(d), na.rm = TRUE) > max_distance) {
      warning("Some observed distances exceed 'max_distance'. Normalized values may exceed 1.\n")
    }
  }
  
  # Compute Gaussian kernel
  res <- exp(- (d^2) / (2 * bandwidth^2))
  
  ## Convert back to dissimilarity
  if(invert == TRUE){
    res <- 1 - res
  }
  
  ## Optional normalization to the value at D_max (done after invert step)
  if(normalize_to_max) {
    s_Dmax <- exp(- (max_distance^2) / (2 * bandwidth^2))
    denom <- 1 - s_Dmax
    if(denom <= 0) {
      stop("Normalization denominator is non-positive; check 'max_distance' and 'bandwidth'.\n")
    }
    if(denom < 1e-12) {
      warning("Normalization denominator is very small; results may be numerically unstable.\n")
    }
    if(invert == TRUE) {
      res <- res / denom
    } else {
      res <- (res - s_Dmax) / denom
    }
  }
  
  return(res)
}


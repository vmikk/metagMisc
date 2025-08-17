
#' Convert between Bray-Curtis and Jaccard/Ružička dissimilarity measures
#'
#' @description 
#' Exact elementwise conversion between Bray–Curtis (BC) and Jaccard/Ružička (JR) 
#' dissimilarities. The mapping is a one-to-one monotonic transform:
#' \deqn{BC = JR / (2 - JR), \qquad JR = (2 * BC) / (1 + BC).}
#'
#' @param x Numeric vector, matrix, or dist object of dissimilarity values (0-1 range)
#' @param to Character string specifying conversion direction. Either "bray" (convert 
#'   from Jaccard to Bray-Curtis) or "jaccard" (convert from Bray-Curtis to Jaccard)
#'
#' @details
#' Let \eqn{W = \sum_i \min(x_i, y_i)} (shared abundance) and 
#' \eqn{M = \sum_i \max(x_i, y_i)} (union abundance). Then:
#' \itemize{
#'   \item Bray–Curtis: \eqn{BC = (M - W) / (M + W)}.
#'   \item Jaccard/Ružička (quantitative Jaccard): \eqn{JR = 1 - W/M = (M - W)/M}.
#' }
#' Eliminating \eqn{W} and \eqn{M} gives the conversion formulas above.  
#'
#' This identity holds for both presence–absence (binary) and abundance data:
#' with binary data, \eqn{JR} reduces to the usual (binary) Jaccard dissimilarity; 
#' with quantitative data, \eqn{JR} is the Ružička (a.k.a. quantitative Jaccard) 
#' dissimilarity. In both cases joint absences are ignored.
#'
#' For presence–absence Jaccard, all species contribute equally (no abundance weighting). 
#' For quantitative JR and Bray–Curtis, both are abundance-based; 
#' they differ only by normalization (JR uses \eqn{M}, 
#' BC uses \eqn{M+W}). 
#' The conversion is strictly increasing, so ranks are preserved, but 
#' metric properties may differ (BC is a semimetric; JR is a metric).
#'
#' NB! This function expects \emph{dissimilarities}. If you have a
#' Jaccard \emph{similarity} \eqn{J}, first convert to dissimilarity \eqn{JR = 1 - J}.
#'
#'
#' @return An object of the same shape/class as \code{x} with converted values.
#'
#' @examples
#' # Numeric vector
#' bc <- c(0.1, 0.5, 1.0)
#' jr <- bray_jaccard_convert(bc, to = "jaccard")   # -> (2*bc)/(1+bc)
#' all.equal(jr, (2*bc)/(1+bc))
#'
#' # With vegan distances (presence–absence Jaccard and Bray–Curtis on binary data)
#' \dontrun{
#'   library(vegan); data(dune)
#'   bc_bin  <- vegdist(dune > 0, method = "bray")
#'   jac_bin <- vegdist(dune > 0, method = "jaccard")
#'   all.equal(bray_jaccard_convert(bc_bin, to = "jaccard"), jac_bin)
#'
#'   # Quantitative Jaccard/Ružička
#'   bc_abund <- vegdist(dune, method = "bray")
#'   jr_abund <- vegdist(dune, method = "ruzicka")
#'   all.equal(bray_jaccard_convert(bc_abund, to = "jaccard"), jr_abund)
#' }
#'
#' @export
#' 
bray_jaccard_convert <- function(x, to = "bray"){

  if(length(to) > 1){
    stop("Argument 'to' must be a single character string.\n")
  }
  to <- match.arg(to)

  if(any(x < 0 | x > 1, na.rm = TRUE)){
    stop("All dissimilarities must be in [0, 1].\n")
  }

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

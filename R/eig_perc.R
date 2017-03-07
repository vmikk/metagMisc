
#' @title Estimate the proportion of explained variance from the eigenvalues.
#'
#' @param eig Numeric vector of eigenvalues (e.g., from PCA or PCoA)
#' @param positive Logical, preserve only positive eigenvalues (if TRUE negative eigenvalues will be removed)
#' @param plot Logical, show a scree plot
#' @param ... Additinal arguments will be passed to a plot function
#'
#' @details
#' Negative eigenvalue could be generated in the analysis of semi- or non-metric dissimilarity measures.
#' @return Numeric vector of the amounts of total variance captured by principal components or principal coordinates.
#' @export
#'
#' @examples
#' pco <- cmdscale(d = eurodist, k = 10, eig = TRUE)
#' eig_perc(pco$eig, positive = T, plot = T)
#'
eig_perc <- function(eig, positive = T, plot = F, ...){

  if(positive == T) { eig <- eig[eig > 0] }

	res <- round(eig / sum(eig), 2)

	if(plot == T) {
	  plot(res, type="b", pch=16, ylab = "Variance explained", las = 1, ...)
	}

	return(res)
}

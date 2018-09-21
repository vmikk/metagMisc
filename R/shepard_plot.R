
#' @title Plot Shepard diagram
#' @description This function helps to estimate the representativeness of ordinations obtained using any reduced-space ordination method.
#' @param dis Original dissimilarity matrix (class "dist")
#' @param ord Ordination scores (class "matrix" or "data.frame") or (class "dist")
#' @param k Number of dimensions (default, 2) that will be taken into account (only if the ordination scores were provided in 'ord')
#'
#' @details This function compares original dissimilarities among the objects with Euclidean distances in reduced space obtained with ordination (object scores or coordinates of the objects in the reduced space).
#' If the points are close to diagonal then the projection in reduced space accounts for a high fraction of the variance.
#'
#' @return Invisible returns ggplot-object.
#' @export
#' @references Legendre P. & Legendre L. Numerical Ecology. 2012. 3rd Ed. - Chapter 9 (Ordination in reduced space).
#' @seealso \code{\link{vegan::stressplot}}, \code{\link{MAAS::Shepard}}
#' @examples
#' library(vegan)
#' data(dune)
#' ds <- decostand(dune, method = "total", MARGIN = 1)  # standardize data to species relative abundance
#' dd <- vegdist(ds, method = "bray")                   # compute dissimilarity between samples
#' ord <- monoMDS(dd, k = 5)                            # NMDS (non-metric multidimensional scaling)
#' dd.ord <- dist( scores(ord) )                        # Euclidean distance in reduced space
#'
#' shepard_plot(dis = dd, ord = dd.ord)               # Bray-Curtis vs distance in reduced space (k-dimensional)
#' shepard_plot(dis = dd, ord = scores(ord), k = 2)   # Bray-Curtis vs distance in reduced space (only first 2 dimensions)
#' shepard_plot(dis = dd, ord = scores(ord), k = 3)   # Bray-Curtis vs distance in reduced space (first 3 dimensions)
#'
#' # Example with phyloseq wrappers
#' data("GlobalPatterns")
#' dd <- distance(GlobalPatterns, method = "unifrac", type = "samples")
#' ord <- ordinate(GlobalPatterns, method = "NMDS", distance = dd)
#' shepard_plot(dis = dd, ord = scores(ord))
#'
shepard_plot <- function(dis, ord, k=2){

  require(ggplot2)

  ## If the ordination scores were provided
  if("matrix" %in% class(ord) || "data.frame" %in% class(ord)){
    if(is.null(k)){ dd_ord <- dist(ord) }       # take all axes
    if(!is.null(k) & ncol(ord) > 1){            # take only a few ordination axes
     dd_ord <- dist(ord[,1:k])
    }
  }

  ## If the distance matrix from ordination was provided
  if("dist" %in% class(ord)){
    dd_ord <- ord
  }

  ## Prepare data.frame for ggplot
  dtt <- data.frame(Observed = as.vector(dis), Ordination = as.vector(dd_ord))

  ## Plot
  pp <- ggplot(data = dtt, aes(x = Observed, y = Ordination)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=Inf, color="grey", linetype = "longdash") +
    geom_point() +
    theme_classic() +
    labs(x = "Distance in multidimensional space", y = "Distance in reduced space")

  print(pp)
  invisible(pp)
}

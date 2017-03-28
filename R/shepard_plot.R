
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


## Example
library(vegan)
data(dune)
ds <- decostand(dune, method = "total", MARGIN = 1)  # standardize data to species relative abundance
dd <- vegdist(ds, method = "bray")                   # compute dissimilarity between samples
ord <- monoMDS(dd, k = 5)                            # NMDS (non-metric multidimensional scaling)
dd.ord <- dist( scores(ord) )                        # Euclidean distance in reduced space

shepard_plot(dis = dd, ord = dd.ord)
shepard_plot(dis = dd, ord = scores(ord), k = 2)
shepard_plot(dis = dd, ord = scores(ord), k = 3)


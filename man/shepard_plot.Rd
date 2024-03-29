% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shepard_plot.R
\name{shepard_plot}
\alias{shepard_plot}
\title{Plot Shepard diagram}
\usage{
shepard_plot(dis, ord, k = 2, show = TRUE)
}
\arguments{
\item{dis}{Original dissimilarity matrix (class "dist")}

\item{ord}{Ordination scores (class "matrix" or "data.frame") or (class "dist")}

\item{k}{Number of dimensions (default, 2) that will be taken into account (only if the ordination scores were provided in 'ord')}

\item{show}{Logical; show the resulting plot on screen}
}
\value{
Invisible returns ggplot-object.
}
\description{
This function helps to estimate the representativeness of ordinations obtained using any reduced-space ordination method.
}
\details{
This function compares original dissimilarities among the objects with Euclidean distances in reduced space obtained with ordination (object scores or coordinates of the objects in the reduced space).
If the points are close to diagonal then the projection in reduced space accounts for a high fraction of the variance.
}
\examples{
library(vegan)
data(dune)
ds <- decostand(dune, method = "total", MARGIN = 1)  # standardize data to species relative abundance
dd <- vegdist(ds, method = "bray")                   # compute dissimilarity between samples
ord <- monoMDS(dd, k = 5)                            # NMDS (non-metric multidimensional scaling)
dd.ord <- dist( scores(ord) )                        # Euclidean distance in reduced space

shepard_plot(dis = dd, ord = dd.ord)               # Bray-Curtis vs distance in reduced space (k-dimensional)
shepard_plot(dis = dd, ord = scores(ord), k = 2)   # Bray-Curtis vs distance in reduced space (only first 2 dimensions)
shepard_plot(dis = dd, ord = scores(ord), k = 3)   # Bray-Curtis vs distance in reduced space (first 3 dimensions)

# Example with phyloseq wrappers
data("GlobalPatterns")
dd <- distance(GlobalPatterns, method = "unifrac", type = "samples")
ord <- ordinate(GlobalPatterns, method = "NMDS", distance = dd)
shepard_plot(dis = dd, ord = scores(ord))

}
\references{
Legendre P. & Legendre L. Numerical Ecology. 2012. 3rd Ed. - Chapter 9 (Ordination in reduced space).
}
\seealso{
\code{\link{vegan::stressplot}}, \code{\link{MAAS::Shepard}}
}

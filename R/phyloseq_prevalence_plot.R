
#' @title Prevalence plot (total OTU abundance vs OTU prevalence)
#' @description This function will plot total OTU abundance vs the fraction of samples in which an OTU is observed.
#' @param physeq A phyloseq-class object
#' @param prev.trh Add horizontal line with prevalence threshold (default is NULL, but  0.05 = 5\% of samples)
#' @param taxcolor Taxonomy rank for coloring the points (e.g. "Phylum")
#' @param facet Logical, split to separate panels by taxonomy rank used for coloring the points
#' @param point_alpha Point transparency value
#' @param showplot Logical, show plot on screen
#'
#' @return Plot of class 'ggplot'.
#'
#' @seealso \code{\link{phyloseq_filter_prevalence}}
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' # Subset data
#' GP <- subset_taxa(GlobalPatterns, Phylum %in% c("Acidobacteria", "Actinobacteria", "Firmicutes", "Verrucomicrobia"))
#' phyloseq_prevalence_plot(GP, taxcolor = "Phylum", facet = TRUE, point_alpha = 0.5, prev.trh = 0.05)
#'
phyloseq_prevalence_plot <- function(physeq, prev.trh = NULL, taxcolor = NULL, facet = FALSE, point_alpha = 0.7, showplot = T){

  ## Compute prevalence of each species
  prevdf <- prevalence(physeq)
  prevdf$PrevFrac <- with(prevdf, Prevalence / phyloseq::nsamples(physeq))

  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot2::ggplot(prevdf, ggplot2::aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot2::ggplot(prevdf, ggplot2::aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }

  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + ggplot2::geom_hline(yintercept = prev.trh, alpha = 0.5, linetype = 2) }

  pp <- pp +
      ggplot2::geom_point(size = 2, alpha = point_alpha) +
      ggplot2::scale_x_log10() +
      ggplot2::xlab("Total Abundance") +
      ggplot2::ylab("Prevalence [Frac. Samples]") +
      ggplot2::theme(legend.position="none")

  if(facet == TRUE){ pp <- pp + ggplot2::facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}

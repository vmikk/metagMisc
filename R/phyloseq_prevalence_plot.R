
#' @title Prevalence plot (total OTU abundance vs OTU prevalence).
#' @description This function will plot total OTU abundance vs the fraction of samples in which an OTU is observed.
#' @param physeq Phyloseq object
#' @param prev.trh Add horizontal line with prevalence threshold (default is NULL, but  0.05 = 5\% of samples)
#' @param taxcolor Taxonomy rank for coloring the points (e.g. "Phylum")
#' @param facet Logical, split to separate panels by taxonomy rank used for coloring the points
#' @param point_alpha Point transparency value
#' @param showplot Logical, show plot on screen
#'
#' @return Plot of class 'ggplot'.
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

  require(ggplot2)

  ## Compute prevalence of each species
  prevdf <- prevalence(physeq)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(physeq))

  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }

  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.5, linetype = 2) }

  pp <- pp +
      geom_point(size = 2, alpha = point_alpha) +
      scale_x_log10() +
      xlab("Total Abundance") +
      ylab("Prevalence [Frac. Samples]") +
      theme(legend.position="none")

  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}


#' @title Compute prevalence of each species.
#'
#' @param physeq Phyloseq object
#' @param add_tax Logical, add taxonomy to the results
#'
#' @return Data frame
#'
#' @examples
#' data(esophagus)
#' prevalence(esophagus)
#'
#' data(GlobalPatterns)
#' head( prevalence(GlobalPatterns, add_tax = T) )
#'
prevalence <- function(physeq, add_tax = TRUE){
  prevdf <- apply(X = otu_table(physeq),
                  MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})

  ## Add taxonomy and total read counts to this data.frame
  prevdf <- data.frame(Prevalence = prevdf,
                       TotalAbundance = taxa_sums(physeq))

  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(physeq))
  }
  return(prevdf)
}


phyloseq_prevalence_plot <- function(physeq, prev.trh = NULL, taxcolor = NULL, facet = FALSE, point_alpha = 0.7, showplot = T){

  require(ggplot2)
  
  ## Compute prevalence of each species, store as data.frame
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


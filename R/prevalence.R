
#' @title Compute prevalence and abundance summary of each species/OTU.
#'
#' @param physeq A phyloseq-class object
#' @param add_tax Logical, add taxonomy to the results
#' @param package Which package to use for data processing ("base" or "data.table")
#'
#' @details
#' This function can be executed using functions from base R
#' (specified by setting the package argument to "base")
#' or the data.table package
#' (specified by setting the package argument to "data.table").
#' By default, it uses the data.table package, which is much faster and more efficient;
#' for instance, on the GlobalPatterns dataset, the implementation using data.table
#' is 7 times faster than that using base R.
#' E.g.:
#' # microbenchmark( \cr
#' #   dt = prevalence(GlobalPatterns, package = "data.table"), \cr
#' #   bb = prevalence(GlobalPatterns, package = "base"), \cr
#' #   times = 10) \cr
#'
#' @return
#' Data frame or data.table with the following columns:
#' \itemize{
#'  \item{"Taxa"}          {Taxa names (species, OTU, or ASV)}
#'  \item{"Prevalence"}    {The number of occurrences (the count of samples where the particular taxon is found)}
#'  \item{"TotalAbundance"}{The total abundance of a taxon across all samples}
#'  \item{"MeanAbundance"} {The average taxon abundance}
#'  \item{"MedianAbundance"}{The median abundance of a taxon}
#'  \item{"Taxonomy"}       {An optional set of columns containing taxonomic information present in phyloseq object}
#' }
#'
#' @export
#'
#' @examples
#' data(esophagus)
#' prevalence(esophagus)
#'
#' data(GlobalPatterns)
#' head( prevalence(GlobalPatterns, add_tax = T) )
#'
prevalence <- function(physeq, add_tax = TRUE, package = "data.table"){

  ## Check if taxa are rows
  trows <- taxa_are_rows(physeq)

  ## Extract OTU table
  otutab <- as.data.frame(otu_table(physeq))

  ## Process data with base R
  if(package %in% "base"){

    ## Transpose OTU table (species should be arranged by rows)
    if(trows == FALSE){
      otutab <- t(otutab)
    }

    ## Estimate prevalence (number of samples with OTU present)
    prevdf <- apply(X = otutab,
                    # MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
                    MARGIN = 1,
                    FUN = function(x){sum(x > 0)})

    ## Add total and average read counts per OTU
    prevdf <- data.frame(
      Taxa            = rownames(otutab),
      Prevalence      = prevdf,
      TotalAbundance  = taxa_sums(physeq),
      MeanAbundance   = rowMeans(otutab),
      MedianAbundance = apply(otutab, 1, median))

    ## Add taxonomy table
    if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
      prevdf <- cbind(prevdf, tax_table(physeq))
    }

  } # end of `base` package


  ## Process data using `data.table` package
  if(package %in% "data.table"){

    setDT(otutab)

    ## Transpose OTU table (species should be arranged by rows)
    if(trows == FALSE){
      otutab <- t(otutab)
    }

    ## Row-wise medians
    ## Based on https://stackoverflow.com/a/48885574  by Jaap Walhout
    meds <- melt(otutab, measure.vars = names(otutab))[, r := 1:.N, variable][, median(value), by = r]$V1

    ## Add total and average read counts per OTU
    prevdf <- data.table(
      Taxa            = taxa_names(physeq),
      Prevalence      = rowSums(otutab > 0, na.rm = TRUE),
      TotalAbundance  = rowSums(otutab, na.rm = TRUE),
      MeanAbundance   = rowMeans(otutab, na.rm = TRUE),
      MedianAbundance = meds)

    ## Add taxonomy table
    if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){

      prevdf <- cbind(prevdf, as.data.frame(tax_table(physeq)))
    }

  }

  return(prevdf)
}
